use std::{
    io::{BufRead, Error, ErrorKind},
    rc::Rc,
};

use clap::{Parser, Subcommand};
use datafusion::{
    arrow::{
        array::DictionaryArray,
        datatypes::{DataType, Schema, UInt8Type},
    },
    common::JoinType,
    functions_aggregate::{
        count::count,
        expr_fn::{array_agg, bit_or},
        min_max::min,
        string_agg::string_agg,
    },
    functions_window::expr_fn::row_number,
    prelude::{
        DataFrame, SessionContext, abs, array_concat, array_distinct, case, cast, coalesce, col,
        concat, lit, make_array, sha256, unnest,
    },
};
use noodles::vcf::{self, Header};
use svelt::{chroms::ChromSet, tables::load_vcf_core, vcf_reader::VcfReader};

/// Structuaral Variant (SV) VCF merging
#[derive(Debug, Parser)]
#[command(name = "svelt")]
#[command(about = "Merge structural variants, aligning similar ones.", long_about = None)]
struct Cli {
    #[command(subcommand)]
    command: Commands,
}

#[derive(Debug, Subcommand)]
enum Commands {
    /// Apply annotations from one VCF to another
    #[command(arg_required_else_help = true)]
    Merge {
        /// The output filename
        #[arg(short, long)]
        out: String,

        /// SV VCF files to merge
        #[arg(num_args(1..))]
        vcf: Vec<String>,
    },
}

fn load_chroms(path: &str) -> std::io::Result<ChromSet> {
    let reader = autocompress::autodetect_open(path)?;
    let mut reader: vcf::io::Reader<Box<dyn BufRead>> =
        vcf::io::reader::Builder::default().build_from_reader(reader)?;
    let header: Header = reader.read_header()?;

    let mut names: Vec<&str> = Vec::new();
    for contig in header.contigs() {
        let chrom = contig.0;
        names.push(chrom);
    }
    let chroms = ChromSet::from(names.as_ref());
    Ok(chroms)
}

#[tokio::main]
async fn main() -> std::io::Result<()> {
    env_logger::init();

    let cli = Cli::parse();

    match cli.command {
        Commands::Merge { out, vcf } => {
            let chroms = load_chroms(&vcf[0])?;
            let chroms = Rc::new(chroms);
            let mut readers = Vec::new();
            for vcf in vcf.iter() {
                let reader = VcfReader::new(vcf, chroms.clone())?;
                readers.push(reader);
            }
            let n = readers.len();

            let ctx = SessionContext::new();

            let mut acc: Option<DataFrame> = None;
            for vix in 0..readers.len() {
                log::info!("reading {}", readers[vix].path);
                let reader: &mut VcfReader = &mut readers[vix];
                let records = load_vcf_core(reader)?;
                let df = ctx
                    .read_batch(records)
                    .map_err(|e| Error::new(ErrorKind::Other, e))?;
                let df = df
                    .with_column("vix", lit(1u32 << vix))?
                    .with_column("row_id", lit(vix as u32) + (col("row_num") * lit(100)))?;

                if let Some(df0) = acc {
                    let df = df0.union(df)?;
                    acc = Some(df);
                } else {
                    acc = Some(df)
                }
            }
            let orig = acc.unwrap();

            let exact = orig
                .clone()
                .aggregate(
                    vec![
                        col("chrom_id"),
                        col("start"),
                        col("end"),
                        col("kind"),
                        col("length"),
                        col("chrom2_id"),
                        col("end2"),
                        col("seq_hash"),
                    ],
                    vec![
                        min(col("row_id")).alias("row_key"),
                        count(col("vix")).alias("vix_count"),
                        bit_or(col("vix")).alias("vix_set"),
                    ],
                )?
                .select(vec![
                    col("chrom_id").alias("exact_chrom_id"),
                    col("start").alias("exact_start"),
                    col("end").alias("exact_end"),
                    col("kind").alias("exact_kind"),
                    col("length").alias("exact_length"),
                    col("chrom2_id").alias("exact_chrom2_id"),
                    col("end2").alias("exact_end2"),
                    col("seq_hash").alias("exact_seq_hash"),
                    col("row_key").alias("exact_row_key"),
                    col("vix_count").alias("exact_vix_count"),
                    col("vix_set").alias("exact_vix_set"),
                ])?
                .filter(col("exact_vix_count").eq(lit(n as u32)))?;

            //orig.clone().show().await?;
            //exact.clone().show().await?;

            let exact_ins = exact.clone().filter(col("exact_kind").eq(lit("INS")))?;
            let exact_del = exact.clone().filter(col("exact_kind").eq(lit("DEL")))?;
            let exact_dup = exact.clone().filter(col("exact_kind").eq(lit("DUP")))?;
            let exact_inv = exact.clone().filter(col("exact_kind").eq(lit("INV")))?;
            let exact_bnd = exact.clone().filter(col("exact_kind").eq(lit("BND")))?;

            exact_bnd.clone().show().await?;

            let results = orig
                .join(
                    exact_ins,
                    JoinType::Left,
                    &["chrom_id", "start", "end", "kind", "length", "seq_hash"],
                    &[
                        "exact_chrom_id",
                        "exact_start",
                        "exact_end",
                        "exact_kind",
                        "exact_length",
                        "exact_seq_hash",
                    ],
                    None,
                )?
                .drop_columns(&[
                    "exact_chrom_id",
                    "exact_start",
                    "exact_end",
                    "exact_kind",
                    "exact_length",
                    "exact_chrom2_id",
                    "exact_end2",
                    "exact_seq_hash",
                ])?
                .with_column_renamed("exact_row_key", "row_key")?
                .with_column_renamed("exact_vix_count", "vix_count")?
                .with_column_renamed("exact_vix_set", "vix_set")?
                .join(
                    exact_del,
                    JoinType::Left,
                    &["chrom_id", "start", "end", "kind", "length"],
                    &[
                        "exact_chrom_id",
                        "exact_start",
                        "exact_end",
                        "exact_kind",
                        "exact_length",
                    ],
                    None,
                )?
                .drop_columns(&[
                    "exact_chrom_id",
                    "exact_start",
                    "exact_end",
                    "exact_kind",
                    "exact_length",
                    "exact_chrom2_id",
                    "exact_end2",
                    "exact_seq_hash",
                ])?
                .with_column(
                    "row_key",
                    coalesce(vec![col("row_key"), col("exact_row_key")]),
                )?
                .with_column(
                    "vix_count",
                    coalesce(vec![col("vix_count"), col("exact_vix_count")]),
                )?
                .with_column(
                    "vix_set",
                    coalesce(vec![col("vix_set"), col("exact_vix_set")]),
                )?
                .drop_columns(&["exact_row_key", "exact_vix_count", "exact_vix_set"])?
                .join(
                    exact_dup,
                    JoinType::Left,
                    &["chrom_id", "start", "end", "kind", "length"],
                    &[
                        "exact_chrom_id",
                        "exact_start",
                        "exact_end",
                        "exact_kind",
                        "exact_length",
                    ],
                    None,
                )?
                .drop_columns(&[
                    "exact_chrom_id",
                    "exact_start",
                    "exact_end",
                    "exact_kind",
                    "exact_length",
                    "exact_chrom2_id",
                    "exact_end2",
                    "exact_seq_hash",
                ])?
                .with_column(
                    "row_key",
                    coalesce(vec![col("row_key"), col("exact_row_key")]),
                )?
                .with_column(
                    "vix_count",
                    coalesce(vec![col("vix_count"), col("exact_vix_count")]),
                )?
                .with_column(
                    "vix_set",
                    coalesce(vec![col("vix_set"), col("exact_vix_set")]),
                )?
                .drop_columns(&["exact_row_key", "exact_vix_count", "exact_vix_set"])?
                .join(
                    exact_inv,
                    JoinType::Left,
                    &["chrom_id", "start", "end", "kind", "length"],
                    &[
                        "exact_chrom_id",
                        "exact_start",
                        "exact_end",
                        "exact_kind",
                        "exact_length",
                    ],
                    None,
                )?
                .drop_columns(&[
                    "exact_chrom_id",
                    "exact_start",
                    "exact_end",
                    "exact_kind",
                    "exact_length",
                    "exact_chrom2_id",
                    "exact_end2",
                    "exact_seq_hash",
                ])?
                .with_column(
                    "row_key",
                    coalesce(vec![col("row_key"), col("exact_row_key")]),
                )?
                .with_column(
                    "vix_count",
                    coalesce(vec![col("vix_count"), col("exact_vix_count")]),
                )?
                .with_column(
                    "vix_set",
                    coalesce(vec![col("vix_set"), col("exact_vix_set")]),
                )?
                .drop_columns(&["exact_row_key", "exact_vix_count", "exact_vix_set"])?
.join(
                    exact_bnd,
                    JoinType::Left,
                    &["chrom_id", "start", "end", "kind", "chrom2_id", "end2"],
                    &[
                        "exact_chrom_id",
                        "exact_start",
                        "exact_end",
                        "exact_kind",
                        "exact_chrom2_id",
                        "exact_end2"
                    ],
                    None,
                )?
                .drop_columns(&[
                    "exact_chrom_id",
                    "exact_start",
                    "exact_end",
                    "exact_kind",
                    "exact_length",
                    "exact_chrom2_id",
                    "exact_end2",
                    "exact_seq_hash",
                ])?
                .with_column(
                    "row_key",
                    coalesce(vec![col("row_key"), col("exact_row_key")]),
                )?
                .with_column(
                    "vix_count",
                    coalesce(vec![col("vix_count"), col("exact_vix_count")]),
                )?
                .with_column(
                    "vix_set",
                    coalesce(vec![col("vix_set"), col("exact_vix_set")]),
                )?
                .drop_columns(&["exact_row_key", "exact_vix_count", "exact_vix_set"])?                ;

            results
                .clone()
                .sort_by(vec![col("chrom_id"), col("start"), col("end"), col("row_id")])?
                .show()
                .await?;

            /*


            .join(
                exact_del,
                JoinType::Left,
                &["chrom_id", "start", "end", "kind", "length"],
                &[
                    "exact_chrom_id",
                    "exact_start",
                    "exact_end",
                    "exact_kind",
                    "exact_length",
                ],
                None,
            )?
            .select(vec![
                col("chrom_id"),
                col("start"),
                col("end"),
                col("kind"),
                col("length"),
                col("chrom2_id"),
                col("end2"),
                col("seq_hash"),
                coalesce(vec![col("row_key"), col("exact_row_key")]).alias("row_key"),
                col("vix_count"),
                col("vix_set"),
            ])?

            let exact = df
                .clone()
                .filter(col("vix_count").eq(lit(n as i32)))?
                .select(vec![
                    col("chrom_id"),
                    col("start"),
                    col("end"),
                    col("kind"),
                    col("rows"),
                    col("vix_set"),
                ])?;

            let df = df.filter(col("vix_count").not_eq(lit(n as i32)))?;

            let starts = df
                .clone()
                .select(vec![
                    col("chrom_id"),
                    col("start"),
                    col("end"),
                    col("kind"),
                    col("length"),
                    col("rows"),
                    col("vix_set"),
                ])?
                .sort(vec![
                    col("chrom_id").sort(true, true),
                    col("start").sort(true, true),
                ])?;
            //.collect()
            //.await?;

            let lhs = starts
                .clone()
                .with_column_renamed("chrom_id", "lhs_chrom_id")?
                .with_column_renamed("start", "lhs_start")?
                .with_column_renamed("end", "lhs_end")?
                .with_column_renamed("kind", "lhs_kind")?
                .with_column_renamed("length", "lhs_length")?
                .with_column_renamed("rows", "lhs_rows")?
                .with_column_renamed("vix_set", "lhs_vix_set")?;
            let rhs = starts
                .clone()
                .with_column_renamed("chrom_id", "rhs_chrom_id")?
                .with_column_renamed("start", "rhs_start")?
                .with_column_renamed("end", "rhs_end")?
                .with_column_renamed("kind", "rhs_kind")?
                .with_column_renamed("length", "rhs_length")?
                .with_column_renamed("rows", "rhs_rows")?
                .with_column_renamed("vix_set", "rhs_vix_set")?;

            let almost_exact = lhs
                .join(
                    rhs,
                    datafusion::common::JoinType::Inner,
                    &["lhs_chrom_id", "lhs_kind"],
                    &["rhs_chrom_id", "rhs_kind"],
                    Some(
                        col("lhs_start")
                            .lt_eq(col("rhs_start"))
                            .and(col("lhs_vix_set").lt(col("rhs_vix_set")))
                            .and(abs(col("lhs_start") - col("rhs_start")).lt(lit(25)))
                            .and(abs(col("lhs_end") - col("rhs_end")).lt(lit(25)))
                            .and((col("lhs_vix_set") & col("rhs_vix_set")).eq(lit(0))),
                    ),
                )?
                .with_column(
                    "max_length",
                    case(col("lhs_length").gt_eq(col("rhs_length")))
                        .when(lit(true), col("lhs_length"))
                        .when(lit(false), col("rhs_length"))
                        .end()?,
                )?
                .with_column(
                    "min_length",
                    case(col("lhs_length").lt_eq(col("rhs_length")))
                        .when(lit(true), col("lhs_length"))
                        .when(lit(false), col("rhs_length"))
                        .end()?,
                )?
                .with_column(
                    "length_ratio",
                    cast(col("min_length"), DataType::Float64)
                        / cast(col("max_length"), DataType::Float64),
                )?
                .filter(col("length_ratio").gt(lit(0.9)))?
                .select(vec![
                    col("lhs_chrom_id").alias("chrom_id"),
                    col("lhs_start").alias("start"),
                    col("lhs_end").alias("end"),
                    col("lhs_kind").alias("kind"),
                    array_distinct(array_concat(vec![col("lhs_rows"), col("rhs_rows")]))
                        .alias("rows"),
                    (col("lhs_vix_set") | col("rhs_vix_set")).alias("vix_set"),
                ])?;

            let results = almost_exact.collect().await?; //.union(almost_exact)?.collect().await?;
            let df = ctx.read_batches(results)?;

            let res = df
                .sort(vec![
                    col("chrom_id").sort(true, true),
                    col("start").sort(true, true),
                    //col("end").sort(true, true),
                    //col("kind").sort(true, true),
                ])?
                .show()
                .await;
            match res {
                Ok(_) => {}
                Err(err) => {
                    log::error!("{}", err.message());
                },
            }
             */
        }
    }
    Ok(())
}
