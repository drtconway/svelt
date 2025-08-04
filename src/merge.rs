use std::{
    io::{BufRead, Error, ErrorKind},
    rc::Rc,
};

use datafusion::{
    arrow::{
        array::{
            Array, BooleanArray, GenericStringArray, Int64Array, RecordBatch, StringArray,
            UInt32Array,
        },
        datatypes::DataType,
    },
    common::JoinType,
    functions_aggregate::expr_fn::first_value,
    prelude::{DataFrame, cast, col, concat, length, lit, nullif, to_hex},
};
use noodles::{
    fasta::{self, repository::adapters::IndexedReader},
    vcf::{self, Header, Record, header::SampleNames},
};

use crate::{
    breakends::unpaired_breakend_check,
    chroms::ChromSet,
    construct::{MergeBuilder, add_svelt_header_fields},
    errors::{Context, FileContext, SveltError, as_io_error},
    merge::{
        approx::{approx_bnd_here_there_join, approx_bnd_there_here_join, approx_near_join},
        exact::{full_exact_bnd, full_exact_indel_join, full_exact_locus_ins_join},
        report::produce_reporting_table,
        union::merge_with,
        variant_id::construct_variant_ids,
    },
    options::{CommonOptions, MergeOptions, make_session_context},
    record_seeker::RecordSeeker,
    row_key::RowKey,
    tables::load_vcf_core,
    vcf_reader::VcfReader,
};

mod approx;
mod classify;
mod exact;
mod report;
mod union;
mod variant_id;

pub async fn merge_vcfs(
    out: &str,
    vcf: &Vec<String>,
    options: Rc<MergeOptions>,
    common: &CommonOptions,
) -> std::io::Result<()> {
    options.check().map_err(as_io_error)?;

    if vcf.len() > 64 {
        log::error!(
            "svelt can only merge up to 64 VCF files at a time ({} given)",
            vcf.len()
        );
        return Err(as_io_error(SveltError::TooManyVcfs(vcf.len())));
    }

    let chroms = load_chroms(&vcf[0])?;
    let chroms = Rc::new(chroms);
    let mut readers = Vec::new();
    for vcf in vcf.iter() {
        let reader = VcfReader::new(vcf, chroms.clone())?;
        readers.push(reader);
    }
    let n = readers.len();

    let ctx = make_session_context(common);

    let mut acc: Option<DataFrame> = None;
    for vix in 0..readers.len() {
        log::info!("reading {}", readers[vix].path);
        let reader: &mut VcfReader = &mut readers[vix];
        let records = load_vcf_core(reader)?;
        let df = ctx
            .read_batch(records)
            .map_err(|e| Error::new(ErrorKind::Other, e))?;
        let df = df
            .with_column("vix", lit(1u64 << vix))?
            .with_column("row_id", RowKey::make(col("row_num"), vix as u32))?;

        if let Some(df0) = acc {
            let df = df0.union(df)?;
            acc = Some(df);
        } else {
            acc = Some(df)
        }
    }
    let orig = acc.unwrap();
    let orig = orig
        .with_column("row_key", cast(col("row_id"), DataType::UInt32))?
        .with_column("vix_count", lit(1))?
        .with_column("vix_set", col("vix"))?
        .with_column("criteria", nullif(lit(""), lit("")))?;

    let mut results = orig.clone();

    if true {
        log::info!("looking for exact matches on indel type variants");
        let join = full_exact_indel_join(results.clone(), n)?;
        results = merge_with(results, join, &ctx, "exact").await?;
    }
    if true {
        log::info!("looking for almost exact matches on insertions");
        let join = full_exact_locus_ins_join(results.clone(), n, &options)?;
        results = merge_with(results, join, &ctx, "locus").await?;
    }
    if true {
        log::info!("looking for exact matches on breakends");
        let join = full_exact_bnd(results.clone(), n)?;
        results = merge_with(results, join, &ctx, "exact").await?;
    }
    if true {
        log::info!("looking for approximate matches on breakends (here-there)");
        let join = approx_bnd_here_there_join(results.clone(), n, &options)?;
        results = merge_with(results, join, &ctx, "here").await?;
    }
    if true {
        log::info!("looking for approximate matches on breakends (there-here)");
        let join = approx_bnd_there_here_join(results.clone(), n, &options)?;
        results = merge_with(results, join, &ctx, "there").await?;
    }
    if true {
        log::info!("looking for nearby matches on indel type variants");
        let join = approx_near_join(results.clone(), n, &options, &ctx).await?;
        results = merge_with(results, join, &ctx, "near").await?;
    }

    results = unpaired_breakend_check(results).await?;

    let mut reference = None;
    if let Some(reference_filename) = &options.reference {
        let reference_reader =
            fasta::io::indexed_reader::Builder::default().build_from_path(reference_filename)?;
        let adapter = IndexedReader::new(reference_reader);
        reference = Some(Rc::new(fasta::Repository::new(adapter)));
    }

    // Annotation needs the seq_hash as a string.
    results = results.with_column("seq_hash", to_hex(col("seq_hash")))?;

    let mut annot = false;
    if let Some(features) = &options.annotate_insertions {
        let ins = results
            .clone()
            .select(vec![
                col("seq_hash"),
                length(col("alt_seq")).alias("length"),
                col("alt_seq"),
            ])?
            .filter(length(col("alt_seq")).gt(lit(0)))?
            .distinct()?;

        let ins = ins.collect().await?;
        let classifications = classify::find_classifications(ins, features, &ctx).await?;
        if false {
            classifications
                .clone()
                .sort_by(vec![col("distance")])?
                .show()
                .await?;
        }

        results = results
            .join(
                classifications,
                JoinType::Left,
                &["seq_hash"],
                &["query_name"],
                None,
            )?
            .drop_columns(&["query_name"])?
            .with_column("class", concat(vec![lit(""), col("class")]))?
            .with_column("strand", concat(vec![lit(""), col("strand")]))?;
        annot = true;
    }

    results = construct_variant_ids(results, &ctx).await?;

    results = add_primary_cols(results)?;

    if let Some(table_out) = &options.write_merge_table {
        produce_reporting_table(results.clone(), &table_out).await?;
    }

    let table = results
        .sort_by(vec![
            col("chrom_id"),
            col("primary_start"),
            col("primary_end"),
            col("row_key"),
            col("row_id"),
        ])?
        .collect()
        .await?;

    let mut sample_names: Vec<String> = Vec::new();
    let mut vix_samples = Vec::new();

    for h in readers.iter().map(|r| &r.header) {
        for s in h.sample_names() {
            sample_names.push(s.clone());
        }
        vix_samples.push(h.sample_names().len());
    }

    let mut seekers = Vec::new();
    for path in vcf.iter() {
        let seeker = RecordSeeker::new(path, chroms.clone())?;
        seekers.push(seeker);
    }

    let mut header = readers[0].header.clone();
    *header.sample_names_mut() = SampleNames::from_iter(sample_names.iter().map(|s| s.clone()));
    add_svelt_header_fields(&mut header, &options.unwanted_info)?;

    let mut builder = MergeBuilder::new(out, options, header, reference)?;

    let mut current_row_key = u32::MAX;
    let mut current_row: Vec<Option<u32>> = (0..n).into_iter().map(|_| None).collect();
    let mut current_row_ids: Vec<String> = (0..n).into_iter().map(|_| String::new()).collect();
    let mut current_row_alts: Vec<Option<String>> = (0..n).into_iter().map(|_| None).collect();
    let mut current_row_paired_bnd = false;
    let mut current_row_criteria = String::new();
    let mut current_row_classification = None;

    for recs in table.into_iter() {
        if false {
            for field in recs.schema_ref().fields().iter() {
                log::info!("field: {:?}", field);
            }
        }
        let variant_ids = get_array::<GenericStringArray<i32>>(&recs, "variant_id");
        let row_ids = get_array::<Int64Array>(&recs, "row_id");
        let row_keys = get_array::<UInt32Array>(&recs, "row_key");
        let alt_seqs = get_array::<GenericStringArray<i32>>(&recs, "alt_seq");
        let paired_bnds = get_array::<BooleanArray>(&recs, "paired_bnd");
        let criteria = get_array::<GenericStringArray<i32>>(&recs, "criteria");
        let classifications = if annot {
            let class = get_array::<StringArray>(&recs, "class");
            let strand = get_array::<StringArray>(&recs, "strand");
            Some((class, strand))
        } else {
            None
        };

        for i in 0..row_ids.len() {
            let row_id = row_ids.value(i) as u32;
            let row_key = row_keys.value(i) as u32;

            if row_key != current_row_key {
                log::debug!("flushing group: {} {:?}", current_row_key, current_row);
                let mut is_empty = true;
                let mut recs: Vec<Option<(Rc<Header>, Record)>> =
                    (0..n).into_iter().map(|_| None).collect();
                for vix in 0..n {
                    if let Some(rn) = current_row[vix] {
                        is_empty = false;
                        let hnr = seekers[vix].take(rn)?.unwrap();
                        recs[vix] = Some(hnr)
                    }
                }
                if !is_empty {
                    let feat = if let Some((class, strand)) = &current_row_classification {
                        format!("{}{}", class, strand)
                    } else {
                        String::new()
                    };

                    builder.construct(
                        recs,
                        &vix_samples,
                        &current_row_ids,
                        &current_row_alts,
                        current_row_paired_bnd,
                        &current_row_criteria,
                        &feat,
                    )?;
                }

                current_row = (0..n).into_iter().map(|_| None).collect();
                current_row_key = row_key;
                current_row_ids = (0..n).into_iter().map(|_| String::new()).collect();
                current_row_alts = (0..n).into_iter().map(|_| None).collect();
                current_row_paired_bnd = false;
                current_row_criteria = String::new();
                current_row_classification = None;
            }

            let (vix, rn) = RowKey::decode(row_id);
            current_row[vix as usize] = Some(rn);

            let variant_id = String::from(variant_ids.value(i));
            current_row_ids[vix as usize] = variant_id;

            let alt_seq = alt_seqs.value(i);
            if alt_seq.len() > 0 {
                current_row_alts[vix as usize] = Some(String::from(alt_seq));
            }

            current_row_paired_bnd |= paired_bnds.value(i);

            let crit = criteria.value(i);
            if crit.len() > current_row_criteria.len() {
                current_row_criteria = String::from(crit);
            }

            if let Some((class, strand)) = &classifications {
                let cls = class.value(i);
                if cls.len() > 0 {
                    current_row_classification =
                        Some((String::from(class.value(i)), String::from(strand.value(i))));
                }
            }
        }
    }
    log::debug!("flushing group: {} {:?}", current_row_key, current_row);
    let mut is_empty = true;
    let mut recs: Vec<Option<(Rc<Header>, Record)>> = (0..n).into_iter().map(|_| None).collect();
    for vix in 0..n {
        if let Some(rn) = current_row[vix] {
            is_empty = false;
            let hnr = seekers[vix].take(rn)?.unwrap();
            recs[vix] = Some(hnr)
        }
    }
    if !is_empty {
        let feat = if let Some((class, strand)) = &current_row_classification {
            format!("{}{}", class, strand)
        } else {
            String::new()
        };

        builder.construct(
            recs,
            &vix_samples,
            &current_row_ids,
            &current_row_alts,
            current_row_paired_bnd,
            &current_row_criteria,
            &feat,
        )?;
    }

    Ok(())
}

fn load_chroms(path: &str) -> std::io::Result<ChromSet> {
    FileContext::new(path).with(|| {
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
    })
}

fn get_array<'a, Type: 'static>(recs: &'a RecordBatch, name: &str) -> &'a Type {
    if false {
        log::info!("getting {}", name);
    }
    recs.column_by_name(name)
        .unwrap()
        .as_any()
        .downcast_ref::<Type>()
        .unwrap()
}

fn add_primary_cols(tbl: DataFrame) -> std::io::Result<DataFrame> {
    let rhs = tbl
        .clone()
        .aggregate(
            vec![col("row_key")],
            vec![
                first_value(col("row_id"), vec![col("vix").sort(true, false)])
                    .alias("left_row_id"),
            ],
        )?
        .drop_columns(&["row_key"])?;

    let primary = tbl
        .clone()
        .join(rhs, JoinType::Inner, &["row_id"], &["left_row_id"], None)?
        .select(vec![
            col("row_key").alias("primary_row_key"),
            col("start").alias("primary_start"),
            col("end").alias("primary_end"),
            col("end2").alias("primary_end2"),
            col("length").alias("primary_length"),
        ])?;

    let tbl = tbl
        .clone()
        .join(
            primary,
            JoinType::Left,
            &["row_key"],
            &["primary_row_key"],
            None,
        )?
        .drop_columns(&["primary_row_key"])?;

    Ok(tbl)
}
