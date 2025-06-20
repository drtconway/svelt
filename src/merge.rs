use std::{
    io::{BufRead, Error, ErrorKind},
    rc::Rc,
};

use datafusion::{
    arrow::{
        array::{BooleanArray, GenericStringArray, Int64Array, RecordBatch},
        datatypes::DataType,
    },
    config::CsvOptions,
    dataframe::DataFrameWriteOptions,
    prelude::{DataFrame, cast, col, lit, nullif, to_hex},
};
use noodles::{
    fasta::{self, repository::adapters::IndexedReader},
    vcf::{self, Header, Record, header::SampleNames},
};

use crate::{
    almost::find_almost_exact,
    backward::find_backwards_bnds,
    chroms::ChromSet,
    construct::{MergeBuilder, add_svelt_header_fields},
    errors::as_io_error,
    exact::find_exact,
    options::{CommonOptions, MergeOptions, make_session_context},
    record_seeker::RecordSeeker,
    row_key::RowKey,
    tables::load_vcf_core,
    vcf_reader::VcfReader,
};

pub async fn merge_vcfs(
    out: &str,
    vcf: &Vec<String>,
    options: Rc<MergeOptions>,
    common: &CommonOptions,
) -> std::io::Result<()> {
    options.check().map_err(as_io_error)?;

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
    let orig = orig
        .with_column("row_key", col("row_id"))?
        .with_column("vix_count", lit(1))?
        .with_column("vix_set", cast(col("vix"), DataType::Int64))?
        .with_column("flip", lit(false))?
        .with_column("criteria", nullif(lit(1), lit(1)))?;

    let results = find_exact(&ctx, orig, n as u32, options.as_ref()).await?;

    let results = find_almost_exact(&ctx, results, n as u32, options.as_ref()).await?;

    let mut results = results;
    if options.allow_breakend_flipping {
        results = find_backwards_bnds(&ctx, results, n as u32, &options).await?;
    }

    let mut reference = None;
    if let Some(reference_filename) = &options.reference {
        let reference_reader =
            fasta::io::indexed_reader::Builder::default().build_from_path(reference_filename)?;
        let adapter = IndexedReader::new(reference_reader);
        reference = Some(Rc::new(fasta::Repository::new(adapter)));
    }

    if let Some(table_out) = &options.write_merge_table {
        let opts = DataFrameWriteOptions::default();
        let csv_opts = CsvOptions::default().with_delimiter(b'\t');
        let csv_opts = Some(csv_opts);
        results
            .clone()
            .with_column("seq_hash", to_hex(col("seq_hash")))?
            .sort_by(vec![
                col("chrom_id"),
                col("start"),
                col("end"),
                col("row_key"),
                col("row_id"),
            ])?
            .drop_columns(&["row_num", "chrom_id", "chrom2_id"])?
            .write_csv(table_out, opts, csv_opts)
            .await?;
    }

    let table = results
        .sort_by(vec![
            col("chrom_id"),
            col("start"),
            col("end"),
            col("row_key"),
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
    let mut current_row_alts: Vec<Option<String>> = (0..n).into_iter().map(|_| None).collect();
    let mut current_row_flip = false;
    let mut current_row_criteria = String::new();

    for recs in table.into_iter() {
        for field in recs.schema_ref().fields().iter() {
            log::debug!("field: {:?}", field);
        }
        let row_ids = get_array::<Int64Array>(&recs, "row_id");
        let row_keys = get_array::<Int64Array>(&recs, "row_key");
        let alt_seqs = get_array::<GenericStringArray<i32>>(&recs, "alt_seq");
        let flips = get_array::<BooleanArray>(&recs, "flip");
        let criteria = get_array::<GenericStringArray<i32>>(&recs, "criteria");

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
                    builder.construct(
                        recs,
                        &vix_samples,
                        &current_row_alts,
                        current_row_flip,
                        &current_row_criteria,
                    )?;
                }

                current_row = (0..n).into_iter().map(|_| None).collect();
                current_row_key = row_key;
                current_row_alts = (0..n).into_iter().map(|_| None).collect();
                current_row_flip = false;
                current_row_criteria = String::new();
            }

            let (vix, rn) = RowKey::decode(row_id);
            current_row[vix as usize] = Some(rn);

            let alt_seq = alt_seqs.value(i);
            if alt_seq.len() > 0 {
                current_row_alts[vix as usize] = Some(String::from(alt_seq));
            }

            current_row_flip |= flips.value(i);

            let crit = criteria.value(i);
            if crit.len() > current_row_criteria.len() {
                current_row_criteria = String::from(crit);
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
        builder.construct(
            recs,
            &vix_samples,
            &current_row_alts,
            current_row_flip,
            &current_row_criteria,
        )?;
    }

    Ok(())
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

fn get_array<'a, Type: 'static>(recs: &'a RecordBatch, name: &str) -> &'a Type {
    recs.column_by_name(name)
        .unwrap()
        .as_any()
        .downcast_ref::<Type>()
        .unwrap()
}
