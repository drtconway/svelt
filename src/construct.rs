use std::collections::HashSet;
use std::io::{Error, ErrorKind};
use std::rc::Rc;

use noodles::core::Position;
use noodles::fasta::Repository;
use noodles::vcf;
use noodles::vcf::header::record::value::map::Builder;
use noodles::vcf::header::record::value::map::info::{Number, Type};
use noodles::vcf::variant::record::{
    AlternateBases as AlternateBases_, Filters as Filters_, Ids as Ids_,
};
use noodles::vcf::variant::record_buf::info::field::Value as InfoValue;
use noodles::vcf::variant::record_buf::samples::Keys;
use noodles::vcf::variant::record_buf::samples::sample::Value;
use noodles::vcf::variant::record_buf::samples::sample::value::genotype::Allele;
use noodles::vcf::variant::record_buf::samples::sample::value::{Array, Genotype};
use noodles::vcf::variant::record_buf::{AlternateBases, Filters, Ids, Info, Samples};
use noodles::vcf::{Header, Record, variant::RecordBuf};
use vcf::variant::record_buf::info::field::value::Array as InfoArray;

use crate::breakends::BreakEnd;
use crate::tables::is_seq;

pub fn add_svelt_header_fields(
    header: &mut Header,
    unwanted_info: &Vec<String>,
) -> std::io::Result<()> {
    let infos = header.infos_mut();
    infos.insert(
        String::from("SVELT_CRITERIA"),
        Builder::default()
            .set_number(Number::Unknown)
            .set_type(Type::String)
            .set_description("The list of criteria that resulted in the merging of these variants.")
            .build()
            .map_err(|e| Error::new(ErrorKind::Other, e))?,
    );
    infos.insert(
        String::from("SVELT_ALT_SEQ"),
        Builder::default()
            .set_number(Number::Unknown)
            .set_type(Type::String)
            .set_description("The list of alt sequences that were replaced with the ALT tag.")
            .build()
            .map_err(|e| Error::new(ErrorKind::Other, e))?,
    );

    for unwanted in unwanted_info.iter() {
        log::info!("removing INFO tag '{}'", unwanted);
        infos.shift_remove(unwanted);
    }

    Ok(())
}

pub fn construct_record(
    header: &Header,
    recs: Vec<Option<(Rc<Header>, Record)>>,
    vix_samples: &Vec<usize>,
    force_alt_tags: bool,
    unwanted_info: &Vec<String>,
    alts: &Vec<Option<String>>,
    flip: bool,
    criteria: &str,
    repo: &Option<Repository>
) -> std::io::Result<RecordBuf> {
    let _ = header;
    let mut the_record = None;
    for vix in 0..recs.len() {
        if let Some(hnr) = &recs[vix] {
            the_record = Some(hnr);
            break;
        }
    }
    let (the_header, the_record) = the_record.unwrap();

    let mut chrom = String::from(the_record.reference_sequence_name());

    let mut variant_start = the_record.variant_start().unwrap()?;

    let mut ids = Vec::new();
    for id in the_record.ids().iter() {
        ids.push(String::from(id));
    }
    let ids = Ids::from_iter(ids.into_iter());

    let (reference_bases, alternate_bases) = make_ref_and_alt(&the_record, force_alt_tags)?;
    let mut reference_bases = reference_bases;
    let mut alternate_bases = alternate_bases;

    if flip {
        assert_eq!(alternate_bases.len(), 1);
        let orig = BreakEnd::new(&chrom, variant_start.get(), &alternate_bases[0])
            .map_err(|e| Error::new(ErrorKind::Other, e))?;
        log::info!("unflipped: {:?}", orig);
        let flipped = orig.flip();
        let flipped = flipped.format(repo.clone().unwrap())?;
        log::info!("flipped: {:?}", flipped);
        chrom = flipped.0;
        variant_start = Position::try_from(flipped.1).unwrap();
        reference_bases = String::from(flipped.2);
        alternate_bases[0] = flipped.3;
    }

    let alternate_bases = AlternateBases::from(alternate_bases);

    let mut quality_score: f32 = 0.0;
    for vix in 0..recs.len() {
        if let Some(hnr) = &recs[vix] {
            if let Some(this_quality_score) = hnr.1.quality_score() {
                let this_quality_score = this_quality_score?;
                if this_quality_score > quality_score {
                    quality_score = this_quality_score;
                }
            }
        }
    }

    let mut filters = HashSet::new();
    for vix in 0..recs.len() {
        if let Some(hnr) = &recs[vix] {
            for filter in hnr.1.filters().iter(hnr.0.as_ref()) {
                let filter = filter?;
                filters.insert(String::from(filter));
            }
        }
    }
    let mut filters: Vec<String> = filters.into_iter().collect();
    filters.sort();
    let filters = Filters::from_iter(filters.into_iter());

    let mut info: Vec<(String, Option<InfoValue>)> = Vec::new();
    for item in the_record.info().iter(&the_header) {
        let (name, value) = item?;
        let name = String::from(name);
        let value = value.map(make_info_value);
        info.push((name, value));
    }
    if criteria.len() > 0 {
        let criteria: Vec<Option<String>> =
            criteria.split(',').map(|s| Some(String::from(s))).collect();
        info.push((
            String::from("SVELT_CRITERIA"),
            Some(InfoValue::Array(InfoArray::String(criteria))),
        ));
    }
    let mut alt_sequences = Vec::new();
    for vix in 0..alts.len() {
        if let Some(alt) = &alts[vix] {
            alt_sequences.push(Some(String::from(alt)));
        }
    }
    if alt_sequences.len() > 0 {
        info.push((
            String::from("SVELT_ALT_SEQ"),
            Some(InfoValue::Array(InfoArray::String(alt_sequences))),
        ));
    }
    let info: Vec<(String, Option<InfoValue>)> = info
        .into_iter()
        .filter(|item| unwanted_info.iter().all(|unwanted| &item.0 != unwanted))
        .collect();
    let info = Info::from_iter(info.into_iter());

    let keys: Vec<String> = the_record
        .samples()
        .keys()
        .iter()
        .map(|k| String::from(k))
        .collect();
    let mut samples = Vec::new();
    for vix in 0..recs.len() {
        match &recs[vix] {
            Some((header, record)) => {
                for sample in record.samples().iter() {
                    let mut fields: Vec<Option<Value>> = Vec::new();
                    for value in sample.values(header) {
                        if let Some(value) = value {
                            let value = value?;
                            let value = make_sample_value(value);
                            fields.push(Some(value));
                        } else {
                            fields.push(None);
                        }
                    }
                    samples.push(fields);
                }
            }
            None => {
                for _ in 0..vix_samples[vix] {
                    let fields: Vec<Option<Value>> = (0..keys.len()).map(|_| None).collect();
                    samples.push(fields);
                }
            }
        }
    }
    let keys = Keys::from_iter(keys.into_iter());
    let samples = Samples::new(keys, samples);

    let bldr = RecordBuf::builder();
    let res = bldr
        .set_reference_sequence_name(chrom)
        .set_variant_start(variant_start)
        .set_ids(ids)
        .set_reference_bases(reference_bases)
        .set_alternate_bases(alternate_bases)
        .set_quality_score(quality_score)
        .set_filters(filters)
        .set_info(info)
        .set_samples(samples)
        .build();

    log::debug!("record: {:?}", res.clone());

    Ok(res)
}

fn make_ref_and_alt(rec: &Record, force_alt_tags: bool) -> std::io::Result<(String, Vec<String>)> {
    let mut reference_bases = String::from(rec.reference_bases());

    let mut alternate_bases = Vec::new();
    for alt in rec.alternate_bases().iter() {
        let alt = alt?;
        alternate_bases.push(String::from(alt));
    }
    assert_eq!(alternate_bases.len(), 1);

    if force_alt_tags && is_seq(&reference_bases) && alternate_bases.iter().all(|s| is_seq(s)) {
        if reference_bases.len() > 1 {
            reference_bases.truncate(1);
            alternate_bases[0] = String::from("<DEL>");
        } else {
            alternate_bases[0] = String::from("<INS>");
        }
    }

    Ok((reference_bases, alternate_bases))
}

fn make_info_value(
    orig: vcf::variant::record::info::field::value::Value,
) -> vcf::variant::record_buf::info::field::value::Value {
    match orig {
        vcf::variant::record::info::field::Value::Integer(x) => {
            vcf::variant::record_buf::info::field::value::Value::Integer(x)
        }
        vcf::variant::record::info::field::Value::Float(x) => {
            vcf::variant::record_buf::info::field::value::Value::Float(x)
        }
        vcf::variant::record::info::field::Value::Flag => {
            vcf::variant::record_buf::info::field::value::Value::Flag
        }
        vcf::variant::record::info::field::Value::Character(x) => {
            vcf::variant::record_buf::info::field::value::Value::Character(x)
        }
        vcf::variant::record::info::field::Value::String(x) => {
            vcf::variant::record_buf::info::field::value::Value::String(String::from(x))
        }
        vcf::variant::record::info::field::Value::Array(array) => {
            let array = match array {
                vcf::variant::record::info::field::value::Array::Integer(values) => {
                    let values: Vec<Option<i32>> = values.iter().map(|v| v.unwrap()).collect();
                    vcf::variant::record_buf::info::field::value::Array::Integer(values)
                }
                vcf::variant::record::info::field::value::Array::Float(values) => {
                    let values: Vec<Option<f32>> = values.iter().map(|v| v.unwrap()).collect();
                    vcf::variant::record_buf::info::field::value::Array::Float(values)
                }
                vcf::variant::record::info::field::value::Array::Character(values) => {
                    let values: Vec<Option<char>> = values.iter().map(|v| v.unwrap()).collect();
                    vcf::variant::record_buf::info::field::value::Array::Character(values)
                }
                vcf::variant::record::info::field::value::Array::String(values) => {
                    let values: Vec<Option<String>> = values
                        .iter()
                        .map(|v| v.unwrap().map(|s| String::from(s)))
                        .collect();
                    vcf::variant::record_buf::info::field::value::Array::String(values)
                }
            };
            vcf::variant::record_buf::info::field::value::Value::Array(array)
        }
    }
}

fn make_sample_value(orig: vcf::variant::record::samples::series::value::Value) -> Value {
    match orig {
        vcf::variant::record::samples::series::Value::Integer(x) => Value::Integer(x),
        vcf::variant::record::samples::series::Value::Float(x) => Value::Float(x),
        vcf::variant::record::samples::series::Value::Character(x) => Value::Character(x),
        vcf::variant::record::samples::series::Value::String(x) => Value::String(String::from(x)),
        vcf::variant::record::samples::series::Value::Genotype(genotype) => {
            Value::Genotype(Genotype::from_iter(genotype.iter().map(|g| {
                let g = g.unwrap();
                Allele::new(g.0, g.1)
            })))
        }
        vcf::variant::record::samples::series::Value::Array(array) => match array {
            vcf::variant::record::samples::series::value::Array::Integer(values) => {
                let values: Vec<Option<i32>> = values.iter().map(|v| v.unwrap()).collect();
                Value::Array(Array::Integer(values))
            }
            vcf::variant::record::samples::series::value::Array::Float(values) => {
                let values: Vec<Option<f32>> = values.iter().map(|v| v.unwrap()).collect();
                Value::Array(Array::Float(values))
            }
            vcf::variant::record::samples::series::value::Array::Character(values) => {
                let values: Vec<Option<char>> = values.iter().map(|v| v.unwrap()).collect();
                Value::Array(Array::Character(values))
            }
            vcf::variant::record::samples::series::value::Array::String(values) => {
                let values: Vec<Option<String>> = values
                    .iter()
                    .map(|v| v.unwrap().map(|s| String::from(s)))
                    .collect();
                Value::Array(Array::String(values))
            }
        },
    }
}
