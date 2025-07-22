use std::collections::HashSet;
use std::io::{BufWriter, Error, ErrorKind};
use std::rc::Rc;
use std::str::FromStr;

use autocompress::io::ProcessorWriter;
use autocompress::{CompressionLevel, Processor, autodetect_create};
use noodles::core::Position;
use noodles::fasta::Repository;
use noodles::vcf;
use noodles::vcf::header::record::value::Map;
use noodles::vcf::header::record::value::map::info::{Number, Type};
use noodles::vcf::header::record::value::map::{Builder, Filter};
use noodles::vcf::variant::io::Write;
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

use crate::breakends::{BreakEnd, parse_breakend};
use crate::options::MergeOptions;
use crate::tables::is_seq;

pub type InnerWriter =
    BufWriter<ProcessorWriter<Box<dyn Processor + Send + Unpin + 'static>, std::fs::File>>;
pub struct MergeBuilder {
    writer: vcf::io::Writer<InnerWriter>,
    options: Rc<MergeOptions>,
    header: Header,
    reference: Option<Rc<Repository>>,
    current_chrom: String,
}

impl MergeBuilder {
    pub fn new(
        out: &str,
        options: Rc<MergeOptions>,
        header: Header,
        reference: Option<Rc<Repository>>,
    ) -> std::io::Result<MergeBuilder> {
        let writer = autodetect_create(out, CompressionLevel::Default)?;
        let writer = BufWriter::new(writer);
        let mut writer = vcf::io::Writer::new(writer);
        writer.write_header(&header)?;

        Ok(MergeBuilder {
            writer,
            options,
            header,
            reference,
            current_chrom: String::new(),
        })
    }

    pub fn construct(
        &mut self,
        recs: Vec<Option<(Rc<Header>, Record)>>,
        vix_samples: &Vec<usize>,
        vids: &Vec<String>,
        alts: &Vec<Option<String>>,
        paired_bnd: bool,
        criteria: &str,
        feature: &str,
    ) -> std::io::Result<()> {
        let rec = construct_record(
            &self.header,
            recs,
            &vix_samples,
            vids,
            &alts,
            paired_bnd,
            &criteria,
            feature,
            self.options.as_ref(),
            &self.reference,
        )?;
        self.writer.write_variant_record(&self.header, &rec)?;
        if rec.reference_sequence_name() != &self.current_chrom {
            self.current_chrom = String::from(rec.reference_sequence_name());
            log::info!("writing variants for {}", self.current_chrom);
        }
        Ok(())
    }
}

pub fn add_svelt_header_fields(
    header: &mut Header,
    unwanted_info: &Vec<String>,
) -> std::io::Result<()> {
    let filters = header.filters_mut();

    if filters.get("UNPAIRED_BND").is_none() {
        filters.insert(
            String::from("UNPAIRED_BND"),
            Map::<Filter>::new("Breakend variant does not have a symmetric pair"),
        );
    }

    let infos = header.infos_mut();

    if infos.get("CHR2").is_none() {
        infos.insert(
            String::from("CHR2"),
            Builder::default()
                .set_number(Number::Unknown)
                .set_type(Type::String)
                .set_description("Chromosome for end coordinate in case of a translocation")
                .build()
                .map_err(|e| Error::new(ErrorKind::Other, e))?,
        );
    }

    if infos.get("END2").is_none() {
        infos.insert(
            String::from("END2"),
            Builder::default()
                .set_number(Number::Unknown)
                .set_type(Type::Integer)
                .set_description("End position of the structural variant on CHR2")
                .build()
                .map_err(|e| Error::new(ErrorKind::Other, e))?,
        );
    }

    infos.insert(
        String::from("ORIGINAL_IDS"),
        Builder::default()
            .set_number(Number::Unknown)
            .set_type(Type::String)
            .set_description("The variant IDs from the original VCFs")
            .build()
            .map_err(|e| Error::new(ErrorKind::Other, e))?,
    );

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

    infos.insert(
        String::from("SVELT_ALT_CLASS"),
        Builder::default()
            .set_number(Number::Unknown)
            .set_type(Type::String)
            .set_description("Classification of the inserted sequence.")
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
    vids: &Vec<String>,
    alts: &Vec<Option<String>>,
    paired_bnd: bool,
    criteria: &str,
    feature: &str,
    options: &MergeOptions,
    reference: &Option<Rc<Repository>>,
) -> std::io::Result<RecordBuf> {
    let _ = header;
    let mut the_variant_id = String::new();
    let mut the_record = None;
    for vix in 0..recs.len() {
        if let Some(hnr) = &recs[vix] {
            the_variant_id = vids[vix].clone();
            the_record = Some(hnr);
            break;
        }
    }
    let (the_header, the_record) = the_record.unwrap();

    let chrom = String::from(the_record.reference_sequence_name());

    let variant_start = if let Some(start) = the_record.variant_start() {
        start?.get()
    } else {
        0
    };

    let ids = vec![the_variant_id];
    let ids = Ids::from_iter(ids.into_iter());

    let (reference_bases, alternate_bases) = make_ref_and_alt(&the_record, options.force_alt_tags)?;
    let mut reference_bases = reference_bases;
    let mut alternate_bases = alternate_bases;

    let variant_start = if variant_start > 0 {
        Position::try_from(variant_start).unwrap()
    } else {
        Position::try_from(1).unwrap()
    };

    let mut chrom2 = None;
    let mut end2 = None;
    if let Ok(bnd) = BreakEnd::new(&chrom, variant_start.get(), &alternate_bases[0]) {
        chrom2 = Some(bnd.chrom2.clone());
        end2 = Some(bnd.end2);
    }

    if options.fill_in_refs {
        if reference_bases == "N" || reference_bases == "n" {
            let seq = reference.clone().unwrap().get(chrom.as_ref()).unwrap()?;
            let pos = Position::try_from(variant_start).unwrap();
            let b: &u8 = if pos.get() < seq.len() {
                seq.get(pos).unwrap()
            } else {
                &b'N'
            };
            let b = *b as char;
            reference_bases = String::from(b);
        }
        for i in 0..alternate_bases.len() {
            let alt = &alternate_bases[i] as &str;
            if is_seq(alt) {
                continue;
            }
            if let Ok(_) = parse_breakend(alt) {
                if alt.starts_with('N') || alt.starts_with('n') {
                    let seq = reference.clone().unwrap().get(chrom.as_ref()).unwrap()?;
                    let pos = Position::try_from(variant_start).unwrap();
                    let b: &u8 = if pos.get() < seq.len() {
                        seq.get(pos).unwrap()
                    } else {
                        &b'N'
                    };
                    let b = *b as char;
                    alternate_bases[i] = format!("{}{}", b, &alt[1..]);
                } else if alt.ends_with('N') || alt.starts_with('n') {
                    let seq = reference.clone().unwrap().get(chrom.as_ref()).unwrap()?;
                    let pos = Position::try_from(variant_start).unwrap();
                    let b: &u8 = if pos.get() < seq.len() {
                        seq.get(pos).unwrap()
                    } else {
                        &b'N'
                    };
                    let b = *b as char;
                    alternate_bases[i] = format!("{}{}", &alt[0..(alt.len() - 1)], b);
                }
            }
        }
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
    if chrom2.is_some() && !paired_bnd {
        filters.insert(String::from("UNPAIRED_BND"));
    }
    let mut filters: Vec<String> = filters.into_iter().collect();
    filters.sort();
    let filters = Filters::from_iter(filters.into_iter());

    // Gather up the original IDs
    let mut original_ids = Vec::new();
    for vix in 0..recs.len() {
        if let Some(hnr) = &recs[vix] {
            for id in hnr.1.ids().iter() {
                original_ids.push(Some(String::from(id)));
            }
        }
    }

    let mut info: Vec<(String, Option<InfoValue>)> = Vec::new();
    for item in the_record.info().iter(&the_header) {
        let (name, value) = item.unwrap();
        let name = String::from(name);
        let value = value.map(make_info_value);
        info.push((name, value));
    }
    if original_ids.len() > 0 {
        info.push((
            String::from("ORIGINAL_IDS"),
            Some(InfoValue::Array(InfoArray::String(original_ids))),
        ));
    }
    if let (Some(chrom2), Some(end2)) = (&chrom2, &end2) {
        info = info
            .into_iter()
            .filter(|item| item.0 != "CHR2" && item.0 != "END2")
            .collect();
        info.push((
            String::from("CHR2"),
            Some(InfoValue::String(chrom2.clone())),
        ));
        info.push((String::from("END2"), Some(InfoValue::Integer(*end2 as i32))));
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
    if feature.len() > 0 {
        info.push((
            String::from("SVELT_ALT_CLASS"),
            Some(InfoValue::String(String::from(feature))),
        ));
    }
    let info: Vec<(String, Option<InfoValue>)> = info
        .into_iter()
        .filter(|item| {
            options
                .unwanted_info
                .iter()
                .all(|unwanted| &item.0 != unwanted)
        })
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
                    let fields: Vec<Option<Value>> = keys
                        .iter()
                        .map(|k| make_empty_fmt_value(options, k))
                        .collect();
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

fn make_empty_fmt_value(options: &MergeOptions, key: &str) -> Option<Value> {
    if key == "GT" {
        let gt = if options.use_ref_alleles {
            Genotype::from_str("0/0").unwrap()
        } else {
            Genotype::from_str("./.").unwrap()
        };
        Some(Value::Genotype(gt))
    } else {
        None
    }
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
