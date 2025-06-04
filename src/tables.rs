use std::sync::Arc;

use datafusion::arrow::{
    array::{PrimitiveBuilder, RecordBatch, StringDictionaryBuilder},
    datatypes::{DataType, Field, Int32Type, Schema, UInt16Type, UInt32Type},
};
use noodles::vcf::{
    Header,
    variant::record::{AlternateBases, Record, info::field::Value},
};

use crate::{
    breakends::parse_breakend,
    chroms::ChromSet,
    errors::{SveltError, as_io_error},
    vcf_reader::VcfReader,
};

pub fn vcf_core_schema() -> Arc<Schema> {
    Arc::new(Schema::new(vec![
        Field::new("chrom_id", DataType::UInt16, false),
        Field::new_dictionary("chrom", DataType::UInt16, DataType::Utf8, false),
        Field::new("start", DataType::UInt32, false),
        Field::new("end", DataType::UInt32, false),
        Field::new_dictionary("kind", DataType::UInt8, DataType::Utf8, false),
        Field::new("length", DataType::Int32, false),
        Field::new("chrom2_id", DataType::UInt16, true),
        Field::new_dictionary("chrom2", DataType::UInt16, DataType::Utf8, false),
        Field::new("end2", DataType::UInt32, true),
        Field::new("seq_hash", DataType::UInt64, true),
    ]))
}

pub fn load_vcf_core(reader: &mut VcfReader) -> std::io::Result<RecordBatch> {
    let header: &Header = &reader.header;
    let chroms: &ChromSet = reader.chroms.as_ref();

    let mut chrom_id_builder = PrimitiveBuilder::<UInt16Type>::new();
    let mut chrom_builder = StringDictionaryBuilder::<UInt16Type>::new();
    let mut start_builder = PrimitiveBuilder::<UInt32Type>::new();
    let mut end_builder = PrimitiveBuilder::<UInt32Type>::new();
    let mut kind_builder = StringDictionaryBuilder::<UInt16Type>::new();
    let mut length_builder = PrimitiveBuilder::<Int32Type>::new();
    let mut chrom2_id_builder = PrimitiveBuilder::<UInt16Type>::new();
    let mut chrom2_builder = StringDictionaryBuilder::<UInt16Type>::new();
    let mut end2_builder = PrimitiveBuilder::<UInt32Type>::new();
    let mut seq_hash_builder = PrimitiveBuilder::<UInt32Type>::new();

    for rec in reader.reader.records() {
        let rec = rec?;

        let chrom_id = chroms.index(rec.reference_sequence_name()).unwrap();
        let chrom = String::from(rec.reference_sequence_name());
        let start = rec.variant_start().unwrap()?.get();
        let end = rec.variant_end(header)?.get();
        let kind = if let Some(Value::String(value)) = rec.info().get(header, "SVTYPE").unwrap()? {
            Ok(String::from(value))
        } else {
            Err(SveltError::MissingType(
                String::from(rec.reference_sequence_name()),
                start,
            ))
        }
        .map_err(as_io_error)?;
        let length =
            if let Some(Value::Integer(value)) = rec.info().get(header, "SVLEN").unwrap()? {
                Some(value)
            } else {
                None
            };
        let chrom2 = if let Some(Value::String(value)) = rec.info().get(header, "CHR2").unwrap()? {
            Ok(Some(String::from(value)))
        } else {
            if kind == "BND" {
                Err(SveltError::MissingChr2(
                    String::from(rec.reference_sequence_name()),
                    start,
                ))
            } else {
                Ok(None)
            }
        }
        .map_err(as_io_error)?;
        let chrom2_id = if let Some(name2) = chrom2.as_ref() {
            chroms.index(name2).map(|x| x as u16)
        } else {
            None
        };
        let end2: Option<u32> = if kind == "BND" {
            if let Some(alt) = rec.alternate_bases().iter().next() {
                let alt = alt?;
                let (chr2, pos2, here, there) = parse_breakend(alt).map_err(as_io_error)?;
                if &chr2 != chrom2.as_ref().unwrap() {
                    Err(SveltError::BadChr2(
                        String::from(rec.reference_sequence_name()),
                        start,
                        chrom2.as_ref().unwrap().clone(),
                        chr2.clone(),
                    ))
                } else {
                    Ok(Some(pos2 as u32))
                }
            } else {
                Err(SveltError::MissingAlt(
                    String::from(rec.reference_sequence_name()),
                    start,
                ))
            }
        } else {
            Ok(None)
        }
        .map_err(as_io_error)?;

        chrom_id_builder.append_value(chrom_id as u16);
        chrom_builder.append_value(chrom);
        start_builder.append_value(start as u32);
        end_builder.append_value(end as u32);
        kind_builder.append_value(kind);
        length_builder.append_option(length);
        chrom2_id_builder.append_option(chrom2_id);
        chrom2_builder.append_option(chrom2);
        end2_builder.append_option(end2);
    }

    todo!()
}
