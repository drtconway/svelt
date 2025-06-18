use std::{
    io::{Error, ErrorKind},
    sync::Arc,
};

use blake2::{Blake2b512, Digest};
use datafusion::arrow::{
    array::{GenericStringBuilder, PrimitiveBuilder, RecordBatch, StringDictionaryBuilder},
    datatypes::{DataType, Field, Int32Type, Int64Type, Schema, UInt8Type, UInt16Type, UInt32Type},
};
use noodles::vcf::{
    Header,
    variant::record::{AlternateBases, Record},
};

use crate::{
    chroms::ChromSet,
    errors::{SveltError, as_io_error},
    inputs::{get_breakend, get_svtype},
    vcf_reader::VcfReader,
};

/// Return a schema for the core SV VCF fields we use.
pub fn vcf_core_schema() -> Arc<Schema> {
    Arc::new(Schema::new(vec![
        Field::new("row_num", DataType::UInt32, false),
        Field::new("chrom_id", DataType::UInt16, false),
        Field::new_dictionary("chrom", DataType::UInt16, DataType::Utf8, false),
        Field::new("start", DataType::Int32, false),
        Field::new("end", DataType::Int32, false),
        Field::new_dictionary("kind", DataType::UInt8, DataType::Utf8, false),
        Field::new("length", DataType::Int32, true),
        Field::new("chrom2_id", DataType::UInt16, true),
        Field::new_dictionary("chrom2", DataType::UInt16, DataType::Utf8, true),
        Field::new("end2", DataType::Int32, true),
        Field::new("alt_seq", DataType::Utf8, true),
        Field::new("seq_hash", DataType::Int64, true),
    ]))
}

/// Read all the records in the VCF are return them as a `RecordBatch`
pub fn load_vcf_core(reader: &mut VcfReader) -> std::io::Result<RecordBatch> {
    let header: &Header = &reader.header;
    let chroms: &ChromSet = reader.chroms.as_ref();

    let mut row_num_builder = PrimitiveBuilder::<UInt32Type>::new();
    let mut chrom_id_builder = PrimitiveBuilder::<UInt16Type>::new();
    let mut chrom_builder = StringDictionaryBuilder::<UInt16Type>::new();
    let mut start_builder = PrimitiveBuilder::<Int32Type>::new();
    let mut end_builder = PrimitiveBuilder::<Int32Type>::new();
    let mut kind_builder = StringDictionaryBuilder::<UInt8Type>::new();
    let mut length_builder = PrimitiveBuilder::<Int32Type>::new();
    let mut chrom2_id_builder = PrimitiveBuilder::<UInt16Type>::new();
    let mut chrom2_builder = StringDictionaryBuilder::<UInt16Type>::new();
    let mut end2_builder = PrimitiveBuilder::<Int32Type>::new();
    let mut alt_seq_builder = GenericStringBuilder::<i32>::new();
    let mut seq_hash_builder = PrimitiveBuilder::<Int64Type>::new();

    for (rn, rec) in reader.reader.records().enumerate() {
        log::debug!("processing record {}", rn);
        let rec = rec?;

        let chrom_id = chroms.index(rec.reference_sequence_name()).unwrap();
        let chrom = String::from(rec.reference_sequence_name());

        let start = if let Some(start) = rec.variant_start() {
            start?.get()
        } else {
            0
        };

        let mut end = if start > 0 {
            rec.variant_end(header)?.get()
        } else {
            if let Some(value) = VcfReader::info_as_int(&rec, header, "END")? {
                value as usize
            } else if let Some(value) = VcfReader::info_as_int(&rec, header, "SVLEN")? {
                start + (value as usize)
            } else {
                0
            }
        };

        let kind = get_svtype(&rec, header)?;

        let length = if let Some(value) = VcfReader::info_as_int(&rec, header, "SVLEN")? {
            Some(value)
        } else {
            None
        };

        if let Some(l) = &length {
            if kind == "DEL" && start + l.abs() as usize != end {
                let d = (start as i32) + l.abs() - (end as i32);
                log::warn!("at {}:{}, end ({}) and length ({}) are inconsistent (off by {}).", chrom, start, end, l, d);
                end = start + l.abs() as usize;
            }
        }

        let bnd = get_breakend(&rec)?;

        let chrom2 = if let Some(value) = VcfReader::info_as_str(&rec, header, "CHR2")? {
            if let Some(bnd) = &bnd {
                if value != bnd.0 {
                    Err(SveltError::BadChr2(
                        chrom.clone(),
                        start,
                        bnd.0.clone(),
                        value.clone(),
                    ))
                } else {
                    Ok(Some(String::from(value)))
                }
            } else {
                Ok(Some(String::from(value)))
            }
        } else {
            if let Some(bnd) = &bnd {
                Ok(Some(bnd.0.clone()))
            } else if kind == "BND" {
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
        let end2: Option<i32> = if kind == "BND" {
            if let Some(bnd) = &bnd {
                let (chr2, pos2, _here, _there) = bnd;
                if chr2 != chrom2.as_ref().unwrap() {
                    Err(SveltError::BadChr2(
                        String::from(rec.reference_sequence_name()),
                        start,
                        chrom2.as_ref().unwrap().clone(),
                        chr2.clone(),
                    ))
                } else {
                    Ok(Some(*pos2 as i32))
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
        let seq: Option<String> = if kind == "INS" {
            if let Some(alt) = rec.alternate_bases().iter().next() {
                let alt = alt?;
                if is_seq(alt) {
                    Ok(Some(String::from(&alt[1..])))
                } else {
                    Ok(None)
                }
            } else {
                Ok(None)
            }
        } else {
            Ok(None)
        }
        .map_err(as_io_error)?;
        let seq_hash = if kind == "INS" {
            if let Some(alt) = rec.alternate_bases().iter().next() {
                let alt = alt?;
                if is_seq(alt) {
                    Ok(Some(digest(&alt[1..]) as i64))
                } else {
                    // hash the tag
                    Ok(Some(digest(alt) as i64))
                }
            } else {
                Ok(None)
            }
        } else {
            Ok(None)
        }
        .map_err(as_io_error)?;

        row_num_builder.append_value(rn as u32);
        chrom_id_builder.append_value(chrom_id as u16);
        chrom_builder.append_value(chrom);
        start_builder.append_value(start as i32);
        end_builder.append_value(end as i32);
        kind_builder.append_value(kind);
        length_builder.append_option(length);
        chrom2_id_builder.append_option(chrom2_id);
        chrom2_builder.append_option(chrom2);
        end2_builder.append_option(end2);
        alt_seq_builder.append_option(seq);
        seq_hash_builder.append_option(seq_hash);
    }

    let row_num_array = row_num_builder.finish();
    let chrom_id_array = chrom_id_builder.finish();
    let chrom_array = chrom_builder.finish();
    let start_array = start_builder.finish();
    let end_array = end_builder.finish();
    let kind_array = kind_builder.finish();
    let length_array = length_builder.finish();
    let chrom2_id_array = chrom2_id_builder.finish();
    let chrom2_array = chrom2_builder.finish();
    let end2_array = end2_builder.finish();
    let alt_seq_array = alt_seq_builder.finish();
    let seq_hash_array = seq_hash_builder.finish();

    let res = RecordBatch::try_new(
        vcf_core_schema(),
        vec![
            Arc::new(row_num_array),
            Arc::new(chrom_id_array),
            Arc::new(chrom_array),
            Arc::new(start_array),
            Arc::new(end_array),
            Arc::new(kind_array),
            Arc::new(length_array),
            Arc::new(chrom2_id_array),
            Arc::new(chrom2_array),
            Arc::new(end2_array),
            Arc::new(alt_seq_array),
            Arc::new(seq_hash_array),
        ],
    )
    .map_err(|e| Error::new(ErrorKind::Other, e))?;
    Ok(res)
}

pub fn is_seq(seq: &str) -> bool {
    seq.chars().all(|c| match c {
        'A' | 'C' | 'G' | 'T' | 'N' => true,
        'a' | 'c' | 'g' | 't' | 'n' => true,
        '*' => true,
        _ => false,
    })
}

pub fn digest(seq: &str) -> u64 {
    let mut hasher = Blake2b512::new();
    hasher.update(seq.as_bytes());
    let dig = hasher.finalize();
    let mut h = 0;
    for i in 0..7 {
        h = (h << 8) | (dig[i] as u64);
    }
    h
}
