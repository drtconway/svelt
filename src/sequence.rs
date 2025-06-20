use std::io::{BufRead, BufReader, Error, ErrorKind};

use autocompress::{Processor, autodetect_open, io::ProcessorReader};
use noodles::{
    fasta,
    vcf::{
        self, Header,
        variant::{
            Record as _, RecordBuf,
            record::{AlternateBases, Ids},
            record_buf::info::field::{Value, value::Array},
        },
    },
};

use crate::{
    errors::{SveltError, as_io_error},
    tables::is_seq,
};

pub trait SequenceIterator: Iterator<Item = std::io::Result<(String, String)>> {}

type InnerReader = BufReader<
    ProcessorReader<Box<dyn Processor + Send + Unpin + 'static>, BufReader<std::fs::File>>,
>;

pub struct FastaSequenceIterator {
    reader: fasta::io::Reader<InnerReader>,
}

impl FastaSequenceIterator {
    pub fn new(filename: &str) -> std::io::Result<FastaSequenceIterator> {
        let reader = autodetect_open(filename)?;
        let reader = BufReader::new(reader);
        let reader = fasta::io::reader::Builder::default().build_from_reader(reader)?;
        Ok(FastaSequenceIterator { reader })
    }

    fn read_one(&mut self) -> std::io::Result<Option<(String, String)>> {
        let mut definition = String::new();
        let mut sequence = Vec::new();

        let r1 = self.reader.read_definition(&mut definition)?;
        let definition = String::from(&definition[1..]);

        if r1 == 0 {
            return Ok(None);
        }

        self.reader.read_sequence(&mut sequence)?;
        let sequence = String::from_utf8(sequence).map_err(|e| Error::new(ErrorKind::Other, e))?;

        Ok(Some((definition, sequence)))
    }
}

impl Iterator for FastaSequenceIterator {
    type Item = std::io::Result<(String, String)>;

    fn next(&mut self) -> Option<Self::Item> {
        match self.read_one() {
            Ok(Some(item)) => Some(Ok(item)),
            Ok(None) => None,
            Err(error) => Some(Err(error)),
        }
    }
}

impl SequenceIterator for FastaSequenceIterator {}

pub struct VcfSequenceIterator {
    reader: vcf::io::Reader<Box<dyn BufRead>>,
    header: Header,
}

impl VcfSequenceIterator {
    pub fn new(filename: &str) -> std::io::Result<VcfSequenceIterator> {
        let reader = autodetect_open(&filename)?;
        let mut reader: vcf::io::Reader<Box<dyn BufRead>> =
            vcf::io::reader::Builder::default().build_from_reader(reader)?;
        let header = reader.read_header()?;
        Ok(VcfSequenceIterator { reader, header })
    }

    fn read_one(&mut self) -> std::io::Result<Option<(String, String)>> {
        let mut record = RecordBuf::default();
        loop {
            let r1 = self.reader.read_record_buf(&self.header, &mut record)?;
            if r1 == 0 {
                return Ok(None);
            }

            // First, if it has SVTYPE and it isn't INS, skip.
            let kind = VcfSequenceIterator::get_info_str(&record, "SVTYPE")?;
            if let Some(kind) = kind {
                if kind != "INS" {
                    continue;
                }
            }

            let chrom = record.reference_sequence_name().to_string();
            let pos = record.variant_start().unwrap().get();

            let name = record.ids().iter().next();
            let name = if let Some(name) = name {
                format!("{}:{}:{}", chrom, pos, name)
            } else {
                format!("{}:{}:_", chrom, pos)
            };

            // Now look for a "long" alt.
            for alt in record.alternate_bases().iter() {
                let alt = alt?;
                if alt.len() > 1 && is_seq(alt) {
                    return Ok(Some((name, String::from(alt))));
                }
            }

            // Now look for a "SVELT_ALT_SEQ" tag.
            let alt = VcfSequenceIterator::get_info_str(&record, "SVELT_ALT_SEQ")?;
            if let Some(alt) = alt {
                return Ok(Some((name, alt)));
            }
        }
    }

    fn get_info_str(record: &RecordBuf, tag: &str) -> std::io::Result<Option<String>> {
        if let Some(value) = record.info().get(tag) {
            if let Some(value) = value {
                match value {
                    Value::String(str) => Ok(Some(str.to_string())),
                    Value::Array(Array::String(strs)) => {
                        let str = strs.iter().filter(|s| s.is_some()).next();
                        if let Some(Some(str)) = str {
                            Ok(Some(str.to_string()))
                        } else {
                            let chrom = record.reference_sequence_name().to_string();
                            let pos = record.variant_start().unwrap().get();
                            Err(as_io_error(SveltError::BadInfoType(
                                chrom,
                                pos,
                                String::from(tag),
                                String::from("String"),
                            )))
                        }
                    }
                    _ => {
                        let chrom = record.reference_sequence_name().to_string();
                        let pos = record.variant_start().unwrap().get();
                        Err(as_io_error(SveltError::BadInfoType(
                            chrom,
                            pos,
                            String::from(tag),
                            String::from("String"),
                        )))
                    }
                }
            } else {
                let chrom = record.reference_sequence_name().to_string();
                let pos = record.variant_start().unwrap().get();
                Err(as_io_error(SveltError::MissingInfo(
                    chrom,
                    pos,
                    String::from(tag),
                )))
            }
        } else {
            Ok(None)
        }
    }
}

impl Iterator for VcfSequenceIterator {
    type Item = std::io::Result<(String, String)>;

    fn next(&mut self) -> Option<Self::Item> {
        match self.read_one() {
            Ok(Some(item)) => Some(Ok(item)),
            Ok(None) => None,
            Err(error) => Some(Err(error)),
        }
    }
}

impl SequenceIterator for VcfSequenceIterator {}
