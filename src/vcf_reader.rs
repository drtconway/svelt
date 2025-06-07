use std::{io::BufRead, rc::Rc};

use noodles::vcf::{self, variant::record::info::field::Value, Header, Record};

use crate::{
    chroms::ChromSet,
    errors::{SveltError, as_io_error},
};

pub struct VcfReader {
    pub path: String,
    pub reader: vcf::io::Reader<Box<dyn BufRead>>,
    pub header: Header,
    pub chroms: Rc<ChromSet>,
}

impl VcfReader {
    pub fn new(path: &str, chroms: Rc<ChromSet>) -> std::io::Result<VcfReader> {
        let path = String::from(path);
        let reader = autocompress::autodetect_open(&path)?;
        let mut reader: vcf::io::Reader<Box<dyn BufRead>> =
            vcf::io::reader::Builder::default().build_from_reader(reader)?;
        let header = reader.read_header()?;
        check_chroms(&header, chroms.as_ref()).map_err(as_io_error)?;

        Ok(VcfReader {
            path,
            reader,
            header,
            chroms,
        })
    }

    pub fn rewind(&mut self) -> std::io::Result<()> {
        let reader = autocompress::autodetect_open(&self.path)?;
        let mut reader: vcf::io::Reader<Box<dyn BufRead>> =
            vcf::io::reader::Builder::default().build_from_reader(reader)?;
        let _header = reader.read_header()?;
        self.reader = reader;
        Ok(())
    }

    pub fn info_as_str(rec: &Record, header: &Header, name: &str) -> std::io::Result<Option<String>> {
        if let Some(field) = rec.info().get(header, name) {
            let field = field?;
            match field {
                Some(value) => {
                    if let Value::String(value) = value {
                        Ok(Some(String::from(value)))
                    } else {
                        let chrom = String::from(rec.reference_sequence_name());
                        let pos = rec.variant_start().unwrap()?.get();
                        Err(as_io_error(SveltError::BadInfoType(
                            chrom,
                            pos,
                            String::from(name),
                            String::from("String"),
                        )))
                    }
                }
                None => Ok(None),
            }
        } else {
            Ok(None)
        }
    }

    pub fn info_as_int(rec: &Record, header: &Header, name: &str) -> std::io::Result<Option<i32>> {
        if let Some(field) = rec.info().get(header, name) {
            let field = field?;
            match field {
                Some(value) => {
                    if let Value::Integer(value) = value {
                        Ok(Some(value))
                    } else {
                        let chrom = String::from(rec.reference_sequence_name());
                        let pos = rec.variant_start().unwrap()?.get();
                        Err(as_io_error(SveltError::BadInfoType(
                            chrom,
                            pos,
                            String::from(name),
                            String::from("Integer"),
                        )))
                    }
                }
                None => Ok(None),
            }
        } else {
            Ok(None)
        }
    }
}

pub fn check_chroms(header: &Header, chroms: &ChromSet) -> std::result::Result<(), SveltError> {
    let n = chroms.len();
    let n0 = header.contigs().len();
    if n != n0 {
        return Err(SveltError::Contigs(n, n0));
    }
    for (ix0, contig) in header.contigs().iter().enumerate() {
        let name = contig.0;
        if let Some(ix) = chroms.index(name) {
            if ix != ix0 {
                return Err(SveltError::ContigOrder(name.clone(), ix, ix0));
            }
        } else {
            return Err(SveltError::ContigMissing(name.clone(), ix0));
        }
    }
    Ok(())
}
