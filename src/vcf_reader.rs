use std::{io::BufRead, rc::Rc};

use noodles::vcf::{self, Header};

use crate::{chroms::ChromSet, errors::{as_io_error, SveltError}};

pub struct VcfReader {
    pub path: String,
    pub reader: vcf::io::Reader<Box<dyn BufRead>>,
    pub header: Header,
    pub chroms: Rc<ChromSet>
}

impl VcfReader {
    pub fn new(path: &str, chroms: Rc<ChromSet>) -> std::io::Result<VcfReader> {
        let path = String::from(path);
        let reader = autocompress::autodetect_open(&path)?;
        let mut reader: vcf::io::Reader<Box<dyn BufRead>> =
        vcf::io::reader::Builder::default().build_from_reader(reader)?;
        let header = reader.read_header()?;
        check_chroms(&header, chroms.as_ref()).map_err(as_io_error)?;

        Ok(VcfReader { path, reader, header, chroms })
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