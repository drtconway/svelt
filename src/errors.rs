use std::{
    error::Error,
    fmt::Display,
    io::{Error as IoError, ErrorKind},
};

#[derive(Debug)]
pub enum SveltError {
    BadBreakEnd(String),
    BadChr2(String, String),
    BadChrom(String),
    BadFormatField(String, Box<dyn Error + Send + Sync + 'static>),
    BadInfoField(String, Box<dyn Error + Send + Sync + 'static>),
    BadInfoType(String, String),
    BadKind(String),
    BadSample(String, Box<dyn Error + Send + Sync + 'static>),
    BadVariant(String, usize, Box<dyn Error + Send + Sync + 'static>),
    Contigs(usize, usize),
    ContigMissing(String, usize),
    ContigOrder(String, usize, usize),
    FileError(String, Box<dyn Error + Send + Sync + 'static>),
    MissingAlt,
    MissingChr2,
    MissingInfo(String),
    MissingK(String),
    MissingType,
    NeardexDuplicate(u32),
    OptionReferenceRequired(String),
    TooManyVcfs(usize),
    UnexpectedNull(String),
}

impl Display for SveltError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            SveltError::BadBreakEnd(alt) => {
                write!(f, "Badly formed breakend '{}'", alt)
            }
            SveltError::BadChr2(chrom2, chr2) => {
                write!(
                    f,
                    "Inconsistent BND INFO - BND had {}, CHR2 had {}",
                    chr2, chrom2
                )
            }
            SveltError::BadChrom(chrom) => {
                write!(f, "Use of undeclared contig '{}'", chrom)
            }
            SveltError::BadFormatField(name, _error) => {
                write!(f, "Problem with parsing FORMAT field '{}'", name)
            }
            SveltError::BadInfoField(name, _error) => {
                write!(f, "Problem with parsing INFO field '{}'", name)
            }
            SveltError::BadInfoType(tag, exp) => {
                write!(
                    f,
                    "Unexpected type of INFO - tag was {}, expected type was {}",
                    tag, exp
                )
            }
            SveltError::BadKind(kind) => {
                write!(f, "Unexpected SVTYPE: '{}'", kind)
            }
            SveltError::BadSample(name, _error) => {
                write!(f, "Problem with parsing sample field '{}'", name)
            }
            SveltError::BadVariant(chrom, position, _error) => {
                write!(f, "Problem with variant at {}:{}", chrom, position)
            }
            SveltError::Contigs(exp, got) => write!(
                f,
                "Unexpected number of contigs - expected {}, got {}",
                exp, got
            ),
            SveltError::ContigMissing(chrom, got) => {
                write!(f, "Unexpected contig '{}' (number {}", chrom, got)
            }
            SveltError::ContigOrder(chrom, exp, got) => {
                write!(
                    f,
                    "Contig '{}' out of order - expected {}, got {}",
                    chrom, exp, got
                )
            }
            SveltError::FileError(filename, _error) => {
                write!(f, "Problem processing file '{}'", filename)
            }
            SveltError::MissingAlt => {
                write!(f, "Missing ALT")
            }
            SveltError::MissingChr2 => {
                write!(f, "Missing CHR2 INFO field")
            }
            SveltError::MissingInfo(field) => {
                write!(f, "Expected field '{}' not found", field)
            }
            SveltError::MissingK(name) => {
                write!(f, "index '{}' has missing metadata", name)
            }
            SveltError::MissingType => {
                write!(f, "Missing SVTYPE")
            }
            SveltError::NeardexDuplicate(key) => {
                write!(f, "Cannot construct Neardex with duplicate key {}", key)
            }
            SveltError::OptionReferenceRequired(opt) => {
                write!(
                    f,
                    "Option '{}' requires a reference sequence to be given",
                    opt
                )
            }
            SveltError::TooManyVcfs(n) => {
                write!(f, "Too many VCFs (max 64) ({} given)", n)
            }
            SveltError::UnexpectedNull(src) => {
                write!(f, "unexpected null value at {}", src)
            }
        }
    }
}

impl std::error::Error for SveltError {
    fn source(&self) -> Option<&(dyn Error + 'static)> {
        match self {
            SveltError::BadFormatField(_name, error) => Some(error.as_ref()),
            SveltError::BadInfoField(_name, error) => Some(error.as_ref()),
            SveltError::BadSample(_name, error) => Some(error.as_ref()),
            SveltError::BadVariant(_chrom, _position, error) => Some(error.as_ref()),
            SveltError::FileError(_path, error) => Some(error.as_ref()),
            _ => None,
        }
    }
}

pub fn as_io_error(error: SveltError) -> std::io::Error {
    IoError::new(ErrorKind::Other, error)
}

pub fn wrap_file_error<E: std::error::Error + Send + Sync + 'static>(
    e: E,
    filename: &str,
) -> std::io::Error {
    as_io_error(SveltError::FileError(String::from(filename), Box::new(e)))
}

pub trait Context {
    fn with<R, F: FnOnce() -> std::io::Result<R>>(&self, inner: F) -> std::io::Result<R>;
}

pub struct FileContext<'a> {
    filename: &'a str,
}

impl<'a> FileContext<'a> {
    pub fn new(filename: &'a str) -> Self {
        FileContext { filename }
    }
}

impl<'a> Context for FileContext<'a> {
    fn with<R, F: FnOnce() -> std::io::Result<R>>(&self, inner: F) -> std::io::Result<R> {
        (inner)().map_err(|e| wrap_file_error(e, self.filename))
    }
}

pub struct VariantContext<'a> {
    chrom: &'a str,
    position: usize,
}

impl<'a> VariantContext<'a> {
    pub fn new(chrom: &'a str, position: usize) -> Self {
        VariantContext { chrom, position }
    }
}

impl<'a> Context for VariantContext<'a> {
    fn with<R, F: FnOnce() -> std::io::Result<R>>(&self, inner: F) -> std::io::Result<R> {
        (inner)().map_err(|e| {
            as_io_error(crate::errors::SveltError::BadVariant(
                self.chrom.to_string(),
                self.position,
                Box::new(e),
            ))
        })
    }
}
