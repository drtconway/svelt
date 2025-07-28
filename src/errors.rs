use std::{
    error::Error,
    fmt::Display,
    io::{Error as IoError, ErrorKind},
};

#[derive(Debug)]
pub enum SveltError {
    BadBreakEnd(String),
    BadChr2(String, usize, String, String),
    BadInfoType(String, usize, String, String),
    BadKind(String, usize, String),
    Contigs(usize, usize),
    ContigMissing(String, usize),
    ContigOrder(String, usize, usize),
    FileError(String, Box<dyn Error + Send + Sync + 'static>),
    MissingAlt(String, usize),
    MissingChr2(String, usize),
    MissingInfo(String, usize, String),
    MissingK(String),
    MissingType(String, usize),
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
            SveltError::BadChr2(chrom, pos, chrom2, chr2) => {
                write!(
                    f,
                    "Bad CHR2 for variant at {}:{} - BND had {}, CHR2 had {}",
                    chrom, pos, chr2, chrom2
                )
            }
            SveltError::BadInfoType(chrom, pos, tag, exp) => {
                write!(
                    f,
                    "Unexpected type for variant {}:{} - tag was {}, expected type was {}",
                    chrom, pos, tag, exp
                )
            }
            SveltError::BadKind(chrom, pos, kind) => {
                write!(f, "Unexpected SVTYPE at {}:{}: '{}'", chrom, pos, kind)
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
            SveltError::FileError(filename, error) => {
                write!(
                    f,
                    "With file '{}' the following error occurred: {}",
                    filename,
                    error.as_ref().to_string()
                )
            }
            SveltError::MissingAlt(chrom, pos) => {
                write!(f, "No ALT present at {}:{}", chrom, pos)
            }
            SveltError::MissingChr2(chrom, pos) => {
                write!(f, "Breakend variant without CHR2 at {}:{}", chrom, pos)
            }
            SveltError::MissingInfo(chrom, pos, field) => {
                write!(
                    f,
                    "Expected field '{}' not found at {}:{}",
                    field, chrom, pos
                )
            }
            SveltError::MissingK(name) => {
                write!(f, "index '{}' has missing metadata", name)
            }
            SveltError::MissingType(chrom, pos) => {
                write!(f, "Variant without SVTYPE at {}:{}", chrom, pos)
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

impl std::error::Error for SveltError {}

pub fn as_io_error(error: SveltError) -> std::io::Error {
    IoError::new(ErrorKind::Other, error)
}

pub fn wrap_file_error<E: std::error::Error + Send + Sync + 'static>(
    e: E,
    filename: &str,
) -> std::io::Error {
    eprintln!("error in {}: {:?}", filename, e);
    as_io_error(SveltError::FileError(String::from(filename), Box::new(e)))
}
