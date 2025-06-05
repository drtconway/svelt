use std::{
    fmt::Display,
    io::{Error, ErrorKind},
};

#[derive(Debug)]
pub enum SveltError {
    BadBreakEnd(String),
    BadChr2(String, usize, String, String),
    BadInfoType(String, usize, String, String),
    BadKind(String),
    Contigs(usize, usize),
    ContigMissing(String, usize),
    ContigOrder(String, usize, usize),
    MissingAlt(String, usize),
    MissingChr2(String, usize),
    MissingInfo(String, usize, String),
    MissingType(String, usize),
    NeardexDuplicate(u32),
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
            SveltError::BadKind(kind) => {
                write!(f, "Unexpected SVTYPE '{}'", kind)
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
            SveltError::MissingType(chrom, pos) => {
                write!(f, "Variant without SVTYPE at {}:{}", chrom, pos)
            }
            SveltError::NeardexDuplicate(key) => {
                write!(f, "Cannot construct Neardex with duplicate key {}", key)
            }
        }
    }
}

impl std::error::Error for SveltError {}

pub fn as_io_error(error: SveltError) -> std::io::Error {
    Error::new(ErrorKind::Other, error)
}
