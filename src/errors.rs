use std::{
    fmt::Display,
    io::{Error, ErrorKind},
};

#[derive(Debug)]
pub enum SveltError {
    BadBreakEnd(String),
    BadChr2(String, usize, String, String),
    BadKind(String),
    Contigs(usize, usize),
    ContigMissing(String, usize),
    ContigOrder(String, usize, usize),
    MissingAlt(String, usize),
    MissingChr2(String, usize),
    MissingType(String, usize),
}

impl Display for SveltError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            SveltError::BadBreakEnd(alt) => write!(f, "Badly formed breakend '{}'", alt),
            SveltError::BadChr2(chrom, pos, chrom2, chr2) => {
                write!(f, "Bad CHR2 for variant at {}:{} - BND had {}, CHR2 had {}", chrom, pos, chr2, chrom2)
            }
            SveltError::BadKind(kind) => write!(f, "Unexpected SVTYPE '{}'", kind),
            SveltError::Contigs(exp, got) => write!(
                f,
                "Unexpected number of contigs - expected {}, got {}",
                exp, got
            ),
            SveltError::ContigMissing(chrom, got) => {
                write!(f, "Unexpected contig '{}' (number {}", chrom, got)
            }
            SveltError::ContigOrder(chrom, exp, got) => write!(
                f,
                "Contig '{}' out of order - expected {}, got {}",
                chrom, exp, got
            ),
            SveltError::MissingAlt(chrom, pos) => write!(f, "No ALT present at {}:{}", chrom, pos),
            SveltError::MissingChr2(chrom, pos) => {
                write!(f, "Breakend variant without CHR2 at {}:{}", chrom, pos)
            }
            SveltError::MissingType(chrom, pos) => {
                write!(f, "Variant without SVTYPE at {}:{}", chrom, pos)
            }
        }
    }
}

impl std::error::Error for SveltError {}

pub fn as_io_error(error: SveltError) -> std::io::Error {
    Error::new(ErrorKind::Other, error)
}
