use std::{fmt::Display, io::{Error, ErrorKind}};

#[derive(Debug)]
pub enum SveltError {
    BadKind(String),
    MissingType(String, usize)
}

impl Display for SveltError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            SveltError::BadKind(kind) => write!(f, "Unexpected SVTYPE '{}'", kind),
            SveltError::MissingType(chrom, pos) => write!(f, "Variant without SVTYPE at {}:{}", chrom, pos),
        }
    }
}

impl std::error::Error for SveltError {}

pub fn as_io_error(error: SveltError) -> std::io::Error {
    Error::new(ErrorKind::Other, error)
}