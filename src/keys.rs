

use blake2::{Blake2b512, Digest};
use noodles::vcf::{variant::{record::info::field::Value, record::AlternateBases, Record as VariantRecord}, Header, Record};

use crate::{chroms::ChromSet, errors::{as_io_error, SveltError}};



#[derive(Debug, PartialEq, Eq, PartialOrd, Ord)]
pub enum VariantKind {
    INS, DEL, DUP, INV, BND, CPX
}

impl TryFrom<&str> for VariantKind {
    type Error = SveltError;

    fn try_from(value: &str) -> Result<Self, Self::Error> {
        match value {
            "INS" => Ok(VariantKind::INS),
            "DEL" => Ok(VariantKind::DEL),
            "DUP" => Ok(VariantKind::DUP),
            "INV" => Ok(VariantKind::INV),
            "BND" => Ok(VariantKind::BND),
            "CPX" => Ok(VariantKind::CPX),
            _ => Err(SveltError::BadKind(String::from(value)))
        }
    }
}



#[derive(Debug, PartialEq, Eq, PartialOrd, Ord)]
pub struct VariantKey {
    chrom_id: usize,
    start: usize,
    end: usize,
    kind: VariantKind,
    length: usize,
    seq: u64
}

impl VariantKey {
    pub fn new(rec: &Record, header: &Header, chroms: &ChromSet) -> std::io::Result<VariantKey> {
        let chrom_id = chroms.index(rec.reference_sequence_name()).unwrap();
        let start = rec.variant_start().unwrap()?.get();
        let end = rec.variant_end(header)?.get();
        let kind = if let Some(Value::String(value)) = rec.info().get(header, "SVTYPE").unwrap()? {
            VariantKind::try_from(value.as_ref())
        } else {
            Err(SveltError::MissingType(String::from(rec.reference_sequence_name()), start))
        }.map_err(as_io_error)?;
        let length = rec.variant_span(header)?;
        let seq = if let Some(alt) = rec.alternate_bases().iter().next() { 
            let alt = alt?;
            if VariantKey::is_seq(alt) { VariantKey::digest(alt) } else { 0 }
        } else { 0 };
        Ok(VariantKey { chrom_id, start, end, kind, length, seq })
    }

    pub fn is_seq(seq: &str) -> bool {
        seq.chars().all(|c| match c {
            'A' | 'C' | 'G' | 'T' | 'N' => true,
            'a' | 'c' | 'g' | 't' | 'n' => true,
            '*' => true,
            _ => false
        })
    }

    pub fn digest(seq: &str) -> u64 {
       let mut hasher = Blake2b512::new();
        hasher.update(seq.as_bytes());
        let dig = hasher.finalize();
        let mut h = 0;
        for i in 0..8 {
            h = (h << 8) | (dig[i] as u64);
        }
        h

    }
}