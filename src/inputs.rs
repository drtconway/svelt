use noodles::vcf::{
    Header, Record,
    variant::record::{AlternateBases, info::field::Value},
};

use crate::{
    breakends::{parse_breakend, BreakEndSide}, errors::{as_io_error, SveltError}, tables::is_seq
};

pub fn get_svtype(rec: &Record, header: &Header) -> std::io::Result<String> {
    if let Some(Value::String(value)) = rec.info().get(header, "SVTYPE").unwrap()? {
        return Ok(String::from(value));
    }

    let ref_ = rec.reference_bases();

    if rec.alternate_bases().len() > 1 {
        let start = rec.variant_start().unwrap()?.get();
        return Err(SveltError::MissingType(
            String::from(rec.reference_sequence_name()),
            start,
        ))
        .map_err(as_io_error);
    }
    if let Some(alt) = rec.alternate_bases().iter().next() {
        let alt = alt?;

        if ref_.len() == 1 && is_seq(alt) && alt.len() > 1 {
            return Ok(String::from("INS"));
        }

        if ref_.len() > 1 && is_seq(alt) && alt.len() == 1 {
            return Ok(String::from("DEL"));
        }

        if let Ok(_) = parse_breakend(alt) {
            return Ok(String::from("BND"));
        }
    }

    let start = rec.variant_start().unwrap()?.get();
    Err(SveltError::MissingType(
        String::from(rec.reference_sequence_name()),
        start,
    ))
    .map_err(as_io_error)
}

pub fn get_breakend(rec: &Record) -> std::io::Result<Option<(String, usize, BreakEndSide, BreakEndSide)>> {
    if let Some(alt) = rec.alternate_bases().iter().next() {
        let alt = alt?;

        if let Ok(bnd) = parse_breakend(alt) {
            return Ok(Some(bnd));
        }
    }
    Ok(None)
}