use noodles::vcf::{
    Header, Record,
    variant::record::{AlternateBases, info::field::Value},
};

use crate::{
    breakends::{BreakEndSide, parse_breakend},
    errors::{SveltError, as_io_error},
    tables::is_seq,
};

pub fn get_svtype(rec: &Record, header: &Header) -> std::io::Result<String> {
    match rec.info().get(header, "SVTYPE") {
        Some(value) => {
            let value = value?;
            if let Some(Value::String(value)) = value {
                let item = value.split(':').next();
                match item {
                    Some(kind) => {
                        return Ok(String::from(kind));
                    }
                    None => {
                        return Err(SveltError::BadKind(value.to_string())).map_err(as_io_error);
                    }
                }
            }
        },
        None => {
            // No SVTYPE, so we will have to try and infer it
        },
    }

    let ref_ = rec.reference_bases();

    if rec.alternate_bases().len() > 1 {
        return Err(SveltError::MissingType).map_err(as_io_error);
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

    Err(SveltError::MissingType).map_err(as_io_error)
}

pub fn get_breakend(
    rec: &Record,
) -> std::io::Result<Option<(String, usize, BreakEndSide, BreakEndSide)>> {
    if let Some(alt) = rec.alternate_bases().iter().next() {
        let alt = alt?;

        if let Ok(bnd) = parse_breakend(alt) {
            return Ok(Some(bnd));
        }
    }
    Ok(None)
}
