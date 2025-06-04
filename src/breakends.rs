use regex::Regex;

use crate::errors::SveltError;

#[derive(Debug, PartialEq, Eq, PartialOrd, Ord)]
pub enum BreakEndSide {
    Before,
    After,
}

pub fn parse_breakend(
    alt: &str,
) -> std::result::Result<(String, usize, BreakEndSide, BreakEndSide), SveltError> {
    if alt.starts_with("[") || alt.starts_with("]") || alt.ends_with("[") || alt.ends_with("]") {
        let bnd1 = Regex::new(r"[ACGTNacgtn]\[([^:]+):([0-9]+)\[").unwrap();
        if let Some(caps) = bnd1.captures(&alt) {
            let chrom2 = String::from(&caps[1]);
            let pos2: usize = (&caps[2]).parse().unwrap();
            return Ok((chrom2, pos2, BreakEndSide::After, BreakEndSide::After));
        }
        let bnd2 = Regex::new(r"[ACGTNacgtn]\]([^:]+):([0-9]+)\]").unwrap();
        if let Some(caps) = bnd2.captures(&alt) {
            let chrom2 = String::from(&caps[1]);
            let pos2: usize = (&caps[2]).parse().unwrap();
            return Ok((chrom2, pos2, BreakEndSide::After, BreakEndSide::Before));
        }
        let bnd3 = Regex::new(r"\]([^:]+):([0-9]+)\][ACGTNacgtn]").unwrap();
        if let Some(caps) = bnd3.captures(&alt) {
            let chrom2 = String::from(&caps[1]);
            let pos2: usize = (&caps[2]).parse().unwrap();
            return Ok((chrom2, pos2, BreakEndSide::Before, BreakEndSide::Before));
        }
        let bnd4 = Regex::new(r"\[([^:]+):([0-9]+)\[[ACGTNacgtn]").unwrap();
        if let Some(caps) = bnd4.captures(&alt) {
            let chrom2 = String::from(&caps[1]);
            let pos2: usize = (&caps[2]).parse().unwrap();
            return Ok((chrom2, pos2, BreakEndSide::Before, BreakEndSide::After));
        }
    }
    Err(SveltError::BadBreakEnd(String::from(alt)))
}
