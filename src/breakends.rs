use std::cmp::max;

use noodles::{core::Position, fasta::Repository};
use regex::Regex;

use crate::errors::SveltError;

#[derive(Debug, PartialEq, Eq, PartialOrd, Ord, Clone, Copy)]
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

#[derive(Debug)]
pub struct BreakEnd {
    chrom: String,
    end: usize,
    side: BreakEndSide,
    chrom2: String,
    end2: usize,
    side2: BreakEndSide,
}

impl BreakEnd {
    pub fn new(chrom: &str, end: usize, alt: &str) -> std::result::Result<BreakEnd, SveltError> {
        let (chrom2, end2, side, side2) = parse_breakend(alt)?;
        Ok(BreakEnd {
            chrom: String::from(chrom),
            end,
            side,
            chrom2,
            end2,
            side2,
        })
    }

    pub fn flip(&self) -> BreakEnd {
        let BreakEnd {
            chrom,
            end,
            side,
            chrom2,
            end2,
            side2,
        } = self;
        BreakEnd {
            chrom: chrom2.clone(),
            end: *end2,
            side: *side2,
            chrom2: chrom.clone(),
            end2: *end,
            side2: *side,
        }
    }

    pub fn format(&self, repo: Repository) -> std::io::Result<(String, usize, char, String)> {
        let BreakEnd {
            chrom,
            end,
            side,
            chrom2,
            end2,
            side2,
        } = self;
        let end = max(1, *end); // Urk! Thanks for that Sniffles!
        let pos = Position::try_from(end).unwrap();
        let seq = repo.get(chrom.as_ref()).unwrap()?;
        let b: &u8 = seq.get(pos).unwrap();
        let b = *b as char;
        let alt = match (side, side2) {
            (BreakEndSide::Before, BreakEndSide::Before) => {
            format!("]{}:{}]{}", chrom2, end2, b)
            }
            (BreakEndSide::Before, BreakEndSide::After) => {
                format!("[{}:{}[{}", chrom2, end2, b)
            },
            (BreakEndSide::After, BreakEndSide::Before) => {
                format!("{}]{}:{}]", b, chrom2, end2)
            },
            (BreakEndSide::After, BreakEndSide::After) => {
                format!("{}[{}:{}[", b, chrom2, end2)
            },
        };
        Ok((chrom.clone(), end, b, alt))
    }
}
