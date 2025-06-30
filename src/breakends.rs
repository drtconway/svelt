use std::cmp::max;

use datafusion::{
    common::JoinType,
    prelude::{DataFrame, col, concat_ws, lit},
};
use noodles::{core::Position, fasta::Repository};
use regex::Regex;

use crate::{errors::SveltError, expressions::prefix_cols};

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
    pub chrom: String,
    pub end: usize,
    pub side: BreakEndSide,
    pub chrom2: String,
    pub end2: usize,
    pub side2: BreakEndSide,
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

    pub fn format(&self, repo: &Repository) -> std::io::Result<(String, usize, char, String)> {
        let BreakEnd {
            chrom,
            end,
            side,
            chrom2,
            end2,
            side2,
        } = self;
        let end = max(1, *end); // XXX 0 is actually valid here.
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
            }
            (BreakEndSide::After, BreakEndSide::Before) => {
                format!("{}]{}:{}]", b, chrom2, end2)
            }
            (BreakEndSide::After, BreakEndSide::After) => {
                format!("{}[{}:{}[", b, chrom2, end2)
            }
        };
        Ok((chrom.clone(), end, b, alt))
    }
}

pub(crate) async fn unpaired_breakend_check(tbl: DataFrame) -> std::io::Result<DataFrame> {
    let df = tbl
        .clone()
        .filter(col("kind").eq(lit("BND")))?
        .with_column(
            "here_there",
            concat_ws(
                lit("_"),
                vec![col("chrom"), col("end"), col("chrom2"), col("end2")],
            ),
        )?
        .with_column(
            "there_here",
            concat_ws(
                lit("_"),
                vec![col("chrom2"), col("end2"), col("chrom"), col("end")],
            ),
        )?;

    let lhs = prefix_cols(df.clone(), "lhs")?;
    let rhs = prefix_cols(df.clone(), "rhs")?;

    let paired = lhs
        .join(
            rhs,
            JoinType::Left,
            &["lhs_vix", "lhs_here_there"],
            &["rhs_vix", "rhs_there_here"],
            None,
        )?
        .with_column("paired_bnd", col("rhs_row_id").is_not_null())?
        .select_columns(&["lhs_row_id", "paired_bnd"])?;

    // paired.sort_by(vec![col("lhs_row_id")])?.show().await?;

    let tbl = tbl
        .join(paired, JoinType::Left, &["row_id"], &["lhs_row_id"], None)?
        .drop_columns(&["lhs_row_id"])?;

    Ok(tbl)
}
