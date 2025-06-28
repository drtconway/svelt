use std::{iter::zip, sync::Arc};

use datafusion::{
    arrow::{
        array::{ArrayRef, PrimitiveBuilder},
        datatypes::{DataType, Int32Type},
    },
    common::cast::as_string_array,
    error::DataFusionError,
    logical_expr::{ColumnarValue, ScalarUDF, Volatility},
    prelude::SimpleScalarUDF,
};

pub fn align() -> ScalarUDF {
    ScalarUDF::from(SimpleScalarUDF::new(
        "align",
        vec![DataType::Utf8, DataType::Utf8],
        DataType::Int32,
        Volatility::Immutable,
        Arc::new(|args| do_alignment(args))
    ))
}

fn do_alignment(args: &[ColumnarValue]) -> std::result::Result<ColumnarValue, DataFusionError> {
    let args = ColumnarValue::values_to_arrays(args)?;

    let lhss = as_string_array(&args[0])?;
    let rhss = as_string_array(&args[1])?;

    let a = NeedlemanWunsch::default();

    let mut result_builder = PrimitiveBuilder::<Int32Type>::new();
    for (lhs, rhs) in zip(lhss.iter(), rhss.iter()) {
        let res = match (lhs, rhs) {
            (None, None) => None,
            (None, Some(_)) => None,
            (Some(_), None) => None,
            (Some(lhs), Some(rhs)) => Some(a.align(lhs, rhs)),
        };
        result_builder.append_option(res);
    }

    let result_array = result_builder.finish();

    Ok(ColumnarValue::from(Arc::new(result_array) as ArrayRef))
}

pub struct NeedlemanWunsch {
    symbol_match: i32,
    symbol_mismatch: i32,
    gap_open: i32,
    gap_extend: i32,
}

impl NeedlemanWunsch {
    pub fn new(
        symbol_match: i32,
        symbol_mismatch: i32,
        gap_open: i32,
        gap_extend: i32,
    ) -> NeedlemanWunsch {
        NeedlemanWunsch {
            symbol_match,
            symbol_mismatch,
            gap_open,
            gap_extend,
        }
    }

    pub fn align(&self, lhs: &str, rhs: &str) -> i32 {
        let n = lhs.len();
        let m = rhs.len();

        let mut s: Vec<Vec<i32>> = (0..=n).map(|_| vec![0; m + 1]).collect();
        let mut ix: Vec<Vec<i32>> = (0..=n).map(|_| vec![0; m + 1]).collect();
        let mut iy: Vec<Vec<i32>> = (0..=n).map(|_| vec![0; m + 1]).collect();

        for i in 1..=n {
            s[i][0] = self.gap_open + (i as i32 - 1) * self.gap_extend;
            ix[i][0] = self.gap_open + (i as i32 - 1) * self.gap_extend;
        }
        for j in 1..=m {
            s[0][j] = self.gap_open + (j as i32 - 1) * self.gap_extend;
            iy[0][j] = self.gap_open + (j as i32 - 1) * self.gap_extend;
        }

        for (i, x) in lhs.chars().enumerate() {
            let i = i + 1;
            for (j, y) in rhs.chars().enumerate() {
                let j = j + 1;

                let v = if x == y {
                    self.symbol_match
                } else {
                    self.symbol_mismatch
                };
                let s_m = s[i - 1][j - 1] + v;
                let s_ix = ix[i - 1][j - 1] + v;
                let s_iy = iy[i - 1][j - 1] + v;
                s[i][j] = std::cmp::max(s_m, std::cmp::max(s_ix, s_iy));

                let ix_o = s[i - 1][j] + self.gap_open;
                let ix_e = ix[i - 1][j] + self.gap_extend;
                ix[i][j] = std::cmp::max(ix_o, ix_e);

                let iy_o = s[i][j - 1] + self.gap_open;
                let iy_e = iy[i][j - 1] + self.gap_extend;
                iy[i][j] = std::cmp::max(iy_o, iy_e);
            }
        }

        //for row in s.iter() {
        //    println!("{:?}", row);
        //}

        s[n][m]
    }
}

impl Default for NeedlemanWunsch {
    fn default() -> Self {
        Self {
            symbol_match: 2,
            symbol_mismatch: -3,
            gap_open: -5,
            gap_extend: -2,
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn align_1() {
        let x = "ACCGTTGA";
        let y = "ACCGTTGA";

        let a = NeedlemanWunsch::default();

        assert_eq!(a.align(x, y), 16);
    }

    #[test]
    fn align_2() {
        let x = "ACCGTTGA";
        let y = "ACCATTGA";

        let a = NeedlemanWunsch::default();

        assert_eq!(a.align(x, y), 11);
    }

    #[test]
    fn align_3() {
        let x = "ACCGTTGA";
        let y = "ACCTTGA";

        let a = NeedlemanWunsch::default();

        assert_eq!(a.align(x, y), 9);
    }

    #[test]
    fn align_4() {
        let x = "ACCGTTGA";
        let y = "AGCTTGA";

        let a = NeedlemanWunsch::default();

        assert_eq!(a.align(x, y), 4);
    }

    #[test]
    fn align_5() {
        let x = "ACTGGATGAGCTCCTCAAAGTCTCACTATGTTGCTCAGGCTGGTCTTGAACTCCTGGCCTCAAGCGATCCTCCCACCTTAGCCTCCCAAAGTGTTGGGATTATAGGCATGAGCCACTGCACCTGGCT";
        let y = "TAACTGGATGAGCTCCTCAAAGTCTCACTATGTTGCTCAGGCTGGTCTTGAACTCCTGGCCTCAAGCGATCCTCCCACCTTAGCCTCCCAAAGTGTTGGGATTATAGGCATGAGCCACTGCACCTGGCT";

        let a = NeedlemanWunsch::default();

        assert_eq!(a.align(x, y), 254);
    }
}
