use datafusion::prelude::{lit, Expr};

/// We make a "global" key for each row of each VCF by composing the VCF number (or VCF index `vix`)
/// and the row number. The number of VCFs is constrained to 64 elsewhere, and the number of variants
/// in each VCF is (effectively) unconstrained, so we encode the VCF index in the low order bits,
/// even though it would seem more natural to encode them as the most significant bits. Further, we
/// reserve 0..100 for the `vix` so that when we print them in tables as decimal numbers, it is easy
/// to distinguish the `vix` and the row number.

pub struct RowKey {}

impl RowKey {
    /// Compose a VCF index and a row number into a single identifier.
    pub fn encode(vix: u32, rn: u32) -> u32 {
        vix + 100 * rn
    }

    /// Decompose an identifier into the VCF index and row number
    pub fn decode(key: u32) -> (u32, u32) {
        (key % 100, key / 100)
    }

    pub fn make(row_num: Expr, vix: u32) -> Expr {
        lit(vix) + row_num * lit(100)
    }
}

