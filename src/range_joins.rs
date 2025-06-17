use std::io::{Error, ErrorKind};

use datafusion::{
    common::JoinType,
    prelude::{DataFrame, Expr, lit},
};

pub async fn overlap_join(
    lhs: DataFrame,
    rhs: DataFrame,
    lhs_eq_cols: &[&str],
    rhs_eq_cols: &[&str],
    lhs_range: (Expr, Expr),
    rhs_range: (Expr, Expr),
    filter: Option<Expr>,
) -> std::io::Result<DataFrame> {
    let mut join_filter = lit(true)
        .and(lhs_range.0.lt_eq(rhs_range.1))
        .and(rhs_range.0.lt_eq(lhs_range.1));
    if let Some(filter) = filter {
        join_filter = join_filter.and(filter);
    }
    lhs.join(
        rhs,
        JoinType::Inner,
        lhs_eq_cols,
        rhs_eq_cols,
        Some(join_filter),
    )
    .map_err(|e| Error::new(ErrorKind::Other, e))
}
