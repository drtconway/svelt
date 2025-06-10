use datafusion::prelude::{case, lit, DataFrame, Expr};


pub fn ifelse(cond: Expr, then: Expr, otherwise: Expr) -> datafusion::common::Result<Expr> {
    case(cond)
        .when(lit(true), then)
        .when(lit(false), otherwise)
        .end()
}

pub fn pmin(lhs: Expr, rhs: Expr) -> datafusion::common::Result<Expr> {
    ifelse(lhs.clone().lt_eq(rhs.clone()), lhs, rhs)
}

pub fn pmax(lhs: Expr, rhs: Expr) -> datafusion::common::Result<Expr> {
    ifelse(lhs.clone().gt_eq(rhs.clone()), lhs, rhs)
}

pub fn prefix_cols(df: DataFrame, prefix: &str) -> std::io::Result<DataFrame> {
    let mut df = df;
    for column in df.schema().columns().iter() {
        let old_name = column.name();
        let new_name = format!("{}_{}", prefix, old_name);
        df = df.with_column_renamed(old_name, &new_name)?;
    }
    Ok(df)
}