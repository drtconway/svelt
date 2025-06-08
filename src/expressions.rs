use datafusion::prelude::{case, lit, Expr};


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

