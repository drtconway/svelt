use std::{cmp::Ordering, collections::BTreeMap};

use datafusion::prelude::SessionContext;

use crate::{
    features::FeatureIndex,
    options::{CommonOptions, QueryOptions},
    sequence::{SequenceIterator, fasta::FastaSequenceIterator, vcf::VcfSequenceIterator},
};

pub async fn find_similar(
    features: &str,
    query: &QueryOptions,
    common: &CommonOptions,
    ctx:&SessionContext
) -> std::io::Result<()> {
    let _ = common;
    if let Some(query) = &query.query {
        return find_similar_single(features, query, ctx).await;
    }

    let idx = FeatureIndex::load(features, ctx).await?;

    if let Some(query_file) = &query.query_file {
        let itr = FastaSequenceIterator::new(&query_file)?;
        find_similar_inner(itr, &idx)?;
    }

    if let Some(vcf_file) = &query.vcf {
        let itr = VcfSequenceIterator::new(&vcf_file)?;
        find_similar_inner(itr, &idx)?;
    }

    Ok(())
}

fn find_similar_inner<Itr: SequenceIterator>(itr: Itr, idx: &FeatureIndex) -> std::io::Result<()> {
    for rec in itr {
        let (name, sequence) = rec?;

        let res = find_similar_compile_results(&sequence, idx)?;

        for (nix, (fwd, rev)) in res {
            println!("{}\t{}\t{}\t{}", name, idx.names[nix as usize], fwd, rev);
        }
    }

    Ok(())
}

async fn find_similar_single(features: &str, query: &str, ctx: &SessionContext) -> std::io::Result<()> {
    let idx = FeatureIndex::load(features, ctx).await?;

    let res = find_similar_compile_results(query, &idx)?;

    for (nix, (fwd, rev)) in res {
        println!("{}\t{}\t{}", idx.names[nix as usize], fwd, rev);
    }

    Ok(())
}

fn find_similar_compile_results(sequence: &str, idx: &FeatureIndex) -> std::io::Result<Vec<(u32, (f64, f64))>> {
    let (fwd, rev) = idx.rank(sequence)?;

    let mut res = BTreeMap::new();
    for (nix, score) in fwd {
        if score < 0.5 {
            continue;
        }
        res.entry(nix).or_insert_with(|| (0.0, 0.0)).0 = score;
    }
    for (nix, score) in rev {
        if score < 0.5 {
            continue;
        }
        res.entry(nix).or_insert_with(|| (0.0, 0.0)).1 = score;
    }

    let mut res: Vec<(u32, (f64, f64))> = res.into_iter().collect();
    res.sort_by( |lhs, rhs| cmp_items(rhs, lhs));

    Ok(res)
}

fn cmp_items(lhs: &(u32, (f64, f64)), rhs: &(u32, (f64, f64))) -> Ordering {
    let lhs_max = if lhs.1.0 > lhs.1.1 { lhs.1.0 } else { lhs.1.1 };
    let rhs_max = if rhs.1.0 > rhs.1.1 { rhs.1.0 } else { rhs.1.1 };
    if let Some(res) = lhs_max.partial_cmp(&rhs_max) {
        res
    } else {
        lhs.0.cmp(&rhs.0)
    }
}