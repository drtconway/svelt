use std::collections::BTreeMap;

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

        let (fwd, rev) = idx.rank(&sequence)?;

        let mut res = BTreeMap::new();
        for (nix, score) in fwd {
            res.entry(nix).or_insert_with(|| (0.0, 0.0)).0 = score;
        }
        for (nix, score) in rev {
            res.entry(nix).or_insert_with(|| (0.0, 0.0)).1 = score;
        }

        for (nix, (fwd, rev)) in res {
            println!("{}\t{}\t{}\t{}", name, idx.names[nix as usize], fwd, rev);
        }
    }

    Ok(())
}

async fn find_similar_single(features: &str, query: &str, ctx: &SessionContext) -> std::io::Result<()> {
    let idx = FeatureIndex::load(features, ctx).await?;

    let (fwd, rev) = idx.rank(query)?;

    let mut res = BTreeMap::new();
    for (nix, score) in fwd {
        res.entry(nix).or_insert_with(|| (0.0, 0.0)).0 = score;
    }
    for (nix, score) in rev {
        res.entry(nix).or_insert_with(|| (0.0, 0.0)).1 = score;
    }

    for (nix, (fwd, rev)) in res {
        println!("{}\t{}\t{}", idx.names[nix as usize], fwd, rev);
    }

    Ok(())
}
