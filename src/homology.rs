use datafusion::{
    functions_aggregate::{expr_fn::first_value, min_max::max},
    prelude::{DataFrame, SessionContext, col, lit},
};
use itertools::Itertools;

use crate::{
    features::FeatureIndex,
    kmers::Kmer,
    kmers_table::kmer_frequencies_to_table,
    options::{CommonOptions, QueryOptions, make_session_context},
    sequence::{SequenceIterator, fasta::FastaSequenceIterator, vcf::VcfSequenceIterator},
};

pub async fn find_similar(
    features: &str,
    query: &QueryOptions,
    k: usize,
    common: &CommonOptions,
) -> std::io::Result<()> {
    let ctx = make_session_context(common);

    if let Some(query) = &query.query {
        return find_similar_single(features, query, k, common).await;
    }

    let idx = FeatureIndex::load(features, &ctx).await?;

    if let Some(query_file) = &query.query_file {
        let itr = FastaSequenceIterator::new(&query_file)?;
        find_similar_inner(itr, k, &idx, &ctx).await?;
    }

    if let Some(vcf_file) = &query.vcf {
        let itr = VcfSequenceIterator::new(&vcf_file)?;
        find_similar_inner(itr, k, &idx, &ctx).await?;
    }

    Ok(())
}

async fn find_similar_inner<Itr: SequenceIterator>(
    itr: Itr,
    k: usize,
    idx: &FeatureIndex,
    ctx: &SessionContext,
) -> std::io::Result<()> {
    let mut curr: Option<DataFrame> = None;
    let mut curr_count: usize = 0;
    for rec in itr {
        let (name, sequence) = rec?;
        log::info!("scanning {} ({}bp)", name, sequence.len());

        let mut fwd = Vec::new();
        let mut rev = Vec::new();
        Kmer::with_many_both(k, &sequence, |x, y| {
            fwd.push(x.0);
            rev.push(y.0);
        });

        fwd.sort();
        let fwd: Vec<(u64, usize)> = fwd
            .into_iter()
            .chunk_by(|x| *x)
            .into_iter()
            .map(|(x, xs)| (x, xs.count()))
            .collect();
        curr_count += fwd.len();
        let fwd = kmer_frequencies_to_table(&fwd, &ctx)
            .await?
            .with_column("strand", lit("+"))?;

        rev.sort();
        let rev: Vec<(u64, usize)> = rev
            .into_iter()
            .chunk_by(|x| *x)
            .into_iter()
            .map(|(x, xs)| (x, xs.count()))
            .collect();
        curr_count += rev.len();
        let rev = kmer_frequencies_to_table(&rev, &ctx)
            .await?
            .with_column("strand", lit("-"))?;

        let block = fwd.union(rev)?.with_column("name", lit(name))?;
        curr = if let Some(block0) = curr {
            Some(block0.union(block)?)
        } else {
            Some(block)
        };

        if curr_count > 500000 {
            log::info!("ranking over {} k-mers.", curr_count);
            let tbl = curr.take().unwrap();
            let res = idx.rank_many(tbl).await?;
            res.aggregate(
                vec![col("name")],
                vec![
                    max(col("dot")).alias("dot"),
                    first_value(col("idx_name"), Some(vec![col("dot").sort(false, false)]))
                        .alias("feature"),
                    first_value(col("class"), Some(vec![col("dot").sort(false, false)]))
                        .alias("class"),
                ],
            )?
            .show()
            .await?;
            curr_count = 0;
        }
    }

    if curr.is_none() {
        return Ok(());
    }

    log::info!("ranking over {} k-mers.", curr_count);
    let all = curr.unwrap();
    let res = idx.rank_many(all).await?;

    res.aggregate(
        vec![col("name")],
        vec![
            first_value(col("idx_name"), Some(vec![col("dot").sort(false, false)]))
                .alias("feature"),
            first_value(col("class"), Some(vec![col("dot").sort(false, false)])).alias("class"),
        ],
    )?
    .show()
    .await?;

    Ok(())
}

async fn find_similar_single(
    features: &str,
    query: &str,
    k: usize,
    common: &CommonOptions,
) -> std::io::Result<()> {
    let ctx = make_session_context(common);

    let idx = FeatureIndex::load(features, &ctx).await?;

    let mut fwd = Vec::new();
    let mut rev = Vec::new();
    Kmer::with_many_both(k, &query, |x, y| {
        fwd.push(x.0);
        rev.push(y.0);
    });

    fwd.sort();
    let fwd: Vec<(u64, usize)> = fwd
        .into_iter()
        .chunk_by(|x| *x)
        .into_iter()
        .map(|(x, xs)| (x, xs.count()))
        .collect();
    let fwd = kmer_frequencies_to_table(&fwd, &ctx).await?;
    let fwd = idx.rank(fwd).await?.with_column("strand", lit("+"))?;

    rev.sort();
    let rev: Vec<(u64, usize)> = rev
        .into_iter()
        .chunk_by(|x| *x)
        .into_iter()
        .map(|(x, xs)| (x, xs.count()))
        .collect();
    let rev = kmer_frequencies_to_table(&rev, &ctx).await?;
    let rev = idx.rank(rev).await?.with_column("strand", lit("-"))?;

    fwd.union(rev)?
        .sort(vec![col("distance").sort(true, false)])?
        .show()
        .await?;

    Ok(())
}
