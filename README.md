# Svelt - Structural Variant VCF merging

Structural Variants (SVs) stored in VCFs are not always easy to
merge. With small variants (SNVs), tools such as `bcftools` do a fine
job, because for the most part, variants are called in a consistent way
in different samples.  SVs can be problematic because the same variant
can sometimes be called in a slightly different position in two different
samples (e.g. a child and a parent). This in turn, can make it difficult
to correctly interpret a variant as either de novo, or inherited. It
also makes compiling useful population statistics difficult, since the
differences in how a variant is called, can "smear out" the population
frequency.

Consequently, a tool for merging VCFs into a single VCF has to do
something more than the basic merging of identical variants as done by
`bcftools`, and needs to identify variants that are not identical,
but nonetheless can be judged as representing the same underlying event.

To that end, we present `svelt`. It is built on high performance Rust
components (namely [datafusion](https://datafusion.apache.org/) and
[noodles](https://github.com/zaeleus/noodles)) to provide high performance
and rigorous specification conformance.

## Quick Start

For the impatient:
```bash
git clone https://github.com/drtconway/svelt.git
cd svelt
cargo build --release
./target/release/svelt merge --out output.vcf child.vcf.gz parent1.vcf.gz parent2.vcf.gz
```

## Merging Rules

1. If two variants are the same, merge them. For non-BND variants, this
   entails the chrom, start, end, length and hash of the ALT allele to
   be the same. For BND variants, the chrom, end, chrom2, and end2 must
   be the same (chrom2 and end2 are derived from the ALT tag).
2. If two variants are almost the same, merge them. For non-BND variants,
   of the same type, on the same chromosome, if the starts are within 25bp
   and the ends are within 25bp, we compute the length ratio (shorter/longer),
   and require that to be > 0.9. For BND variants, we require chrom
   to be the same on both, chrom2 to be the same on both, end to be
   within 25bp, and end2 to be within 25bp.

The use of a window of 25bp is not arbitrary. Many SVs are mediated by
mobile element sequences which have short repeat sequences at the end
which are often duplicated in the SV formation process (TSDs in the literature -
terminal sequence duplications). It is common for SVs to be called on
either end of the TSD in a manner that often depends on the data itself.
Accordingly, most times when SVs are called slightly differently, it is
within the window of the TSD length (7-22bp).

The sequence length ratio of 0.9 is chosen as a proxy for constraining the
sequence involved to being highly homologous.

## TODO

## Output Generation

- Currently the INFO column is naively populated from the left-most INFO. For
  standard INFO fields, we should consider merging.
- QUAL is taken from the left-most variant. It should probably use the max for
  the non-null variants
- FILTERS is taken from the left-most variant. It should probably be the union.
- It would be nice to generate INFO tags documenting which rule was used to merge
  the variants.
- It would be nice to generate INFO tags containing the relevant sequence for
  INS/DUP/INV/DEL variants, where this is known.

### Additional Merging Rules

- For BNDs, if the *here* end is the same or close, the *there* end might be
  allowed more latitude (e.g. a couple of hundred bp).
- For BNDs, if the *here* ends are not close, but the *there* ends are, maybe
  we should flip the representation, and merge the flipped representation.
- Sometimes a DUP may be called as an INS because insufficient homology was
  detected, however DUPs at the same location might be considered a prior to
  promote the INS to a DUP, and merge.
