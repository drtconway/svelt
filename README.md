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
   within 25bp, and end2 to be within 150bp.
3. If BND variants have end2 on the same chromosome and within 25bp, and
   the end on the same chromosome but within 150bp, we flip them and merge.
   This can only be done if a genomic reference sequence is supplied.

The use of a window of 25bp is not arbitrary. Many SVs are mediated by
mobile element sequences which have short repeat sequences at the end
which are often duplicated in the SV formation process (TSDs in the literature -
terminal sequence duplications). It is common for SVs to be called on
either end of the TSD in a manner that often depends on the data itself.
Accordingly, most times when SVs are called slightly differently, it is
within the window of the TSD length (7-22bp).

The sequence length ratio of 0.9 is chosen as a proxy for constraining the
sequence involved to being highly homologous.

## Output Details

- The QUAL field is taken as the maximum score across the merged records.
- The FILTER column is taken as the union of the FILTER values across the
  merged recrods.
- If ALT sequences are being replaced with `<ALT>` tags, an INFO field
  `SVELT_ALT_SEQ` is generated with the sequences. NB they do not include
  the context base at the start which is not part of the insertion.
- An INFO tag `SVELT_CRITERIA` is generated which contains the criteria
  used for merging the given alleles.
- If an index of features is supplied, insertion sequences (if present)
  are classified to show the best matching feature, which is included in
  the INFO field `SVELT_ALT_CLASS`.

## TODO

## Output Generation

- Currently the INFO column is naively populated from the left-most INFO. For
  standard INFO fields, we should consider merging.
- It would be nice to compute AN, AC, and AF from the samples we merge and
  include them in the INFO.
- The ID column is currently populated from the left-most VCF with a variant.
  IDs are required to be uniquely identifying within a VCF, but not between
  VCFs, so it's possible (probable!) that we generate non-unique IDs in the
  output. The proposed solution is to put the IDs from the original VCFs into
  an INFO field, and generate a new ID for each row in the output in such a
  way as to ensure it is unique.

### Additional Merging Rules

- Sometimes a DUP may be called as an INS because insufficient homology was
  detected, however DUPs at the same location might be considered a prior to
  promote the INS to a DUP, and merge.
