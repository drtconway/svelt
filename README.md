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