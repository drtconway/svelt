#!/usr/bin/env python3

import argparse
import pysam
import os

def extract_ids_from_vcf(vcf_path):
    """Extracts a set of IDs from a single VCF file."""
    sample = os.path.basename(vcf_path).split('.')[0]
    ids = set()
    with pysam.VariantFile(vcf_path) as vcf:
        for rec in vcf:
            if rec.id and rec.id != ".":
                if rec.id in ids:
                    print(f"Warning: Duplicate ID '{rec.id}' in {vcf_path}")
                ids.add(rec.id)
    return sample, ids

def extract_input_ids(input_files):
    """Extracts IDs for each input file and returns a list preserving sample order."""
    sample_names = []
    sample_to_ids = []
    for path in input_files:
        sample, ids = extract_ids_from_vcf(path)
        sample_names.append(sample)
        sample_to_ids.append(ids)
    return sample_names, sample_to_ids

def extract_svelt_ids_by_sample(svelt_vcf_path, sample_names):
    """
    Extracts ORIGINAL_IDS or IDLIST_EXT from SVelt output,
    assigning them to samples with non-reference genotypes.
    """
    ids_by_sample = [set() for _ in sample_names]
    with pysam.VariantFile(svelt_vcf_path) as vcf:
        for rec in vcf:
            # Extract original input IDs from INFO field
            id_field = rec.info.get("ORIGINAL_IDS") or rec.info.get("IDLIST_EXT")
            if not id_field:
                continue

            if isinstance(id_field, str):
                original_ids = id_field.split(",")
            else:
                original_ids = list(id_field)

            # Assign each original ID to samples with non-reference genotype
            for i, sample in enumerate(sample_names):
                gt = rec.samples[sample].get("GT")
                if gt and any(allele not in (0, None) for allele in gt):
                    for oid in original_ids:
                        if oid != ".":
                            ids_by_sample[i].add(oid)
    return ids_by_sample

def compare_ids_by_sample(sample_names, input_ids_by_sample, svelt_ids_by_sample, fail_on_error=True):
    print("\n=== SVelt Comparison by Sample ===")
    for i, sample in enumerate(sample_names):
        input_ids = input_ids_by_sample[i]
        svelt_ids = svelt_ids_by_sample[i]
        matched = input_ids & svelt_ids
        missing = input_ids - svelt_ids
        unexpected = svelt_ids - input_ids

        print(f"\nSample: {sample}")
        print(f" - Input variant IDs:   {len(input_ids)}")
        print(f" - Matched in SVelt:    {len(matched)}")
        print(f" - Missing from SVelt:  {len(missing)}")
        print(f" - Unexpected in SVelt: {len(unexpected)}")

        if missing:
            print(f"   Example missing: {list(missing)[:5]}")
        if unexpected:
            print(f"   Example unexpected: {list(unexpected)[:5]}")
            if fail_on_error:
                raise ValueError(
                    f"Sample {sample}: Found {len(unexpected)} unexpected ID(s) in SVelt output.\n"
                    f"Example: {list(unexpected)[:5]}"
                )

def main():
    parser = argparse.ArgumentParser(description="Compare SVelt ORIGINAL_IDS or IDLIST_EXT to input VCFs (genotype-aware)")
    parser.add_argument('--inputs', nargs='+', required=True, help='Input VCF files (Sniffles)')
    parser.add_argument('--svelt', required=True, help='SVelt merged VCF file')
    parser.add_argument('--no-error', action='store_true', help='Warn instead of failing on unexpected IDs')
    args = parser.parse_args()

    print("Extracting input IDs...")
    sample_names, input_ids_by_sample = extract_input_ids(args.inputs)

    print("Extracting ORIGINAL_IDS/IDLIST_EXT from SVelt VCF...")
    svelt_ids_by_sample = extract_svelt_ids_by_sample(args.svelt, sample_names)

    print("Comparing input IDs to SVelt output...")
    try:
        compare_ids_by_sample(sample_names, input_ids_by_sample, svelt_ids_by_sample, fail_on_error=not args.no_error)
    except ValueError as e:
        print(str(e))
        exit(1)

if __name__ == "__main__":
    main()
