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

def extract_svelt_ids_by_sample(svelt_vcf_path, num_samples):
    """Extracts ORIGINAL_IDS from SVelt output, preserving per-sample position."""
    ids_by_sample = [set() for _ in range(num_samples)]
    with pysam.VariantFile(svelt_vcf_path) as vcf:
        for rec in vcf:
            original_ids = rec.info.get("ORIGINAL_IDS")
            if not original_ids:
                print(f"Warning: Missing ORIGINAL_IDS at {rec.chrom}:{rec.pos}")
                print(str(rec))
                continue
            if isinstance(original_ids, str):
                ids = original_ids.split(",")
            else:
                ids = list(original_ids)
            if len(ids) != num_samples:
                print(f"Warning: ORIGINAL_IDS has {len(ids)} entries, expected {num_samples} at {rec.chrom}:{rec.pos}")
                continue
            for i, oid in enumerate(ids):
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
    parser = argparse.ArgumentParser(description="Compare SVelt ORIGINAL_IDS to input VCFs (position-aware)")
    parser.add_argument('--inputs', nargs='+', required=True, help='Input VCF files (Sniffles)')
    parser.add_argument('--svelt', required=True, help='SVelt merged VCF file')
    parser.add_argument('--no-error', action='store_true', help='Warn instead of failing on unexpected IDs')
    args = parser.parse_args()

    print("Extracting input IDs...")
    sample_names, input_ids_by_sample = extract_input_ids(args.inputs)

    print("Extracting ORIGINAL_IDS from SVelt VCF...")
    svelt_ids_by_sample = extract_svelt_ids_by_sample(args.svelt, len(args.inputs))

    print("Comparing input IDs to SVelt ORIGINAL_IDS...")
    try:
        compare_ids_by_sample(sample_names, input_ids_by_sample, svelt_ids_by_sample, fail_on_error=not args.no_error)
    except ValueError as e:
        print(str(e))
        exit(1)

if __name__ == "__main__":
    main()
