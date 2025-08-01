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
                ids.add(rec.id)
    return sample, ids

def extract_input_ids(input_files):
    """Extracts per-sample ID sets in order."""
    sample_names = []
    sample_to_ids = []
    for path in input_files:
        sample, ids = extract_ids_from_vcf(path)
        sample_names.append(sample)
        sample_to_ids.append(ids)
    return sample_names, sample_to_ids

def extract_jasmine_ids_by_sample(jasmine_vcf_path, num_samples):
    """Extracts IDs from Jasmine merged VCF, split per sample index."""
    ids_by_sample = [set() for _ in range(num_samples)]
    with pysam.VariantFile(jasmine_vcf_path) as vcf:
        for rec in vcf:
            idlist = rec.info.get("IDLIST_EXT") or rec.info.get("IDLIST")
            if not idlist:
                continue
            if isinstance(idlist, str):
                ids = idlist.split(",")
            else:
                ids = list(idlist)
            if len(ids) != num_samples:
                print(f"Warning: Expected {num_samples} IDs, found {len(ids)} at {rec.chrom}:{rec.pos}")
                continue
            for i, id in enumerate(ids):
                if id != ".":
                    ids_by_sample[i].add(id)
    return ids_by_sample

def compare_ids_by_sample(sample_names, input_ids_by_sample, jasmine_ids_by_sample, fail_on_error=True):
    print("\n=== Jasmine Comparison by Sample ===")
    for i, sample in enumerate(sample_names):
        input_ids = input_ids_by_sample[i]
        jasmine_ids = jasmine_ids_by_sample[i]

        matched = input_ids & jasmine_ids
        missing = input_ids - jasmine_ids
        unexpected = jasmine_ids - input_ids

        print(f"\nSample: {sample}")
        print(f" - Input variant IDs:     {len(input_ids)}")
        print(f" - Matched in Jasmine:    {len(matched)}")
        print(f" - Missing from Jasmine:  {len(missing)}")
        print(f" - Unexpected in Jasmine: {len(unexpected)}")

        if missing:
            print(f"   Example missing: {list(missing)[:5]}")
        if unexpected:
            print(f"   Example unexpected: {list(unexpected)[:5]}")
            if fail_on_error:
                raise ValueError(
                    f"Sample {sample}: {len(unexpected)} unexpected ID(s) in Jasmine output.\n"
                    f"Example: {list(unexpected)[:5]}"
                )

def main():
    parser = argparse.ArgumentParser(description="Compare Jasmine IDLIST[_EXT] to input VCFs (position-aware)")
    parser.add_argument('--inputs', nargs='+', required=True, help='Input VCF files (e.g., Sniffles)')
    parser.add_argument('--jasmine', required=True, help='Jasmine merged VCF file')
    parser.add_argument('--no-error', action='store_true', help='Warn instead of failing on unexpected IDs')
    args = parser.parse_args()

    print("Extracting input IDs...")
    sample_names, input_ids_by_sample = extract_input_ids(args.inputs)

    print("Extracting Jasmine IDs by sample...")
    jasmine_ids_by_sample = extract_jasmine_ids_by_sample(args.jasmine, len(args.inputs))

    print("Comparing input IDs to Jasmine merged output...")
    try:
        compare_ids_by_sample(sample_names, input_ids_by_sample, jasmine_ids_by_sample, fail_on_error=not args.no_error)
    except ValueError as e:
        print(str(e))
        exit(1)

if __name__ == "__main__":
    main()
