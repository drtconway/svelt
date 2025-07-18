#!/usr/bin/env python3

import argparse
import pysam
import os

def extract_ids_from_vcf(vcf_path):
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
    sample_to_ids = {}
    for path in input_files:
        sample, ids = extract_ids_from_vcf(path)
        sample_to_ids[sample] = ids
    return sample_to_ids

def extract_svelt_ids_by_sample(svelt_vcf_path, num_samples):
    ids_by_sample = [set() for _ in range(num_samples)]
    with pysam.VariantFile(svelt_vcf_path) as vcf:
        for rec in vcf:
            original_ids = rec.info.get("ORIGINAL_IDS")
            if not original_ids:
                print(f"Warning: Missing ORIGINAL_IDS for record at {rec.chrom}:{rec.pos}")
                continue
            if isinstance(original_ids, str):
                original_ids = original_ids.split(",")
            for i, oid in enumerate(original_ids):
                if i >= num_samples:
                    print(f"Warning: ORIGINAL_IDS length ({len(original_ids)}) exceeds number of inputs ({num_samples})")
                    break
                if oid != ".":
                    ids_by_sample[i].add(oid)
    return ids_by_sample

def compare_ids(sample_to_ids, svelt_ids_by_sample, fail_on_error=True):
    all_errors = []

    for i, (sample, input_ids) in enumerate(sample_to_ids.items()):
        svelt_ids = svelt_ids_by_sample[i]
        matched = input_ids & svelt_ids
        missing = input_ids - svelt_ids
        unexpected = svelt_ids - input_ids

        print(f"\nSample {i + 1} ({sample}):")
        print(f" - Total input IDs:       {len(input_ids)}")
        print(f" - Matched in SVelt:      {len(matched)}")
        print(f" - Missing from SVelt:    {len(missing)}")
        print(f" - Unexpected in SVelt:   {len(unexpected)}")

        if missing:
            print(f"Example missing IDs: {list(missing)[:5]}")
        if unexpected:
            print(f"Example unexpected IDs: {list(unexpected)[:5]}")
            all_errors.append(
                f"Sample {sample} has {len(unexpected)} unexpected ID(s): {list(unexpected)[:5]}"
            )

    if all_errors and fail_on_error:
        raise ValueError("\nComparison failed due to unexpected IDs:\n" + "\n".join(all_errors))

def main():
    parser = argparse.ArgumentParser(description="Compare SVelt ORIGINAL_IDS to input VCFs")
    parser.add_argument('--inputs', nargs='+', required=True, help='Input VCF files (Sniffles)')
    parser.add_argument('--svelt', required=True, help='SVelt merged VCF file')
    parser.add_argument('--no-error', action='store_true', help='Warn instead of failing on unexpected IDs')
    args = parser.parse_args()

    print("Extracting input IDs...")
    sample_to_ids = extract_input_ids(args.inputs)

    print("Extracting ORIGINAL_IDS from SVelt VCF...")
    svelt_ids_by_sample = extract_svelt_ids_by_sample(args.svelt, len(args.inputs))

    print("\nComparing each input to corresponding ORIGINAL_IDS field in SVelt...")
    try:
        compare_ids(sample_to_ids, svelt_ids_by_sample, fail_on_error=not args.no_error)
    except ValueError as e:
        print(str(e))
        exit(1)

if __name__ == "__main__":
    main()
