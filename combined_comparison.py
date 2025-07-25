#!/usr/bin/env python3

import argparse
import pysam
import os
import sys

def extract_ids_from_vcf(vcf_path):
    """Extracts a set of IDs from a single VCF file and builds a sample map."""
    sample = os.path.basename(vcf_path).split('.')[0]
    ids = set()
    id_to_sample = {}
    with pysam.VariantFile(vcf_path) as vcf:
        for rec in vcf:
            if rec.id and rec.id != ".":
                ids.add(rec.id)
                id_to_sample[rec.id] = sample
    return sample, ids, id_to_sample

def extract_input_ids(input_files, need_sample_map=False):
    """Extracts per-sample ID sets and optionally builds ID-to-sample map."""
    sample_names = []
    sample_to_ids = []
    id_to_sample = {}
    for path in input_files:
        sample, ids, id_map = extract_ids_from_vcf(path)
        sample_names.append(sample)
        sample_to_ids.append(ids)
        if need_sample_map:
            id_to_sample.update(id_map)
    return (sample_names, sample_to_ids, id_to_sample) if need_sample_map else (sample_names, sample_to_ids)

def extract_svelt_ids_by_sample(vcf_path, sample_names, id_to_sample):
    """Extracts ORIGINAL_IDS or IDLIST_EXT from SVelt merged VCF."""
    ids_by_sample = [set() for _ in sample_names]
    with pysam.VariantFile(vcf_path) as vcf:
        for rec in vcf:
            id_field = rec.info.get("ORIGINAL_IDS") or rec.info.get("IDLIST_EXT")
            if not id_field:
                continue
            original_ids = id_field.split(",") if isinstance(id_field, str) else list(id_field)
            for i, sample in enumerate(sample_names):
                gt = rec.samples[sample].get("GT")
                if gt and any(allele not in (0, None) for allele in gt):
                    for oid in original_ids:
                        if oid != "." and id_to_sample.get(oid) == sample:
                            ids_by_sample[i].add(oid)
    return ids_by_sample

def extract_jasmine_ids_by_sample(vcf_path, num_samples):
    """Extracts IDLIST[_EXT] from Jasmine merged VCF, assigning to each sample index."""
    ids_by_sample = [set() for _ in range(num_samples)]
    mismatched_count = 0
    total_records = 0

    with pysam.VariantFile(vcf_path) as vcf:
        for rec in vcf:
            total_records += 1
            idlist = rec.info.get("IDLIST_EXT") or rec.info.get("IDLIST")
            if not idlist:
                continue
            ids = idlist.split(",") if isinstance(idlist, str) else list(idlist)

            if len(ids) != num_samples:
                mismatched_count += 1
                continue

            for i, id_ in enumerate(ids):
                if id_ != ".":
                    ids_by_sample[i].add(id_)

    if mismatched_count > 0:
        print(f"\nJasmine IDLIST mismatch: {mismatched_count} of {total_records} records "
              f"had incorrect number of IDs (expected {num_samples}). These were skipped.\n")

    return ids_by_sample


def compare_ids_by_sample(sample_names, input_ids_by_sample, merged_ids_by_sample, label="Merged", fail_on_error=True):
    print(f"\n=== {label} Comparison by Sample ===")
    for i, sample in enumerate(sample_names):
        input_ids = input_ids_by_sample[i]
        merged_ids = merged_ids_by_sample[i]

        matched = input_ids & merged_ids
        missing = input_ids - merged_ids
        unexpected = merged_ids - input_ids

        print(f"\nSample: {sample}")
        print(f" - Input variant IDs:     {len(input_ids)}")
        print(f" - Matched in {label}:    {len(matched)}")
        print(f" - Missing from {label}:  {len(missing)}")
        print(f" - Unexpected in {label}: {len(unexpected)}")

        if missing:
            print(f"   Example missing: {list(missing)[:5]}")
        if unexpected:
            print(f"   Example unexpected: {list(unexpected)[:5]}")
            if fail_on_error:
                raise ValueError(
                    f"Sample {sample}: {len(unexpected)} unexpected ID(s) in {label} output.\n"
                    f"Example: {list(unexpected)[:5]}"
                )

def main():
    parser = argparse.ArgumentParser(description="Compare merged VCFs (Jasmine or Svelt) to input VCFs.")
    parser.add_argument('--inputs', nargs='+', required=True, help='Input VCF files (e.g., Sniffles)')
    parser.add_argument('--jasmine', help='Jasmine merged VCF file')
    parser.add_argument('--svelt', help='Svelt merged VCF file')
    parser.add_argument('--no-error', action='store_true', help='Warn instead of failing on unexpected IDs')
    args = parser.parse_args()

    if bool(args.jasmine) == bool(args.svelt):
        print("You must provide exactly one of --jasmine or --svelt.")
        sys.exit(1)

    print("Extracting input IDs...")
    if args.svelt:
        sample_names, input_ids_by_sample, id_to_sample = extract_input_ids(args.inputs, need_sample_map=True)
        print("Extracting ORIGINAL_IDS/IDLIST_EXT from Svelt VCF...")
        merged_ids_by_sample = extract_svelt_ids_by_sample(args.svelt, sample_names, id_to_sample)
        label = "Svelt"
    else:
        sample_names, input_ids_by_sample = extract_input_ids(args.inputs, need_sample_map=False)
        print("Extracting Jasmine IDs by sample...")
        merged_ids_by_sample = extract_jasmine_ids_by_sample(args.jasmine, len(sample_names))
        label = "Jasmine"

    print("Comparing input IDs to merged output...")
    try:
        compare_ids_by_sample(sample_names, input_ids_by_sample, merged_ids_by_sample, label=label, fail_on_error=not args.no_error)
    except ValueError as e:
        print(str(e))
        sys.exit(1)

if __name__ == "__main__":
    main()
