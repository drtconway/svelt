#!/usr/bin/env python3

import argparse
import pysam
import os

def extract_ids_from_vcf(vcf_path):
    """Extracts a set of IDs from a single VCF file, along with a sample-to-ID mapping."""
    sample = os.path.basename(vcf_path).split('.')[0]
    ids = set()
    id_to_sample = {}
    with pysam.VariantFile(vcf_path) as vcf:
        for rec in vcf:
            if rec.id and rec.id != ".":
                if rec.id in ids:
                    print(f"Warning: Duplicate ID '{rec.id}' in {vcf_path}")
                ids.add(rec.id)
                id_to_sample[rec.id] = sample
    return sample, ids, id_to_sample

def extract_input_ids(input_files):
    """Extracts IDs for each input file and returns a list preserving sample order."""
    sample_names = []
    sample_to_ids = []
    id_to_sample = {}
    for path in input_files:
        sample, ids, sample_map = extract_ids_from_vcf(path)
        sample_names.append(sample)
        sample_to_ids.append(ids)
        id_to_sample.update(sample_map)
    return sample_names, sample_to_ids, id_to_sample

def extract_svelt_ids_by_sample(svelt_vcf_path, sample_names, id_to_sample):
    """
    Extracts ORIGINAL_IDS or IDLIST_EXT from Svelt output,
    assigning them to samples with non-reference genotypes,
    only if the ID came from that sample’s input.
    """
    ids_by_sample = [set() for _ in sample_names]
    with pysam.VariantFile(svelt_vcf_path) as vcf:
        for rec in vcf:
            id_field = rec.info.get("ORIGINAL_IDS") or rec.info.get("IDLIST_EXT")
            if not id_field:
                continue

            if isinstance(id_field, str):
                original_ids = id_field.split(",")
            else:
                original_ids = list(id_field)

            for i, sample in enumerate(sample_names):
                gt = rec.samples[sample].get("GT")
                if gt and any(allele not in (0, None) for allele in gt):
                    for oid in original_ids:
                        if oid != "." and id_to_sample.get(oid) == sample:
                            ids_by_sample[i].add(oid)
    return ids_by_sample

def compare_ids_by_sample(sample_names, input_ids_by_sample, svelt_ids_by_sample, fail_on_error=True):
    print("\n=== Svelt Comparison by Sample ===")
    for i, sample in enumerate(sample_names):
        input_ids = input_ids_by_sample[i]
        svelt_ids = svelt_ids_by_sample[i]
        matched = input_ids & svelt_ids
        missing = input_ids - svelt_ids
        unexpected = svelt_ids - input_ids

        print(f"\nSample: {sample}")
        print(f" - Input variant IDs:   {len(input_ids)}")
        print(f" - Matched in Svelt:    {len(matched)}")
        print(f" - Missing from Svelt:  {len(missing)}")
        print(f" - Unexpected in Svelt: {len(unexpected)}")

        if missing:
            print(f"   Example missing: {list(missing)[:5]}")
        if unexpected:
            print(f"   Example unexpected: {list(unexpected)[:5]}")
            if fail_on_error:
                raise ValueError(
                    f"Sample {sample}: Found {len(unexpected)} unexpected ID(s) in Svelt output.\n"
                    f"Example: {list(unexpected)[:5]}"
                )

def main():
    parser = argparse.ArgumentParser(description="Compare Svelt ORIGINAL_IDS or IDLIST_EXT to input VCFs (genotype-aware, sample-aware)")
    parser.add_argument('--inputs', nargs='+', required=True, help='Input VCF files (Sniffles)')
    parser.add_argument('--svelt', required=True, help='Svelt merged VCF file')
    parser.add_argument('--no-error', action='store_true', help='Warn instead of failing on unexpected IDs')
    args = parser.parse_args()

    print("Extracting input IDs...")
    sample_names, input_ids_by_sample, id_to_sample = extract_input_ids(args.inputs)

    print("Extracting ORIGINAL_IDS/IDLIST_EXT from Svelt VCF...")
    svelt_ids_by_sample = extract_svelt_ids_by_sample(args.svelt, sample_names, id_to_sample)

    print("Comparing input IDs to Svelt output...")
    try:
        compare_ids_by_sample(sample_names, input_ids_by_sample, svelt_ids_by_sample, fail_on_error=not args.no_error)
    except ValueError as e:
        print(str(e))
        exit(1)

if __name__ == "__main__":
    main()
