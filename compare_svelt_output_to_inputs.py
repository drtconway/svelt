# Usage (sniffles files)
# python3 compare_svelt_output_to_inputs.py \
#   --inputs NA20762.sniffles.vcf NA20809.sniffles.vcf NA20827.sniffles.vcf \
#   --svelt x_svelt.vcf

import argparse
import pysam
import os
from collections import defaultdict

def extract_input_ids(input_files):
    sample_to_ids = {}
    all_ids = set()
    for path in input_files:
        sample = os.path.basename(path).split('.')[0]
        vcf = pysam.VariantFile(path)
        ids = set()
        for rec in vcf:
            if rec.id and rec.id != ".":
                clean_id = rec.id.split('_', 1)[-1] if '_' in rec.id else rec.id
                ids.add(clean_id)
                all_ids.add(clean_id)
        sample_to_ids[sample] = ids
    return sample_to_ids, all_ids

def extract_svelt_ids(svelt_vcf_path):
    svelt_ids = set()
    vcf = pysam.VariantFile(svelt_vcf_path)
    for rec in vcf:
        original_ids = rec.info.get("ORIGINAL_IDS")
        if not original_ids:
            continue
        if isinstance(original_ids, (list, tuple)):
            svelt_ids.update(str(oid) for oid in original_ids if oid != '.')
        else:
            svelt_ids.update(oid for oid in original_ids.split(',') if oid != '.')
    return svelt_ids

def compare_ids(sample_to_ids, all_input_ids, svelt_ids):
    unexpected = svelt_ids - all_input_ids
    if unexpected:
        raise ValueError(
            f"Found {len(unexpected)} ORIGINAL_IDS in SVelt output not present in any input.\n"
            f"Example unexpected IDs: {list(unexpected)[:10]}"
        )

    for sample, ids in sample_to_ids.items():
        matched = ids & svelt_ids
        missing = ids - svelt_ids
        print(f"\nSample: {sample}")
        print(f" - Total input IDs: {len(ids)}")
        print(f" - Matched in SVelt output: {len(matched)}")
        print(f" - Missing from SVelt output: {len(missing)}")
        if missing:
            print(f"Example missing IDs: {list(missing)[:10]}")

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--inputs', nargs='+', required=True, help='List of input VCF files')
    parser.add_argument('--svelt', required=True, help='SVelt merged VCF file')
    args = parser.parse_args()

    print("Extracting input IDs...")
    sample_to_ids, all_input_ids = extract_input_ids(args.inputs)

    print("Extracting ORIGINAL_IDS from SVelt VCF...")
    svelt_ids = extract_svelt_ids(args.svelt)

    print("\nComparing input IDs to SVelt output...")
    compare_ids(sample_to_ids, all_input_ids, svelt_ids)

if __name__ == "__main__":
    main()
