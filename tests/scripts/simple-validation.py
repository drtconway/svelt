import sys

from pysam import VariantFile

def scan_ids(filename):
    with VariantFile(filename) as vcf:
        for rec in vcf:
            if rec.id is not None:
                yield rec.id

def scan_ids_from_info(filename, tagname):
    with VariantFile(filename) as vcf:
        for rec in vcf:
            yield rec.info.get(tagname)

merged = set()
for vids in scan_ids_from_info(sys.argv[1], "IDLIST"):
    for vid in vids:
        if vid in merged:
            print(f'{sys.argv[1]}: variant ID {vid} already seen!')
        merged.add(vid)

seen = set()
for filename in sys.argv[2:]:
    single = set()
    for vid in scan_ids(filename):
        if vid in single:
            print(f'{filename}: variant ID {vid} already seen!')
        single.add(vid)

    here = merged & single
    if here != single:
        missing_from_merged = single - here
        print(f'{filename} had IDs that are missing from the merged VCF: {sorted(missing_from_merged)}')
    merged -= single

if len(merged):
    print(f'Merged had IDs not in singles: {sorted(merged)}')