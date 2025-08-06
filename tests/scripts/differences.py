from itertools import chain
import sys
from pysam import VariantFile
from datafusion import SessionContext, col, functions as f, lit

def index_merged_variants(filename, tag):
    with VariantFile(filename) as vcf:
        idx = {}
        for rec in vcf:
            locus = (rec.chrom, rec.pos, rec.id)
            vids = tuple(sorted(rec.info[tag]))
            idx[vids] = locus
        return idx

def index_original_variants(filename):
    with VariantFile(filename) as vcf:
        idx = {}
        for rec in vcf:
            vid = rec.id
            chrom = rec.chrom
            start = rec.start
            svtype = rec.info['SVTYPE']
            svlen = 0
            if "SVLEN" in rec.info:
                svlen = int(rec.info['SVLEN'])
            end = start
            if svtype != "INS" and svtype != "BND":
                end += abs(svlen)
            idx[vid] = (chrom, start, end, svtype, svlen)
        return idx

def get_rec(vid, orig):
    for idx in orig:
        if vid in idx:
            return idx[vid]
    print(f'vid not found: {vid}')
    return None

def get_variant_items(filename):
    with VariantFile(filename) as vcf:
        for rec in vcf:
            item = {"vid": rec.id,
                    "chrom": rec.chrom,
                   "start": rec.start,
                   "end": rec.start,
                   "svtype": rec.info['SVTYPE'],
                   "svlen": rec.info.get("SVLEN", None)}
            if item["svlen"] is not None and item["svtype"] != "INS" and item["svtype"] != "BND":
                item["end"] = item["start"] + abs(item["svlen"])
            yield item

def get_variant_sets(filename, wanted, tag):
    print(f'loading {filename}')
    with VariantFile(filename) as vcf:
        for rec in vcf:
            id = rec.id
            vids = rec.info.get(tag)
            found = False
            for vid in vids:
                if vid in wanted:
                    found = True
                    break
            if found:
                for vid in vids:
                    yield {"jas_id": id, "jas_vid": vid}

svelt = index_merged_variants(sys.argv[3], "ORIGINAL_IDS")
jasmine = index_merged_variants(sys.argv[4], "IDLIST")
orig = list(map(index_original_variants, sys.argv[5:]))

def write_differences_table(differences, full, name, tag):
    ctx = SessionContext()

    orig_df = ctx.from_pylist(list(chain(*map(get_variant_items, sys.argv[5:]))))
    #orig_df.show()

    jas_df = ctx.from_pylist(list(get_variant_sets(full, set().union(*differences), tag)))
    #jas_df.show()

    wanted_df = jas_df.join(orig_df, left_on="jas_vid", right_on="vid")
    #wanted_df.show()

    summary = wanted_df.aggregate([col("jas_id")], [
        f.min(col("svtype")).alias("svtype"),
        f.min(col("chrom")).alias("chrom"),
        f.min(col("start")).alias("start_min"),
        f.max(col("start")).alias("start_max"),
        f.min(col("end")).alias("end_min"),
        f.max(col("end")).alias("end_max"),
        f.min(f.abs(col("svlen"))).alias("len_min"),
        f.max(f.abs(col("svlen"))).alias("len_max"),
    ])
    summary = summary.with_column("start_delta", col("start_max") - col("start_min")) \
                    .with_column("end_delta", col("end_max") - col("end_min")) \
                    .with_column("length_delta", col("len_max") - col("len_min")) \
                    .with_column("length_ratio", col("len_min") * lit(1.0) / col("len_max"))
    summary.write_csv(name, True)

svelt_idx = {}
for vids in svelt.keys():
    for vid in vids:
        svelt_idx[vid] = vids

jasmine_idx = {}
for vids in jasmine.keys():
    for vid in vids:
        jasmine_idx[vid] = vids

more_svelt = set()
more_jasmine = set()
for vid, svelt_vids in svelt_idx.items():
    jasmine_vids = jasmine_idx[vid]
    if svelt_vids == jasmine_vids:
        continue
    if len(set(jasmine_vids) - set(svelt_vids)) > 0:
        more_jasmine.add(jasmine_vids)
    else:
        more_svelt.add(svelt_vids)

write_differences_table(more_svelt, sys.argv[3], sys.argv[1], "ORIGINAL_IDS")
write_differences_table(more_jasmine, sys.argv[4], sys.argv[2], "IDLIST")