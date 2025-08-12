from itertools import chain
import sys
from pysam import VariantFile
from datafusion import SessionContext, col, functions as f, lit

def index_merged_variants(filename, tag):
    with VariantFile(filename) as vcf:
        idx = {}
        for rec in vcf:
            id = rec.id
            chrom = rec.chrom
            start = rec.start
            end = rec.stop
            length = 0
            svtype = rec.info["SVTYPE"]
            if svtype == "BND":
                continue
            elif svtype == "TRA":
                continue
            elif svtype == "INS":
                length = int(rec.info["SVLEN"])
            else:
                length = int(rec.info["SVLEN"])
                end = start + abs(length)
            locus = (chrom, start, end, length, id)
            vids = tuple(sorted(rec.info[tag]))
            if len(vids) == 1:
                continue
            for vid in vids:
                idx[vid] = locus
        return idx

def index_original_variants(filename):
    with VariantFile(filename) as vcf:
        idx = {}
        for rec in vcf:
            vid = rec.id
            chrom = rec.chrom
            start = rec.start
            end = rec.stop
            length = 0
            svtype = rec.info["SVTYPE"]
            if svtype == "BND":
                continue
            elif svtype == "TRA":
                continue
            elif svtype == "INS":
                length = int(rec.info["SVLEN"])
            else:
                length = int(rec.info["SVLEN"])
                end = start + abs(length)
            locus = (chrom, start, end, length)
            yield (vid, svtype, locus)

svelt = index_merged_variants(sys.argv[2], "ORIGINAL_IDS")
jasmine = index_merged_variants(sys.argv[3], "IDLIST")

rows = []
for orig_fn in sys.argv[4:]:
    for vid, svtype, locus in index_original_variants(orig_fn):
        if vid in svelt:
            (chrom, start, end, length) = locus
            merged_locus = svelt[vid]
            (mchrom, mstart, mend, mlength, mid) = merged_locus
            rows.append({"tool": "svelt", "merged_id": mid, "original_id": vid, "svtype": svtype,
                         "start_offset": (mstart - start), "end_offset": (mend - end), "length_offset": (mlength - length)})
        if vid in jasmine:
            (chrom, start, end, length) = locus
            merged_locus = jasmine[vid]
            (mchrom, mstart, mend, mlength, mid) = merged_locus
            rows.append({"tool": "jasmine", "merged_id": mid, "original_id": vid, "svtype": svtype,
                         "start_offset": (mstart - start), "end_offset": (mend - end), "length_offset": (mlength - length)})
            
ctx = SessionContext()

dt = ctx.from_pylist(rows)
summary = dt.aggregate([col("tool"), col("svtype")], [
    f.mean(f.abs(col("start_offset"))).alias("start_disp_mean"),
    f.stddev(f.abs(col("start_offset"))).alias("start_disp_sd"),
    f.mean(f.abs(col("end_offset"))).alias("end_disp_mean"),
    f.stddev(f.abs(col("end_offset"))).alias("end_disp_sd"),
    f.mean(f.abs(col("length_offset"))).alias("length_disp_mean"),
    f.stddev(f.abs(col("length_offset"))).alias("length_disp_sd"),
]).sort(col("svtype"), col("tool"))
summary.write_csv(sys.argv[1], True)