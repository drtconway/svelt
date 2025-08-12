import hashlib
import random
import sys
import pysam

def sig(seed, x):
    seed = f'{seed}'.encode('utf-8')
    bytes = f'{x}'.encode('utf-8')
    hash = hashlib.sha256()
    hash.update(seed)
    hash.update(bytes)
    hash = hash.digest()
    res = 0.0
    for b in hash:
        res = (res + b) / 256.0
    return res

def select_variants(seed, p, filename):
    with pysam.VariantFile(filename) as vcf:
        for rec in vcf:
            id = rec.id
            if id is None:
                id = f'{rec.chrom}-{rec.pos}-{rec.ref}-{rec.alts}-{rec.info["SVTYPE"]}'
            q = sig(seed, id)
            if q < p:
                yield rec

base_chroms = dict([(f'chr{i}', i) for i in range(1,23)])

def variant_key(v):
    chrom = v.chrom
    chrom_num = 100
    if chrom in base_chroms:
        chrom_num = base_chroms[chrom]
    return (chrom_num, chrom, v.pos, v.ref, tuple(v.alts))

def all_the_variant_ids(variants):
    for md in variants.values():
        for (_, vs) in md.values():
            for ids in vs:
                for id in ids:
                    yield id

def select_n_variants(seed, n, filename, tag):
    random.seed(seed)
    variants = {}
    with pysam.VariantFile(filename) as vcf:
        for rec in vcf:
            svtype = rec.info["SVTYPE"]
            ids = tuple(sorted(rec.info[tag]))
            m = len(ids)
            if svtype not in variants:
                variants[svtype] = {}
            if m not in variants[svtype]:
                variants[svtype][m] = (0, [])
            (j, vs) = variants[svtype][m]
            if len(vs) < n:
                vs.append(ids)
            else:
                u = random.randint(0, j)
                if u < n:
                    vs[u] = ids
            variants[svtype][m] = (j + 1, vs)
    for rec in all_the_variant_ids(variants):
        yield rec

ids = set(select_n_variants(19, 10, sys.argv[1], "IDLIST"))

out_dir = sys.argv[2]

for orig_filename in sys.argv[3:]:
    with pysam.VariantFile(orig_filename) as orig_vcf:
        with pysam.VariantFile(f'{out_dir}/sampled-{orig_filename}', 'w', header=orig_vcf.header) as sampled_vcf:
            for rec in orig_vcf:
                if rec.id in ids:
                    sampled_vcf.write(rec)