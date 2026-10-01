#!/usr/bin/env python3
"""Compare native vs GATK final genotyped VCFs from pipeline_v4."""
import gzip, sys, re

def load(path):
    """Return dict[(chrom,pos,ref,alt)] -> (filter, qual, gt, dp, gq)."""
    out = {}
    with gzip.open(path, 'rt') as f:
        for line in f:
            if line.startswith('#'): continue
            c = line.rstrip('\n').split('\t')
            chrom, pos, _id, ref, alt, qual, filt = c[0], int(c[1]), c[2], c[3], c[4], c[5], c[6]
            gt = dp = gq = None
            if len(c) > 9:
                keys = c[8].split(':')
                vals = c[9].split(':')
                d = dict(zip(keys, vals))
                gt, dp, gq = d.get('GT'), d.get('DP'), d.get('GQ')
            out[(chrom, pos, ref, alt)] = (filt, qual, gt, dp, gq)
    return out

native = load('public_pipeline/pipeline_v4/native/gt/out.vcf.gz')
gatk = load('public_pipeline/pipeline_v4/gatk/full.vcf.gz')

nk, gk = set(native), set(gatk)
common = nk & gk
only_n = nk - gk
only_g = gk - nk

print(f"native records: {len(nk)}")
print(f"gatk records:   {len(gk)}")
print(f"common sites (chrom:pos:ref:alt): {len(common)}")
print(f"native-only: {len(only_n)}   gatk-only: {len(only_g)}")
print(f"site concordance vs union: {len(common)/len(nk|gk)*100:.2f}%")

# genotype concordance on common sites
gt_same = sum(1 for k in common if native[k][2] == gatk[k][2])
print(f"GT identical on common: {gt_same}/{len(common)} ({gt_same/len(common)*100:.2f}%)")

# QUAL correlation-ish summary on common
import statistics
qdiffs = []
for k in common:
    try:
        qn = float(native[k][1]); qg = float(gatk[k][1])
        qdiffs.append(qn - qg)
    except (TypeError, ValueError):
        pass
if qdiffs:
    print(f"QUAL diff (native-gatk): mean={statistics.mean(qdiffs):+.3f} "
          f"median={statistics.median(qdiffs):+.3f} "
          f"stdev={statistics.stdev(qdiffs):.3f} n={len(qdiffs)}")

def show(title, keys, src):
    print(f"\n{title} (up to 10):")
    for k in sorted(keys, key=lambda x: x[1])[:10]:
        filt, qual, gt, dp, gq = src[k]
        print(f"  {k[0]}:{k[1]} {k[2]}>{k[3]} qual={qual} filter={filt} GT={gt} DP={dp} GQ={gq}")

show("native-only", only_n, native)
show("gatk-only", only_g, gatk)

# GT discordant examples
disc = [k for k in common if native[k][2] != gatk[k][2]]
print(f"\nGT discordant (up to 10 of {len(disc)}):")
for k in sorted(disc, key=lambda x: x[1])[:10]:
    print(f"  {k[0]}:{k[1]} {k[2]}>{k[3]} native GT={native[k][2]} GQ={native[k][4]} | gatk GT={gatk[k][2]} GQ={gatk[k][4]}")
