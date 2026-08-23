#!/usr/bin/env python3
"""
merge_pair_vcfs.py
==================
Merge multiple pairwise PASS VCFs for a single tumor sample.

For each variant (CHROM, POS, REF, ALT), count how many pairwise
comparisons called it. Add INFO fields:
  PAIR_COUNT  = number of pairs that called the variant
  PAIR_FREQ   = PAIR_COUNT / total_pairs
  PAIR_NAMES  = comma-separated list of pair names that called it

The output retains the VCF record from the first pair that found
the variant (representative genotype / INFO).

Usage:
  python merge_pair_vcfs.py --vcfs pair1.vcf.gz pair2.vcf.gz ... \
      --output merged.vcf.gz --tumor-id SAMPLE --n-pairs N
"""

import argparse
import gzip
import os
import sys
from collections import OrderedDict


def parse_args():
    parser = argparse.ArgumentParser(description="Merge pairwise PASS VCFs")
    parser.add_argument("--vcfs", nargs="+", required=True, help="Input PASS VCF files")
    parser.add_argument("--output", required=True, help="Output VCF (bgzipped)")
    parser.add_argument("--tumor-id", required=True, help="Tumor sample ID")
    parser.add_argument("--n-pairs", type=int, required=True, help="Total number of pairs")
    return parser.parse_args()


def open_vcf(path):
    """Open a VCF, handling .gz transparently."""
    if path.endswith(".gz"):
        return gzip.open(path, "rt")
    return open(path, "r")


def variant_key(chrom, pos, ref, alt):
    return (chrom, pos, ref, alt)


def main():
    args = parse_args()

    # Collect all variants across pairs: key -> {count, freq, pairs, record}
    variants = OrderedDict()
    header_lines = []
    header_collected = False

    for vcf_path in args.vcfs:
        # Extract pair name from filename: path/PAIR_NAME/PAIR_NAME.mutect2.pass.vcf.gz
        pair_name = os.path.basename(os.path.dirname(vcf_path))
        if not pair_name:
            pair_name = os.path.basename(vcf_path).replace(".mutect2.pass.vcf.gz", "")

        with open_vcf(vcf_path) as fh:
            for line in fh:
                line = line.rstrip("\n")
                if line.startswith("#"):
                    if not header_collected:
                        header_lines.append(line)
                    continue

                fields = line.split("\t")
                if len(fields) < 8:
                    continue

                chrom, pos, _, ref, alt = fields[0], fields[1], fields[2], fields[3], fields[4]
                key = variant_key(chrom, pos, ref, alt)

                if key not in variants:
                    variants[key] = {
                        "count": 0,
                        "pairs": [],
                        "record": line,
                    }
                variants[key]["count"] += 1
                variants[key]["pairs"].append(pair_name)

        header_collected = True

    # Inject custom INFO header lines before #CHROM
    custom_headers = [
        '##INFO=<ID=PAIR_COUNT,Number=1,Type=Integer,Description="Number of pairwise comparisons that called this variant">',
        '##INFO=<ID=PAIR_FREQ,Number=1,Type=Float,Description="Fraction of pairs calling this variant (PAIR_COUNT/total_pairs)">',
        '##INFO=<ID=PAIR_NAMES,Number=.,Type=String,Description="Pair names that called this variant">',
        f'##merge_pair_vcfs_tumor={args.tumor_id}',
        f'##merge_pair_vcfs_n_pairs={args.n_pairs}',
    ]

    # Write output
    with gzip.open(args.output, "wt") as out:
        for hline in header_lines:
            if hline.startswith("#CHROM"):
                for ch in custom_headers:
                    out.write(ch + "\n")
            out.write(hline + "\n")

        # Write variant records with PAIR_COUNT/PAIR_FREQ/PAIR_NAMES added to INFO
        for key, data in variants.items():
            freq = data["count"] / args.n_pairs
            pair_names = ",".join(data["pairs"])

            fields = data["record"].split("\t")
            info = fields[7]
            extra_info = f"PAIR_COUNT={data['count']};PAIR_FREQ={freq:.4f};PAIR_NAMES={pair_names}"
            if info == "." or info == "":
                fields[7] = extra_info
            else:
                fields[7] = f"{info};{extra_info}"

            out.write("\t".join(fields) + "\n")

    n_total = len(variants)
    n_multi = sum(1 for v in variants.values() if v["count"] > 1)
    print(f"Merged {len(args.vcfs)} VCFs for tumor {args.tumor_id}")
    print(f"Total unique variants: {n_total}")
    print(f"Variants called by >1 pair: {n_multi}")
    print(f"Output: {args.output}")


if __name__ == "__main__":
    main()
