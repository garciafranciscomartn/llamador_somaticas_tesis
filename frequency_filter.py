#!/usr/bin/env python3
"""
frequency_filter.py
===================
Filter a merged VCF based on the PAIR_FREQ INFO field.

Keep only variants with PAIR_FREQ >= min_frequency (default 0.75).
True somatic variants should be called in most/all pairwise comparisons,
so high-frequency variants are the real somatic calls.
Variants found in few pairs are likely noise or artifacts.

Usage:
  python frequency_filter.py --input merged.vcf.gz --output filtered.vcf.gz \
      --min-frequency 0.75
"""

import argparse
import gzip
import re
import sys


def parse_args():
    parser = argparse.ArgumentParser(description="Filter VCF by PAIR_FREQ")
    parser.add_argument("--input", required=True, help="Input merged VCF (bgzipped)")
    parser.add_argument("--output", required=True, help="Output filtered VCF (bgzipped)")
    parser.add_argument("--min-frequency", type=float, default=0.75,
                        help="Minimum PAIR_FREQ to keep (inclusive). Default: 0.75")
    return parser.parse_args()


def get_pair_freq(info_field):
    """Extract PAIR_FREQ value from INFO field."""
    match = re.search(r"PAIR_FREQ=([0-9.]+)", info_field)
    if match:
        return float(match.group(1))
    return None


def main():
    args = parse_args()

    n_total = 0
    n_kept = 0
    n_removed = 0

    with gzip.open(args.input, "rt") as fh, gzip.open(args.output, "wt") as out:
        for line in fh:
            if line.startswith("#"):
                out.write(line)
                continue

            n_total += 1
            fields = line.split("\t")
            if len(fields) < 8:
                n_total -= 1
                continue

            freq = get_pair_freq(fields[7])
            if freq is not None and freq >= args.min_frequency:
                out.write(line)
                n_kept += 1
            else:
                n_removed += 1

    print(f"Frequency filter (min_freq >= {args.min_frequency}):")
    print(f"  Total variants:   {n_total}")
    print(f"  Kept (consistent): {n_kept}")
    print(f"  Removed (rare):    {n_removed}")
    print(f"  Output: {args.output}")


if __name__ == "__main__":
    main()
