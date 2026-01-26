#!/usr/bin/env python3
"""
merge_summary.py
Merge depth and Nextclade results using fuzzy matching.
"""

import pandas as pd
import argparse
import sys
import re

def extract_core_id(name):
    name = str(name).split()[0]                 # drop anything after a space
    name = name.split('/')[-1]                  # drop path prefix
    name = name.split('.')[0]                   # drop everything after first dot
    name = name.replace('-', '')                # remove hyphens
    core = re.sub(r'[^A-Za-z0-9]', '', name).upper()
    return core

def main():
    parser = argparse.ArgumentParser(
        description="Merge depth summary and Nextclade results (incl. coverage as %)."
    )
    parser.add_argument("--depth", required=True, help="Path to average_depth_summary.tsv")
    parser.add_argument("--nextclade", required=True, help="Path to nextclade.tsv")
    parser.add_argument("--output", required=True, help="Output file (TSV)")

    args = parser.parse_args()

    # ---------- 1. Depth ----------
    depth_df = pd.read_csv(args.depth, sep="\t")
    if depth_df.empty:
        sys.exit("Depth TSV is empty")

    sample_col = next((c for c in depth_df.columns if "sample" in c.lower()), None)
    if not sample_col:
        sys.exit("No column with 'sample' in name found in depth TSV")
    depth_df = depth_df.rename(columns={sample_col: "sample_id"})

    # Standardise depth columns
    rename_map = {}
    for col in depth_df.columns:
        low = col.lower()
        if "mean" in low and "depth" in low:
            rename_map[col] = "mean_depth"
        elif "median" in low and "depth" in low:
            rename_map[col] = "median_depth"
    depth_df = depth_df.rename(columns=rename_map)

    depth_df["_core_id"] = depth_df["sample_id"].apply(extract_core_id)
    print(f"Loaded {len(depth_df)} samples from depth file")
    print("Depth core IDs (first 3):", depth_df["_core_id"].head(3).tolist())

    # ---------- 2. Nextclade ----------
    nc_df = pd.read_csv(args.nextclade, sep="\t")
    if nc_df.empty:
        sys.exit("Nextclade TSV is empty")
    if "seqName" not in nc_df.columns:
        sys.exit("No 'seqName' column in Nextclade TSV")

    nc_df["_core_id"] = nc_df["seqName"].apply(extract_core_id)

    # Pick columns: clade, lineage, and coverage
    clade_col   = next((c for c in nc_df.columns if "clade" in c.lower()), None)
    lineage_col = next((c for c in nc_df.columns if any(x in c.lower() for x in ["lineage", "pango"])), None)
    cov_col     = next((c for c in nc_df.columns if "coverage" in c.lower()), None)

    nc_cols = ["_core_id"]
    rename_nc = {}

    if clade_col:
        nc_cols.append(clade_col)
        rename_nc[clade_col] = "clade"
    if lineage_col:
        nc_cols.append(lineage_col)
        rename_nc[lineage_col] = "lineage"
    if cov_col:
        nc_cols.append(cov_col)
        rename_nc[cov_col] = "coverage(%)"   # final column name with unit

    nc_subset = nc_df[nc_cols].copy()
    nc_subset = nc_subset.rename(columns=rename_nc)

    # Convert coverage to percentage (as float, rounded to 1 decimal)
    if "coverage(%)" in nc_subset.columns:
        nc_subset["coverage(%)"] = (
            pd.to_numeric(nc_subset["coverage(%)"], errors="coerce") * 100
        ).round(1)

    print(f"Loaded {len(nc_subset)} samples from nextclade file")
    print("Nextclade core IDs (first 3):", nc_df["_core_id"].head(3).tolist())

    # ---------- 3. Merge ----------
    merged = depth_df.merge(nc_subset, on="_core_id", how="left")
    merged = merged.drop(columns="_core_id")

    # Fill NAs (except sample_id)
    for col in merged.columns:
        if col != "sample_id":
            merged[col] = merged[col].fillna("N/A").replace("", "N/A")

    # Ensure all expected columns exist
    expected = ["mean_depth", "median_depth", "coverage(%)", "clade", "lineage"]
    for col in expected:
        if col not in merged.columns:
            merged[col] = "N/A"

    final_cols = ["sample_id"] + expected
    merged = merged[[c for c in final_cols if c in merged.columns]]

    # ---------- 4. Write ----------
    merged.to_csv(args.output, sep="\t", index=False)

    # ---------- Diagnostics ----------
    print("\n" + "="*60)
    print(f"Output written to: {args.output}")
    print("="*60)
    print(f"Total rows: {len(merged)}")

    matched = merged[(merged["clade"] != "N/A") | (merged["lineage"] != "N/A")].shape[0]
    print(f"Matched (have Nextclade data): {matched}")
    print(f"Unmatched: {len(merged) - matched}")

    if matched < len(merged):
        print("\nSamples without Nextclade data:")
        for s in merged[(merged["clade"] == "N/A") & (merged["lineage"] == "N/A")]["sample_id"].head(10):
            print(f" - {s}")
        extra = (len(merged) - matched) - 10
        if extra > 0:
            print(f" ... and {extra} more")

    print("="*60 + "\n")

if __name__ == "__main__":
    main()