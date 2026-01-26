#!/usr/bin/env python3
"""
Simplified BAM Depth Plotter

Features:
1. Generates individual depth plots for each BAM file.
2. Generates a single average depth plot across all samples.
3. Outputs a TSV file of average and median depth per sample.
"""

import argparse
import sys
from pathlib import Path
import pysam
import matplotlib.pyplot as plt
import numpy as np
import csv


def find_bam_files(input_path):
    """Find all BAM files in a directory or single BAM file."""
    input_path = Path(input_path)
    if input_path.is_file() and input_path.suffix.lower() == ".bam":
        return [str(input_path)]
    elif input_path.is_dir():
        bam_files = sorted(input_path.glob("*.bam"))
        if not bam_files:
            sys.exit(f"❌ No BAM files found in directory: {input_path}")
        return [str(f) for f in bam_files]
    else:
        sys.exit(f"❌ Invalid path: {input_path}")


def calculate_depth(bam_file):
    """Return positions and depths for a BAM file."""
    try:
        bam = pysam.AlignmentFile(bam_file, "rb")
        reference = bam.references[0]
        ref_length = bam.get_reference_length(reference)

        positions, depths = [], []
        for pileupcolumn in bam.pileup(reference, 0, ref_length):
            positions.append(pileupcolumn.pos)
            depths.append(pileupcolumn.n)

        bam.close()
        return np.array(positions), np.array(depths)
    except Exception as e:
        print(f"Error processing {bam_file}: {e}")
        return np.array([]), np.array([])


def plot_individual_depth(positions, depths, bam_file, output_dir):
    """Plot depth for a single BAM file."""
    sample = Path(bam_file).stem
    plt.figure(figsize=(12, 5))
    plt.plot(positions, depths, color="blue", linewidth=0.8)
    plt.fill_between(positions, depths, color="lightblue", alpha=0.4)
    plt.title(f"Read Depth: {sample}")
    plt.xlabel("Position")
    plt.ylabel("Depth")

    if len(depths) > 0:
        mean_depth = np.mean(depths)
        median_depth = np.median(depths)
        plt.axhline(mean_depth, color="red", linestyle="--", label=f"Mean: {mean_depth:.1f}")
        plt.axhline(median_depth, color="green", linestyle=":", label=f"Median: {median_depth:.1f}")
        plt.legend()

    plt.tight_layout()
    output_file = Path(output_dir) / f"{sample}_depth.png"
    plt.savefig(output_file, dpi=300, bbox_inches="tight")
    plt.close()
    print(f"✅ Saved individual plot: {output_file}")


def plot_average_depth(depth_data, output_dir):
    """Plot average depth curve across all samples."""
    if not depth_data:
        print("⚠️ No data available for average plot.")
        return

    all_positions = [pos for pos, _, _ in depth_data if len(pos) > 0]
    if not all_positions:
        print("⚠️ No valid position data found.")
        return

    min_pos = min(pos.min() for pos in all_positions)
    max_pos = max(pos.max() for pos in all_positions)
    position_grid = np.arange(min_pos, max_pos + 1)

    depth_matrix = np.zeros((len(depth_data), len(position_grid)))
    for i, (positions, depths, _) in enumerate(depth_data):
        if len(positions) > 0:
            depth_matrix[i] = np.interp(position_grid, positions, depths, left=0, right=0)

    avg_depth = np.mean(depth_matrix, axis=0)
    std_depth = np.std(depth_matrix, axis=0)

    plt.figure(figsize=(12, 5))
    plt.plot(position_grid, avg_depth, color="navy", linewidth=1.5, label="Average Depth")
    plt.fill_between(position_grid, avg_depth - std_depth, avg_depth + std_depth,
                     color="lightblue", alpha=0.3, label="±1 SD")
    plt.title("Average Read Depth Across Samples")
    plt.xlabel("Position")
    plt.ylabel("Depth")
    plt.legend()
    plt.tight_layout()

    avg_plot = Path(output_dir) / "average_depth.png"
    plt.savefig(avg_plot, dpi=300, bbox_inches="tight")
    plt.close()
    print(f"✅ Saved average depth plot: {avg_plot}")


def write_depth_summary(depth_data, output_dir):
    """Write average and median depth per sample to TSV."""
    tsv_path = Path(output_dir) / "average_depth_summary.tsv"
    with open(tsv_path, "w", newline="") as tsvfile:
        writer = csv.writer(tsvfile, delimiter="\t")
        writer.writerow(["Sample", "MeanDepth", "MedianDepth"])
        for _, depths, filename in depth_data:
            mean_depth = np.mean(depths) if len(depths) > 0 else 0
            median_depth = np.median(depths) if len(depths) > 0 else 0
            writer.writerow([filename, f"{mean_depth:.2f}", f"{median_depth:.2f}"])
    print(f"✅ Saved depth summary TSV: {tsv_path}")


def main():
    parser = argparse.ArgumentParser(
        description="Plot BAM read depth per file, average plot, and summary TSV."
    )
    parser.add_argument("input_path", help="Path to BAM file or directory of BAMs")
    parser.add_argument("-o", "--output-dir", required=True,
                        help="Output directory for plots and summary TSV")
    args = parser.parse_args()

    # Setup
    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    # Find BAM files
    bam_files = find_bam_files(args.input_path)
    print(f"📂 Found {len(bam_files)} BAM file(s).")

    depth_data = []
    for bam_file in bam_files:
        positions, depths = calculate_depth(bam_file)
        sample = Path(bam_file).stem
        depth_data.append((positions, depths, sample))
        plot_individual_depth(positions, depths, bam_file, output_dir)

    # Average plot + summary
    plot_average_depth(depth_data, output_dir)
    write_depth_summary(depth_data, output_dir)


if __name__ == "__main__":
    main()
