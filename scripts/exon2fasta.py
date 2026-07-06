"""Basic FASTA metrics used by the cLFR consensus workflow."""

import argparse
import os
import tempfile
from collections import defaultdict
from statistics import mean, median

os.environ.setdefault("MPLCONFIGDIR", os.path.join(tempfile.gettempdir(), "matplotlib"))

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import seaborn as sns
from Bio import SeqIO


DEFAULT_BARCODE_LEN = 15
DEFAULT_MAX_PLOT_LEN = 7000


def collect_fasta_lengths(fasta, barcode_len=DEFAULT_BARCODE_LEN):
    """Return per-record lengths and unique barcode counts for a FASTA file."""
    lengths = []
    barcode_counts = defaultdict(int)

    for record in SeqIO.parse(fasta, "fasta"):
        lengths.append(len(record.seq))
        barcode = record.id[:barcode_len]
        barcode_counts[barcode] += 1

    return lengths, barcode_counts


def write_length_summary(fasta, outdir, lengths, barcode_counts):
    summary_path = os.path.join(outdir, "frag_length_distribution.txt")
    with open(summary_path, "a") as handle:
        handle.write(f"{fasta}\n")
        handle.write(f"frag count {len(lengths)}, bc count {len(barcode_counts)}\n")

        if lengths:
            handle.write(f"mean={round(mean(lengths), 4)}, median={median(lengths)}\n")
            handle.write(f"min={min(lengths)}, max={max(lengths)}\n")
        else:
            handle.write("mean=NA, median=NA\n")
            handle.write("min=NA, max=NA\n")


def plot_length_distribution(lengths, outdir, name, cutoff, max_plot_len=DEFAULT_MAX_PLOT_LEN):
    plot_path = os.path.join(outdir, f"{name}_frag_length_distribution.pdf")
    filtered_lengths = [length for length in lengths if length <= max_plot_len]

    sns.set_theme(style="whitegrid", rc={"figure.figsize": (15, 8)})
    plt.figure()

    if filtered_lengths:
        sns.histplot(filtered_lengths, bins=50, kde=True)
    else:
        plt.text(0.5, 0.5, "No FASTA records", ha="center", va="center")

    plt.title(f"{name}_frag_length_distribution_N{cutoff}")
    plt.xlabel("Fragment length")
    plt.ylabel("Count")
    plt.tight_layout()
    plt.savefig(plot_path)
    plt.close()


def metrics_basic(fasta, outdir, cutoff, name, barcode_len=DEFAULT_BARCODE_LEN):
    os.makedirs(outdir, exist_ok=True)
    lengths, barcode_counts = collect_fasta_lengths(fasta, barcode_len)
    write_length_summary(fasta, outdir, lengths, barcode_counts)
    plot_length_distribution(lengths, outdir, name, cutoff)


def parse_args():
    parser = argparse.ArgumentParser(
        description="Generate basic length metrics for consensus FASTA output."
    )
    parser.add_argument("--module", choices=["metrics_basic"], required=True)
    parser.add_argument("--fasta", required=True)
    parser.add_argument("--outdir", required=True)
    parser.add_argument("--minreads_fasta", type=int, required=True)
    parser.add_argument("--name", required=True)
    parser.add_argument("--barcode_len", type=int, default=DEFAULT_BARCODE_LEN)
    return parser.parse_args()


def main():
    args = parse_args()
    metrics_basic(
        fasta=args.fasta,
        outdir=args.outdir,
        cutoff=args.minreads_fasta,
        name=args.name,
        barcode_len=args.barcode_len,
    )


if __name__ == "__main__":
    main()
