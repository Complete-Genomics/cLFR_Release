#!/usr/bin/env python3
"""Correct consensus FASTA direction using minimap2 PAF strand calls.

This is the release-scoped replacement for the legacy denovo_supp.py
`fix_fa_rc` path. It keeps only the consensus FASTA direction-correction
logic needed by cLFR_Release.
"""

import argparse
from pathlib import Path
from Bio import SeqIO


def reverse_complement(seq):
    return seq.translate(str.maketrans("ACGTNacgtn", "TGCANtgcan"))[::-1]


def read_paf_strands(paf, umi_len):
    exact = {}
    prefix = {}
    prefix_len = umi_len + 6
    with open(paf) as handle:
        for line in handle:
            if not line.strip() or line.startswith("#") or line.startswith("@"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 5:
                continue
            query_name = fields[0]
            strand = fields[4]
            exact.setdefault(query_name, strand)
            prefix.setdefault(query_name[:prefix_len], strand)
    return exact, prefix


def adapter_correct(seq, adapter, flank_end):
    if not adapter or flank_end <= 0:
        return seq
    adapter_rc = reverse_complement(adapter)
    if adapter_rc in seq[-flank_end:]:
        return reverse_complement(seq)
    return seq


def choose_strand(record_id, exact_strands, prefix_strands, umi_len):
    if record_id in exact_strands:
        return exact_strands[record_id]
    return prefix_strands.get(record_id[: umi_len + 6])


def correct_fasta(input_fasta, paf, output_fasta, not_mapped_fasta, umi_len, adapter, flank_end):
    exact_strands, prefix_strands = read_paf_strands(paf, umi_len)
    total = 0
    flipped_by_mapping = 0
    not_mapped = 0

    Path(output_fasta).parent.mkdir(parents=True, exist_ok=True)
    Path(not_mapped_fasta).parent.mkdir(parents=True, exist_ok=True)

    with open(output_fasta, "w") as out, open(not_mapped_fasta, "w") as missing:
        for record in SeqIO.parse(input_fasta, "fasta"):
            total += 1
            seq = str(record.seq)
            strand = choose_strand(record.id, exact_strands, prefix_strands, umi_len)
            if strand is None:
                not_mapped += 1
                fixed = adapter_correct(seq, adapter, flank_end)
                missing.write(f">{record.id}\n{fixed}\n")
                out.write(f">{record.id}\n{fixed}\n")
                continue
            if strand != "+":
                flipped_by_mapping += 1
                seq = reverse_complement(seq)
            fixed = adapter_correct(seq, adapter, flank_end)
            out.write(f">{record.id}\n{fixed}\n")

    print(f"records={total}")
    print(f"flipped_by_mapping={flipped_by_mapping}")
    print(f"not_mapped={not_mapped}")


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("--input_fasta", required=True)
    parser.add_argument("--paf", required=True)
    parser.add_argument("--output_fasta", required=True)
    parser.add_argument("--not_mapped_fasta", required=True)
    parser.add_argument("--umi_len", type=int, default=15)
    parser.add_argument("--adapter_seq", default="")
    parser.add_argument("--flank_end", type=int, default=0)
    return parser.parse_args()


def main():
    args = parse_args()
    correct_fasta(
        input_fasta=args.input_fasta,
        paf=args.paf,
        output_fasta=args.output_fasta,
        not_mapped_fasta=args.not_mapped_fasta,
        umi_len=args.umi_len,
        adapter=args.adapter_seq,
        flank_end=args.flank_end,
    )


if __name__ == "__main__":
    main()
