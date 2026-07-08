#!/usr/bin/env python3
"""Legacy-compatible consensus FASTA direction correction.

This is the release-scoped replacement for the LFR_Pipeline
`denovo_supp.py --module fix_fa_rc` consensus path. It intentionally keeps the
legacy output semantics so cLFR_Release can reproduce LFR_Pipeline results:
PAF strand is keyed by the consensus FASTA ID prefix, the first seven PAF lines
are skipped, unmapped fragments are adapter-corrected and appended after mapped
fragments, and adapter reverse-complement correction is applied after mapping
orientation correction.
"""

import argparse
from pathlib import Path
from Bio import SeqIO


def reverse_complement(seq):
    return seq.translate(str.maketrans("ACGTNacgtn", "TGCANtgcan"))[::-1]


def adapter_correct(seq, adapter, flank_end):
    if not adapter or flank_end <= 0:
        return seq
    adapter_rc = reverse_complement(adapter)
    if adapter_rc in seq[-flank_end:]:
        return reverse_complement(seq)
    return seq


def read_fasta_by_prefix(input_fasta, umi_len):
    records = {}
    prefix_len = umi_len + 6
    with open(input_fasta) as handle:
        for record in SeqIO.parse(handle, "fasta"):
            records[record.id[:prefix_len]] = (record.id, str(record.seq))
    return records


def read_legacy_paf_strands(paf, umi_len):
    strands = {}
    prefix_len = umi_len + 6
    with open(paf) as handle:
        for _ in range(7):
            try:
                next(handle)
            except StopIteration:
                return strands
        for line in handle:
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 5:
                continue
            strands[fields[0][:prefix_len]] = fields[4]
    return strands


def write_record(handle, record_id, seq):
    handle.write(f">{record_id}\n{seq}\n")


def correct_fasta(input_fasta, paf, output_fasta, not_mapped_fasta, umi_len, adapter, flank_end):
    records = read_fasta_by_prefix(input_fasta, umi_len)
    strands = read_legacy_paf_strands(paf, umi_len)

    output_path = Path(output_fasta)
    missing_path = Path(not_mapped_fasta)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    missing_path.parent.mkdir(parents=True, exist_ok=True)

    mapped_fixed = []
    missing_fixed = []
    flipped_by_mapping = 0

    for prefix, (record_id, seq) in records.items():
        strand = strands.get(prefix)
        if strand is None:
            missing_fixed.append((record_id, adapter_correct(seq, adapter, flank_end)))
            continue
        if strand != "+":
            flipped_by_mapping += 1
            seq = reverse_complement(seq)
        mapped_fixed.append((record_id, adapter_correct(seq, adapter, flank_end)))

    with open(output_path, "w") as out:
        for record_id, seq in mapped_fixed:
            write_record(out, record_id, seq)
        for record_id, seq in missing_fixed:
            write_record(out, record_id, seq)

    with open(missing_path, "w") as missing:
        for record_id, seq in missing_fixed:
            write_record(missing, record_id, seq)

    print(f"records={len(records)}")
    print(f"flipped_by_mapping={flipped_by_mapping}")
    print(f"not_mapped={len(missing_fixed)}")


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
