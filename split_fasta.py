#!/usr/bin/env python3

import os
import argparse
from Bio import SeqIO

def split_fasta_per_protein(input_dir, output_dir, min_length=80):
    input_dir = os.path.abspath(input_dir)
    output_dir = os.path.abspath(output_dir)
    os.makedirs(output_dir, exist_ok=True)

    refset_dirs = [d for d in os.listdir(input_dir) if d.startswith("RefSet_")]

    for refset in refset_dirs:
        refset_path = os.path.join(input_dir, refset)
        fasta_files = [f for f in os.listdir(refset_path) if f.endswith(('.fa', '.faa', '.fasta'))]

        if not fasta_files:
            print(f"No FASTA files in {refset_path}. Skipping.")
            continue

        out_ref_dir = os.path.join(output_dir, f"seq{refset}")
        os.makedirs(out_ref_dir, exist_ok=True)

        seen_ids = set()
        total, written, short, duplicate = 0, 0, 0, 0

        for fasta in fasta_files:
            fasta_path = os.path.join(refset_path, fasta)
            for record in SeqIO.parse(fasta_path, "fasta"):
                total += 1
                seq = str(record.seq).rstrip("*")

                if len(seq) < min_length:
                    short += 1
                    continue

                if "*" in seq:
                    print(f"Warning: Internal stop codon in {record.id}")

                if record.id in seen_ids:
                    print(f"Duplicate ID found: {record.id}. Skipping.")
                    duplicate += 1
                    continue
                seen_ids.add(record.id)

                output_file = os.path.join(out_ref_dir, f"{record.id}.fa")
                with open(output_file, "w") as out_f:
                    out_f.write(f">{record.id}\n{seq}\n")
                written += 1

        print(f"\n Processed {refset}:")
        print(f"   Total records: {total}")
        print(f"   Written:       {written}")
        print(f"   Skipped (short): {short}")
        print(f"   Skipped (duplicate): {duplicate}")
        print("-" * 40)

def main():
    parser = argparse.ArgumentParser(description="Split protein FASTA into individual files per protein.")
    parser.add_argument("--input_dir", required=True, help="Directory containing RefSet_*/Nr_NIs_PSeqs.faa")
    parser.add_argument("--output_dir", required=True, help="Output directory to store seqRefSet_* folders")
    parser.add_argument("--min_length", type=int, default=80, help="Minimum sequence length (default: 80)")
    args = parser.parse_args()

    split_fasta_per_protein(args.input_dir, args.output_dir, args.min_length)

if __name__ == "__main__":
    main()
