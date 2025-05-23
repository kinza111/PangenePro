#!/usr/bin/env python3

import os
import argparse
import xml.etree.ElementTree as ET
from Bio import SeqIO

IPR_NS = "https://ftp.ebi.ac.uk/pub/software/unix/iprscan/5/schemas"

def extract_domains_from_query_xml(query_xml):
    tree = ET.parse(query_xml)
    root = tree.getroot()
    domains = set()

    for entry in root.findall(f".//{{{IPR_NS}}}entry"):
        if entry.attrib.get('ac'):
            domains.add(entry.attrib['ac'])

    print(f"Extracted {len(domains)} InterPro domain(s) from query XML.")
    return domains

def filter_proteins(domain_file, target_domains):
    filtered = set()

    with open(domain_file) as f:
        next(f)  # skip header
        for line in f:
            parts = line.strip().split("\t")
            if len(parts) != 3:
                continue
            prot_id, domain, desc = parts
            if domain in target_domains:
                filtered.add(prot_id)

    return filtered

def copy_validated_sequences(protein_ids, seq_folder, out_fasta):
    seq_files = [f for f in os.listdir(seq_folder) if f.endswith(".fa")]
    written = 0

    with open(out_fasta, "w") as out:
        for fname in seq_files:
            fa_path = os.path.join(seq_folder, fname)
            for record in SeqIO.parse(fa_path, "fasta"):
                if record.id in protein_ids:
                    SeqIO.write(record, out, "fasta")
                    written += 1
    print(f"Wrote {written} validated protein sequences to {out_fasta}.")

def merge_validated_fastas(validated_dir, output_fasta):
    os.makedirs(os.path.dirname(output_fasta), exist_ok=True)
    merged_records = []
    count = 0

    for fname in os.listdir(validated_dir):
        if fname.startswith("validated_RefSet_") and fname.endswith(".faa"):
            refset = fname.replace("validated_", "").replace(".faa", "")
            fasta_path = os.path.join(validated_dir, fname)

            for record in SeqIO.parse(fasta_path, "fasta"):
                original_id = record.id
                record.id = f"{refset}|{original_id}"
                record.description = ""
                merged_records.append(record)
                count += 1

    with open(output_fasta, "w") as out_f:
        SeqIO.write(merged_records, out_f, "fasta")

    print("\n" + "-"*42)
    print("MODULE 2 - PANGENES ANALYSIS")
    print(f"Merged {count} validated proteins into:")
    print(f"{output_fasta}")
    print("Ready for DIAMOND all-vs-all ortholog clustering.")
    print("-"*42 + "\n")

def main():
    parser = argparse.ArgumentParser(description="Filter proteins with InterPro domains matching query and prepare for pangene clustering.")
    parser.add_argument("--query_xml", required=True, help="InterProScan XML result for the query.fa file.")
    parser.add_argument("--input_dir", required=True, help="Directory containing domains_RefSet_X.tsv files.")
    parser.add_argument("--seqs_dir", required=True, help="Directory containing seqRefSet_X folders with .fa files.")
    parser.add_argument("--output_dir", required=True, help="Output directory to store validated .tsv and .faa files.")
    args = parser.parse_args()

    os.makedirs(args.output_dir, exist_ok=True)
    target_domains = extract_domains_from_query_xml(args.query_xml)

    domain_files = [f for f in os.listdir(args.input_dir) if f.startswith("domains_RefSet_") and f.endswith(".tsv")]

    for domain_file in domain_files:
        refset_num = domain_file.split("_")[-1].replace(".tsv", "")
        domain_path = os.path.join(args.input_dir, domain_file)
        filtered_ids = filter_proteins(domain_path, target_domains)

        out_tsv = os.path.join(args.output_dir, f"validated_RefSet_{refset_num}.tsv")
        with open(out_tsv, "w") as out:
            out.write("Protein ID\n")
            for pid in sorted(filtered_ids):
                out.write(f"{pid}\n")
        print(f"Validated {len(filtered_ids)} proteins in RefSet_{refset_num}")

        seq_folder = os.path.join(args.seqs_dir, f"seqRefSet_{refset_num}")
        out_faa = os.path.join(args.output_dir, f"validated_RefSet_{refset_num}.faa")
        copy_validated_sequences(filtered_ids, seq_folder, out_faa)

    # Final step: merge all .faa files into one for clustering
    merged_fasta = os.path.join("ortholog_clustering", "combined_validated.faa")
    merge_validated_fastas(args.output_dir, merged_fasta)

if __name__ == "__main__":
    main()
