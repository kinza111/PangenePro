#!/usr/bin/env python3

import os
import argparse
import xml.etree.ElementTree as ET

# InterProScan XML uses this namespace
IPR_NS = "https://ftp.ebi.ac.uk/pub/software/unix/iprscan/5/schemas"

def parse_xml_file(xml_file):
    try:
        tree = ET.parse(xml_file)
        root = tree.getroot()

        protein_id = root.find(f".//{{{IPR_NS}}}xref").attrib['id']
        domain_entries = []

        for entry in root.findall(f".//{{{IPR_NS}}}entry"):
            if entry.attrib.get('type') == 'DOMAIN':
                acc = entry.attrib.get('ac')
                desc = entry.attrib.get('desc')
                domain_entries.append((protein_id, acc, desc))

        return domain_entries

    except Exception as e:
        print(f"Error parsing {xml_file}: {e}")
        return []

def process_seqref_folder(folder_path, output_file):
    xml_files = [f for f in os.listdir(folder_path) if f.endswith('.xml')]
    total_proteins, total_domains = 0, 0

    with open(output_file, 'w') as out:
        out.write("Protein ID\tInterPro Domain\tDescription\n")
        for xml in xml_files:
            xml_path = os.path.join(folder_path, xml)
            entries = parse_xml_file(xml_path)
            if entries:
                total_proteins += 1
                total_domains += len(entries)
                for entry in entries:
                    out.write("\t".join(entry) + "\n")

    print(f"Processed {folder_path}: {total_proteins} proteins, {total_domains} domains.")

def main():
    parser = argparse.ArgumentParser(description="Extract InterPro domain annotations from .xml files.")
    parser.add_argument("--input_dir", required=True, help="Directory containing seqRefSet_*/ folders with .xml files")
    parser.add_argument("--output_dir", default="domains", help="Directory to save domain output files (default: domains)")
    args = parser.parse_args()

    input_dir = os.path.abspath(args.input_dir)
    output_dir = os.path.abspath(args.output_dir)
    os.makedirs(output_dir, exist_ok=True)

    refset_dirs = [d for d in os.listdir(input_dir) if d.startswith("seqRefSet_")]
    if not refset_dirs:
        print(f"No seqRefSet_* folders found in {input_dir}. Exiting.")
        return

    for folder in refset_dirs:
        folder_path = os.path.join(input_dir, folder)
        refset_num = folder.split("_")[-1]
        output_file = os.path.join(output_dir, f"domains_RefSet_{refset_num}.tsv")
        process_seqref_folder(folder_path, output_file)

    print("Domain extraction complete.")

if __name__ == "__main__":
    main()
