import os
import sys
import gffutils
import tempfile
from Bio import SeqIO
from Bio.Blast.Applications import NcbiblastpCommandline
import subprocess

######################################
#######Process Annotation files#######
def detect_annotation_format(filepath):
    if filepath.endswith(".gtf"):
        return "gtf"
    elif filepath.endswith(".gff") or filepath.endswith(".gff3"):
        return "gff3"
    else:
        raise ValueError(f"Unsupported annotation format for file: {filepath}")

def run_command(command):
    print(f"Running: {command}")
    process = subprocess.Popen(command, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    out, err = process.communicate()
    if process.returncode != 0:
        raise RuntimeError(f"Command failed:\n{err.decode()}")
    return out.decode()

def build_gff_db(gff_file):
    db_file = tempfile.NamedTemporaryFile(delete=False).name
    db = gffutils.create_db(
        gff_file,
        dbfn=db_file,
        force=True,
        keep_order=True,
        merge_strategy='merge',
        sort_attribute_values=True
    )
    return db

def extract_longest_isoforms(db, proteome_fasta, hits, output_fasta):
    protein_lengths = {}
    best_isoforms = {}
    records = SeqIO.to_dict(SeqIO.parse(proteome_fasta, "fasta"))

    for hit in hits:
        cds_list = [f for f in db.features_of_type('CDS') if hit in f.attributes.get('protein_id', []) or hit in f.attributes.get('ID', [])]
        isoform_lengths = {}

        for cds in cds_list:
            parent = cds.attributes.get('Parent') or cds.attributes.get('transcript_id') or ['unknown']
            parent_id = parent[0]
            length = abs(cds.end - cds.start) + 1

            if parent_id not in isoform_lengths:
                isoform_lengths[parent_id] = 0
            isoform_lengths[parent_id] += length

        if isoform_lengths:
            longest_parent = max(isoform_lengths, key=isoform_lengths.get)
            for cds in cds_list:
                parent = cds.attributes.get('Parent') or cds.attributes.get('transcript_id') or ['unknown']
                if parent[0] == longest_parent:
                    protein_id = cds.attributes.get('protein_id', [cds.id])[0]
                    best_isoforms[protein_id] = True

    with open(output_fasta, "w") as out_f:
        kept = 0
        for prot_id in best_isoforms:
            if prot_id in records:
                SeqIO.write(records[prot_id], out_f, "fasta")
                kept += 1
        print(f"Saved {kept} isoform removed non-redundant protein fasta to: {output_fasta}")

def run_blast_and_filter(query_file, ref_data, out_dir, evalue_threshold=1e-5):
    if not os.path.exists(out_dir):
        os.makedirs(out_dir)

    for i, (genome, proteome, annotation) in enumerate(ref_data, start=1):
        print(f"\n Processing RefSet_{i}")
        ref_out = os.path.join(out_dir, f"RefSet_{i}")
        os.makedirs(ref_out, exist_ok=True)

        # Output filenames
        basename = os.path.splitext(os.path.basename(query_file))[0]
        blast_out_default = os.path.join(ref_out, f"{basename}_RefSet_{i}_default.txt")
        blast_out_fmt6 = os.path.join(ref_out, f"{basename}_RefSet_{i}_format6.txt")

        # Default format
        blast_cmd_default = NcbiblastpCommandline(query=query_file, subject=proteome, out=blast_out_default, evalue=evalue_threshold)
        run_command(str(blast_cmd_default))
        print(f" BLAST (default format) complete → {blast_out_default}")

        # Format 6 (tabular)
        blast_cmd_fmt6 = NcbiblastpCommandline(query=query_file, subject=proteome, out=blast_out_fmt6, outfmt=6, evalue=evalue_threshold)
        run_command(str(blast_cmd_fmt6))
        print(f" BLAST (format 6) complete → {blast_out_fmt6}")

        hits = set()
        with open(blast_out_fmt6) as f:
            for line in f:
                fields = line.strip().split("\t")
                if len(fields) >= 2:
                    hits.add(fields[1])

        hits_file = os.path.join(ref_out, "unique_accessions.txt")
        with open(hits_file, "w") as f:
            f.write("\n".join(hits))
        print(f" Unique protein hits: {len(hits)} → {hits_file}")

        # GFF DB and isoform extraction
        fmt = detect_annotation_format(annotation)
        db = build_gff_db(annotation)
        final_fasta = os.path.join(ref_out, "Nr_NIs_PSeqs.faa")
        extract_longest_isoforms(db, proteome, hits, final_fasta)

def main():
    if len(sys.argv) < 6:
        print("Usage: python get_genes.py <query_file> <genome1 proteome1 annotation1> [<genome2 proteome2 annotation2> ...] <output_directory>")
        sys.exit(1)

    query_file = sys.argv[1]
    out_dir = sys.argv[-1]
    ref_data = [(sys.argv[i], sys.argv[i + 1], sys.argv[i + 2]) for i in range(2, len(sys.argv) - 1, 3)]

    try:
        run_blast_and_filter(query_file, ref_data, out_dir)
    except Exception as e:
        print(f" Error: {e}")

if __name__ == "__main__":
    main()
