from Bio.Blast.Applications import NcbiblastpCommandline
from Bio import SeqIO
import os
import subprocess
import sys


######################################
#######Process Annotation files#######
def process_gff_gtf(file_path, file_format='gff'):
    gene_ids = set()
    with open(file_path, 'r') as file_handle:
        for line in file_handle:
            fields = line.strip().split('\t')
            if len(fields) >= 9:
                attributes = fields[8]
                attr_list = attributes.split(';') if file_format == 'gff' else attributes.split()
                for attr in attr_list:
                    if 'ID=' in attr:
                        gene_id = attr.replace('ID=', '')
                        gene_ids.add(gene_id)
    return gene_ids

def run_command(command, output_file=None, error_file=None):
    stdout_redirect = f"> {output_file}" if output_file else ""
    stderr_redirect = f"2> {error_file}" if error_file else ""
    modified_command = f"{command} {stdout_redirect} {stderr_redirect}"
    process = subprocess.Popen(modified_command, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    out, err = process.communicate()
    if output_file:
        with open(output_file, 'w') as out_file:
            out_file.write(out.decode())
    if error_file:
        with open(error_file, 'w') as err_file:
            err_file.write(err.decode())
    return out.decode(), err.decode()

#######################
#######Run BLAST#######
def run_blast(query_file, ref_data, out_dir, blast_type, evalue_threshold=1e-5):
    if not os.path.isfile(query_file):
        raise FileNotFoundError(f"Query file '{query_file}' does not exist")
    for ref_set in ref_data:
        genome_file, proteome_file, annotation_file = ref_set
        if not all(os.path.isfile(file) for file in ref_set):
            raise FileNotFoundError(f"One or more files in reference set {ref_set} do not exist")

    if not os.path.exists(out_dir):
        os.makedirs(out_dir)

    for i, ref_set in enumerate(ref_data, start=1):
        genome_file, proteome_file, annotation_file = ref_set
        ref_set_out_dir = os.path.join(out_dir, f"RefSet_{i}")
        os.makedirs(ref_set_out_dir, exist_ok=True)

###################################
######Running BLAST commands#######

        default_out_file = os.path.join(ref_set_out_dir, f"{os.path.splitext(os.path.basename(query_file))[0]}_RefSet_{i}_default.txt")
        format6_out_file = os.path.join(ref_set_out_dir, f"{os.path.splitext(os.path.basename(query_file))[0]}_RefSet_{i}_format6.txt")

        if blast_type == "blastp":
            blast_cline_default = NcbiblastpCommandline(query=query_file, subject=proteome_file, out=default_out_file)
            blast_cline_format6 = NcbiblastpCommandline(query=query_file, subject=proteome_file, out=format6_out_file, outfmt=6, evalue=evalue_threshold)

            print(f"Running {blast_type} for Reference Set {i} (format 6 output) with E-value filter ({evalue_threshold})...")
            run_command(str(blast_cline_default), output_file=os.path.join(ref_set_out_dir, "blast_default_output.txt"), error_file=os.path.join(ref_set_out_dir, "blast_default_error.txt"))
            print(f"{blast_type} (default output) complete. Results saved to {default_out_file}")

            print(f"Running {blast_type} for Reference Set {i} (format 6 output)...")
            run_command(str(blast_cline_format6), output_file=os.path.join(ref_set_out_dir, "blast_format6_output.txt"), error_file=os.path.join(ref_set_out_dir, "blast_format6_error.txt"))
            print(f"{blast_type} (format 6 output) complete. Results saved to {format6_out_file}")

#######################################
######Extracting unique Accessions#####

            unique_accessions = set()
            with open(format6_out_file, 'r') as format6_file:
                for line in format6_file:
                    fields = line.strip().split('\t')
                    if len(fields) >= 2:
                        unique_accessions.add(fields[1])

            unique_accessions_file = os.path.join(ref_set_out_dir, f"unique_accessions_RefSet_{i}.txt")
            with open(unique_accessions_file, 'w') as unique_accessions_out:
                for accession in unique_accessions:
                    unique_accessions_out.write(f"{accession}\n")

            print(f"Non-redundant accession numbers saved to {unique_accessions_file}")

############################
######Isoform removal#######

            gene_ids = process_gff_gtf(annotation_file, file_format='gff' if 'gff' in annotation_file else 'gtf')

            command1 = f"grep -Fwf {unique_accessions_file} {annotation_file} | awk '$3 == \"CDS\" {{print}}' > {os.path.join(ref_set_out_dir, 'Extracted_cds.gff')}"
            run_command(command1)

            command2 = f'''
            awk -F '[\t ;=]' '
              BEGIN {{OFS="\t"; print "Chromosome", "GeneID", "rna", "protein_id", "Start", "End";}}
              $3 == "CDS" {{
                chromosome = $1;
                start_position = $4;
                end_position = $5;

                match($0, /gene=([^;]+)/, gene_match);
                gene_id = gene_match[1];

                match($0, /Parent=rna-([^;]+)/, rna_match);
                rna_id = rna_match[1];  

                match($0, /protein_id=([^;]+)/, protein_match);
                protein_id = protein_match[1];

                print chromosome, gene_id, rna_id, protein_id, start_position, end_position;
              }}
            ' {os.path.join(ref_set_out_dir, 'Extracted_cds.gff')} > {os.path.join(ref_set_out_dir, 'Table1.txt')}
            '''
            run_command(command2)

            command3 = f"cut -f 4 {os.path.join(ref_set_out_dir, 'Table1.txt')} | seqtk subseq {proteome_file} - | awk '/^>/ {{if (seqlen) print id, seqlen; id = substr($1, 2); seqlen = 0; next;}} {{seqlen += length($0)}} END {{if (seqlen) print id, seqlen;}}' | sort -k1,1 > {os.path.join(ref_set_out_dir, 'Proteins_length.txt')} && awk 'FNR == NR {{lengths[$1]=$2; next}} {{print $0, (FNR==1) ? \"protein length\" : lengths[$4]}}' {os.path.join(ref_set_out_dir, 'Proteins_length.txt')} {os.path.join(ref_set_out_dir, 'Table1.txt')} > {os.path.join(ref_set_out_dir, 'Table2.txt')}"
            run_command(command3)

            command4 = f"awk 'NR==1 {{print $0; next}} {{key = $1 FS $2 FS $3 FS $4; if (lengths[key] == \"\" || lengths[key] < $6) lengths[key] = $6; data[key] = $0}} END {{for (key in data) print data[key]}}' {os.path.join(ref_set_out_dir, 'Table2.txt')} > {os.path.join(ref_set_out_dir, 'Table_Nr_NIs.txt')}"
            run_command(command4)

            command5 = f"awk 'NR > 1 {{print $4}}' {os.path.join(ref_set_out_dir, 'Table_Nr_NIs.txt')} > {os.path.join(ref_set_out_dir, 'Nr_NIs_Accessions.txt')}"
            run_command(command5)

            command6 = f"seqtk subseq {proteome_file} {os.path.join(ref_set_out_dir, 'Nr_NIs_Accessions.txt')} > {os.path.join(ref_set_out_dir, 'Nr_NIs_PSeqs.faa')}"
            run_command(command6)

        else:
            print(f"Invalid BLAST type: {blast_type}. Skipping Reference Set {i}.")

def main():
    if len(sys.argv) < 6:
        print("Usage: python get_genes.py <query_file> <genome_file1 proteome_file1 annotation_file1> [<genome_file2 proteome_file2 annotation_file2> ...] <output_directory>")
        sys.exit(1)

    query_file = sys.argv[1]
    ref_data = [(sys.argv[i], sys.argv[i + 1], sys.argv[i + 2]) for i in range(2, len(sys.argv) - 1, 3)]
    out_dir = sys.argv[-1]
    blast_type = "blastp"

    try:
        run_blast(query_file, ref_data, out_dir, blast_type)
    except Exception as e:
        print(f"Error: {e}")

if __name__ == "__main__":
    main()
