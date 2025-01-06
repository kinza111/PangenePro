import os
import shutil

def process_file(input_file):
    data = {}
    with open(input_file, 'r') as file:
        next(file)  # Skip the header
        for line in file:
            protein_id, interpro_domain, _ = line.strip().split('\t')
            if protein_id not in data:
                data[protein_id] = set()
            data[protein_id].add(interpro_domain)
    return data

def filter_proteins(data):
    common_domains = set.intersection(*data.values())
    filtered_proteins = {}
    for protein_id, domains in data.items():
        if domains & common_domains:
            filtered_proteins[protein_id] = domains
    return filtered_proteins

def write_output(filtered_proteins, input_file):
    output_prefix = os.path.splitext(input_file)[0]
    filtered_proteins_file = output_prefix + '_filtered_proteins.txt'
    final_pids_file = output_prefix + '_Final_PIDs.txt'
    
    with open(filtered_proteins_file, 'w') as file:
        file.write('Protein ID\tInterPro Domains\n')
        for protein_id, domains in filtered_proteins.items():
            file.write('{}\t{}\n'.format(protein_id, ', '.join(domains)))
    print(f"Filtered proteins saved to '{filtered_proteins_file}'.")

    protein_ids = list(filtered_proteins.keys())
    with open(final_pids_file, 'w') as file:
        file.write('\n'.join(protein_ids))
    print(f"Protein IDs saved to '{final_pids_file}'.")

def extract_sequences(filtered_proteins):
    seqs_directory = 'seqs'
    domanals_directory = '/PangenePro/DomAnals'
    for subdir in os.listdir(seqs_directory):
        if subdir.startswith('seqRefSet_') and os.path.isdir(os.path.join(seqs_directory, subdir)):
            seq_refset_directory = os.path.join(seqs_directory, subdir)
            output_sequences_directory = os.path.join(domanals_directory, subdir)
            os.makedirs(output_sequences_directory, exist_ok=True)
            
            for protein_id in filtered_proteins.keys():
                fasta_file = os.path.join(seq_refset_directory, f"{protein_id}.fa")
                if os.path.exists(fasta_file):
                    shutil.copy(fasta_file, output_sequences_directory)
                    print(f"Sequence file '{protein_id}.fa' copied to '{output_sequences_directory}'.")

def main():
    domanals_directory = '/PangenePro/DomAnals'
    for filename in os.listdir(domanals_directory):
        if filename.startswith('protein_domains_all') and filename.endswith('.txt'):
            input_file = os.path.join(domanals_directory, filename)
            data = process_file(input_file)
            filtered_proteins = filter_proteins(data)
            write_output(filtered_proteins, input_file)
            extract_sequences(filtered_proteins)

if __name__ == "__main__":
    main()
