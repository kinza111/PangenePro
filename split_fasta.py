#!/usr/bin/env python3

import os
import re
from Bio import SeqIO

def explode_fasta():
    # Configuration
    minSeqLength = 80

    # Get the current working directory
    current_dir = os.getcwd()

    # Define input and output directories
    input_dir = os.path.join(current_dir, "OP")  
    output_dir = os.path.join(current_dir, "seqs") 

    # Create the output directory if it doesn't exist
    os.makedirs(output_dir, exist_ok=True)

    # Iterate through each RefSet folder in the input directory
    for refset_folder in os.listdir(input_dir):
        refset_path = os.path.join(input_dir, refset_folder)
        if not os.path.isdir(refset_path):
            continue

        # Find all FASTA files in the RefSet folder
        fasta_files = [file for file in os.listdir(refset_path) if file.endswith(('.fasta', '.fa', '.faa'))]
        # If no FASTA files are found, skip this RefSet
        if not fasta_files:
            print(f"No FASTA files found in {refset_path}. Skipping.")
            continue

        # Iterate through each FASTA file
        for fasta_file in fasta_files:
            count = 0

            # Output directory for domain analysis results
            output_dir_refset = os.path.join(output_dir, f"seq{refset_folder}")
            os.makedirs(output_dir_refset, exist_ok=True)

            # Open the current FASTA file for reading
            with open(os.path.join(refset_path, fasta_file), 'r') as fasta:
                # Iterate through all records
                for record in SeqIO.parse(fasta, 'fasta'):
                    count += 1

                    # Only process if sequence is longer than 80 amino acids
                    if len(record.seq) >= minSeqLength:
                        sequence = str(record.seq)

                        if re.match(r'(.*)\*$', sequence):
                            print(f"Warning: Entry {record.id} contains non-trailing stop codon.")

                        # Write sequence to individual fasta file
                        with open(os.path.join(output_dir_refset, f"{record.id}.fa"), 'w') as entry:
                            entry.write(f">{record.id}\n")
                            entry.write(sequence + '\n')

                        if count % 1000 == 0:
                            print('Records processed: ' + str(count))
                    else:
                        print(f"Warning: Entry {record.id} is shorter than minimum sequence length of {minSeqLength}")

            print('====================================')
            print(f'End of file: {fasta_file}')
            print('====================================')
            print('Total records processed: ' + str(count))
            print('====================================') 

if __name__ == "__main__":
    explode_fasta()
