#!/bin/bash

# Variables
waitevery=30
seqs_dir="seqs"
IPS_out_dir="IPS_out"

# Create IPS_out directory if it doesn't exist
mkdir -p "${IPS_out_dir}"

# Count the number of seqRefSet directories
num_seq_ref_set=$(ls -d "${seqs_dir}"/seqRefSet* | wc -l)

# Create output directories based on the number of seqRefSet directories
for ((i=1; i<=$num_seq_ref_set; i++)); do
    out_dir="${IPS_out_dir}/out${i}"
    mkdir -p "${out_dir}"
done

for seq_ref_set_dir in "${seqs_dir}"/*; do
    if [ -d "$seq_ref_set_dir" ]; then
        echo "Processing seqRefSet directory: $seq_ref_set_dir"
        seq_ref_set_name=$(basename "${seq_ref_set_dir}")
        out_dir="${IPS_out_dir}/out${seq_ref_set_name}"
        
        
        for fa_file in "${seq_ref_set_dir}"/*.fa; do
            if [ -f "$fa_file" ]; then
                echo "Processing file: $fa_file"
                genome_name=$(basename "${fa_file}" .fa)
                
                python3 iprscan5_urllib3.py \
                    --goterms \
                    --pathways \
                    --email=kinzaamjad456@gmail.com \
                    --outfile="${out_dir}/${genome_name}.xml" \
                    --outformat=xml \
                    "${fa_file}"
                
                
                ((i++ % waitevery == 0)) && wait
            fi
        done
    fi
done


wait

exit 0
