#!/bin/bash

set -e

# ----------------------
# USER INPUTS
# ----------------------
QUERY_FA="AT_CRK.faa"  ## REPLACE THE NAME OF THE FILE HERE WITH YOUR QUERY GENE FAMILY PROTEIN SEQUENCE FILE WITH, provide same name with .xml extension in the next line. 
QUERY_XML="AT_CRK.xml" 
GENOME_FILES=(
  "x.fna x.faa x.gff.gff" ## Replace these files with your genome files, also specify if these are in different directory
)
EMAIL="kinzaamjad456@gmail.com"

# ----------------------
# OUTPUT DIRECTORIES
# ----------------------
OUTDIR="output"
SPLITDIR="split_seqs"
DOMAINDIR="domains"
VALIDATEDDIR="validated_proteins"
CLUSTERDIR="ortholog_clustering"

echo "========== PangenePro Pipeline =========="

# ----------------------
# TIME REQUIRED AT EACH STEP 
# ----------------------
time_step () {
  local START=$(date +%s)
  echo -e "\n $1"
  eval "$2"
  local END=$(date +%s)
  echo "Completed in $((END - START)) seconds"
}

# ----------------------
# 1. CANDIDATE GENE IDENTIFICATION
# ----------------------
INPUTS=""
for trio in "${GENOME_FILES[@]}"; do
  INPUTS+=" $trio"
done
time_step "Step 1: Identifying candidate genes..." \
"python get_genes.py $QUERY_FA $INPUTS $OUTDIR"

### INTERPROSCAN REST API ###

# ----------------------
# 2. Split FASTA
# ----------------------
time_step "Step 2: Splitting protein FASTAs..." \
"python split_fasta.py --input_dir $OUTDIR --output_dir $SPLITDIR"

# ----------------------
# 3. Run InterProScan (batch)
# ----------------------
time_step "Step 3: Submitting to InterProScan..." \
"bash interpro_batch_submit.sh $EMAIL $SPLITDIR"

# ----------------------
# 4. Parse XML OUTPUT & GENERATE DOMAIN TABLE
# ----------------------
time_step "Step 4: Extracting InterPro domains..." \
"python get_domains1.py --input_dir $SPLITDIR --output_dir $DOMAINDIR"

# ----------------------
# 5. Filter Final Validated Proteins
# ----------------------
time_step "Step 5: Domain-based filtering..." \
"python get_domains2.py --query_xml $QUERY_XML --input_dir $DOMAINDIR --seqs_dir $SPLITDIR --output_dir $VALIDATEDDIR"

# ----------------------
# 6. Cluster + Classify Pangenes
# ----------------------
time_step "Step 6: Ortholog clustering and pangene classification..." \
"python ortholog_clustering/pangene_clustering.py --input_fasta $CLUSTERDIR/combined_validated.faa --validated_dir $VALIDATEDDIR --output_dir $CLUSTERDIR"

# ----------------------
# 7. UpSet Plot
# ----------------------
time_step "Step 7: Generating UpSet plot..." \
"Rscript Graph.R $CLUSTERDIR/pangene_matrix.tsv"

echo -e "\n All steps completed successfully!"
echo "Final results:"
echo "- Pangene Matrix:       $CLUSTERDIR/pangene_matrix.tsv"
echo "- Pangene Summary:      $CLUSTERDIR/pangene_summary.tsv"
echo "- UpSet Plot (PNG):     $CLUSTERDIR/pangene_matrix_upset_plot.png"
