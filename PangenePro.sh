#!/bin/bash

# Define directory paths
PangenePro_dir="PangenePro"

#Specify the query protein sequence, proteome, genome, and anotation files 
python get_genes.py <query> <genome_1 proteome_1 annotation_1> <genome_2 proteome_2 annotation_2> <genome_x proteome_x annotation_x>
if [ $? -ne 0 ]; then
    echo "Error running get_genes.py"
    exit 1
fi

#InterProScan REST API
python split_fasta.py
if [ $? -ne 0 ]; then
    echo "Error running split_fasta.py"
    exit 1
fi

./interpro_batch_submit.sh
if [ $? -ne 0 ]; then
    echo "Error running interpro_batch_submit.sh"
    exit 1
fi

#Retrieve domains from .xml files and provide final validated members
python get_domains1.py
if [ $? -ne 0 ]; then
    echo "Error running get_domains1.py"
    exit 1
fi

python get_domains2.py
if [ $? -ne 0 ]; then
    echo "Error running get_domains2.py"
    exit 1
fi


#-i the path of protein fasta files folder
#-w the path of output folder
#-b path to ibn folder
python ./OrthoVenn.py -i ./DomAnals -w ./ -b ./bin
if [ $? -ne 0 ]; then
    echo "Error running OrthoVenn2.py"
    exit 1
fi

python Graph.py
if [ $? -ne 0 ]; then
    echo "Error running Graph.py"
    exit 1
fi
echo "Check Current PangenePro directory for Pangenes sets and Graphs"

