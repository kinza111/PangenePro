#!/usr/python/bin

import os

fasta_folder = "/PangenePro/OrthoVenn/fasta"
gff_folder = "/PangenePro/OrthoVenn/GFF"

fasta_names_array = os.listdir(fasta_folder)
gff_names_array = os.listdir(gff_folder)

for i in gff_names_array:
    id_table = {}
    with open(os.path.join(gff_folder, i)) as gff_file:
        for line in gff_file:
            line = line.rstrip("\n")
            array = line.split("\t")
            id_table[array[1]] = ""

    id_num = 0
    with open(os.path.join(fasta_folder, i.replace(".gff", ".fasta"))) as fasta_file:
        for line in fasta_file:
            line = line.rstrip("\n")
            if(line.startswith(">")):
                id = line.replace(">", "")
                if(id in id_table):
                    id_num += 1

    #print(i + " " + str(len(id_table)) + " " + str(id_num))
    if(len(id_table) != id_num):
        print(i + " " + str(len(id_table)) + " " + str(id_num))
