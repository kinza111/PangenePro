#!/usr/bin/env python

import os

input_folder = "/home/wzy/Documents/PLANTS/test"

file_names = os.listdir(input_folder)
print(file_names)
for i in file_names:
    if (i.endswith(".gff3")):  # 该文件以.gff3结尾
        name_array = i.split(".")  # 文件名按.分割，将列表第一个值赋给species_name
        species_name = name_array[0]
        species_name_array = species_name.split("_")  # species_name按_分割，赋给species_name_array
        print(species_name)
        output_file = open(os.path.join(input_folder, species_name + ".gff"), "w")  # 输出文件
        with open(os.path.join(input_folder, i)) as f:  # 打开以.gff3结尾的文件，命名为f
            gene_mRNA = {}
            for line in f:  # 依次读取f的每一行
                line = line.rstrip("\n")  # 删除行尾的换行符
                array = line.split("\t")  # 每一行以\t键分割
                if (len(array) > 3 and array[2] == "mRNA" and array[0] != "Mt" and array[ 0] != "Pt"):
                    sub_array = array[8].split(";")  # 以;分割第9个值
                    gene = sub_array[1].replace("Parent=gene:", "")
                    mRNA = sub_array[0].replace("ID=transcript:", "")
                    gene_mRNA.setdefault(gene, []).append(mRNA)
            for key in gene_mRNA:
                gene_mRNA[key].sort()
                first = gene_mRNA[key][0]
                gene_mRNA[key] = first

        with open(os.path.join(input_folder, i)) as f:
            for line in f:
                line = line.rstrip("\n")  # 删除行尾的换行符
                array = line.split("\t")  # 每一行以\t键分割
                if (len(array) > 3 and array[2] == "mRNA"):  # 如果该行有超过4个值以及第三个值是gene
                    sub_array = array[8].split(";")  # 以;分割第9个值
                    mRNA = sub_array[0].replace("ID=transcript:", "")
                    if (mRNA in gene_mRNA.values()):
                        output_file.write(species_name_array[0][0] + species_name_array[1][0] + array[0] + "\t" +
                                          species_name + "|" +
                                          sub_array[1].replace("Parent=gene:", "").replace("\"", "")
                                          + "\t" + array[3] + "\t" + array[4] + "\n")
        output_file.close()

        for j in file_names:
            if (j.startswith(species_name) and j.endswith(".fa")):
                output_file = open(os.path.join(input_folder, species_name + ".fasta"), "w")  # 输出文件
                seq = ""
                seq_id = ""
                index = 0

                with open(os.path.join(input_folder, j)) as f:
                    for line in f:
                        line = line.rstrip("\n")
                        if (line.startswith(">")):
                            if (not seq_id == ""):
                                if (seq_id in gene_mRNA.values()):
                                    gene = list(gene_mRNA.keys())[list(gene_mRNA.values()).index(seq_id)]
                                    output_file.write(">" + gene + "\n" + seq + "\n")
                                    array = line.split(" ")
                                    seq_id = array[4].replace("transcript:", "")
                                    seq = ""
                                    index += 1
                                array = line.split(" ")
                                seq_id = array[4].replace("transcript:", "")
                                seq = ""
                            array = line.split(" ")
                            seq_id = array[4].replace("transcript:", "")
                        else:
                            seq += line

                if (seq_id in gene_mRNA.values()):
                    gene = list(gene_mRNA.keys())[list(gene_mRNA.values()).index(seq_id)]
                    output_file.write(">" + gene + "\n" + seq + "\n")
                    index += 1
                output_file.close()
