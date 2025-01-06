#!/usr/bin/python

import os
import subprocess
import getopt
import sys

def readSpeciesList(index_file_path):
    species_list = []
    with open(index_file_path) as index_file:
        for line in index_file:
            line = line.rstrip("\n")
            array = line.split("\t")
            species_list.append(array[1])
    return species_list

def parseMCLClusterFile(mcl_file_path, species_list, overlap_folder):
    each_cluster_protein_number = {}
    cluster_num = 0
    with open(mcl_file_path) as mcl_file:
        for line in mcl_file:
            cluster_num += 1
            line = line.rstrip("\n")
            array = line.split("\t")
            each_cluster_protein_number["cluster" + str(cluster_num)] = float(len(array))

            index_list = [] #the index number of species name
            format_index_list = [] #the index number of species name by remove duplication
            for i in array:
                sub_array = i.split("|")
                index_num = 0
                for k in species_list:
                    index_num += 1
                    if(sub_array[0] == k):
                        index_list.append(index_num)
                        break
            format_index_list = list(set(index_list))

            bin_str = ""
            for i in range(1, len(species_list) + 1):
                flag = 0
                for j in format_index_list:
                    if(i == j):
                        flag = 1
                        break
                if(flag == 1):
                    bin_str += "1"
                else:
                    bin_str += "0"

            bin_overlap_folder = os.path.join(overlap_folder, bin_str)
            cluster_file = open(os.path.join(bin_overlap_folder, "cluster"), "a+")
            cluster_file.write("cluster" + str(cluster_num) + "\t" + line.replace("\t", ";") + "\t" + array[0] + "\n")
    return each_cluster_protein_number

def formatFasta(input_fasta_folder, format_fasta_floder):
    species_protein_count = {}
    fasta_name_list = os.listdir(input_fasta_folder)
    for i in fasta_name_list:
        species_name = i.replace(".fasta", "").replace(".fa", "")
        count = 0
        format_fasta_file = open(os.path.join(format_fasta_floder, i), "w")
        with open(os.path.join(input_fasta_folder, i)) as file:
            for line in file:
                line = line.rstrip("\n")
                if(">" in line):
                    count += 1
                    format_fasta_file.write(">" + species_name + "|" + line.replace(">", "") + "\n")
                else:
                    format_fasta_file.write(line + "\n")
        species_protein_count[species_name] = count
        format_fasta_file.close()
    return species_protein_count

def doStatistc(mcl_file_path, species_list, species_protein_count, work_dir):
    species_cluster_count = {}
    species_protein_in_cluster_count = {}
    for i in species_list:
        species_protein_in_cluster_count[i] = 0

    with open(mcl_file_path) as mcl_file:
        for line in mcl_file:
            line = line.rstrip("\n")
            array = line.split("\t")
            for i in species_list:
                if(i in line):
                    if(i in species_cluster_count):
                        count = species_cluster_count[i] + 1
                        species_cluster_count[i] = count
                    else:
                        species_cluster_count[i] = 1
            for j in array:
                sub_array = j.split("|")
                if(sub_array[0] in species_protein_in_cluster_count):
                    count = species_protein_in_cluster_count[sub_array[0]] + 1
                    species_protein_in_cluster_count[sub_array[0]] = count

    print(len(species_protein_count))
    for key in species_protein_count:
        print(key + " " + str(species_protein_count[key]))

    statistic_file = open(os.path.join(work_dir, "speciesStatistic"), "w")
    for i in species_list:
        statistic_file.write(i + "\t" + str(species_protein_count[i]) + "\t" + str(species_cluster_count[i]) + "\t" + str(species_protein_count[i] - species_protein_in_cluster_count[i]) + "\n")
    statistic_file.close()

def generateCytoscapeXML(ids_str, all_interactions_table, species_list):
    xml = ""
    xml += "<graphml>"
    xml += "<key id=\"label\" for=\"node\" attr.name=\"label\" attr.type=\"string\"/>"
    xml += "<key id=\"species\" for=\"node\" attr.name=\"species\" attr.type=\"string\"/>"
    xml += "<key id=\"label\" for=\"edge\" attr.name=\"label\" attr.type=\"string\"/>"
    xml += "<key id=\"weight\" for=\"edge\" attr.name=\"weight\" attr.type=\"double\"/>"
    xml += "<graph edgedefault=\"undirected\">"

    ids_array = ids_str.split("\t")
    for i in ids_array:
        temp_array = i.split("|")
        xml += "<node id=\"" + temp_array[1] + "\">"
        xml += "<data key=\"label\">" + temp_array[1] + "</data>"
        index = 0
        flag = 0
        for j in species_list:
            index += 1
            if(temp_array[0] == j):
                xml += "<data key=\"species\">" + "species" + str(index) + "</data>"
                flag = 1
                break
        if(flag == 0):
            xml += "<data key=\"species\">" + temp_array[0] + "</data>"
            xml += "</node>"

    for i in ids_array:
        for j in ids_array:
            if((i + "-vsvsvs-" + j) in all_interactions_table):
                name_array1 = i.split("|")
                name_array2 = j.split("|")
                xml += "<edge source=\"" + name_array1[1] + "\" target=\"" + name_array2[1] + "\">"
                xml += "<data key=\"weight\">" + all_interactions_table[i + "-vsvsvs-" + j] + "</data>"
                xml += "</edge>"

    xml += "</graph></graphml>"
    return xml

def generateClusterProteinRelationship(orthagogue_result_folder, mcl_file_path, all_interaction_file_path, species_list):
    all_interactions_table = {}
    with open(all_interaction_file_path) as all_interaction_file:
        for line in all_interaction_file:
            line = line.rstrip("\n")
            array = line.split("\t")
            if((array[0] + "-vsvsvs-" + array[1]) in all_interactions_table or (array[1] + "-vsvsvs-" + array[0]) in all_interactions_table):
                continue
            else:
                all_interactions_table[array[0] + "-vsvsvs-" + array[1]] = array[2]

    all_cluster_interactions_file = open(os.path.join(orthagogue_result_folder, "all.cluster.interactions"), "w")
    line_num = 0
    with open(mcl_file_path) as mcl_file:
        for line in mcl_file:
            line_num += 1
            line = line.rstrip("\n")
            array = line.split("\t")
            if(len(array) > 1):
                xml = generateCytoscapeXML(line, all_interactions_table, species_list)
                all_cluster_interactions_file.write("cluster" + str(line_num) + "\t" + xml + "\n")
    all_cluster_interactions_file.close()

def generateClusterRelationship(orthagogue_result_folder, mcl_file_path, all_interaction_file_path):
    all_interactions_table = {}
    with open(all_interaction_file_path) as all_interaction_file:
        for line in all_interaction_file:
            line = line.rstrip("\n")
            array = line.split("\t")
            if ((array[0] + "-vsvsvs-" + array[1]) in all_interactions_table or (
                    array[1] + "-vsvsvs-" + array[0]) in all_interactions_table):
                continue
            else:
                all_interactions_table[array[0] + "-vsvsvs-" + array[1]] = array[2]

    id_cluster_table = {}
    line_num = 0
    with open(mcl_file_path) as mcl_file:
        for line in mcl_file:
            line_num += 1
            line = line.rstrip("\n")
            array = line.split("\t")
            if (len(array) > 1):
                for i in array:
                    id_cluster_table[i] = "cluster" + str(line_num)

    cluster_relationships = {}
    for key in all_interactions_table:
        array = key.split("-vsvsvs-")
        if(array[0] in id_cluster_table and array[1] in id_cluster_table):
            if(id_cluster_table[array[0]] != id_cluster_table[array[1]]):
                c_interaction = id_cluster_table[array[0]] + "-vsvsvs-" + id_cluster_table[array[1]]
                if(c_interaction in cluster_relationships):
                    weight = cluster_relationships[c_interaction] + 0.1
                    cluster_relationships[c_interaction] = weight
                else:
                    cluster_relationships[c_interaction] = 0.1

    all_cluster_relationships_file = open(os.path.join(orthagogue_result_folder, "all.cluster.relationship"), "w")
    for key in cluster_relationships:
        all_cluster_relationships_file.write(key.replace("-vsvsvs-", "\t") + "\t" + str(cluster_relationships[key]) + "\n")
    all_cluster_relationships_file.close()

def clusterrelationship2xml(orthagogue_result_folder, all_cluster_relationships_file_path, each_cluster_protein_number):
    cluster_list = []
    all_interaction_table = {}
    with open(all_cluster_relationships_file_path) as all_cluster_relationships_file:
        for line in all_cluster_relationships_file:
            line = line.rstrip("\n")
            array = line.split("\t")
            if(not array[0] in cluster_list):
                cluster_list.append(array[0])
            if (not array[1] in cluster_list):
                cluster_list.append(array[1])
            all_interaction_table[array[0] + "-vsvsvs-" + array[1]] = array[2]

    all_cluster_relationships_xml = open(os.path.join(orthagogue_result_folder, "all.cluster.relationship_table"), "w")
    for cluster in cluster_list:
        node_list = []
        edge_table = {}
        for key in all_interaction_table:
            if(cluster in key):
                node_list.append(cluster)
                edge_table[key] = all_interaction_table[key]

        xml = ""
        xml += "<graphml>"
        xml += "<key id=\"label\" for=\"node\" attr.name=\"label\" attr.type=\"string\"/>"
        xml += "<key id=\"size\" for=\"node\" attr.name=\"size\" attr.type=\"double\"/>"
        xml += "<key id=\"label\" for=\"edge\" attr.name=\"label\" attr.type=\"string\"/>"
        xml += "<key id=\"weight\" for=\"edge\" attr.name=\"weight\" attr.type=\"double\"/>"
        xml += "<graph edgedefault=\"undirected\">"

        for node in node_list:
            xml += "<node id=\"" + node + "\">"
            if(node in each_cluster_protein_number):
                xml += "<data key=\"size\">" + str(each_cluster_protein_number[node]) + "</data>"
            else:
                xml += "<data key=\"size\">" + "2.0" + "</data>"
            xml += "<data key=\"label\">" + node + "</data>"
            xml += "</node>"

        for edge in edge_table:
            temp_array = edge.split("-vsvsvs-")
            xml += "<edge source=\"" + temp_array[0] + "\" target=\"" + temp_array[1] + "\">"
            xml += "<data key=\"weight\">" + edge_table[edge] + "</data>"
            xml += "<data key=\"label\">" + str(10 * edge_table[edge]) + "</data>"
            xml += "</edge>"

        xml += "</graph></graphml>"
        all_cluster_relationships_xml.write(cluster + "\t" + xml + "\n")
    all_cluster_relationships_xml.close()

input_fasta_folder = "/home/yiwang/data/Public/OrthoVenn2_data/fasta_1"
index_file_path = "/home/yiwang/data/var/www/OrthoVenn/temp/b67bd2255e2bfa711b04ea9d6606b395/index_file"
mcl_file_path = "/home/yiwang/data/var/www/OrthoVenn/temp/b67bd2255e2bfa711b04ea9d6606b395/orthagogue/out.all.abc.I15"
all_interaction_file_path = "/home/yiwang/data/Public/OrthoVenn2_data/work_dir_1/orthagogue/all.abc"
overlap_folder = "/home/yiwang/data/Public/OrthoVenn2_data/work_dir_1/overlap"

work_dir = "/home/yiwang/data/var/www/OrthoVenn/temp/b67bd2255e2bfa711b04ea9d6606b395"
format_fasta_folder = os.path.join(work_dir, "format_fasta")
orthagogue_result_folder = os.path.join(work_dir, "orthagogue", "aaa")

#species_list = readSpeciesList(index_file_path)
#species_protein_count = formatFasta(input_fasta_folder, format_fasta_folder)
#doStatistc(mcl_file_path, species_list, species_protein_count, work_dir)
#generateClusterProteinRelationship(orthagogue_result_folder, mcl_file_path, all_interaction_file_path, species_list)
#generateClusterRelationship(orthagogue_result_folder, mcl_file_path, all_interaction_file_path)
#clusterrelationship2xml(orthagogue_result_folder, os.path.join(orthagogue_result_folder, "all.cluster.relationship"))

raw_input_folder = ""
work_dir = ""
orthovenn_data_folder = ""
uniprot_data_folder = ""
species_group = ""
evalue = "1e-5"
ivalue = 1.5
threads = 6
opts,args = getopt.getopt(sys.argv[1:],'-h-v-i:-w:-r:-u:-s:-e:-f:-t:',['help','version','input_folder_name'])
for opt_name,opt_value in opts:
    if opt_name in ('-h','--help'):
        print("[*] Help info")
        exit()
    if opt_name in ('-v','--version'):
        print("[*] Version is 0.60 ")
        exit()
    if opt_name in ('-i','--input_folder_name'):
        raw_input_folder = opt_value
    if opt_name in ('-w'):
        work_dir = opt_value
    if opt_name in ('-r'):
        orthovenn_data_folder = opt_value
    if opt_name in ('-u'):
        uniprot_data_folder = opt_value
    if opt_name in ('-s'):
        species_group = opt_value
    if evalue in ('-e'):
        evalue =opt_value
    if opt_name in ('-f'):
        opt_name = opt_value
    if threads in ('-w'):
        threads = opt_value
