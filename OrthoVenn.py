#!/usr/bin/python

import os
import shutil
import subprocess
import re
import OrthoVennGSEA
import VennOverlapping
import sys
import getopt

def checkProteinFastaFileName(input_fasta_folder):
    flag = 1
    fasta_name_list = os.listdir(input_fasta_folder)
    for i in fasta_name_list:
        if(i.endswith(".fasta") or i.endswith(".fa") or i.endswith(".fas") or i.endswith(".fasta.tar.gz") or i.endswith(".asta.zip")):
            continue
        else:
            print("The name of file (" + i + ") is not right!")
            flag = 0
    if(flag == 0):
        return False
    else:
        return True

def checkProteinSequence(input_fasta_folder):
    flag = 1
    fasta_name_list = os.listdir(input_fasta_folder)
    for i in fasta_name_list:
        id_table = {}
        seq_num = 0
        seq_id = ""
        seq = ""
        with open(os.path.join(input_fasta_folder, i)) as file:
            for line in file:
                line = line.rstrip("\n")
                if(line.startswith(">")):
                    if(seq_id != ""):
                        seq_id_array_1 = seq_id.split(" ")
                        seq_id_array_2 = seq_id_array_1[0].split("|")
                        if(seq_id_array_2[0] in id_table):
                            print("The file (" + i + ") contains duplicate id!")
                            flag = 0
                            break
                    seq_id = line.replace(">", "")
                    seq = ""
                else:
                    seq += line
    if(flag == 0):
        return False
    else:
        return True

def mergeFile(index_file_path, database_protein_folder, merge_file_path):
    merge_file = open(merge_file_path, "w")
    with open(index_file_path) as index_file:
        for line in index_file:
            line = line.rstrip("\n")
            array = line.split("\t")
            if(array[2] == "database"):
                with open(os.path.join(database_protein_folder, array[1] + ".fasta")) as file:
                    for line in file:
                        merge_file.write(line)
            else:
                with open(array[2]) as file:
                    for line in file:
                        merge_file.write(line)
    merge_file.close()

def mergeGFFInput(gff_folder, mcscanx_folder):
    all_gff = open(os.path.join(mcscanx_folder, "all.gff"), "w")
    file_names = os.listdir(gff_folder)
    for i in file_names:
        with open(os.path.join(gff_folder, i)) as file:
            for line in file:
                all_gff.write(line)
    all_gff.close()

def getSpeciesProteinCount(merge_file_path):
    species_protein_count = {}
    with open(merge_file_path) as merge_file:
        for line in merge_file:
            line = line.rstrip("\n")
            if(">" in line):
                array = line.replace(">", "").split("|")
                if(array[0] in species_protein_count):
                    count = species_protein_count[array[0]] + 1
                    species_protein_count[array[0]] = count
                else:
                    species_protein_count[array[0]] = 1
    return species_protein_count

def readAlignmentTable(alignment_file):
    query_target_table = {}
    with open(alignment_file) as file:
        for line in file:
            line = line.rstrip("\n")
            array = line.split("\t")
            query_target_table[array[0]] = array[1]
    return query_target_table

def readAnnotationTable(annotation_file):
    annotation_table = {}
    with open(annotation_file) as file:
        for line in file:
            line = line.rstrip("\n")
            array = line.split("\t")
            annotation_table[array[0]] = line
    return annotation_table

def formatFasta(input_fasta_folder, format_fasta_floder):
    fasta_name_list = os.listdir(input_fasta_folder)
    for i in fasta_name_list:
        species_name = i.replace(".fasta", "").replace(".fa", "").replace(".tar.gz", "").replace(".zip", "")
        format_fasta_file = open(os.path.join(format_fasta_floder, species_name + ".fasta"), "w")
        with open(os.path.join(input_fasta_folder, i.replace(".tar.gz", "").replace(".zip", ""))) as file:
            for line in file:
                line = line.rstrip("\n")
                if(">" in line):
                    array = line.split(" ")
                    format_fasta_file.write(">" + species_name + "|" + array[0].replace(">", "") + "\n")
                else:
                    format_fasta_file.write(line + "\n")
        format_fasta_file.close()

def readIndexTable(index_file_path):
    index_dict = {}
    with open(index_file_path) as index_file:
        for line in index_file:
            line = line.rstrip("\n")
            array = line.split("\t")
            index_dict[array[0]] = array[1]
    return index_dict

def readIndexFastaPath(index_file_path):
    index_dict = {}
    with open(index_file_path) as index_file:
        for line in index_file:
            line = line.rstrip("\n")
            array = line.split("\t")
            index_dict[array[0]] = array[2]
    return index_dict

def readSpeciesList(input_fasta_folder, species_id_file_path):
    species_list = []

    if(os.path.exists(species_id_file_path)):
        with open(species_id_file_path) as species_id_file:
            for line in species_id_file:
                line = line.rstrip("\n")
                if(line == "no_select"):
                    break
                species_list.append(line)
        return species_list
    else:
        fasta_file_names = os.listdir(input_fasta_folder)
        for i in fasta_file_names:
            species_list.append(i.replace(".fasta", "").replace(".fa", ""))
        return species_list

def generateIndexTableFile(input_fasta_folder, species_list):
    index_file_path = os.path.join(work_dir, "index_file")
    index_file = open(index_file_path, "w")
    files_name = os.listdir(input_fasta_folder)
    index = 0
    for i in species_list:
        index += 1
        flag = 0
        for j in files_name:
            if(i == j.replace(".fasta", "").replace(".fa", "")):
                flag = 1
                index_file.write(str(index) + "\t" + i + "\t" + os.path.join(input_fasta_folder, j) + "\n")
                break
        if(flag == 0):
            index_file.write(str(index) + "\t" + i + "\t" + "database" + "\n")
    index_file.close()
    return index_file_path

def generateSpeicesIDFile(species_list, work_dir):
    species_id_file = open(os.path.join(work_dir, "species_id"), "w")
    for i in species_list:
        species_id_file.write(i + "\n")
    species_id_file.close()

def generateBinSpeciesTable(index_file_path):
    bin_species = {}
    species_list = []
    with open(index_file_path) as index_file:
        for line in index_file:
            line = line.rstrip("\n")
            array = line.split("\t")
            species_list.append(array[2])
    for i in range(1, 2**len(species_list) + 1):
        raw_bin = str(bin(i)).replace("0b", "")

def generateBinSpeciesTable(index_file_path):
    bin_species = {}
    species_list = []
    with open(index_file_path) as index_file:
        for line in index_file:
            line = line.rstrip("\n")
            array = line.split("\t")
            species_list.append(array[1])
    for i in range(1, 2 ** len(species_list)):
        raw_bin = str(bin(i)).replace("0b", "")
        bin_str = ""
        for j in range(0, len(species_list) - len(raw_bin)):
            bin_str += "0"
        bin_str += raw_bin
        species_str = ""
        num = 0;
        for k in bin_str:
            num += 1
            if("1" == k):
                species_str += species_list[num - 1] + "-vs-"
        bin_species[bin_str] = species_str[:-4]
    return bin_species

def overlappingAnnotation(overlap_folder, mcl_file_path, orthagogue_result_folder, index_file_path, species_list, work_dir, output_folder, all_fasta_file_path, database_uniprot_folder, species_group, tmp_folder, diamond_old_version_path, usearch_path, go_slim_path, obo_file_path, species_protein_count):
    # Ensure the overlap folder exists
    if not os.path.exists(overlap_folder):
        os.makedirs(overlap_folder)

    # Ensure the output folder exists
    if not os.path.exists(output_folder):
        os.makedirs(output_folder)

    # File path for overlap_id_table
    overlap_id_table_path = os.path.join(overlap_folder, "overlap_id_table.txt")

    try:
        # Code to generate the overlap_id_table file
        # Replace this with your actual implementation

        # For example:
        with open(overlap_id_table_path, 'w') as overlap_id_table:
            overlap_id_table.write("Example content")

        print("Overlap ID table generated successfully.")
    except Exception as e:
        print("Error generating overlap ID table:", e)

def parseMCLClusterFile(mcl_file_path, species_list, work_dir, output_folder, bin_species_table):
    each_cluster_protein_number = {}
    bin_count_table = {}
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
                        output_file = open(os.path.join(output_folder, k), "a+")
                        output_file.write("cluster" + str(cluster_num) + "\n")
                        output_file.close()
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
            if(bin_str in bin_count_table):
                count = bin_count_table[bin_str] + 1
                bin_count_table[bin_str] = count
            else:
                bin_count_table[bin_str] = 1

            bin_overlap_folder = os.path.join(os.path.join(work_dir, "overlap"), bin_str)
            cluster_file = open(os.path.join(bin_overlap_folder, "cluster"), "a+")
            cluster_file.write("cluster" + str(cluster_num) + "\t" + str(len(array)) + "\t" + line.replace("\t", ";") + "\t" + array[0] + "\n")
            cluster_file.close()

    overlap_id_file = open(os.path.join(work_dir, "overlap_id_table"), "w")
    for key in bin_species_table:
        if(key in bin_count_table):
            overlap_id_file.write(bin_species_table[key] + "\t" + key + "\t" + str(bin_count_table[key]) + "\n")
    overlap_id_file.close()
    return each_cluster_protein_number

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
                    weight = float(cluster_relationships[c_interaction]) + 0.1
                    cluster_relationships[c_interaction] = weight
                else:
                    cluster_relationships[c_interaction] = 0.1

    all_cluster_relationships_file = open(os.path.join(orthagogue_result_folder, "all.cluster.relationship"), "w")
    for key in cluster_relationships:
        all_cluster_relationships_file.write(key.replace("-vsvsvs-", "\t") + "\t" + str(round(cluster_relationships[key],1)) + "\n")
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
            all_interaction_table[array[0] + "-vsvsvs-" + array[1]] = float(array[2])

    all_cluster_relationships_xml = open(os.path.join(orthagogue_result_folder, "all.cluster.relationship_table"), "w")
    for cluster in cluster_list:
        node_list = []
        edge_table = {}
        for key in all_interaction_table:
            temp_array = key.split("-vsvsvs-")
            if(cluster == temp_array[0] or cluster == temp_array[1]):
                node_list.append(temp_array[0])
                node_list.append(temp_array[1])
                id1 = temp_array[0] + "-vsvsvs-" + temp_array[1]
                id2 = temp_array[1] + "-vsvsvs-" + temp_array[0]
                if(id1 in edge_table):
                    temp_w = edge_table[id1]
                    edge_table[id1] = temp_w + float(all_interaction_table[key])
                elif(id2 in edge_table):
                    temp_w = edge_table[id2]
                    edge_table[id2] = temp_w + float(all_interaction_table[key])
                else:
                    edge_table[key] = float(all_interaction_table[key])

        node_list = list(set(node_list))

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
            xml += "<data key=\"weight\">" + str(edge_table[edge]) + "</data>"
            temp_w= float(10 * float(edge_table[edge]))
            xml += "<data key=\"label\">" + str(temp_w) + "</data>"
            xml += "</edge>"

        xml += "</graph></graphml>"
        all_cluster_relationships_xml.write(cluster + "\t" + xml + "\n")
    all_cluster_relationships_xml.close()

def getAllCluster(cluster_file_path):
    cluster_table = {}
    line_num = 0
    with open(cluster_file_path) as file:
        for line in file:
            line_num += 1
            line = line.rstrip("\n")
            cluster_table["cluster" + str(line_num)] = line
    return cluster_table

def generatePairwiseTable(species_list, cluster_table, pair_file_path):
    output_file = open(pair_file_path, "w")
    output_file.write("Species" + "\t")
    for i in species_list:
        output_file.write("\t" + i)
    output_file.write("\n")
    for i in species_list:
        output_file.write(i)
        for j in species_list:
            pair_num = 0
            for key in cluster_table:
                if(i in cluster_table[key] and j in cluster_table[key]):
                    pair_num += 1
            output_file.write("\t" + str(pair_num))
        output_file.write("\n")
    output_file.close()

def parseMCScanxOutput(mcscanx_collinearity_file_path, all_collinearity_table_path):
    output_file = open(all_collinearity_table_path, "w")
    with open(mcscanx_collinearity_file_path) as file:
        for line in file:
            line = line.rstrip("\n")
            if(line.startswith("#")):
                continue
            array = line.split("\t")
            if(len(array) == 4):
                output_file.write(array[1] + "\t" + array[2] + "\t" + array[3].strip() + "\n")
    output_file.close()


raw_input_folder = ""
work_dir = ""
bin_path = ""
orthovenn_data_folder = ""
uniprot_data_folder = ""
species_group = ""
evalue = "1e-5"
ivalue = "1.5"
threads = "6"

opts,args = getopt.getopt(sys.argv[1:],'-h-v-i:-w:-b:-r:-u:-s:-e:-f:-t:',['help','version','input_folder_name'])
if len(opts) == 0:
    print("Usage:")
    print("-i the path of protein fasta files folder")
    print("-w the path of output folder")
    print("-b the path of bin folder")
    exit()
for opt_name,opt_value in opts:
    if opt_name in ('-h','--help'):
        print("[*] Help info")
        exit()
    if opt_name in ('-v','--version'):
        print("[*] Version is 0.90 ")
        exit()
    if opt_name in ('-i','--input_folder_name'):
        raw_input_folder = opt_value
    if opt_name in ('-w'):
        work_dir = opt_value
    if opt_name in ('-b'):
        bin_path = opt_value
    if opt_name in ('-r'):
        orthovenn_data_folder = opt_value
    if opt_name in ('-u'):
        uniprot_data_folder = opt_value
    if opt_name in ('-s'):
        species_group = opt_value
    if evalue in ('-e'):
        evalue =opt_value
    if opt_name in ('-f'):
        ivalue = opt_value
    if threads in ('-t'):
        threads = opt_value

gff_input_folder = ""

"""
raw_input_folder = "/home/yiwang/data/Public/OrthoVenn2_data/fasta_5"
gff_input_folder = "/home/yiwang/data/Public/OrthoVenn2_data/gff_5"
work_dir = "/home/yiwang/data/Public/OrthoVenn2_data/work_dir_5"
orthovenn_data_folder = "/home/yiwang/data/database/OrthoVenn/data"
uniprot_data_folder = "/home/yiwang/data/database/OrthoVenn/uniprot"
species_group = "bacteria"
evalue = "1e-5"
ivalue = "1.5"
threads = "6"
"""

input_fasta_folder = ""
if(raw_input_folder == ""):
    input_fasta_folder = os.path.join(work_dir, "input_fasta")
else:
    input_fasta_folder = raw_input_folder
species_id_file_path = os.path.join(work_dir, "species_id")

database_protein_folder = os.path.join(orthovenn_data_folder, species_group, "protein_fasta")
database_alignment_folder = os.path.join(orthovenn_data_folder, species_group, "ublast_result")
database_uniprot_folder = uniprot_data_folder

obo_file_path = os.path.join(bin_path, "go.obo")

diamond_path = os.path.join(bin_path, "diamond")
#diamond_old_version_path = os.path.join(bin_path, "diamond")
orthagogue_path = os.path.join(bin_path, "orthAgogue")
mcl_path = "mcl"
go_slim_path = os.path.join(bin_path, "goslimviewer_standalone/goslimviewer_standalone.pl")
usearch_path = os.path.join(bin_path, "usearch7.0.1090_i86linux32")
mcscanx_path = os.path.join(bin_path, "MCScanX")

circle_num = 0
while(True):
    str_output = subprocess.getstatusoutput("top -b -n 1")
    flag = 0
    for i in str_output:
        i = str(i)
        if (i.endswith("diamond") or i.endswith("orthAgogue")):
            flag = 1
            break

    if (circle_num > 300):
        busy_file = open(os.path.join(work_dir, "busy"), "w")
        busy_file.write("The server is busy, please try again later!")
        busy_file.close()
        exit()
    if (flag == 1):
        circle_num += 1
        time.sleep(60)
    else:
        break

def callCMD(cmd):
    FNULL = open(os.devnull, 'w')
    subprocess.call(cmd, shell=True, stdout=FNULL, stderr=subprocess.STDOUT)
    FNULL.close()
    #subprocess.call(cmd, shell=True)

format_fasta_folder = os.path.join(work_dir, "format_fasta")
all_fasta_folder = os.path.join(work_dir, "all_fasta")
diamond_db_folder = os.path.join(work_dir, "diamond_db")
diamond_result_folder = os.path.join(work_dir, "diamond_result")
orthagogue_result_folder = os.path.join(work_dir, "orthagogue")
overlap_folder = os.path.join(work_dir, "overlap")
output_folder = os.path.join(work_dir, "output_folder")
tmp_folder = os.path.join(work_dir, "tmp")
daa_folder = os.path.join(work_dir, "daa")
mcscanx_folder = os.path.join(work_dir, "mcscanx")
mcscanx_overlap_folder = os.path.join(work_dir, "mcscanx_overlap")
mcscanx_output_folder = os.path.join(work_dir, "mcscanx_output_folder")

if(os.path.exists(work_dir)):
    print("Output folder exists, Please remove it!")
    exit()
else:
    os.makedirs(work_dir)

os.makedirs(format_fasta_folder)
os.makedirs(all_fasta_folder)
os.makedirs(diamond_db_folder)
os.makedirs(diamond_result_folder)
os.makedirs(orthagogue_result_folder)
os.makedirs(overlap_folder)
os.makedirs(output_folder)
#os.makedirs(tmp_folder)
#os.makedirs(daa_folder)
#os.makedirs(mcscanx_folder)
#os.makedirs(mcscanx_overlap_folder)
#os.makedirs(mcscanx_output_folder)

input_fasta_name_list = os.listdir(input_fasta_folder)
for i in input_fasta_name_list:
    if (i.endswith(".tar.gz")):
        tar_cmd = "tar -zxvf " + os.path.join(input_fasta_folder, i) + " -C " + input_fasta_folder
        callCMD(tar_cmd)
        rm_cmd = "rm " + os.path.join(input_fasta_folder, i)
        callCMD(rm_cmd)
    if (i.endswith(".zip")):
        unzip_cmd = "unzip " + os.path.join(input_fasta_folder, i) + " -d " + input_fasta_folder
        callCMD(unzip_cmd)
        rm_cmd = "rm " + os.path.join(input_fasta_folder, i)
        callCMD(rm_cmd)

if(checkProteinFastaFileName(input_fasta_folder) == False):
    exit()
if(checkProteinSequence(input_fasta_folder) == False):
    exit()

species_list = readSpeciesList(input_fasta_folder, species_id_file_path)
if(not os.path.exists(species_id_file_path)):
    species_id_file = open(species_id_file_path, "w")
    for i in species_list:
        species_id_file.write(i + "\n")
    species_id_file.close()
formatFasta(input_fasta_folder, format_fasta_folder)
generateSpeicesIDFile(species_list, work_dir)

index_file_path = generateIndexTableFile(format_fasta_folder, species_list)
index_name_table = readIndexTable(index_file_path)
index_fasta_table = readIndexFastaPath(index_file_path)

all_fasta_file_path = os.path.join(all_fasta_folder, "all.fasta")
mergeFile(index_file_path, database_protein_folder, all_fasta_file_path)
species_protein_count = getSpeciesProteinCount(all_fasta_file_path)

for key in index_fasta_table:
    if(index_fasta_table[key] != "database"):
        diamond_makedb_cmd = diamond_path + " makedb " + "--in " + index_fasta_table[key] + " -d " + os.path.join(diamond_db_folder, index_name_table[key]) + " --threads " + str(threads)
        #print(diamond_makedb_cmd)
        callCMD(diamond_makedb_cmd)

for i in index_fasta_table:
    for j in index_fasta_table:
        if(index_fasta_table[i] == "database" and index_fasta_table[j] != "database"):
            diamond_blastp_cmd = diamond_path + " blastp " + "-d " + os.path.join(diamond_db_folder, index_name_table[j]) + " -q " + os.path.join(database_protein_folder, index_name_table[i] + ".fasta") + " -o " + os.path.join(diamond_result_folder, index_name_table[i] + "-vs-" + index_name_table[j]) + " --max-target-seqs 2000  --evalue 0.01" + " --threads " + str(threads)
            callCMD(diamond_blastp_cmd)
            #diamond_cmd_1 = diamond_old_version_path + " blastp -d " + os.path.join(diamond_db_folder, index_name_table[j]) + ".dmnd" + " -q " + os.path.join(database_protein_folder, index_name_table[i] + ".fasta") + " -a " + os.path.join(daa_folder, index_name_table[i] + "-vs-" + index_name_table[j]) + ".daa" + " -t " + tmp_folder + " --max-target-seqs 1000 --evalue 0.01 --threads 4"
            #callCMD(diamond_cmd_1)
            #diamond_cmd_2 = diamond_old_version_path + " view -a " + os.path.join(daa_folder, index_name_table[i] + "-vs-" + index_name_table[j]) + ".daa" + " -o " + os.path.join(diamond_result_folder, index_name_table[i] + "-vs-" + index_name_table[j])
            #callCMD(diamond_cmd_2)
        if(index_fasta_table[i] != "database" and index_fasta_table[j] == "database"):
            diamond_blastp_cmd = diamond_path + " blastp " + "-d " + os.path.join(database_protein_folder, index_name_table[j]) + " -q " + index_fasta_table[i] + " -o " + os.path.join(diamond_result_folder, index_name_table[i] + "-vs-" + index_name_table[j]) + " --max-target-seqs 2000  --evalue 0.01" + " --threads " + str(threads)
            callCMD(diamond_blastp_cmd)
            #diamond_cmd_1 = diamond_old_version_path + " blastp -d " + os.path.join(database_protein_folder, index_name_table[j]) + ".dmnd" + " -q " + index_fasta_table[i] + " -a " + os.path.join(daa_folder, index_name_table[i] + "-vs-" + index_name_table[j]) + ".daa" + " -t " + tmp_folder + " --max-target-seqs 1000 --evalue 0.01 --threads 4"
            #callCMD(diamond_cmd_1)
            #diamond_cmd_2 = diamond_old_version_path + " view -a " + os.path.join(daa_folder, index_name_table[i] + "-vs-" + index_name_table[j]) + ".daa" + " -o " + os.path.join(diamond_result_folder, index_name_table[i] + "-vs-" + index_name_table[j])
            #callCMD(diamond_cmd_2)
        if(index_fasta_table[i] == "database" and index_fasta_table[j] == "database"):
            cp_cmd = "cp " + os.path.join(database_alignment_folder, index_name_table[i] + "-vs-" + index_name_table[j] + ".gz") + " " + os.path.join(diamond_result_folder, index_name_table[i] + "-vs-" + index_name_table[j] + ".gz")
            callCMD(cp_cmd)
            gzip_cmd = "gzip -d -c " + os.path.join(diamond_result_folder, index_name_table[i] + "-vs-" + index_name_table[j] + ".gz") + " > " + os.path.join(diamond_result_folder, index_name_table[i] + "-vs-" + index_name_table[j])
            callCMD(gzip_cmd)
            rm_cmd = "rm " + os.path.join(diamond_result_folder, index_name_table[i] + "-vs-" + index_name_table[j] + ".gz")
            callCMD(rm_cmd)
        if(index_fasta_table[i] != "database" and index_fasta_table[j] != "database"):
            diamond_blastp_cmd = diamond_path + " blastp " + "-d " + os.path.join(diamond_db_folder, index_name_table[j]) + " -q " + index_fasta_table[i] + " -o " + os.path.join(diamond_result_folder, index_name_table[i] + "-vs-" + index_name_table[j]) + " --max-target-seqs 2000  --evalue 0.01" + " --threads " + str(threads)
            callCMD(diamond_blastp_cmd)
            #diamond_cmd_1 = diamond_old_version_path + " blastp -d " + os.path.join(diamond_db_folder, index_name_table[j]) + ".dmnd" + " -q " + index_fasta_table[i] + " -a " + os.path.join(daa_folder, index_name_table[i] + "-vs-" + index_name_table[j]) + ".daa" + " -t " + tmp_folder + " --max-target-seqs 1000 --evalue 0.01 --threads 4"
            #callCMD(diamond_cmd_1)
            #diamond_cmd_2 = diamond_old_version_path + " view -a " + os.path.join(daa_folder, index_name_table[i] + "-vs-" + index_name_table[j]) + ".daa" + " -o " + os.path.join(diamond_result_folder, index_name_table[i] + "-vs-" + index_name_table[j])
            #callCMD(diamond_cmd_2)

cat_diamond_result_cmd = "cat " + os.path.join(diamond_result_folder, "*") + " > " + os.path.join(diamond_result_folder, "all_alignment_result_raw")
subprocess.call(cat_diamond_result_cmd, shell=True)

all_alignment_result_path = os.path.join(diamond_result_folder, "all_alignment_result")
sort_diamond_result_cmd = "sort " + os.path.join(diamond_result_folder, "all_alignment_result_raw") + " -t $\t -k 1,1 -k 2,2 -k 12rn,12 -o " + all_alignment_result_path
subprocess.call(sort_diamond_result_cmd, shell=True)

orthagogue_cmd = orthagogue_path + " -i " + all_alignment_result_path + " -e " + "5" + " -O " + orthagogue_result_folder + " -c 4"
FNULL = open(os.devnull, 'w')
subprocess.call(orthagogue_cmd, shell=True, stdout=FNULL, stderr=subprocess.STDOUT)

mcl_file_path = os.path.join(orthagogue_result_folder, "out.all.abc.I" + str(ivalue).replace(".", ""))
mcl_cmd = mcl_path + " " + os.path.join(orthagogue_result_folder, "all.abc") + " --abc -I " + str(ivalue) + " -te 4 -show-log n -o " + mcl_file_path
FNULL = open(os.devnull, 'w')
subprocess.call(mcl_cmd, shell=True, stdout=FNULL, stderr=subprocess.STDOUT)

#VennOverlapping.overlappingAnnotation(overlap_folder, mcl_file_path, orthagogue_result_folder, index_file_path, species_list, work_dir, output_folder, all_fasta_file_path, database_uniprot_folder, species_group, tmp_folder, diamond_old_version_path, usearch_path, go_slim_path, obo_file_path, species_protein_count)

if(gff_input_folder == ""):
    exit()

"""
exists_name = {}
for i in species_list:
    for j in species_list:
        name = i + "-vs-" + j
        name_1 = j + "-vs-" + i
        if(name in exists_name or name_1 in exists_name):
            continue
        os.makedirs(os.path.join(mcscanx_folder, name))
        if(name == name_1):
            ln_alignment_cmd = "ln -s " + os.path.join(diamond_result_folder, name) + " " + os.path.join(mcscanx_folder, name, "all.blast")
            callCMD(ln_alignment_cmd)
            ln_gff_cmd = "ln -s " + os.path.join(gff_input_folder, i + ".gff") + " " + os.path.join(mcscanx_folder, name, "all.gff")
            callCMD(ln_gff_cmd)
            mcscanx_cmd = mcscanx_path + " " + os.path.join(mcscanx_folder, name) + "/" + "all" + " -a"
            callCMD(mcscanx_cmd)
        else:
            ln_alignment_cmd = "ln -s " + os.path.join(diamond_result_folder, name) + " " + os.path.join(mcscanx_folder, name, name)
            callCMD(ln_alignment_cmd)
            ln_alignment_cmd_1 = "ln -s " + os.path.join(diamond_result_folder, name_1) + " " + os.path.join(mcscanx_folder, name, name_1)
            callCMD(ln_alignment_cmd_1)
            cat_alignment_cmd = "cat " + os.path.join(mcscanx_folder, name, "*") + " > " + os.path.join(mcscanx_folder, name, "all.blast")
            callCMD(cat_alignment_cmd)
            ln_gff_cmd_1 = "ln -s " + os.path.join(gff_input_folder, i + ".gff") + " " + os.path.join(mcscanx_folder, name, i + ".gff")
            callCMD(ln_gff_cmd_1)
            ln_gff_cmd_2 = "ln -s " + os.path.join(gff_input_folder, j + ".gff") + " " + os.path.join(mcscanx_folder, name, j + ".gff")
            callCMD(ln_gff_cmd_2)
            cat_gff_cmd = "cat " + os.path.join(mcscanx_folder, name, "*.gff") + " > " + os.path.join(mcscanx_folder, name, "all.gff")
            callCMD(cat_gff_cmd)
            mcscanx_cmd = mcscanx_path + " " + os.path.join(mcscanx_folder, name) + "/" + "all" + " -a"
            callCMD(mcscanx_cmd)

file_names_in_mcscanx = os.listdir(mcscanx_folder)
all_collinearity_table_path = os.path.join(mcscanx_folder, "all.abc")
all_collinearity_table_file = open(all_collinearity_table_path, "w")
for i in file_names_in_mcscanx:
    collinearity_output_path = os.path.join(mcscanx_folder, i, "all.collinearity")
    with open(collinearity_output_path) as file:
        for line in file:
            line = line.rstrip("\n")
            if(line.startswith("#")):
                continue
            array = line.split("\t")
            if(len(array) == 4):
                all_collinearity_table_file.write(array[1] + "\t" + array[2] + "\t" + array[3].strip() + "\n")
all_collinearity_table_file.close()

mcl_collinearity_file_path = os.path.join(mcscanx_folder, "out.all.collinearity.abc.I15")
mcl_collinearity_cmd = mcl_path + " " + all_collinearity_table_path + " --abc -I 1.5 -te 4 -show-log n -o " + mcl_collinearity_file_path
callCMD(mcl_collinearity_cmd)

format_collinearity_file_path = os.path.join(mcscanx_folder, "out.all.abc.I15")
format_collinearity_file = open(format_collinearity_file_path, "w")
with open(mcl_collinearity_file_path) as mcl_collinearity_file:
    for line in mcl_collinearity_file:
        line = line.rstrip("\n")
        array = line.split("\t")
        if(len(array) == 1):
            continue
        format_collinearity_file.write(line + "\n")
format_collinearity_file.close()

VennOverlapping.overlappingAnnotation(mcscanx_overlap_folder, format_collinearity_file_path, mcscanx_folder, index_file_path, species_list, work_dir, mcscanx_output_folder, all_fasta_file_path, database_uniprot_folder, species_group, tmp_folder, diamond_old_version_path, usearch_path, go_slim_path, obo_file_path, species_protein_count)
"""

#cp_alignment_cmd = "cp " + all_alignment_result_path + " " + os.path.join(mcscanx_folder, "all.blast")
#callCMD(cp_alignment_cmd)
#mergeGFFInput(gff_input_folder, mcscanx_folder)
#mcscanx_cmd = mcscanx_path + " " + mcscanx_folder + "/" + "all"
#callCMD(mcscanx_cmd)



#all_collinearity_table_path = os.path.join(mcscanx_folder, "all.collinearity.abc")
#mcl_collinearity_file_path = os.path.join(mcscanx_folder, "out.all.collinearity.abc.I15")
#parseMCScanxOutput(os.path.join(mcscanx_folder, "all.collinearity"), all_collinearity_table_path)
#mcl_collinearity_cmd = mcl_path + " " + all_collinearity_table_path + " --abc -I 1.5 -te 4 -show-log n -o " + mcl_collinearity_file_path
#callCMD(mcl_collinearity_cmd)

#format_collinearity_file_path = os.path.join(mcscanx_folder, "format.all.collinearity.abc.I15")
#format_collinearity_file = open(format_collinearity_file_path, "w")
#with open(mcl_collinearity_file_path) as mcl_collinearity_file:
#    for line in mcl_collinearity_file:
#        line = line.rstrip("\n")
#        array = line.split("\t")
#        if(len(array) == 1):
#            continue
#        format_collinearity_file.write(line + "\n")
#format_collinearity_file.close()
