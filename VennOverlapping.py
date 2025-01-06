#!/usr/bin/python

import os
import subprocess
import re
import OrthoVennGSEA
import OrthoVennNetwork

def test1():
    print("___")

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

def parseMCLClusterFile(mcl_file_path, species_list, work_dir, overlap_folder, output_folder, bin_species_table):
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

            bin_overlap_folder = os.path.join(overlap_folder, bin_str)
            cluster_file = open(os.path.join(bin_overlap_folder, "cluster"), "a+")
            cluster_file.write("cluster" + str(cluster_num) + "\t" + str(len(array)) + "\t" + line.replace("\t", ";") + "\t" + array[0] + "\n")
            cluster_file.close()

    overlap_id_file = open(os.path.join(work_dir, "overlap_id_table"), "w")
    for key in bin_species_table:
        if(key in bin_count_table):
            overlap_id_file.write(bin_species_table[key] + "\t" + key + "\t" + str(bin_count_table[key]) + "\n")
    overlap_id_file.close()
    return each_cluster_protein_number

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

def callCMD(cmd):
    FNULL = open(os.devnull, 'w')
    subprocess.call(cmd, shell=True, stdout=FNULL, stderr=subprocess.STDOUT)
    FNULL.close()
    subprocess.call(cmd, shell=True)

def overlappingAnnotation(overlap_folder, mcl_file_path, orthagogue_or_mcscanx_result_folder, index_file_path, species_list, work_dir, output_folder, all_fasta_file_path, database_uniprot_folder, species_group, tmp_folder, diamond_old_version_path, usearch_path, go_slim_path, obo_file_path, species_protein_count):
    bin_species_table = generateBinSpeciesTable(index_file_path)
    for key in bin_species_table:
        os.makedirs(os.path.join(overlap_folder, key))

    each_cluster_protein_number = parseMCLClusterFile(mcl_file_path, species_list, work_dir, overlap_folder, output_folder, bin_species_table)
    all_interaction_file_path = os.path.join(orthagogue_or_mcscanx_result_folder, "all.abc")

    overlap_folder_list = os.listdir(overlap_folder)
    for i in overlap_folder_list:
        overlap_bin_folder = os.path.join(overlap_folder, i)
        cluster_fasta_id_file = open(os.path.join(overlap_bin_folder, "cluster_fasta_id"), "w")
        if (not os.path.exists(os.path.join(overlap_bin_folder, "cluster"))):
            continue
        selected_id_num = 0
        with open(os.path.join(overlap_bin_folder, "cluster")) as cluster_file:
            for line in cluster_file:
                line = line.rstrip("\n")
                array = line.split("\t")
                cluster_fasta_id_file.write(array[3] + "\n")
                selected_id_num += 1
        cluster_fasta_id_file.close()
        seqtk_cmd = "seqtk subseq -l 50 " + all_fasta_file_path + " " + os.path.join(overlap_bin_folder, "cluster_fasta_id") + " > " + os.path.join(overlap_bin_folder, "cluster_fasta")
        subprocess.call(seqtk_cmd, shell=True)
        if (selected_id_num > 300):
            # diamond_blastp_cmd = diamond_path + " blastp " + "-d " + os.path.join(database_uniprot_folder, species_group + "_annotation") + " -q " + os.path.join(overlap_bin_folder, "cluster_fasta") + " -o " + os.path.join(overlap_bin_folder, "annotation_alignment") + " --max-target-seqs 1  --evalue 0.01" + " --threads " + str(threads)
            # FNULL = open(os.devnull, 'w')
            # subprocess.call(diamond_blastp_cmd, shell=True, stdout=FNULL, stderr=subprocess.STDOUT)
            diamond_cmd_1 = diamond_old_version_path + " blastp -d " + os.path.join(database_uniprot_folder, species_group + "_annotation") + " -q " + os.path.join(overlap_bin_folder, "cluster_fasta") + " -a " + os.path.join(overlap_bin_folder, "annotation_alignment") + ".daa" + " -t " + tmp_folder + " --max-target-seqs 1000 --evalue 0.01 --threads 4"
            callCMD(diamond_cmd_1)
            diamond_cmd_2 = diamond_old_version_path + " view -a " + os.path.join(overlap_bin_folder, "annotation_alignment") + ".daa" + " -o " + os.path.join(overlap_bin_folder, "annotation_alignment")
            callCMD(diamond_cmd_2)
        else:
            usearch_cmd = usearch_path + " -ublast " + os.path.join(overlap_bin_folder, "cluster_fasta") + " -db " + os.path.join(database_uniprot_folder, species_group + "_annotation.fasta") + " -evalue 1e-1 -top_hits_only -blast6out " + os.path.join(overlap_bin_folder, "annotation_alignment")
            FNULL = open(os.devnull, 'w')
            subprocess.call(usearch_cmd, shell=True, stdout=FNULL, stderr=subprocess.STDOUT)

        cluster_table_file = open(os.path.join(overlap_bin_folder, "cluster_table"), "w")
        goslim_input_file = open(os.path.join(overlap_bin_folder, "goslim_input"), "w")
        comp = re.compile("GO:[0-9]+")
        query_target_table = readAlignmentTable(os.path.join(overlap_bin_folder, "annotation_alignment"))
        annotation_table = readAnnotationTable(os.path.join(database_uniprot_folder, species_group + "_annotation"))
        with open(os.path.join(overlap_bin_folder, "cluster")) as cluster_file:
            for line in cluster_file:
                line = line.rstrip("\n")
                array = line.split("\t")
                if (array[3] in query_target_table):
                    cluster_table_file.write(line + "\t" + annotation_table[query_target_table[array[3]]] + "\n")
                    go_array = comp.findall(annotation_table[query_target_table[array[3]]])
                    for g in go_array:
                        goslim_input_file.write(query_target_table[array[3]] + "\t" + g + "\n")
                else:
                    cluster_table_file.write(line + "\t" + "N/A" + "\n")
        cluster_table_file.close()
        goslim_input_file.close()

    overlap_folder_list = os.listdir(overlap_folder)
    for i in overlap_folder_list:
        overlap_bin_folder = os.path.join(overlap_folder, i)
        goslim_input_path = os.path.join(overlap_bin_folder, "goslim_input")
        goslim_cmd = "perl " + go_slim_path + " -i " + goslim_input_path + " -s pir -o " + os.path.join(overlap_bin_folder, "goslim_result")
        FNULL = open(os.devnull, 'w')
        subprocess.call(goslim_cmd, shell=True, stdout=FNULL, stderr=subprocess.STDOUT)

    doStatistc(mcl_file_path, species_list, species_protein_count, work_dir)

    OrthoVennGSEA.doGSEA(os.path.join(work_dir, "overlap"), obo_file_path)
    cluster_table = getAllCluster(mcl_file_path)
    generatePairwiseTable(species_list, cluster_table, os.path.join(orthagogue_or_mcscanx_result_folder, "pairwise.csv"))

    OrthoVennNetwork.generateClusterProteinRelationship(orthagogue_or_mcscanx_result_folder, mcl_file_path, all_interaction_file_path, species_list)
    OrthoVennNetwork.generateClusterRelationship(orthagogue_or_mcscanx_result_folder, mcl_file_path, all_interaction_file_path)
    OrthoVennNetwork.clusterrelationship2JSON(orthagogue_or_mcscanx_result_folder, os.path.join(orthagogue_or_mcscanx_result_folder, "all.cluster.relationship"), each_cluster_protein_number)

    #generateClusterProteinRelationship(orthagogue_or_mcscanx_result_folder, mcl_file_path, all_interaction_file_path, species_list)
    #generateClusterRelationship(orthagogue_or_mcscanx_result_folder, mcl_file_path, all_interaction_file_path)
    #clusterrelationship2xml(orthagogue_or_mcscanx_result_folder, os.path.join(orthagogue_or_mcscanx_result_folder, "all.cluster.relationship"), each_cluster_protein_number)
