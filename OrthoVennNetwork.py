#!/usr/bin/python

import os
import json

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

def generateCytoscapeJSON(ids_str, all_interactions_table, species_list):
    elements_array = []
    ids_array = ids_str.split("\t")

    for i in ids_array:
        temp_array = i.split("|")
        node = {"data": {"id": temp_array[1], "species": temp_array[0]}}
        elements_array.append(node)

    for i in ids_array:
        for j in ids_array:
            if((i + "-vsvsvs-" + j) in all_interactions_table):
                name_array1 = i.split("|")
                name_array2 = j.split("|")
                edge = {"data": {"id": name_array1[1] + "-" + name_array2[1], "source": name_array1[1], "target": name_array2[1], "weight": all_interactions_table[i + "-vsvsvs-" + j]}}
                elements_array.append(edge)

    elements_json = json.dumps(elements_array)
    return elements_json

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
                network_json = generateCytoscapeJSON(line, all_interactions_table, species_list)
                all_cluster_interactions_file.write("cluster" + str(line_num) + "\t" + network_json + "\n")
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

def clusterrelationship2JSON(orthagogue_result_folder, all_cluster_relationships_file_path, each_cluster_protein_number):
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

    all_cluster_relationships_json = open(os.path.join(orthagogue_result_folder, "all.cluster.relationship_table"), "w")

    for cluster in cluster_list:
        node_list = []
        for key in all_interaction_table:
            temp_array = key.split("-vsvsvs-")
            if(cluster == temp_array[0] or cluster == temp_array[1]):
                node_list.append(temp_array[0])
                node_list.append(temp_array[1])

        node_list = list(set(node_list))

        elements_array = []
        for i in node_list:
            node = {"data": {"id": i, "weight": each_cluster_protein_number[i]}}
            elements_array.append(node)

        for i in node_list:
            for j in node_list:
                if ((i + "-vsvsvs-" + j) in all_interaction_table):
                    edge = {"data": {"id": i + "-" + j, "source": i, "target": j, "weight": all_interaction_table[i + "-vsvsvs-" + j]}}
                    elements_array.append(edge)

        elements_json = json.dumps(elements_array)
        all_cluster_relationships_json.write(cluster + "\t" + elements_json + "\n")
    all_cluster_relationships_json.close()


'''
input_fasta_folder = "/home/yiwang/data/Public/OrthoVenn2_data/fasta_6"
index_file_path = "/home/yiwang/data/Public/OrthoVenn2_data/work_dir_6/index_file"
mcl_file_path = "/home/yiwang/data/Public/OrthoVenn2_data/work_dir_6/orthagogue/out.all.abc.I15"
all_interaction_file_path = "/home/yiwang/data/Public/OrthoVenn2_data/work_dir_6/orthagogue/all.abc"
all_cluster_relationships_file_path = "/home/yiwang/data/Public/OrthoVenn2_data/work_dir_6/orthagogue/all.cluster.relationship"
overlap_folder = "/home/yiwang/data/Public/OrthoVenn2_data/work_dir_6/overlap"

work_dir = "/home/yiwang/data/Public/OrthoVenn2_data/work_dir_6"
format_fasta_folder = os.path.join(work_dir, "format_fasta")
orthagogue_result_folder = os.path.join(work_dir, "orthagogue")
species_id_file_path = os.path.join(work_dir, "species_id")

species_list = readSpeciesList(input_fasta_folder, species_id_file_path)
generateClusterProteinRelationship(orthagogue_result_folder, mcl_file_path, all_interaction_file_path, species_list)
clusterrelationship2JSON(orthagogue_result_folder, all_cluster_relationships_file_path)
'''
