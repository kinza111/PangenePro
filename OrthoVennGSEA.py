#!/usr/bin/python

import os
import scipy.stats as ss

def getGOTable(obo_file_path):
    go_table = {}
    id = ""
    alt_id_array = []
    name = ""
    namespace = ""
    with open(obo_file_path) as obo_file:
        for line in obo_file:
            line = line.rstrip("\n")
            if("[Term]" in line):
                if(id != ""):
                    go_table[id] = name + "|" + namespace
                    for i in alt_id_array:
                        go_table[i] = name + "|" + namespace
                alt_id_array = []
            elif(line.startswith("id:")):
                id = line.replace("id: ", "")
            elif(line.startswith("name:")):
                name = line.replace("name: ", "")
            elif(line.startswith("namespace:")):
                namespace = line.replace("namespace: ", "")
            elif(line.startswith("alt_id:")):
                alt_id_array.append(line.replace("alt_id: ", ""))
    go_table[id] = name + "|" + namespace
    return go_table

def getReferenceTable(overlap_folder):
    reference_table = {}
    overlap_folder_list = os.listdir(overlap_folder)
    for i in overlap_folder_list:
        overlap_bin_folder = os.path.join(overlap_folder, i)
        if(not os.path.exists(os.path.join(overlap_bin_folder, "goslim_input"))):
            continue
        with open(os.path.join(overlap_bin_folder, "goslim_input")) as goslim_input_file:
            for line in goslim_input_file:
                line = line.rstrip("\n")
                array = line.split("\t")
                if(array[1] in reference_table):
                    count = reference_table[array[1]] + 1
                    reference_table[array[1]] = count
                else:
                    reference_table[array[1]] = 1
    return reference_table

def getListTable(overlap_bin_folder):
    list_table = {}
    if (os.path.exists(os.path.join(overlap_bin_folder, "goslim_input"))):
        with open(os.path.join(overlap_bin_folder, "goslim_input")) as goslim_input_file:
            for line in goslim_input_file:
                line = line.rstrip("\n")
                array = line.split("\t")
                if (array[1] in list_table):
                    count = list_table[array[1]] + 1
                    list_table[array[1]] = count
                else:
                    list_table[array[1]] = 1
    return list_table

def doGSEA(overlap_folder, obo_file_path):
    go_table = getGOTable(obo_file_path)
    reference_table = getReferenceTable(overlap_folder)
    overlap_folder_list = os.listdir(overlap_folder)
    for i in overlap_folder_list:
        p_value_table = {}
        overlap_bin_folder = os.path.join(overlap_folder, i)
        list_table = getListTable(overlap_bin_folder)
        for key in list_table:
            if(key in reference_table):
                hpd = ss.hypergeom(len(go_table), len(reference_table), reference_table[key])
                p = hpd.pmf(list_table[key])
                p_value_table[key] = p
        sorted_p_value_list = sorted(p_value_table.items(), key=lambda x: x[1])

        enrichment_file = open(os.path.join(overlap_bin_folder, "enrichment.txt"), "w")
        index = 0
        for i in sorted_p_value_list:
            index += 1
            value = float(i[1]) * float(len(sorted_p_value_list)) / float(index)
            if(value <= 0.05 and list_table[i[0]] >= 5):
                enrichment_file.write(i[0] + "\t" + str(list_table[i[0]]) + "\t" + go_table[i[0]].replace("|", "\t") + "\t" + str(i[1]) + "\n")
                if(i[1] > 0.05):
                    break
        enrichment_file.close()
