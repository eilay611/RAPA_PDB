import os.path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scipy.stats
from matplotlib.lines import Line2D
from matplotlib.patches import Patch

from pdb_guide_Eilay_jj import *
from rapa_vikas_prediction import protein

username = os.getlogin()

import sys
peogect_path = os.path.dirname(os.path.abspath("contacts_between_spike_to_antibodies.py"))
if  username!="eilay":
    sys.path.insert(0, "C:\\Users\\Eilay Koren\\PycharmProjects\\lab_projects\\LielTools_4")
else:
    sys.path.insert(0, "/home/eilay/PycharmProjects/LielTools")
import PlotTools
import importlib
importlib.reload(PlotTools)

"""
------------------------------------------------------------------------------------------------------------------------
parameters
------------------------------------------------------------------------------------------------------------------------
"""
lab_project_path = os.path.dirname(os.path.abspath("change_file_type.py"))
H1N1_parameters = (0,20,15)#first_index_of_peptides,len_of_peptides,len_of_overlap
H3N2_parameters = (0,20,15)#first_index_of_peptides,len_of_peptides,len_of_overlap
H5N1_parameters = 20#len_of_peptides
significant_folder = os.path.join(lab_project_path,"rapa","seq")
"""
------------------------------------------------------------------------------------------------------------------------
function
------------------------------------------------------------------------------------------------------------------------
"""
def start_end_of_peptide_by_id_number(peptide_id,first_index_of_peptides=0,len_of_peptides=20,len_of_overlap=15,last_index_of_peptides=0):
    """
    :param peptide_id: the id number of the peptide you desire to know the start,end index tupel
    :param first_index_of_peptides: the first index you start the peptides
    :param len_of_peptides: the len of each peptide
    :param len_of_overlap: how many amino acids are overlap
    :param last_index_of_peptides: the last index you start the peptides
    :return: (first_index,last_index) of the desire peptide id in respect to the protein numbering
    """
    first_index = first_index_of_peptides + (len_of_peptides-len_of_overlap) * (peptide_id-1)
    last_index = first_index_of_peptides + len_of_peptides - 1 + (len_of_peptides-len_of_overlap) * (peptide_id-1)
    if last_index_of_peptides != 0:
        if last_index > last_index_of_peptides:
            last_index = last_index_of_peptides
    return (first_index,last_index)

def start_end_of_peptide_by_first_position(first_position,len_of_peptides=20,last_index_of_peptides=0):
    """
    :param peptide_id: the id number of the peptide you desire to know the start,end index tupel
    :param first_index_of_peptides: the first index you start the peptides
    :param len_of_peptides: the len of each peptide
    :param len_of_overlap: how many amino acids are overlap
    :param last_index_of_peptides: the last index you start the peptides
    :return: (first_index,last_index) of the desire peptide id in respect to the protein numbering
    """
    first_index = first_position-1
    last_index = first_index + (len_of_peptides-1)
    if last_index_of_peptides != 0:
        if last_index > last_index_of_peptides:
            last_index = last_index_of_peptides
    return (first_index,last_index)

def the_seq_from_peptide(H1N1_xlxs_path = os.path.join(significant_folder,"H1N1_Cal09.xlsx") ,H3N2_xlxs_path = os.path.join(significant_folder,"H3N2_X31.xlsx"),H5N1_xlxs_path = os.path.join(significant_folder,"H5N1_vietnam.xlsx")):
    H1N1 = pd.read_excel(H1N1_xlxs_path)
    H3N2 = pd.read_excel(H3N2_xlxs_path)
    H5N1 = pd.read_excel(H5N1_xlxs_path)

    def rstrip_exact_kk(a):
        if a.endswith("KK"):
            return a[:-len("KK")]
        return a
    def lstrip_exact_kk(a):
        if a.startswith("KK"):
            return a[len("KK"):]
        return a

    H3N2["seq_no_k"] = H3N2["Sequence"].apply(lstrip_exact_kk)
    H5N1["seq_no_k"] = H5N1["Sequence"].apply(rstrip_exact_kk)
    H1N1["seq_no_k"] = H1N1["Sequence"]#did it manually
    H3N2["seq_no_k_len"] = H3N2["seq_no_k"].apply(lambda x: len(x))
    H5N1["seq_no_k_len"] = H5N1["seq_no_k"].apply(lambda x: len(x))
    H1N1["seq_no_k_len"] = H1N1["seq_no_k"].apply(lambda x: len(x))

    def get_the_bigger_from_list_of_overlap(lst):
        """
        :param lst: list_of_string_with_over_lap:
        example
        for lst = ["AAANNV","NNVSDGSD","bad"]
        :return:  'AAANNVSDGSDbad'
        """
        def overlap(a, b):
            return max(i for i in range(len(b) + 1) if a.endswith(b[:i]))
        res = lst[0]
        for s in lst[1:]:
            o = overlap(res, s)
            res += s[o:]
        return res

    sequences = []
    def add_to_seq(position_of_end_HA,name_HA,name_NA,sequences,data):

        record = SeqRecord(
            Seq(get_the_bigger_from_list_of_overlap(list(data.iloc[:position_of_end_HA, :]["seq_no_k"]))),
            id="",
            name="",
            description=name_HA)
        sequences.append(record)
        record = SeqRecord(
            Seq(get_the_bigger_from_list_of_overlap(list(data.iloc[position_of_end_HA:, :]["seq_no_k"]))),
            id="",
            name="",
            description=name_NA)
        sequences.append(record)

    add_to_seq(111, "H1 A/California/07/2009 | A / H1N1 | H1", "Cal_N1 A/California/07/2009 | A / H1N1 | N1", sequences, H1N1)
    add_to_seq(110, "H3 H3(A/aichi/2/68)", "N2 N2(A/aichi/2/68)", sequences, H3N2)
    add_to_seq(110, "H5 H5(Vn1203)", "Vie_N1 N1(Vn1203)", sequences, H5N1)

    with open(os.path.join(significant_folder, "the_seq_from_pep.fasta"), "w") as output_handle:
        SeqIO.write(sequences, output_handle, "fasta")

def make_significent_db(data_path,is_4_groups=True,output_path=os.path.join(lab_project_path, "rapa", "seq", "significant", "figs")):
    """
    :param data_path:
    :param is_4_groups:
    :param output_path:
    :return: values are index
    """
    if is_4_groups==True:
        ending_for_save = ""
        cols=["IgM_rapalive_rapadead", "IgG_rapalive_rapadead", "IgM_pbslive_pbsdead", "IgG_pbslive_pbsdead",
                 "IgM_rapalive_pbslive", "IgG_rapalive_pbslive", "IgM_rapadead_pbsdead", "IgG_rapadead_pbsdead"]
    else:
        ending_for_save = "_live_dead"
        cols =["IgM_live_dead", "IgG_live_dead"]
    significnt_db_result = pd.DataFrame(
        columns=cols,
        index=["H1", "H3", "H5", "Cal_N1", "N2", "Vie_N1"])  # value are index!!
    significnt_db_result_simple = pd.DataFrame(
        columns=cols,
        index=["H1", "H3", "H5", "Cal_N1", "N2", "Vie_N1"])  # value are index!!
    # dict_of_a_b_c_d = {"a":"rapadead","b":"rapalive","c":"pbsdead","d":"pbslive"}
    # dict_of_a_b_c_d_invert = {v: k for k, v in dict_of_a_b_c_d.items()}

    for path, subdirs, files in os.walk(data_path):
        for name in files:
            file_path = os.path.join(path, name)
            file_name = name
            ab_type = file_name.lstrip("Rapa4")[:3]  # IgM or IgG
            if is_4_groups==True:
                splited = file_name.lower().split("_")
                group = splited[4] + splited[5] + "_" + splited[7] + splited[8]  # like rapalive_rapadead or pbslive_pbsdead
            else:
                group = "live_dead"
            subgroup = ab_type + "_" + group
            xl = pd.ExcelFile(file_path)
            for sheet in xl.sheet_names:
                significant_db = pd.read_excel(file_path, sheet_name=sheet, index_col=0)
                if significant_db.empty:
                    continue
                def f(x):
                    splited_fisher = x.split(", ")  # ['a=1', 'b=9', 'c=5', 'd=3']
                    splited_fisher = [int(x[2:]) for x in splited_fisher]  # [1, 9, 5, 3]
                    tow_groups = group.split("_")
                    if (splited_fisher[0] + 1) / (splited_fisher[1] + 2) > (splited_fisher[2] + 1) / (
                            splited_fisher[3] + 2):
                        return tow_groups[0]
                    else:
                        return tow_groups[1]

                significant_db["Fisher_wining"] = significant_db["Fisher Counts"].apply(f)
                if is_4_groups==True:
                    if sheet == "N2" or sheet == "H3":
                        list_of_list_significant_residius = list(significant_db["Antigen"].map(lambda x: list(
                            range(start_end_of_peptide_by_id_number(int(x.split("_")[-1]), *H3N2_parameters)[0],
                                  start_end_of_peptide_by_id_number(int(x.split("_")[-1]), *H3N2_parameters)[1] + 1))))
                        simple = list(significant_db["Antigen"].map(
                            lambda x: list(start_end_of_peptide_by_id_number(int(x.split("_")[-1]), *H3N2_parameters))))
                        fisher_win = list(significant_db["Fisher_wining"])
                        counter = 0
                        for pep in simple:
                            pep.append(fisher_win[counter])
                            counter += 1
                    else:
                        list_of_list_significant_residius = list(significant_db["Antigen"].map(lambda x: list(
                            range(start_end_of_peptide_by_first_position(int(x.split("_")[-1]), H5N1_parameters)[0],
                                  start_end_of_peptide_by_first_position(int(x.split("_")[-1]), H5N1_parameters)[1] + 1))))
                        simple = list(significant_db["Antigen"].map(
                            lambda x: list(start_end_of_peptide_by_first_position(int(x.split("_")[-1]), H5N1_parameters))))
                        fisher_win = list(significant_db["Fisher_wining"])
                        counter = 0
                        for pep in simple:
                            pep.append(fisher_win[counter])
                            counter += 1
                else:
                    list_of_list_significant_residius = list(significant_db["Antigen"].map(lambda x: list(
                        range(start_end_of_peptide_by_first_position(int(x.split("_")[-1]), H5N1_parameters)[0],
                              start_end_of_peptide_by_first_position(int(x.split("_")[-1]), H5N1_parameters)[1] + 1))))
                    simple = list(significant_db["Antigen"].map(
                        lambda x: list(start_end_of_peptide_by_first_position(int(x.split("_")[-1]), H5N1_parameters))))
                    fisher_win = list(significant_db["Fisher_wining"])
                    counter = 0
                    for pep in simple:
                        pep.append(fisher_win[counter])
                        counter += 1

                flat_list_significant_residius = list(
                    set([item for sublist in list_of_list_significant_residius for item in sublist]))
                significnt_db_result.loc[sheet, subgroup] = flat_list_significant_residius
                significnt_db_result_simple.loc[sheet, subgroup] = simple
                # flat_list_of_tupel_significant_residius = [(i,50) for i in flat_list_significant_residius]
    significnt_db_result.to_csv(
        os.path.join(output_path, "significnt_db_result{}.tsv".format(ending_for_save)), "\t")
    significnt_db_result_simple.to_csv(
        os.path.join(output_path, "significnt_db_result_simple{}.tsv".format(ending_for_save)), "\t")
    return significnt_db_result,significnt_db_result_simple

def make_seq_of_protein(the_seq_fasta_file_path=os.path.join(lab_project_path, "rapa", "seq", "the_seq_from_pep.fasta")):
    seq_of_proteins = {}
    for record in SeqIO.parse(the_seq_fasta_file_path, format="fasta"):
        eilay_name = record.description.split(" ")[0]
        seq_of_proteins[eilay_name] = record.seq
    lenght_of_proteins = {}
    for key, value in seq_of_proteins.items():
        lenght_of_proteins[key] = len(value)
    return seq_of_proteins,lenght_of_proteins

def make_protein_dict(significnt_db_result,significnt_db_result_simple,lenght_of_proteins):
    protein_dict = {}
    protein_dict_simple = {}
    for protein in significnt_db_result.index:
        temp_db = pd.DataFrame(index=significnt_db_result.columns, dtype=float,
                               columns=list(range(1, lenght_of_proteins[protein] + 1)))
        temp_db_simple = pd.DataFrame(index=significnt_db_result.columns, dtype=float,
                                      columns=list(range(1, lenght_of_proteins[protein] + 1)))
        counter = 1
        for subgroup in significnt_db_result.columns:
            if type(significnt_db_result.loc[protein, subgroup]) is list:
                for index_of_protein in significnt_db_result.loc[protein, subgroup]:
                    try:
                        a = temp_db.loc[subgroup, index_of_protein + 1]
                        temp_db.loc[subgroup, index_of_protein + 1] = counter
                    except KeyError:
                        continue
                for peptide_tupel in significnt_db_result_simple.loc[protein, subgroup]:
                    for pep_index in range(peptide_tupel[0], peptide_tupel[1] + 1):
                        try:
                            if pd.isna(temp_db_simple.loc[subgroup, pep_index + 1]):
                                temp_db_simple.loc[subgroup, pep_index + 1] = peptide_tupel[2]
                            elif temp_db_simple.loc[subgroup, pep_index + 1] == peptide_tupel[2]:
                                temp_db_simple.loc[subgroup, pep_index + 1] = peptide_tupel[2]
                            else:
                                temp_db_simple.loc[subgroup, pep_index + 1] = "ambivalent"
                        except KeyError:
                            continue
            counter = counter + 1
        protein_dict[protein] = temp_db.copy()
        protein_dict_simple[protein] = temp_db_simple.copy()
    return protein_dict,protein_dict_simple

def make_align_protein_dict(protein_dict,protein_dict_simple,HA_align_path=os.path.join(lab_project_path, "rapa", "seq", "the_seq_from_pep_HA_align_muscel.fas"),NA_align_path=os.path.join(lab_project_path, "rapa", "seq", "the_seq_from_pep_NA_align_muscel.fas")):
    align_protein_dict = {}
    align_protein_dict_simple = {}
    HA_align_db = pd.DataFrame()
    NA_align_db = pd.DataFrame()
    for record in SeqIO.parse(HA_align_path,"fasta"):
        HA_align_db[record.description.split(" ")[0]] = list(record.seq)
    for record in SeqIO.parse(NA_align_path, "fasta"):
        NA_align_db[record.description.split(" ")[0]] = list(record.seq)
    HA_align_db.index = HA_align_db.index + 1
    NA_align_db.index = NA_align_db.index + 1
    align_db = HA_align_db.merge(NA_align_db,how="outer",left_index=True, right_index=True)
    dict_of_maps_of_coloumns={}
    for protein in protein_dict:
        matrix = protein_dict[protein]
        simple_matrix = protein_dict_simple[protein]
        map_of_colomns = {}
        for col in matrix.columns:
            map_of_colomns[col] = align_db.loc[align_db[protein]!="-",protein].iloc[[col-1]].index[0]
        dict_of_maps_of_coloumns[protein] = map_of_colomns
        align_matrix = matrix.rename(columns=map_of_colomns)
        align_simple_matrix = simple_matrix.rename(columns=map_of_colomns)
        align_protein_dict[protein] = fill_col_aligment_with_nan(align_matrix)
        align_protein_dict_simple[protein] = fill_col_aligment_with_nan(align_simple_matrix)
    return align_protein_dict,align_protein_dict_simple,align_db,dict_of_maps_of_coloumns

def fill_col_aligment_with_nan(protein):
    """
    :param protein: protein from protein dictionary
    :return: the same but insure that the cols are in OKEV one by one
    """
    last_col = protein.columns[0]
    for col in protein.columns[1:]:
        if col - last_col != 1:  # their is a jump
            first_insert_pos = last_col + 1
            last_insert_pos = col - 1
            added = pd.DataFrame(np.nan, index=protein.index, columns=range(first_insert_pos, last_insert_pos + 1))
            protein = pd.concat([protein.loc[:, :last_col], added, protein.loc[:, col:]], axis=1)
        last_col = col
    return protein


def make_structral_compare(significnt_db_result_simple,is_4_groups=True,endresult_per_peptide_file_path=os.path.join(lab_project_path,"rapa","endresulte_per_peptide.tsv"),output_path=os.path.join(lab_project_path, "rapa", "seq", "significant", "figs")):
    if is_4_groups==True:
        end_name = ""
    elif is_4_groups == "rearengment_morin":
        end_name = "_rearengment"
    else:
        end_name = "_live_dead"
    structrul_compare = pd.DataFrame()
    structrul_compare_only_sig = pd.DataFrame()
    endresult_per_peptide_file_db = pd.read_csv(endresult_per_peptide_file_path,sep="\t",index_col=0)
    for protein in significnt_db_result_simple.index:#"H1"...
        dict_groups = {}#{"rapadead":[1,6,..]#first positions,...}
        first_positions_in_use = set()
        for comparison in significnt_db_result_simple.columns:
            if not significnt_db_result_simple.loc[protein,comparison] is np.nan:
                for sig_peptide in significnt_db_result_simple.loc[protein,comparison]:#[35, 54, rapadead] [first index,last index,sig high in]
                    group_name = sig_peptide[2] #rapadead or rapalive or live ...
                    dict_groups.setdefault(group_name, [])
                    dict_groups[group_name].append(int(sig_peptide[0])+1)
                    first_positions_in_use.add(int(sig_peptide[0])+1)
        dict_groups["all_other"] = list(set(endresult_per_peptide_file_db.loc[endresult_per_peptide_file_db["conservation_nonbinery_"+protein]!=np.nan,"first_position"]).difference(first_positions_in_use))
        for col in ["sasa_"+protein,"conservation_nonbinery_"+protein]:
            if col == 'sasa_N2':
                continue
            all_sifnifcent = []
            for group_name in dict_groups:#rapadead or rapalive or live ...
                stat, p_val = scipy.stats.ranksums(endresult_per_peptide_file_db.loc[endresult_per_peptide_file_db["first_position"].isin(dict_groups[group_name]),col],endresult_per_peptide_file_db.loc[~endresult_per_peptide_file_db["first_position"].isin(dict_groups[group_name]),col])
                stat, p_val = (round(stat, 3),round(p_val, 3))
                structrul_compare.loc[col.split("_")[0]+"_"+group_name+"_vs_all",protein] = str(stat)+"_"+str(p_val)
                if p_val<0.05:
                    structrul_compare_only_sig.loc[col.split("_")[0] + "_" + group_name + "_vs_all", protein] = str(
                        stat) + "_" + str(p_val)
                all_sifnifcent = all_sifnifcent + dict_groups[group_name]

            stat, p_val = scipy.stats.ranksums(endresult_per_peptide_file_db.loc[
                                                   endresult_per_peptide_file_db["first_position"].isin(
                                                       all_sifnifcent), col],
                                               endresult_per_peptide_file_db.loc[
                                                   ~endresult_per_peptide_file_db["first_position"].isin(
                                                       dict_groups["all_other"]), col])
            stat, p_val = (round(stat, 3), round(p_val, 3))
            structrul_compare.loc[col.split("_")[0] + "_" + "all_significant" + "_vs_all_other", protein] = str(stat) + "_" + str(
                p_val)
            if p_val < 0.05:
                structrul_compare_only_sig.loc[col.split("_")[0] + "_" + "all_significant" + "_vs_all_other", protein] = str(
                    stat) + "_" + str(p_val)
    structrul_compare.to_csv(os.path.join(output_path,"structrul_compare{}.tsv".format(end_name)),"\t")
    structrul_compare_only_sig.to_csv(os.path.join(output_path,"structrul_compare_only_sig{}.tsv".format(end_name)),"\t")
    return structrul_compare,structrul_compare_only_sig

def make_structral_compare_per_position(significnt_db_result_simple,to_plot=True,is_4_groups=True,endresult_per_position_file_path=os.path.join(lab_project_path,"rapa","endresulte_per_position.tsv"),endresult_per_peptide_file_path=os.path.join(lab_project_path,"rapa","endresulte_per_peptide.tsv"),output_path=os.path.join(lab_project_path, "rapa", "seq", "significant", "figs")):
    if is_4_groups==True:
        end_name = ""
    elif is_4_groups == "rearengment_morin":
        end_name = "_rearengment"
    else:
        end_name = "_live_dead"
    structrul_compare = pd.DataFrame()
    structrul_compare_only_sig = pd.DataFrame()
    endresult_per_peptide_file_db = pd.read_csv(endresult_per_peptide_file_path,sep="\t",index_col=0)
    endresult_per_position_file_db = pd.read_csv(endresult_per_position_file_path,sep="\t",index_col=0)
    #endresult_per_position_file_db.index=endresult_per_position_file_db["pep_count"]
    for protein in significnt_db_result_simple.index:#"H1"...
        dict_groups = {}#{"rapadead":[1,6,..]#first positions,...}
        first_positions_in_use = set()
        for comparison in significnt_db_result_simple.columns:
            if not significnt_db_result_simple.loc[protein,comparison] is np.nan:
                for sig_peptide in significnt_db_result_simple.loc[protein,comparison]:#[35, 54, rapadead] [first index,last index,sig high in]
                    group_name = sig_peptide[2] #rapadead or rapalive or live ...
                    dict_groups.setdefault(group_name, [])
                    if (group_name in ["live","dead"] and protein in ["H3","N2"]):#TODO figure out why
                        dict_groups[group_name].append(int(sig_peptide[0]) )
                        first_positions_in_use.add(int(sig_peptide[0]) )
                    else:
                        dict_groups[group_name].append(int(sig_peptide[0])+1)
                        first_positions_in_use.add(int(sig_peptide[0])+1)
        #dict_groups["all_other"] = list(set(endresult_per_peptide_file_db.loc[endresult_per_peptide_file_db["conservation_nonbinery_"+protein]!=np.nan,"first_position"]).difference(first_positions_in_use))

        for col in ["sasa_"+protein,"conservation_nonbinery_"+protein]:
            if col == 'sasa_N2':
                continue
            #all_sifnifcent = []
            for group_name in dict_groups:#rapadead or rapalive or live ...
                the_significent_position_align = list(endresult_per_peptide_file_db.loc[endresult_per_peptide_file_db["first_position"].isin(dict_groups[group_name]),["first_position","last_position"]].apply(lambda x:list(range(x[0],x[1]+1)),axis=1))
                the_non_significent_position_align = list(endresult_per_peptide_file_db.loc[(~endresult_per_peptide_file_db["first_position"].isin(dict_groups[group_name])) & (~endresult_per_peptide_file_db[col].isna()),["first_position","last_position"]].apply(lambda x:list(range(x[0],x[1]+1)),axis=1))#can be added some final fictiv position
                the_significent_position_align = pd.Series([item for sublist in the_significent_position_align for item in sublist])
                the_non_significent_position_align = pd.Series([item for sublist in the_non_significent_position_align for item in sublist])
                the_last_index_of_the_protein = endresult_per_position_file_db["conservation_nonbinery_"+protein].dropna().index[-1]
                the_non_significent_position_align.loc[the_non_significent_position_align > the_last_index_of_the_protein] = np.nan #get rid of fictiv end
                the_significent_position_align.loc[the_non_significent_position_align > the_last_index_of_the_protein] = np.nan #get rid of fictiv end
                the_non_significent_position_align.dropna(inplace=True) #get rid of fictiv end
                the_significent_position_align.dropna(inplace=True) #get rid of fictiv end
                the_significent_score = endresult_per_position_file_db.loc[the_significent_position_align,col].dropna()
                the_non_significent_score = endresult_per_position_file_db.loc[the_non_significent_position_align,col].dropna()
                stat, p_val = scipy.stats.ranksums(the_significent_score,the_non_significent_score)
                if to_plot==True:
                    fig, ax = plt.subplots(1, 1, figsize=(8, 6))
                    colors = ['#aacfcf', '#d291bc']
                    sns.set(font_scale=1)

                    df =pd.DataFrame({col: pd.concat([the_significent_score, the_non_significent_score], axis=0),
                                     'group of positions': np.repeat([group_name, "not "+group_name], [len(the_significent_score),
                                                                        len(the_non_significent_score)])})
                    df.rename(columns={col:col.replace("_", " ").replace("conservation nonbinery","variability")},inplace=True)
                    sns.violinplot(x='group of positions', y=col.replace("_", " ").replace("conservation nonbinery","variability"),
                                data=df,
                                ax=ax, palette=colors, width=0.5
                                )
                    # sns.stripplot(x="group of positions",
                    #               y=col.replace("_", " "),
                    #               color='black',
                    #               alpha=0.3,
                    #               data=df)
                    if p_val<0.05:
                        plt.savefig(os.path.join(output_path,"structrul_compare_sig_{}_{}_{}.png".format(col,group_name,end_name)),dpi=600)
                    else:
                        plt.savefig(os.path.join(output_path,"structrul_compare_{}_{}_{}.png".format(col,group_name,end_name)),dpi=600)
                    plt.close('all')

                stat, p_val = (round(stat, 3),round(p_val, 6))
                structrul_compare.loc[col.split("_")[0]+"_"+group_name+"_vs_all",protein] = str(stat)+"_"+str(p_val)
                if p_val<0.05:
                    structrul_compare_only_sig.loc[col.split("_")[0] + "_" + group_name + "_vs_all", protein] = str(
                        stat) + "_" + str(p_val)
                #all_sifnifcent = all_sifnifcent + dict_groups[group_name]

            # stat, p_val = scipy.stats.ranksums(endresult_per_peptide_file_db.loc[
            #                                        endresult_per_peptide_file_db["first_position"].isin(
            #                                            all_sifnifcent), col],
            #                                    endresult_per_peptide_file_db.loc[
            #                                        ~endresult_per_peptide_file_db["first_position"].isin(
            #                                            dict_groups["all_other"]), col])
            # stat, p_val = (round(stat, 3), round(p_val, 3))
            # structrul_compare.loc[col.split("_")[0] + "_" + "all_significant" + "_vs_all_other", protein] = str(stat) + "_" + str(
            #     p_val)
            # if p_val < 0.05:
            #     structrul_compare_only_sig.loc[col.split("_")[0] + "_" + "all_significant" + "_vs_all_other", protein] = str(
            #         stat) + "_" + str(p_val)
    structrul_compare.to_csv(os.path.join(output_path,"structrul_compare{}.tsv".format(end_name)),"\t")
    structrul_compare_only_sig.to_csv(os.path.join(output_path,"structrul_compare_only_sig{}.tsv".format(end_name)),"\t")
    # minus say that the_non_significent_score is greater
    # plus say that the_significent_score is greater
    if to_plot:
        def split_the_name(x):
            if x == 'nan_nan':
                return np.nan
            return x.split("_")[0]
        def split_the_name2(x):
            if x == 'nan_nan':
                return np.nan
            return x.split("_")[1]
        statist = structrul_compare.applymap(lambda x: float(split_the_name(x)) if (np.all(pd.notnull(x))) else x)
        statist.loc[["conservation" in x for x in statist.index],:] = statist.loc[["conservation" in x for x in statist.index],:]*-1# to make that the higher the score the higher the conservation
        pvals =  structrul_compare.applymap(lambda x:  split_the_name2(x) if (np.all(pd.notnull(x))) else x)
        pvals =  pvals.applymap(lambda x:  "*" *int(-1*math.log10(float(x)*2)) if (np.all(pd.notnull(x))) else x )
        statist.sort_index(inplace=True,axis=0)
        pvals.sort_index(inplace=True,axis=0)

        ax_heatmap = PlotTools.plot_heatmap(statist.rename(lambda x:x.rstrip("_vs_all").replace("_"," ").replace("conservation","variability"),axis='index'), cmap="coolwarm",center=0, figsize=(20, 20),
                 title='', title_fontsize=13, ax=None,
                 font_scale=3, snsStyle='ticks', xRotation=0,
                 yRotation=0,
                 xlabel='protein', ylabel='comparison', colormap_label='rank sum statist',
                 vmin=None, vmax=None, supress_ticks=True,
                 annotate_text=pvals, annotate_fontsize=20,
                 annotation_format="",
                 mask=None, colorbar_ticks=None,
                 hide_colorbar=False,
                 xy_labels_fontsize=None,
                 grid_linewidths=0, grid_linecolor='white')
        #now i will add an horizontal line to separate the sasa from the rest

        for xc,index in enumerate(statist.rename(lambda x:x.rstrip("_vs_all").replace("_"," ").replace("conservation","variability"),axis='index').index):
            if index.startswith("sasa"):
                ax_heatmap.axhline(y=xc, linewidth=4, color='black')
                break

        plt.savefig(os.path.join(output_path,"structrul_compare_plus_is_the_group_higher{}.png".format(end_name)),dpi=600)
    return structrul_compare,structrul_compare_only_sig

def plots(protein_dict,protein_dict_simple,is_4_groups=True,
          is_align=False,morin_2=True,morin3=True,
          endresult_per_position_file_path = os.path.join(lab_project_path, "rapa", "endresulte_per_position.tsv"),dict_of_maps_of_coloumns=None,
          colors_map={"ambivalent": [255 / 255, 193 / 255, 7 / 255],
                      "live": [128 / 255, 115 / 255, 172 / 255],
                      "dead": [216 / 255, 218 / 255, 235 / 255],
                      "diff": [0 / 255, 141 / 255, 255 / 255]}
          ):
    plt.close("all")
    vmax = protein_dict["H1"].shape[0]
    labels = []
    for index in protein_dict["H1"].index:
        labels.append(index.split("_")[1])
        labels.append(index.split("_")[2])
    labels = list(set(labels))
    labels.sort()
    labels.append('ambivalent')
    number_of_lables = len(labels)
    value_to_int = {k:j for k,j in zip( labels,range(1,len(labels)+1))}
    if is_4_groups==True:
        #number_of_lables=5
        #vmax=8
        #value_to_int = {'ambivalent': 1, 'rapadead': 2, 'rapalive': 3, 'pbsdead': 4, 'pbslive': 5}
        #labels = ['ambivalent', 'rapadead', 'rapalive', 'pbsdead', 'pbslive']
        ending_for_save = ""
    elif is_4_groups == "rearengment_morin":
        ending_for_save = "_rearengment"
    else:
        #number_of_lables=3
        #vmax = 2
        #value_to_int = {'ambivalent': 1, 'dead': 2, 'live': 3}
        #labels = ['ambivalent', 'dead', 'live']
        ending_for_save = "_live_dead"
    if is_align:
        align_text = "_align"
        shere_x = "col"
    else:
        align_text = ""
        shere_x = False
    for protein_name, matrix in protein_dict.items():
        grid = PlotTools.plot_clustermap(matrix,
                                         xlabel='Position'+align_text.replace("_"," "), ylabel='compration',
                                         adjRight=0.75, adjBottom=0.1,
                                         hide_cbar=True, vmin=1, vmax=vmax,
                                         title=protein_name, title_fontsize=40,
                                         xy_labels_fontsize=35,
                                         xticklabels=20, cmap=sns.color_palette("Paired")[:vmax],
                                         row_clustering=False, col_clustering=False,
                                         figsize=(20, 15), font_scale=2, xRotation=90)
        ax = grid.ax_heatmap
        for xc in range(0, len(matrix.columns), 5):
            ax.axvline(x=xc, linewidth=0.3, color='b')
        matrix.to_csv(os.path.join(lab_project_path, "rapa", "seq", "significant", "figs", protein_name +ending_for_save+align_text+ ".tsv"), "\t")
        plt.savefig(os.path.join(lab_project_path, "rapa", "seq", "significant", "figs", protein_name +ending_for_save+align_text+ ".png"),dpi=600)
        #plt.show()
        plt.close("all")

    fig, axs = plt.subplots(3, 2, sharey=True,sharex=shere_x, figsize=(65, 35))  # subplots(row, columns)
    for protein_name, matrix in protein_dict.items():
        dictonery = {"H1": (0, 0), "Cal_N1": (0, 1),
                     "H3": (1, 0), "N2": (1, 1),
                     "H5": (2, 0), "Vie_N1": (2, 1)}
        counter_x, counter_y = dictonery[protein_name]
        PlotTools.plot_heatmap(matrix, cmap=sns.color_palette("Paired")[:vmax], figsize=(20, 15),
                               title=protein_name, title_fontsize=40, ax=axs[counter_x, counter_y],
                               font_scale=2, snsStyle='ticks', xRotation=90,
                               yRotation=0,
                               xlabel='Position'+align_text.replace("_"," "), ylabel='compration', colormap_label='',
                               vmin=None, vmax=None, supress_ticks=True,
                               annotate_text=False, annotate_fontsize=8,
                               annotation_format=".2f",
                               mask=None, colorbar_ticks=None,
                               hide_colorbar=True,
                               xy_labels_fontsize=45,
                               grid_linewidths=0, grid_linecolor='white')
    plt.savefig(os.path.join(lab_project_path, "rapa", "seq", "significant", "figs", "all_proteins" +ending_for_save+align_text+ ".png"),dpi=600)
    #plt.show()
    plt.close("all")


    fig, axs = plt.subplots(3, 2, sharey=True,sharex=shere_x, figsize=(65, 35))  # subplots(row, columns)
    for protein_name, matrix in protein_dict_simple.items():
        dictonery = {"H1": (0, 0), "Cal_N1": (0, 1),
                     "H3": (1, 0), "N2": (1, 1),
                     "H5": (2, 0), "Vie_N1": (2, 1)}
        counter_x, counter_y = dictonery[protein_name]
        ax = PlotTools.plot_heatmap(matrix.replace(value_to_int),
                                    cmap=sns.color_palette("Paired", number_of_lables), figsize=(20, 15),
                                    title=protein_name, title_fontsize=40, ax=axs[counter_x, counter_y],
                                    font_scale=3, snsStyle='ticks', xRotation=90,
                                    yRotation=0,
                                    xlabel='Position'+align_text.replace("_"," "), ylabel='compration', colormap_label='',
                                    vmin=1, vmax=number_of_lables, supress_ticks=True,
                                    annotate_text=False, annotate_fontsize=8,
                                    annotation_format=".2f",
                                    mask=None, colorbar_ticks=None,
                                    hide_colorbar=False,
                                    xy_labels_fontsize=45,
                                    grid_linewidths=0, grid_linecolor='white')
        colorbar = ax.collections[0].colorbar
        r = colorbar.vmax - colorbar.vmin
        colorbar.set_ticks([colorbar.vmin + r / number_of_lables * (0.5 + i) for i in range(number_of_lables)])
        # colorbar.set_ticklabels(list(value_to_int.keys()))
        colorbar.set_ticklabels(labels)
    plt.savefig(os.path.join(lab_project_path, "rapa", "seq", "significant", "figs",
                             "all_proteins_with_cmap" + ending_for_save +align_text+ ".png"),dpi=600)
    #plt.show()
    plt.close("all")

    if is_4_groups == "rearengment_morin":
        fig, axs = plt.subplots(3, 2, sharey=True,sharex=shere_x, figsize=(65, 35))  # subplots(row, columns)
        for protein_name, matrix in protein_dict_simple.items():
            dictonery = {"H1": (0, 0), "Cal_N1": (0, 1),
                         "H3": (1, 0), "N2": (1, 1),
                         "H5": (2, 0), "Vie_N1": (2, 1)}
            counter_x, counter_y = dictonery[protein_name]
            value_to_int = {'ambivalent': 1,
                            'rapadead': 2,
                            'rapalive': 3,
                            'dead': 2,
                            'live': 3,
                            'pbsdead': 2,
                            'pbslive': 3}
            the_indexs =['IgM_live_dead',
                         'IgM_rapalive_rapadead',
                         'IgM_pbslive_pbsdead',
                         'IgG_live_dead',
                         'IgG_rapalive_rapadead',
                         'IgG_pbslive_pbsdead']
            index_to_change = {'IgM_live_dead': "",
                               'IgM_rapalive_rapadead': "",
                               'IgM_pbslive_pbsdead': "",
                               'IgG_live_dead': "",
                               'IgG_rapalive_rapadead': "",
                               'IgG_pbslive_pbsdead': ""}
            cmap = [[255 / 255, 193 / 255, 7 / 255], "#CC2D3A", "#4AA518"]  # color,green,purple
            ax = PlotTools.plot_heatmap(matrix.replace(value_to_int).loc[the_indexs,:].rename(index_to_change),
                                        cmap=cmap, figsize=(60, 15),
                                        title="", title_fontsize=40, ax=axs[counter_x, counter_y],
                                        font_scale=4, snsStyle='ticks', xRotation=90,
                                        yRotation=0,
                                        xlabel='', ylabel='', colormap_label='',
                                        vmin=1, vmax=3, supress_ticks=True,
                                        annotate_text=False, annotate_fontsize=8,
                                        annotation_format=".2f",
                                        mask=None, colorbar_ticks=None,
                                        hide_colorbar=True,
                                        xy_labels_fontsize=45,
                                        grid_linewidths=0, grid_linecolor='white')
            for xc in range(0, len(matrix.columns)+1, 5):
                ax.axvline(x=xc, linewidth=0.3, color='b')
            for xc in range(0, len(matrix.index)+1, 1):
                if xc == 3:
                    ax.axhline(y=xc, linewidth=0.6, color='b')
                else:
                    ax.axhline(y=xc, linewidth=0.3, color='b')
        plt.savefig(os.path.join(lab_project_path, "rapa", "seq", "significant", "figs",
                                 "all_proteins_with_cmap" + ending_for_save +align_text+ ".png"),dpi=600)
        #plt.show()
        plt.close("all")

    if morin_2 ==True and is_4_groups == "rearengment_morin" and is_align==True:
        HA_df = pd.DataFrame()
        NA_df = pd.DataFrame()
        for protein_name, matrix in protein_dict_simple.items():
            index_to_change = {'IgM_live_dead': protein_name+" - IgM - All",
                               'IgM_rapalive_rapadead': protein_name+" - IgM - RAP",
                               'IgM_pbslive_pbsdead': protein_name+" - IgM - PBS",
                               'IgG_live_dead': protein_name+" - IgG - All",
                               'IgG_rapalive_rapadead': protein_name+" - IgG - RAP",
                               'IgG_pbslive_pbsdead': protein_name+" - IgG - PBS"}
            matrix=matrix.rename(index=index_to_change )
            #matrix = main_dict['rearengment_morin']['protein_dict_simple']['H1']
            columns_of_protein = protein_dict_simple[protein_name].columns
            #spaser1 = pd.DataFrame(np.nan, [protein_name + " - IgM"], columns_of_protein)
            #spaser = pd.DataFrame(np.nan, [protein_name + " - IgG"], columns_of_protein)
            HANA_df = matrix.copy()#pd.concat([spaser,spaser1,matrix])
            HANA_df["protein"]=protein_name
            HANA_df["type"] = HANA_df.index.map(lambda x:x.split(" - ")[1])
            if protein_name.startswith("H"):
                HA_df = pd.concat([HA_df,HANA_df])
            else:
                NA_df = pd.concat([NA_df,HANA_df])
        index_order_HA =[ 'H3 - IgG - All', 'H3 - IgG - RAP', 'H3 - IgG - PBS',
                       'H1 - IgG - All', 'H1 - IgG - RAP', 'H1 - IgG - PBS',
                      'H5 - IgG - All', 'H5 - IgG - RAP', 'H5 - IgG - PBS',
                       'H3 - IgM - All', 'H3 - IgM - RAP', 'H3 - IgM - PBS',
                       'H1 - IgM - All', 'H1 - IgM - RAP', 'H1 - IgM - PBS',
                       'H5 - IgM - All', 'H5 - IgM - RAP', 'H5 - IgM - PBS']
        index_order_NA = [ 'N2 - IgG - All', 'N2 - IgG - RAP',
         'N2 - IgG - PBS', 'Cal_N1 - IgG - All', 'Cal_N1 - IgG - RAP',
         'Cal_N1 - IgG - PBS', 'Vie_N1 - IgG - All',
         'Vie_N1 - IgG - RAP', 'Vie_N1 - IgG - PBS', 'N2 - IgM - All', 'N2 - IgM - RAP',
         'N2 - IgM - PBS',
         'Cal_N1 - IgM - All', 'Cal_N1 - IgM - RAP', 'Cal_N1 - IgM - PBS',
          'Vie_N1 - IgM - All', 'Vie_N1 - IgM - RAP',
         'Vie_N1 - IgM - PBS']
        def fix_index_name(x):
            splited = x.split(" - ")
            if len(splited) == 3:
                return splited[2]
            return splited[0]
        HA_df =HA_df.loc[index_order_HA,:]
        NA_df =NA_df.loc[index_order_NA,:]

        HA_df.rename(index=fix_index_name,inplace=True)
        NA_df.rename(index=fix_index_name,inplace=True)
        HA_df.to_csv(os.path.join(lab_project_path,"rapa","HA_df.csv"))
        HA_df.to_csv(os.path.join(lab_project_path,"rapa","NA_df.csv"))
        def make_IGMdeath_IGGlive(HNA_df):
            HNA_df_last_num =HNA_df.loc[:,:"protein"].columns[-2]
            HNA_df = HNA_df.copy()
            deadth_life = pd.DataFrame()
            for protein in HNA_df["protein"].unique():
                small_df = HNA_df.loc[HNA_df["protein"]==protein].copy()
                for compresionm in HNA_df.index.unique():
                    all_all = small_df.loc[compresionm,:HNA_df_last_num].apply(lambda x: x[0]== 2 and x[1]==3)
                    deadth_life.loc[protein,compresionm+"_"+compresionm] = ",".join(all_all.loc[all_all==True].index.astype(str))
                #small_df["ched"] = small_df.index +"_"+ small_df.loc[:,"type"]
                small_df["ched"] = small_df.index +"_"+ small_df["type"]
                all_all = small_df.loc[small_df["ched"].isin(["RAP_IgM","PBS_IgG"]), :HNA_df_last_num].apply(lambda x: x[0] == 2 and x[1] == 3)
                deadth_life.loc[protein,"RAP_IgM" + "__" + "PBS_IgG"] = ",".join(all_all.loc[all_all == True].index.astype(str))
            return deadth_life
        HA_deadth_life = make_IGMdeath_IGGlive(HA_df.replace(value_to_int))
        HA_deadth_life.to_csv(os.path.join(lab_project_path,"rapa","HA_death_life_df.csv"))
        NA_deadth_life = make_IGMdeath_IGGlive(NA_df.replace(value_to_int))
        NA_deadth_life.to_csv(os.path.join(lab_project_path,"rapa","NA_death_life_df.csv"))
        def make_it_unalign(HNA_deadth_life):
            HNA_deadth_life_copy = HNA_deadth_life.copy()
            for protein in HNA_deadth_life_copy.index:
                dict_of_maps_of_coloumns_reverse = {v: k for k, v in dict_of_maps_of_coloumns[protein].items()}
                for column in HNA_deadth_life_copy.columns:
                    if not pd.isna(HNA_deadth_life_copy.loc[protein,column]) and HNA_deadth_life_copy.loc[protein,column]!="":#Check if there is a value
                        HNA_deadth_life_copy.loc[protein,column] = ",".join([str(dict_of_maps_of_coloumns_reverse[int(x)]) for x in HNA_deadth_life_copy.loc[protein,column].split(",")])
            return HNA_deadth_life_copy
        unalign_HA_deadth_life = make_it_unalign(HA_deadth_life)
        unalign_HA_deadth_life.to_csv(os.path.join(lab_project_path,"rapa","HA_death_life_df_unalign.csv"))
        unalign_NA_deadth_life = make_it_unalign(NA_deadth_life)
        unalign_NA_deadth_life.to_csv(os.path.join(lab_project_path,"rapa","NA_death_life_df_unalign.csv"))

        HA_df =HA_df.drop(["protein","type"],axis=1)
        NA_df =NA_df.drop(["protein","type"],axis=1)
        value_to_int = {'ambivalent': 1,
                        'rapadead': 2,
                        'rapalive': 3,
                        'dead': 2,
                        'live': 3,
                        'pbsdead': 2,
                        'pbslive': 3}
        # cmap = [(160 / 255, 126 / 255, 91 / 255),
        #         (91 / 255, 160 / 255, 154 / 255),
        #         (108 / 255, 88 / 255, 150 / 255)]  # color,green,purple
        # cmap = [[0.352941176, 0.68627451, 0.37254902],
        #         [0.058823529, 0.352941176, 0.137254902],
        #         [0.607843137, 0.431372549, 0.666666667]
        # ]
        # cmap = {"ambivalent": [255 / 255, 193 / 255, 7 / 255],
        #         "dead": [216 / 218, 218 / 255, 235 / 255],
        #         "live": [128 / 255, 115 / 255, 172 / 255],
        #         "diff": [228 / 255, 130 / 255, 217 / 255]
        #               }  # morin
        # cmap = {"ambivalent":[255 / 255, 193 / 255, 7 / 255],
        #         "dead": [179 / 255, 88 / 255, 6 / 255],
        #         "live": [84 / 255, 39 / 255, 136 / 255] ,
        #         "diff": [228 / 255, 130 / 255, 217 / 255]
        #               }#eilay after morin

        cmap = [colors_map["ambivalent"],colors_map["dead"],colors_map["live"]]
        #cmap = [[255 / 255, 193 / 255, 7 / 255], "#CC2D3A", "#4AA518"]  # color,green,purple

        endresult_per_position_file_db = pd.read_csv(endresult_per_position_file_path, sep="\t", index_col=0)
        dict_of_maps_of_coloumns = main_dict['live_dead']['dict_of_maps_of_coloumns']
        def make_endresult_per_pos_align(endresult_per_position_file_db,dict_of_maps_of_coloumns):
            align_sasa_conser = pd.DataFrame(columns=["pep_count"])
            for proteinn_name,mapper in dict_of_maps_of_coloumns.items():
                relevant_cols =[col for col in endresult_per_position_file_db.columns if col.endswith(proteinn_name)]
                relevant_cols.append("pep_count")
                protein_db=endresult_per_position_file_db.loc[:,relevant_cols]
                protein_db["pep_count"] = protein_db["pep_count"].map(mapper)
                align_sasa_conser = align_sasa_conser.merge(protein_db,on="pep_count",how="outer")
            align_sasa_conser = align_sasa_conser.dropna(axis=0,how="all")
            align_sasa_conser = align_sasa_conser.sort_values(by="pep_count")
            align_sasa_conser.index =  align_sasa_conser["pep_count"]
            cols_cons = [col for col in align_sasa_conser.columns if col.startswith("conservation_nonbinery")]
            cols_sasa = [col for col in align_sasa_conser.columns if col.startswith("sasa_")]
            fig_HA, axes = plt.subplots(1,2, sharey=True, sharex=False, figsize=(20, 50))  # subplots(row, columns)
            PlotTools.plot_heatmap(align_sasa_conser.loc[:,cols_cons],
                                   cmap="vlag",
                                   title="conservation", title_fontsize=40,ax=axes[0],
                                   font_scale=4, snsStyle='ticks', xRotation=90,
                                   yRotation=0,
                                   xlabel='', ylabel='', colormap_label='',
                                   supress_ticks=True,
                                   annotate_text=False, annotate_fontsize=8,
                                   annotation_format=".2f",
                                   mask=None, colorbar_ticks=None,
                                   hide_colorbar=False,
                                   xy_labels_fontsize=45,
                                   grid_linewidths=0, grid_linecolor='white')
            PlotTools.plot_heatmap(align_sasa_conser.loc[:,cols_sasa],
                                   cmap="Blues",
                                   title="sasa", title_fontsize=40,ax=axes[1],
                                   font_scale=4, snsStyle='ticks', xRotation=90,
                                   yRotation=0,
                                   xlabel='', ylabel='', colormap_label='',
                                    supress_ticks=True,
                                   annotate_text=False, annotate_fontsize=8,
                                   annotation_format=".2f",
                                   mask=None, colorbar_ticks=None,
                                   hide_colorbar=False,
                                   xy_labels_fontsize=45,
                                   grid_linewidths=0, grid_linecolor='white')
            fig_HA.suptitle("sasa and conservation score comparison", fontsize=60)
            plt.tight_layout()
            plt.savefig(os.path.join(lab_project_path, "rapa", "seq", "significant", "figs",
                                     "sasa and conservation score comparison".replace(" ","_") + ending_for_save + align_text + ".png"),dpi=600)
            #plt.show()
            PlotTools.plot_heatmap(align_sasa_conser.loc[:,cols_sasa].corr("spearman"),
                                   cmap="Blues",figsize=(20,20),
                                   title="spearman correlation between conservations score", title_fontsize=40,
                                   font_scale=4, snsStyle='ticks', xRotation=90,
                                   yRotation=0,
                                   xlabel='', ylabel='', colormap_label='spearman',
                                    supress_ticks=True,
                                   annotate_text=True, annotate_fontsize=25,
                                   annotation_format=".2f",
                                   mask=None, colorbar_ticks=None,
                                   hide_colorbar=False,
                                   xy_labels_fontsize=45,
                                   grid_linewidths=0.1, grid_linecolor='white')
            plt.tight_layout()
            plt.savefig(os.path.join(lab_project_path, "rapa", "seq", "significant", "figs",
                                     "spearman correlation between conservations score".replace(" ",
                                                                                      "_") + ending_for_save + align_text + ".png"),dpi=600)
            #plt.show()

            PlotTools.plot_heatmap(align_sasa_conser.loc[:,cols_cons].corr("spearman"),
                                   cmap="Blues",figsize=(20,20),
                                   title="spearman correlation between sasa score", title_fontsize=40,
                                   font_scale=4, snsStyle='ticks', xRotation=90,
                                   yRotation=0,
                                   xlabel='', ylabel='', colormap_label='spearman',
                                    supress_ticks=True,
                                   annotate_text=True, annotate_fontsize=25,
                                   annotation_format=".2f",
                                   mask=None, colorbar_ticks=None,
                                   hide_colorbar=False,
                                   xy_labels_fontsize=45,
                                   grid_linewidths=0.1, grid_linecolor='white')
            plt.tight_layout()
            plt.savefig(os.path.join(lab_project_path, "rapa", "seq", "significant", "figs",
                                     "spearman correlation between sasa score".replace(" ",
                                                                                      "_") + ending_for_save + align_text + ".png"),dpi=600)
            #plt.show()
            plt.close("all")
            means_sasa_conser=pd.DataFrame([align_sasa_conser.loc[:,[col for col in cols_cons if col.split("_")[-1].startswith("H")]].mean(axis = 1),
                                           align_sasa_conser.loc[:,[col for col in cols_cons if not col.split("_")[-1].startswith("H")]].mean(axis = 1),
                                           align_sasa_conser.loc[:,[col for col in cols_sasa if col.split("_")[-1].startswith("H")]].mean(axis = 1),
                                           align_sasa_conser.loc[:,[col for col in cols_sasa if not col.split("_")[-1].startswith("H")]].mean(axis = 1)],
                                           index=["HA_conservation","NA_conservation","HA_sasa","NA_sasa"])
            return align_sasa_conser,means_sasa_conser
        align_sasa_conser,means_sasa_conser=make_endresult_per_pos_align(endresult_per_position_file_db,dict_of_maps_of_coloumns)

        #fig_HA, (axs_HA,down,downest) = plt.subplots(3,1,gridspec_kw={'height_ratios':[0.9,0.05,0.05]}, sharey=False,sharex=True, figsize=(40, 30))  # subplots(row, columns)
        fig_HA, axs_HA = plt.subplots(1,1, figsize=(40, 30))  # subplots(row, columns)
        HA_df.columns = HA_df.columns.astype(float).astype(int)
        HA_df.replace(value_to_int).to_csv(os.path.join(lab_project_path,"rapa","HA_df_ploted.csv"))
        PlotTools.plot_heatmap(HA_df.replace(value_to_int),
                               cmap=cmap,
                               title="", title_fontsize=40,ax=axs_HA,
                               font_scale=4, snsStyle='ticks', xRotation=90,
                               yRotation=0,
                               xlabel='', ylabel='', colormap_label='',
                               xticklabels=20,
                               vmin=1, vmax=3, supress_ticks=False,
                               annotate_text=False, annotate_fontsize=8,
                               annotation_format=".2f",
                               mask=None, colorbar_ticks=None,
                               hide_colorbar=True,
                               xy_labels_fontsize=45,
                               grid_linewidths=0, grid_linecolor='white')
        for xc in range(0, len(HA_df.columns) + 1, 5):
            axs_HA.axvline(x=xc, linewidth=0.25, color='black')

        for xc in range(0, len(HA_df.index) + 1, 1):
            if xc==9:
                axs_HA.axhline(y=xc, linewidth=4, color='black')
            elif xc%3==0 and xc!=0:
                axs_HA.axhline(y=xc, linewidth=1, color='black')
            else:
                axs_HA.axhline(y=xc, linewidth=0.25, color='black')
        axs_HA.text(x=-70 ,y =1.5 ,s="H3")
        axs_HA.text(x= -70,y =4.5 ,s="H1")
        axs_HA.text(x= -70,y =7.5 ,s="H5")
        axs_HA.text(x= -70,y =10.5 ,s="H3")
        axs_HA.text(x= -70,y = 13.5,s="H1")
        axs_HA.text(x= -70,y = 16.5,s="H5")
        axs_HA.text(x= -120,y = 4.5,s="IgG")
        axs_HA.text(x= -120,y = 13.5,s="IgM")

        #axs_HA.set_yticklabels(labels=axs_HA.get_yticklabels(),ha="center")
        # plt.yticks(np.arange(len(axs_HA.get_yticklabels())) + 0.5, axs_HA.get_yticklabels(),
        #            rotation=0, ha="center")
        axs_HA.vlines(x=-3, ymin=0, ymax=10, linewidth=4, color='pink')
        plt.legend(loc="upper left")
        x, y = np.array([[-40, -40], [0.3,2.7]])
        line = Line2D(x, y, lw=5., color='black', alpha=0.7)#from matplotlib.lines import Line2D
        line.set_clip_on(False)
        x, y = np.array([[-40, -40], [3.3,5.7]])
        axs_HA.add_line(line)
        line = Line2D(x, y, lw=5., color='black', alpha=0.7)#from matplotlib.lines import Line2D
        line.set_clip_on(False)
        axs_HA.add_line(line)
        x, y = np.array([[-40, -40], [6.3,8.7]])
        line = Line2D(x, y, lw=5., color='black', alpha=0.7)#from matplotlib.lines import Line2D
        line.set_clip_on(False)
        axs_HA.add_line(line)
        x, y = np.array([[-40, -40], [9.3,11.7]])
        line = Line2D(x, y, lw=5., color='black', alpha=0.7)#from matplotlib.lines import Line2D
        line.set_clip_on(False)
        axs_HA.add_line(line)
        x, y = np.array([[-40, -40], [12.3,14.7]])
        line = Line2D(x, y, lw=5., color='black', alpha=0.7)#from matplotlib.lines import Line2D
        line.set_clip_on(False)
        axs_HA.add_line(line)
        x, y = np.array([[-40, -40], [15.3,17.7]])
        line = Line2D(x, y, lw=5., color='black', alpha=0.7)#from matplotlib.lines import Line2D
        line.set_clip_on(False)
        axs_HA.add_line(line)
        x, y = np.array([[-80, -80], [0.3,8.7]])
        line = Line2D(x, y, lw=5., color='black', alpha=0.7)#from matplotlib.lines import Line2D
        line.set_clip_on(False)
        axs_HA.add_line(line)
        x, y = np.array([[-80, -80], [9.3,17.7]])
        line = Line2D(x, y, lw=5., color='black', alpha=0.7)#from matplotlib.lines import Line2D
        line.set_clip_on(False)
        axs_HA.add_line(line)
        axs_HA.legend(handles=[Patch(label='ambivalent', color=cmap[0]),
                            Patch(label='dead', color=cmap[1]),
                            Patch(label='life', color=cmap[2])],
                   loc='upper left', bbox_to_anchor=(1.02, 1), frameon=False, title='Significant Abs response')

        #
        # PlotTools.plot_heatmap(means_sasa_conser.loc[["HA_conservation"],:],
        #                        cmap="vlag",
        #                        title="", title_fontsize=40,ax=down,
        #                        font_scale=4, snsStyle='ticks', xRotation=90,
        #                        yRotation=0,center=0,
        #                        xlabel='', ylabel='', colormap_label='',
        #                        supress_ticks=True,
        #                        annotate_text=False, annotate_fontsize=8,
        #                        annotation_format=".2f",
        #                        mask=None, colorbar_ticks=None,
        #                        hide_colorbar=True,
        #                        xy_labels_fontsize=45,
        #                        grid_linewidths=0, grid_linecolor='white')
        #
        # PlotTools.plot_heatmap(means_sasa_conser.loc[["HA_sasa"],:],
        #                        cmap="Blues",
        #                        title="", title_fontsize=40,ax=downest,
        #                        font_scale=4, snsStyle='ticks', xRotation=90,
        #                        yRotation=0,
        #                        xlabel='', ylabel='', colormap_label='',
        #                        supress_ticks=True,
        #                        annotate_text=False, annotate_fontsize=8,
        #                        annotation_format=".2f",
        #                        mask=None, colorbar_ticks=None,
        #                        hide_colorbar=True,
        #                        xy_labels_fontsize=45,
        #                        grid_linewidths=0, grid_linecolor='white')
        #
        # downest.remove()
        # down.remove()
        plt.tight_layout()
        plt.savefig(os.path.join(lab_project_path, "rapa", "seq", "significant", "figs",
                                 "all_proteins_with_cmap_morin2_HA"+ ".png"),dpi=600)
        #plt.show()


        #fig_NA, (axs_NA,down1,downest1) = plt.subplots(3,1,gridspec_kw={'height_ratios':[0.9,0.05,0.05]}, sharey=False,sharex=True, figsize=(40, 30))  # subplots(row, columns)
        fig_NA, axs_NA = plt.subplots(1,1, figsize=(40, 30))  # subplots(row, columns)
        NA_df.columns=NA_df.columns.astype(float).astype(int)
        NA_df.replace(value_to_int).to_csv(os.path.join(lab_project_path,"rapa","NA_df_ploted.csv"))
        PlotTools.plot_heatmap(NA_df.replace(value_to_int),
                               cmap=cmap,
                               title="", title_fontsize=40,ax=axs_NA,
                               font_scale=4, snsStyle='ticks', xRotation=90,
                               yRotation=0,
                               xlabel='', ylabel='', colormap_label='',
                               vmin=1, vmax=3, supress_ticks=False,
                               xticklabels=20,
                               annotate_text=False, annotate_fontsize=8,
                               annotation_format=".2f",
                               mask=None, colorbar_ticks=None,
                               hide_colorbar=True,
                               xy_labels_fontsize=45,
                               grid_linewidths=0, grid_linecolor='white')
        for xc in range(0, len(NA_df.columns) + 1, 5):
            axs_NA.axvline(x=xc, linewidth=0.25, color='black')

        for xc in range(0, len(NA_df.index) + 1, 1):
            if xc==9:
                axs_NA.axhline(y=xc, linewidth=4, color='black')
            elif xc%3==0 and xc!=0:
                axs_NA.axhline(y=xc, linewidth=1, color='black')
            else:
                axs_NA.axhline(y=xc, linewidth=0.25, color='black')
        axs_NA.text(x=-105 ,y =1.5 ,s="   N2")
        axs_NA.text(x= -105,y =4.5 ,s="Cal_N1")
        axs_NA.text(x= -105,y =7.5 ,s="Vie_N1")
        axs_NA.text(x= -105,y =10.5 ,s="   N2")
        axs_NA.text(x= -105,y = 13.5,s="Cal N1")
        axs_NA.text(x= -105,y = 16.5,s="Vie N1")
        axs_NA.text(x= -150,y = 4.5,s="IgG")
        axs_NA.text(x= -150,y = 13.5,s="IgM")

        axs_NA.vlines(x=-3, ymin=0, ymax=10, linewidth=4, color='pink')
        #axs_NA.set_yticklabels(labels=axs_NA.get_yticklabels(),ha="center")


        x, y = np.array([[-40, -40], [0.3,2.7]])
        line = Line2D(x, y, lw=5., color='black', alpha=0.7)#from matplotlib.lines import Line2D
        line.set_clip_on(False)
        x, y = np.array([[-40, -40], [3.3,5.7]])
        axs_NA.add_line(line)
        line = Line2D(x, y, lw=5., color='black', alpha=0.7)#from matplotlib.lines import Line2D
        line.set_clip_on(False)
        axs_NA.add_line(line)
        x, y = np.array([[-40, -40], [6.3,8.7]])
        line = Line2D(x, y, lw=5., color='black', alpha=0.7)#from matplotlib.lines import Line2D
        line.set_clip_on(False)
        axs_NA.add_line(line)
        x, y = np.array([[-40, -40], [9.3,11.7]])
        line = Line2D(x, y, lw=5., color='black', alpha=0.7)#from matplotlib.lines import Line2D
        line.set_clip_on(False)
        axs_NA.add_line(line)
        x, y = np.array([[-40, -40], [12.3,14.7]])
        line = Line2D(x, y, lw=5., color='black', alpha=0.7)#from matplotlib.lines import Line2D
        line.set_clip_on(False)
        axs_NA.add_line(line)
        x, y = np.array([[-40, -40], [15.3,17.7]])
        line = Line2D(x, y, lw=5., color='black', alpha=0.7)#from matplotlib.lines import Line2D
        line.set_clip_on(False)
        axs_NA.add_line(line)
        x, y = np.array([[-110, -110], [0.3,8.7]])
        line = Line2D(x, y, lw=5., color='black', alpha=0.7)#from matplotlib.lines import Line2D
        line.set_clip_on(False)
        axs_NA.add_line(line)
        x, y = np.array([[-110, -110], [9.3,17.7]])
        line = Line2D(x, y, lw=5., color='black', alpha=0.7)#from matplotlib.lines import Line2D
        line.set_clip_on(False)
        axs_NA.add_line(line)
        axs_NA.legend(handles=[Patch(label='ambivalent', color=cmap[0]),
                            Patch(label='dead', color=cmap[1]),
                            Patch(label='life', color=cmap[2])],
                   loc='upper left', bbox_to_anchor=(1.02, 1), frameon=False, title='Significant Abs response')

        #
        # PlotTools.plot_heatmap(means_sasa_conser.loc[["NA_conservation"],:481],
        #                        cmap="vlag",
        #                        title="", title_fontsize=40,ax=down1,
        #                        font_scale=4, snsStyle='ticks', xRotation=90,
        #                        yRotation=0,
        #                        xlabel='', ylabel='', colormap_label='',
        #                        supress_ticks=True,
        #                        annotate_text=False, annotate_fontsize=8,
        #                        annotation_format=".2f",
        #                        mask=None, colorbar_ticks=None,
        #                        hide_colorbar=True,
        #                        xy_labels_fontsize=45,
        #                        grid_linewidths=0, grid_linecolor='white')
        #
        # PlotTools.plot_heatmap(means_sasa_conser.loc[["NA_sasa"],:481],
        #                        cmap="Blues",
        #                        title="", title_fontsize=40,ax=downest1,
        #                        font_scale=4, snsStyle='ticks', xRotation=90,
        #                        yRotation=0,
        #                        xlabel='', ylabel='', colormap_label='',
        #                        supress_ticks=True,
        #                        annotate_text=False, annotate_fontsize=8,
        #                        annotation_format=".2f",
        #                        mask=None, colorbar_ticks=None,
        #                        hide_colorbar=True,
        #                        xy_labels_fontsize=45,
        #                        grid_linewidths=0, grid_linecolor='white')
        # downest1.remove()
        # down1.remove()
        plt.tight_layout()
        plt.savefig(os.path.join(lab_project_path, "rapa", "seq", "significant", "figs",
                                 "all_proteins_with_cmap_morin2_NA"+ ".png"),dpi=600)
        #plt.show()

        if morin3 == True:
            NA_df_ploted = NA_df.replace(value_to_int)
            HA_df_ploted = HA_df.replace(value_to_int)

            # NA_df_ploted = pd.read_csv(os.path.join(lab_project_path, "rapa", "NA_df_ploted.csv"),index_col=0)
            # HA_df_ploted = pd.read_csv(os.path.join(lab_project_path, "rapa", "HA_df_ploted.csv"),index_col=0)
            # colors_map = {"ambivalent": [255 / 255, 193 / 255, 7 / 255],
            #               "live": [128 / 255, 115 / 255, 172 / 255],
            #               "dead": [216 / 255, 218 / 255, 235 / 255],
            #               "diff": [0 / 255, 141 / 255, 255 / 255]}
            # colors_map = {"ambivalent": [255 / 255, 255 / 255, 255 / 255],
            #               "live": [128 / 255, 115 / 255, 172 / 255],
            #               "dead": [216 / 255, 218 / 255, 235 / 255],
            #               "diff": [0 / 255, 141 / 255, 255 / 255]}
            # cmap = [colors_map["ambivalent"], colors_map["dead"], colors_map["live"]]
            colors_map_rapa_pbs ={"ambivalent": [255 / 255, 255 / 255, 255 / 255],
                              'PBS_Live': [33 / 255, 102 / 255, 172 / 255],
             'PBS_Dead':  [209 / 255, 229 / 255, 240 / 255],
             'Rapa_Live': [178 / 255, 24 / 255, 43 / 255],
             'Rapa_Dead': [253 / 255, 219 / 255, 199 / 255]}
            cmap2 = [colors_map_rapa_pbs["ambivalent"], colors_map_rapa_pbs["Rapa_Dead"], colors_map_rapa_pbs["Rapa_Live"],
                        colors_map_rapa_pbs["ambivalent"],colors_map_rapa_pbs["PBS_Dead"], colors_map_rapa_pbs["PBS_Live"]]


            def morin3(HA_df_ploted,HA=True):
                #HA_df_ploted = HA_df_ploted.replace(value_to_int)
                fig_6ab, (ax_6a,ax_6b) = plt.subplots(2, 1, figsize=(40, 15),sharey=True)  # subplots(row, columns)
                if HA==True:
                    indexx= ["H3","H1","H5"]
                    HA_s ="HA"
                else:
                    indexx= ["N2","(H1)N1","(H5)N1"]
                    HA_s ="NA"

                df_6a = HA_df_ploted.iloc[[0,3,6],:]
                df_6a.index = indexx
                PlotTools.plot_heatmap(df_6a,
                                       cmap=cmap,
                                       title="IgG", title_fontsize=70, ax=ax_6a,
                                       font_scale=5, snsStyle='ticks', xRotation=90,
                                       yRotation=0,
                                       xlabel='', ylabel='', colormap_label='',
                                       vmin=1, vmax=3, supress_ticks=False,
                                       xticklabels=20,
                                       annotate_text=False, annotate_fontsize=8,
                                       annotation_format=".2f",
                                       mask=None, colorbar_ticks=None,
                                       hide_colorbar=True,
                                       xy_labels_fontsize=70,
                                       grid_linewidths=0, grid_linecolor='black')
                for xc in range(0, len(HA_df_ploted.columns) + 1, 5):
                    ax_6a.axvline(x=xc, linewidth=0.25, color='black')

                df_6b = HA_df_ploted.iloc[[9,12,15],:]
                df_6b.index = indexx
                PlotTools.plot_heatmap(df_6b,
                                       cmap=cmap,
                                       title="IgM", title_fontsize=70, ax=ax_6b,
                                       font_scale=5, snsStyle='ticks', xRotation=90,
                                       yRotation=0,
                                       xlabel='', ylabel='', colormap_label='',
                                       vmin=1, vmax=3, supress_ticks=False,
                                       xticklabels=20,
                                       annotate_text=False, annotate_fontsize=8,
                                       annotation_format=".2f",
                                       mask=None, colorbar_ticks=None,
                                       hide_colorbar=True,
                                       xy_labels_fontsize=70,
                                       grid_linewidths=0, grid_linecolor='black')
                for xc in range(0, len(HA_df_ploted.columns) + 1, 5):
                    ax_6b.axvline(x=xc, linewidth=0.25, color='black')

                plt.savefig(os.path.join(lab_project_path, "rapa", "seq", "significant", "figs",
                                         "6ab_morin3_{}".format(HA_s) + ".png"),dpi=600)
                HA_df_ploted.loc["PBS",:] = HA_df_ploted.loc["PBS",:]+3
                fig_6c, (ax_6c1, ax_6c2,ax_6c3) = plt.subplots(3, 1, figsize=(40, 18), sharey=True)  # subplots(row, columns)
                PlotTools.plot_heatmap(HA_df_ploted.iloc[[2,1], :],
                                       cmap=cmap2,
                                       title="H3", title_fontsize=70, ax=ax_6c1,
                                       font_scale=5, snsStyle='ticks', xRotation=90,
                                       yRotation=0,
                                       xlabel='', ylabel='', colormap_label='',
                                       vmin=1, vmax=6, supress_ticks=False,
                                       xticklabels=20,
                                       annotate_text=False, annotate_fontsize=8,
                                       annotation_format=".2f",
                                       mask=None, colorbar_ticks=None,
                                       hide_colorbar=True,
                                       xy_labels_fontsize=70,
                                       grid_linewidths=0, grid_linecolor='black')

                PlotTools.plot_heatmap(HA_df_ploted.iloc[[5,4], :],
                                       cmap=cmap2,
                                       title="H1", title_fontsize=70, ax=ax_6c2,
                                       font_scale=5, snsStyle='ticks', xRotation=90,
                                       yRotation=0,
                                       xlabel='', ylabel='', colormap_label='',
                                       vmin=1, vmax=6, supress_ticks=False,
                                       xticklabels=20,
                                       annotate_text=False, annotate_fontsize=8,
                                       annotation_format=".2f",
                                       mask=None, colorbar_ticks=None,
                                       hide_colorbar=True,
                                       xy_labels_fontsize=70,
                                       grid_linewidths=0, grid_linecolor='black')

                PlotTools.plot_heatmap(HA_df_ploted.iloc[[8,7], :],
                                       cmap=cmap2,
                                       title="H5", title_fontsize=70, ax=ax_6c3,
                                       font_scale=5, snsStyle='ticks', xRotation=90,
                                       yRotation=0,
                                       xlabel='', ylabel='', colormap_label='',
                                       vmin=1, vmax=6, supress_ticks=False,
                                       xticklabels=20,
                                       annotate_text=False, annotate_fontsize=8,
                                       annotation_format=".2f",
                                       mask=None, colorbar_ticks=None,
                                       hide_colorbar=True,
                                       xy_labels_fontsize=70,
                                       grid_linewidths=0, grid_linecolor='black')
                for ax in [ax_6c1,ax_6c2,ax_6c3]:
                    for xc in range(0, len(HA_df_ploted.columns) + 1, 5):
                        ax.axvline(x=xc, linewidth=0.25, color='black')

                plt.savefig(os.path.join(lab_project_path, "rapa", "seq", "significant", "figs",
                                     "6c_morin3_{}".format(HA_s) + ".png"), dpi=600)

                fig_6d, (ax_6d1, ax_6d2,ax_6d3) = plt.subplots(3, 1, figsize=(40, 18), sharey=True)  # subplots(row, columns)
                PlotTools.plot_heatmap(HA_df_ploted.iloc[[11,10], :],
                                       cmap=cmap2,
                                       title="H3", title_fontsize=70, ax=ax_6d1,
                                       font_scale=5, snsStyle='ticks', xRotation=90,
                                       yRotation=0,
                                       xlabel='', ylabel='', colormap_label='',
                                       vmin=1, vmax=6, supress_ticks=False,
                                       xticklabels=20,
                                       annotate_text=False, annotate_fontsize=8,
                                       annotation_format=".2f",
                                       mask=None, colorbar_ticks=None,
                                       hide_colorbar=True,
                                       xy_labels_fontsize=70,
                                       grid_linewidths=0, grid_linecolor='black')

                PlotTools.plot_heatmap(HA_df_ploted.iloc[[14,13], :],
                                       cmap=cmap2,
                                       title="H1", title_fontsize=70, ax=ax_6d2,
                                       font_scale=5, snsStyle='ticks', xRotation=90,
                                       yRotation=0,
                                       xlabel='', ylabel='', colormap_label='',
                                       vmin=1, vmax=6, supress_ticks=False,
                                       xticklabels=20,
                                       annotate_text=False, annotate_fontsize=8,
                                       annotation_format=".2f",
                                       mask=None, colorbar_ticks=None,
                                       hide_colorbar=True,
                                       xy_labels_fontsize=70,
                                       grid_linewidths=0, grid_linecolor='black')

                PlotTools.plot_heatmap(HA_df_ploted.iloc[[17,16], :],
                                       cmap=cmap2,
                                       title="H5", title_fontsize=70, ax=ax_6d3,
                                       font_scale=5, snsStyle='ticks', xRotation=90,
                                       yRotation=0,
                                       xlabel='', ylabel='', colormap_label='',
                                       vmin=1, vmax=6, supress_ticks=False,
                                       xticklabels=20,
                                       annotate_text=False, annotate_fontsize=8,
                                       annotation_format=".2f",
                                       mask=None, colorbar_ticks=None,
                                       hide_colorbar=True,
                                       xy_labels_fontsize=70,
                                       grid_linewidths=0, grid_linecolor='black')
                for ax in [ax_6d1,ax_6d2,ax_6d3]:
                    for xc in range(0, len(HA_df_ploted.columns) + 1, 5):
                        ax.axvline(x=xc, linewidth=0.25, color='black')

                plt.savefig(os.path.join(lab_project_path, "rapa", "seq", "significant", "figs",
                                     "6d_morin3_{}".format(HA_s) + ".png"), dpi=600)

                fig_6e, (ax_6e1, ax_6e2,ax_6e3) = plt.subplots(3, 1, figsize=(40, 18), sharey=True)  # subplots(row, columns)
                PlotTools.plot_heatmap(HA_df_ploted.iloc[[2,10], :].rename(index={"PBS":"IgG-PBS", "RAP":"IgM-RAP"}),
                                       cmap=cmap2,
                                       title="H3", title_fontsize=70, ax=ax_6e1,
                                       font_scale=5, snsStyle='ticks', xRotation=90,
                                       yRotation=0,
                                       xlabel='', ylabel='', colormap_label='',
                                       vmin=1, vmax=6, supress_ticks=False,
                                       xticklabels=20,
                                       annotate_text=False, annotate_fontsize=8,
                                       annotation_format=".2f",
                                       mask=None, colorbar_ticks=None,
                                       hide_colorbar=True,
                                       xy_labels_fontsize=70,
                                       grid_linewidths=0, grid_linecolor='black')

                PlotTools.plot_heatmap(HA_df_ploted.iloc[[5,13], :].rename(index={"PBS":"IgG-PBS", "RAP":"IgM-RAP"}),
                                       cmap=cmap2,
                                       title="H1", title_fontsize=70, ax=ax_6e2,
                                       font_scale=5, snsStyle='ticks', xRotation=90,
                                       yRotation=0,
                                       xlabel='', ylabel='', colormap_label='',
                                       vmin=1, vmax=6, supress_ticks=False,
                                       xticklabels=20,
                                       annotate_text=False, annotate_fontsize=8,
                                       annotation_format=".2f",
                                       mask=None, colorbar_ticks=None,
                                       hide_colorbar=True,
                                       xy_labels_fontsize=70,
                                       grid_linewidths=0, grid_linecolor='black')

                PlotTools.plot_heatmap(HA_df_ploted.iloc[[8,16], :].rename(index={"PBS":"IgG-PBS", "RAP":"IgM-RAP"}),
                                       cmap=cmap2,
                                       title="H5", title_fontsize=70, ax=ax_6e3,
                                       font_scale=5, snsStyle='ticks', xRotation=90,
                                       yRotation=0,
                                       xlabel='', ylabel='', colormap_label='',
                                       vmin=1, vmax=6, supress_ticks=False,
                                       xticklabels=20,
                                       annotate_text=False, annotate_fontsize=8,
                                       annotation_format=".2f",
                                       mask=None, colorbar_ticks=None,
                                       hide_colorbar=True,
                                       xy_labels_fontsize=70,
                                       grid_linewidths=0, grid_linecolor='black')
                for ax in [ax_6e1,ax_6e2,ax_6e3]:
                    for xc in range(0, len(HA_df_ploted.columns) + 1, 5):
                        ax.axvline(x=xc, linewidth=0.25, color='black')
                plt.savefig(os.path.join(lab_project_path, "rapa", "seq", "significant", "figs",
                                     "6e_morin3_{}".format(HA_s) + ".png"), dpi=600)

            morin3(HA_df_ploted, HA=True)
            morin3(NA_df_ploted, HA=False)
        plt.close("all")



def rearengment_of_main_dict(main_dict,output_path=os.path.join(lab_project_path, "rapa", "seq", "significant", "figs"),endresult_per_peptide_file_path=os.path.join(lab_project_path,"rapa","endresulte_per_peptide.tsv"),pdb_folder_path = os.path.join(lab_project_path,"rapa","pdb","pdb_fix_insertions"),to_plot=False,colors_map = {"ambivalent": [255 / 255, 193 / 255, 7 / 255],
              "live": [128 / 255, 115 / 255, 172 / 255],
              "dead": [216 / 255, 218 / 255, 235 / 255],
              "diff": [0 / 255, 141 / 255, 255 / 255]}):
    name = "rearengment_morin"
    temp_dict = {}
    relevant_col = ["IgM_rapalive_rapadead","IgG_rapalive_rapadead","IgM_pbslive_pbsdead","IgG_pbslive_pbsdead"]
    temp_dict["significnt_db_result"] = main_dict["live_dead"]["significnt_db_result"].merge(main_dict["4_groups"]["significnt_db_result"].loc[:,relevant_col],left_index=True, right_index=True)
    temp_dict["significnt_db_result"].to_csv(os.path.join(output_path,"significnt_db_result{}.tsv".format(name)),sep="\t")
    temp_dict["significnt_db_result_simple"] = main_dict["live_dead"]["significnt_db_result_simple"].merge(main_dict["4_groups"]["significnt_db_result_simple"].loc[:,relevant_col],left_index=True, right_index=True)
    temp_dict["significnt_db_result_simple"].to_csv(os.path.join(output_path,"significnt_db_result_simple{}.tsv".format(name)),sep="\t")
    temp_dict["seq_of_proteins"] = main_dict["live_dead"]["seq_of_proteins"]
    temp_dict["lenght_of_proteins"] = main_dict["live_dead"]["lenght_of_proteins"]
    temp_dict["protein_dict"] = {protein:matrix_live_dead.T.merge(matrix.T.loc[:,relevant_col],left_index=True, right_index=True).T for protein,matrix_live_dead,matrix in zip(main_dict["live_dead"]["protein_dict"].keys(),main_dict["live_dead"]["protein_dict"].values(),main_dict["4_groups"]["protein_dict"].values()) }
    temp_dict["protein_dict_simple"] = {protein:matrix_live_dead.T.merge(matrix.T.loc[:,relevant_col],left_index=True, right_index=True).T for protein,matrix_live_dead,matrix in zip(main_dict["live_dead"]["protein_dict_simple"].keys(),main_dict["live_dead"]["protein_dict_simple"].values(),main_dict["4_groups"]["protein_dict_simple"].values()) }
    temp_dict["align_protein_dict"] = {protein:matrix_live_dead.T.merge(matrix.T.loc[:,relevant_col],left_index=True, right_index=True).T for protein,matrix_live_dead,matrix in zip(main_dict["live_dead"]["align_protein_dict"].keys(),main_dict["live_dead"]["align_protein_dict"].values(),main_dict["4_groups"]["align_protein_dict"].values()) }
    temp_dict["align_protein_dict_simple"] = {protein:matrix_live_dead.T.merge(matrix.T.loc[:,relevant_col],left_index=True, right_index=True).T for protein,matrix_live_dead,matrix in zip(main_dict["live_dead"]["align_protein_dict_simple"].keys(),main_dict["live_dead"]["align_protein_dict_simple"].values(),main_dict["4_groups"]["align_protein_dict_simple"].values()) }
    structrul_compare,structrul_compare_only_sig = make_structral_compare_per_position(temp_dict["significnt_db_result_simple"] ,to_plot=to_plot, is_4_groups=name,
                               endresult_per_peptide_file_path=endresult_per_peptide_file_path,
                               output_path=output_path)
    temp_dict["structrul_compare"] = structrul_compare
    temp_dict["structrul_compare_only_sig"] = structrul_compare_only_sig
    temp_dict["signature_df"] = make_signature_df_to_rearngment_morin(temp_dict["significnt_db_result_simple"])
    if to_plot == True:
        plots(temp_dict["protein_dict"], temp_dict["protein_dict_simple"], is_4_groups=name, is_align=False,dict_of_maps_of_coloumns=main_dict["live_dead"]["dict_of_maps_of_coloumns"])
        plots(temp_dict["align_protein_dict"], temp_dict["align_protein_dict_simple"], is_4_groups=name, is_align=True,dict_of_maps_of_coloumns=main_dict["live_dead"]["dict_of_maps_of_coloumns"],colors_map=colors_map)
    temp_dict["dict_of_comands"] = write_function_of_coloring_rearngment_morin(temp_dict["signature_df"], temp_dict["seq_of_proteins"], output_path=output_path,
                                                 pdb_folder_path=pdb_folder_path,colors_map=colors_map)
    temp_dict["align_db"] = main_dict["live_dead"]["align_db"]
    temp_dict["dict_of_maps_of_coloumns"] = main_dict["live_dead"]["dict_of_maps_of_coloumns"]
    main_dict[name] = temp_dict
    return main_dict

def make_signature_df(significnt_db_result_simple):
    comparisons = set()
    for col in significnt_db_result_simple.columns:
         comparisons.add(col[4:])#as live_dead or rapalive_rapadead...
    groups=[]
    for comparison in comparisons:
        groups.append(comparison.split("_")[0])
        groups.append(comparison.split("_")[1])
    signature_df = pd.DataFrame(dtype=object,index=significnt_db_result_simple.index,columns=groups,data=set())

    for protein in significnt_db_result_simple.index:
        for comparison in significnt_db_result_simple.columns:
            if type(significnt_db_result_simple.loc[protein,comparison]) is type(list()):
                for first_index,last_index,group in significnt_db_result_simple.loc[protein,comparison]:#[[10, 29, 'dead'], [50, 69, 'dead']]
                    # try:
                    #     signature_df.loc[protein,group]
                    # except KeyError:
                    #     signature_df[group] = signature_df[group].astype(object)
                    #     signature_df.at[protein, group] = set(range(first_index+1,last_index+2))
                    if not type(signature_df.loc[protein, group]) is type(set()):
                        signature_df.loc[protein, group]=set(range(first_index+1,last_index+2))
                    else:
                        signature_df.loc[protein, group] = signature_df.loc[protein, group].union(set(range(first_index+1,last_index+2)))
    def f(x):
        if type(x) is type(set()):
            a = list(x)
            a.sort()
            return a
        else:
            return np.nan
    signature_df = signature_df.applymap(f)
    return signature_df
def make_signature_df_to_rearngment_morin(significnt_db_result_simple):
    comparisons = set()
    for col in significnt_db_result_simple.columns:
         comparisons.add(col[4:])#as live_dead or rapalive_rapadead...
    groups=[]
    for comparison in comparisons:
        groups.append("IgM_"+comparison.split("_")[0])
        groups.append("IgG_"+comparison.split("_")[1])
        groups.append("IgG_"+comparison.split("_")[0])
        groups.append("IgM_"+comparison.split("_")[1])
    signature_df = pd.DataFrame(dtype=object,index=significnt_db_result_simple.index,columns=groups,data=set())
    for comparison in significnt_db_result_simple.columns:
        for protein in significnt_db_result_simple.index:
            if type(significnt_db_result_simple.loc[protein,comparison]) is type(list()):
                for first_index,last_index,group in significnt_db_result_simple.loc[protein,comparison]:#[[10, 29, 'dead'], [50, 69, 'dead']]
                    # try:
                    #     signature_df.loc[protein,group]
                    # except KeyError:
                    #     signature_df[group] = signature_df[group].astype(object)
                    #     signature_df.at[protein, group] = set(range(first_index+1,last_index+2))
                    group = comparison.split("_")[0]+"_"+group
                    if not type(signature_df.loc[protein,group]) is type(set()):
                        signature_df.loc[protein, group]=set(range(first_index+1,last_index+2))
                    else:
                        signature_df.loc[protein, group] = signature_df.loc[protein, group].union(set(range(first_index+1,last_index+2)))
    temp_df = []
    counter=0
    for i in range(0,len(signature_df.columns),2):
        def f(x):
            """
            'make tow sets or npnan to the intersection between sets only if they exists else npnna'
            :param set1:
            :param set2:
            :return:
            """
            set1 = x[0]
            set2 = x[1]

            if type(set1) is type(set()) and type(set2) is type(set()):
                return set1.intersection(set2) if len(set1.intersection(set2))>0 else np.nan
            else:
                return np.nan

        temp_df.append([i+2+counter,"diff_"+signature_df.columns[i]+"__"+signature_df.columns[i+1],signature_df.loc[:,[signature_df.columns[i],signature_df.columns[i+1]]].apply(f,axis=1)])
        counter+=1
        #signature_df[signature_df.columns[i].replace("IgM","dif")] = signature_df.loc[:,[signature_df.columns[i],signature_df.columns[i+1]]].apply(f)
    for i in temp_df:
        signature_df.insert(i[0],i[1],i[2])
    def f(x):
        if type(x) is type(set()):
            a = list(x)
            a.sort()
            return a
        else:
            return np.nan
    signature_df = signature_df.applymap(f)
    return signature_df

def write_function_of_coloring_rearngment_morin(signature_df, seq_of_proteins, output_path=os.path.join(lab_project_path, "rapa",
                                                                                       "seq", "significant", "figs", "0602_structchral_figs", "live_vs_dead_non_regarding_igg_igm"),
                                                pdb_folder_path = os.path.join(lab_project_path, "rapa","pdb", "pdb_fix_insertions"),
                                                colors_map={"ambivalent": [255 / 255, 193 / 255, 7 / 255],
                                                            "live": [128 / 255, 115 / 255, 172 / 255],
                                                            "dead": [216 / 255, 218 / 255, 235 / 255],
                                                            "diff": [0 / 255, 141 / 255, 255 / 255]}
                                                ):
    dictonery_of_small_name_to_big_name = {"H5": "4kdm_H5N1_vietnam_H5.pdb",
                                           "Vie_N1": "2hty_H5N1_vietnam_N1.pdb",
                                           "N2":"3nss_H1N1_california2009_N1.pdb",
                                           "H1":"3lzg_H1N1_california2009_H1.pdb",
                                           "Cal_N1":"3nss_H1N1_california2009_N1.pdb",
                                           "H3": "3vun_H3N2_x31_H3.pdb"}
    #colors_map = {"live": [74 / 255, 165 / 255, 24 / 255] , "dead": [204 / 255, 45 / 255, 58 / 255],"diff":[0, 141/255, 255/255]}#orig
    # cmd.set_color("differ", [0, 141 / 255, 255 / 255])
    # cmd.set_color("live", [0.2901960784313726, 0.6470588235294118, 0.09411764705882353])
    # cmd.set_color("dead", [0.8, 0.17647058823529413, 0.22745098039215686])
    #colors_map = {"live": [128 / 255, 115 / 255, 172 / 255] , "dead": [216 / 255, 218 / 255, 235 / 255],"diff":[0/255, 141/255, 255/255]}#morin
    colors_map=colors_map.copy()
    del colors_map["ambivalent"]
    #colors_map = {"live": [84 / 255, 39 / 255, 136 / 255] , "dead": [179 / 255, 88 / 255, 6 / 255],"diff":[228/255, 130/255, 217/255]}#eilay after morin

    dict_of_comands = {}
    for protein in signature_df.index:

        colors = [colors_map[x[-4:]] if (not "diff" in x) else colors_map["diff"] for x in signature_df.loc[protein,:].dropna().index]
        colors = tuple(colors)

        if protein in ["H1","H3","H5"]:
            chains = ["A","B"]
        else:
            chains=["A"]
        tuple_of_resi=[]
        tuple_of_selection=[]

        for group in signature_df.columns:
            if type(signature_df.loc[protein,group]) is type(list()):
                tuple_of_resi.append(signature_df.loc[protein,group])
                tuple_of_selection.append(group)
        tuple_of_selection = [protein+"_"+selection for selection in tuple_of_selection]
        if protein=="N2":
            tuple_of_selection2 = []
            for x in tuple_of_selection:# 'N2_IgG_dead'
                splited_by__ = x.split("_")
                #insert "onCalN1" to the second place
                splited_by__.insert(1,"onCalN1")
                tuple_of_selection2.append("_".join(splited_by__))
            tuple_of_selection = tuple_of_selection2
        tuple_of_resi = tuple(tuple_of_resi)
        tuple_of_selection = tuple(tuple_of_selection)
        if protein.startswith("H"):
            view="cmd.set_view((-0.632156432,0.241126940,-0.736363709,-0.753249407,0.031504747,0.656967521,0.181611776,0.969977736,0.161710218,-0.000299945,-0.000327468,-394.290985107,-58.018157959,-0.366636276,22.139713287,-41151.570312500,41940.089843750,-20.000000000))"
        else:
            view="cmd.set_view((-0.596710622,0.698415458,-0.395157874,0.416881382,0.690569997,0.591022372,0.685664713,0.187936321,-0.703230739,-0.000304744,-0.000427932,-195.570053101,16.614908218,35.565357208,10.529327393,-4724.140136719,5115.273437500,-20.000000000))"
        dict_of_comands[protein] = """run /home/eilay/PycharmProjects/color_peptides/coloring.py;{};cmd.set_color("differ", {});cmd.set_color("live", {});cmd.set_color("dead", {});color_pep_on_pdb("{}","{}",{},{},{},{},"{}","/home/eilay/musel",ray=1,flag_eilay=True,base_color="grey90") """.format(view,colors_map["diff"],colors_map["live"],colors_map["dead"],os.path.join(pdb_folder_path,dictonery_of_small_name_to_big_name[protein]),str(seq_of_proteins[protein]),str(chains),str(tuple_of_resi),str(colors),str(tuple_of_selection),output_path)
    return dict_of_comands


def write_function_of_coloring(signature_df,seq_of_proteins,output_path=os.path.join(lab_project_path, "rapa", "seq", "significant", "figs","0602_structchral_figs","live_vs_dead_non_regarding_igg_igm"),pdb_folder_path = os.path.join(lab_project_path,"rapa","pdb","pdb_fix_insertions") ):
    # dictonery_of_big_name_to_small_name = {"4kdm_H5N1_vietnam_H5.pdb": "H5",
    #                                        "2hty_H5N1_vietnam_N1.pdb": "Vie_N1",
    #                                        "2viu_H3N2_x31_H3.pdb": "H3",
    #                                        "":"N2",
    #                                        "3lzg_H1N1_california2009_H1.pdb":"H1",
    #                                        "3nss_H1N1_california2009_N1.pdb":"Cal_N1",
    #                                        "3vun_H3N2_x31_H3.pdb": "H3"}
    # dictonery_of_small_name_to_big_name = {v: k for k, v in dictonery_of_big_name_to_small_name.items()}
    dictonery_of_small_name_to_big_name = {
        'H5': '4kdm_H5N1_vietnam_H5.pdb',
         'Vie_N1': '2hty_H5N1_vietnam_N1.pdb',
         'H3': '3vun_H3N2_x31_H3.pdb',
         'N2': '',
         'H1': '3lzg_H1N1_california2009_H1.pdb',
         'Cal_N1': '3nss_H1N1_california2009_N1.pdb'}
    colors=[[204 / 255, 45 / 255, 58 / 255],[74 / 255, 165 / 255, 24 / 255],"magenta","orange","green","red","blue","wheat"]
    colors = colors[:signature_df.shape[1]-1]
    colors.reverse()
    colors = tuple(colors)
    dict_of_comands = {}
    for protein in signature_df.index:
        if protein in ["H1","H3","H5"]:
            chains = ["A","B"]
        else:
            chains=["A"]
        tuple_of_resi=[]
        tuple_of_selection=[]
        for group in signature_df.columns:
            if type(signature_df.loc[protein,group]) is type(list()):
                tuple_of_resi.append(signature_df.loc[protein,group])
                tuple_of_selection.append(group)
        tuple_of_selection = [protein+"_"+selection for selection in tuple_of_selection]
        tuple_of_resi = tuple(tuple_of_resi)
        tuple_of_selection = tuple(tuple_of_selection)
        dict_of_comands[protein] = """color_pep_on_pdb("{}","{}",{},{},{},{},"{}","/home/eilay/musel",ray=1) """.format(os.path.join(pdb_folder_path,dictonery_of_small_name_to_big_name[protein]),str(seq_of_proteins[protein]),str(chains),str(tuple_of_resi),str(colors),str(tuple_of_selection),output_path)
    return dict_of_comands

def main_function(significant_path=os.path.join(lab_project_path, "rapa", "seq", "significant"),
                  output_path=os.path.join(lab_project_path, "rapa", "seq", "significant", "figs","0602_structchral_figs","live_vs_dead_non_regarding_igg_igm"),
                  the_seq_fasta_file_path=os.path.join(lab_project_path, "rapa", "seq", "the_seq_from_pep.fasta"),
                  HA_align_path=os.path.join(lab_project_path, "rapa", "seq", "the_seq_from_pep_HA_align_muscel.fas"),
                  NA_align_path=os.path.join(lab_project_path, "rapa", "seq", "the_seq_from_pep_NA_align_muscel.fas"),
                  endresult_per_peptide_file_path=os.path.join(lab_project_path,"rapa","endresulte_per_peptide.tsv"),
                  pdb_folder_path = os.path.join(lab_project_path,"rapa","pdb","pdb_fix_insertions"),
                  to_plot=False,
                  colors_map={"ambivalent": [255 / 255, 193 / 255, 7 / 255],
              "live": [128 / 255, 115 / 255, 172 / 255],
              "dead": [216 / 255, 218 / 255, 235 / 255],
              "diff": [0 / 255, 141 / 255, 255 / 255]}):
    main_dict = {}
    for is_4_groups in [False,True]:
        if is_4_groups==True:
            name = "4_groups"
            data_path = os.path.join(significant_path,"data")
        else:
            name = "live_dead"
            data_path = os.path.join(significant_path, "data_live_vs_dead")
        significnt_db_result,significnt_db_result_simple = make_significent_db(data_path,is_4_groups=is_4_groups,output_path=output_path)
        signature_df = make_signature_df(significnt_db_result_simple)
        structrul_compare,structrul_compare_only_sig = make_structral_compare_per_position(significnt_db_result_simple,to_plot=to_plot,endresult_per_peptide_file_path=endresult_per_peptide_file_path,is_4_groups=is_4_groups)
        seq_of_proteins,lenght_of_proteins = make_seq_of_protein(the_seq_fasta_file_path=the_seq_fasta_file_path)
        protein_dict,protein_dict_simple = make_protein_dict(significnt_db_result, significnt_db_result_simple, lenght_of_proteins)
        align_protein_dict, align_protein_dict_simple, align_db,dict_of_maps_of_coloumns = make_align_protein_dict(protein_dict,protein_dict_simple,HA_align_path=HA_align_path,NA_align_path=NA_align_path)
        if to_plot==True:
            plots(protein_dict, protein_dict_simple, is_4_groups=is_4_groups,is_align=False)
            plots(align_protein_dict, align_protein_dict_simple, is_4_groups=is_4_groups, is_align=True)
        dict_of_comands = write_function_of_coloring(signature_df,seq_of_proteins,output_path=output_path,pdb_folder_path=pdb_folder_path)
        main_dict[name] = {"significnt_db_result":significnt_db_result,"significnt_db_result_simple":significnt_db_result_simple,
                           "signature_df":signature_df,
                           "seq_of_proteins":seq_of_proteins,"lenght_of_proteins":lenght_of_proteins,
                           "protein_dict":protein_dict,"protein_dict_simple":protein_dict_simple,
                           "align_protein_dict":align_protein_dict,"align_protein_dict_simple":align_protein_dict_simple,"align_db":align_db,"dict_of_maps_of_coloumns":dict_of_maps_of_coloumns,
                           "structrul_compare":structrul_compare,"structrul_compare_only_sig":structrul_compare_only_sig,
                           "dict_of_comands":dict_of_comands}
    main_dict = rearengment_of_main_dict(main_dict=main_dict,output_path=output_path,endresult_per_peptide_file_path= endresult_per_peptide_file_path,pdb_folder_path=pdb_folder_path,to_plot=to_plot,colors_map=colors_map)
    return main_dict
def take_dict_of_maps_of_coloumns_to_df(dict_of_maps_of_coloumns):
    df = pd.DataFrame(dtype=int)
    for protein,protein_dict_of_map in dict_of_maps_of_coloumns.items():
        for original_position,alighn_position in protein_dict_of_map.items():
            df.loc[original_position,"align_pos_of_"+protein] = alighn_position
    df.index.name = "original_position"
    #change any values to int unless it np.nan if so keep it np.nan
    df = df.fillna(-1)
    df = df.applymap(lambda x: int(x) if not np.isnan(x) else x)
    df.to_excel(os.path.join(lab_project_path, "rapa","alignment.xlsx"))
    return df

with open(os.path.join(lab_project_path, "rapa","main_dict_colorschem1_basecolordark.pkl"), 'rb') as f:
    main_dict = pickle.load(f)
    def apply_function_recursively(d, func):
        for key in d:
            if type(d[key]) == dict:
                d[key] = apply_function_recursively(d[key], func)
            elif type(d[key]) == str:
                d[key] = func(d[key])
        return d
    def replace_oem_with_eilay(s):
        return s.replace("oem", "eilay")
    main_dict_old = apply_function_recursively(main_dict, replace_oem_with_eilay)

main_dict4 = main_function(output_path=os.path.join(lab_project_path,"rapa","seq","significant","figs","color_scheme4_figs","3d"),to_plot=False,
                           colors_map={"ambivalent": [255 / 255, 193 / 255, 7 / 255],
                                       "live": [84 / 255, 39 / 255, 136 / 255],
                                       "dead": [179 / 255, 88 / 255, 6 / 255],
                                       "diff": [0 / 255, 141 / 255, 255 / 255]}
                           )
with open(os.path.join(lab_project_path, "rapa","main_dict_colorschem4.pkl"), 'wb') as f:
    pickle.dump(main_dict4,f,protocol=pickle.HIGHEST_PROTOCOL)
main_dict5 = main_function(output_path=os.path.join(lab_project_path,"rapa","seq","significant","figs","color_scheme5_figs","3d"),to_plot=False,
                           colors_map={"ambivalent": [255 / 255, 193 / 255, 7 / 255],
                                       "live": [84 / 255, 39 / 255, 136 / 255],
                                       "dead": [179 / 255, 88 / 255, 6 / 255],
                                       "diff": [0 / 255, 141 / 255, 255 / 255]}
                           )
with open(os.path.join(lab_project_path, "rapa","main_dict_colorschem5.pkl"), 'wb') as f:
    pickle.dump(main_dict5,f,protocol=pickle.HIGHEST_PROTOCOL)
main_dict6 = main_function(output_path=os.path.join(lab_project_path,"rapa","seq","significant","figs","color_scheme6_figs","3d"),to_plot=False,
                           colors_map={"ambivalent": [255 / 255, 193 / 255, 7 / 255],
                                       "live": [84 / 255, 39 / 255, 136 / 255],
                                       "dead": [179 / 255, 88 / 255, 6 / 255],
                                       "diff": [0 / 255, 141 / 255, 255 / 255]}
                           )
with open(os.path.join(lab_project_path, "rapa","main_dict_colorschem6.pkl"), 'wb') as f:
    pickle.dump(main_dict6,f,protocol=pickle.HIGHEST_PROTOCOL)

with open(os.path.join(lab_project_path, "rapa","main_dict_colorschem6.pkl"), 'rb') as f:
    main_dict = pickle.load(f)
    def apply_function_recursively(d, func):
        for key in d:
            if type(d[key]) == dict:
                d[key] = apply_function_recursively(d[key], func)
            elif type(d[key]) == str:
                d[key] = func(d[key])
        return d
    def replace_oem_with_eilay(s):
        return s.replace("oem", "eilay")
    main_dict_old = apply_function_recursively(main_dict, replace_oem_with_eilay)


def generate_marwa_aligmnet_table_fix(marwa_xl_file_path=os.path.join(lab_project_path,"rapa","sup_table_456_corection.xlsx"),
                                      aligment_filepath=os.path.join(lab_project_path,"rapa","alignment.xlsx")):
    alignmet_db = pd.read_excel(aligment_filepath,index_col=0)
    marwa_xl_db = pd.read_excel(marwa_xl_file_path)
    def func(x):
        if x =="Vn1023":
            return "Vn1203"
        else:
            return x
    marwa_xl_db["Strain"]= marwa_xl_db["Strain"].apply(func)
    def func(x):
        peptide_num = x["Peptide Number"]
        protein = x["Protein"]
        strain = x["Strain"]
        dictonery = {'Vn1023': 'Vie',
                     'CA09': 'Cal',
                     'Vn1203': 'Vie'}
        if protein == "N1":
            name = dictonery[strain]+"_N1"
        else:
            name = protein
        return alignmet_db.loc[peptide_num,"align_pos_of_" + name]
    marwa_xl_db["Adjusted Peptide_fix"] = marwa_xl_db.apply(func,axis=1)
    marwa_xl_db["differ"] = marwa_xl_db["Adjusted Peptide_fix"] == marwa_xl_db["Adjusted Peptide"]
    marwa_xl_db.to_excel(os.path.join(lab_project_path,"rapa","sup_table_456_corection_fixxxx.xlsx"))

# #main_dict2 = main_function(output_path="/home/eilay/PycharmProjects/lab_projects/rapa/seq/significant/figs",to_plot=False)
# main_dict1 = main_function(output_path=os.path.join(lab_project_path,"rapa","seq","significant","figs","color_scheme1_figs","3d"),to_plot=False,
#                            colors_map={"ambivalent": [255 / 255, 193 / 255, 7 / 255],
#                                        "live": [128 / 255, 115 / 255, 172 / 255],
#                                        "dead": [216 / 255, 218 / 255, 235 / 255],
#                                        "diff": [0 / 255, 141 / 255, 255 / 255]}
#                            )
# with open(os.path.join(lab_project_path, "rapa","main_dict_colorschem1.pkl"), 'wb') as f:
#     pickle.dump(main_dict1,f,protocol=pickle.HIGHEST_PROTOCOL)
# for file in os .listdir(os.path.join(lab_project_path,"rapa","seq","significant","figs")):
#     if os.path.isfile(os.path.join(lab_project_path,"rapa","seq","significant","figs",file)):
#         import shutil
#         shutil.move(os.path.join(lab_project_path,"rapa","seq","significant","figs",file),os.path.join(lab_project_path,"rapa","seq","significant","figs","color_scheme1_figs",file))
# main_dict2 = main_function(output_path=os.path.join(lab_project_path,"rapa","seq","significant","figs","color_scheme2_figs","3d"),to_plot=False,
#                            colors_map={"ambivalent": [255 / 255, 193 / 255, 7 / 255],
#                                        "live": [84 / 255, 39 / 255, 136 / 255],
#                                        "dead": [179 / 255, 88 / 255, 6 / 255],
#                                        "diff": [0 / 255, 141 / 255, 255 / 255]}
#                            )
# with open(os.path.join(lab_project_path, "rapa","main_dict_colorschem2.pkl"), 'wb') as f:
#     pickle.dump(main_dict1,f,protocol=pickle.HIGHEST_PROTOCOL)
#
#
# for file in os .listdir(os.path.join(lab_project_path,"rapa","seq","significant","figs")):
#     if os.path.isfile(os.path.join(lab_project_path,"rapa","seq","significant","figs",file)):
#         import shutil
#         shutil.move(os.path.join(lab_project_path,"rapa","seq","significant","figs",file),os.path.join(lab_project_path,"rapa","seq","significant","figs","color_scheme2_figs",file))
#
# main_dict3 = main_function(output_path=os.path.join(lab_project_path,"rapa","seq","significant","figs","color_scheme3_figs","3d"),to_plot=False,
#                            colors_map={"ambivalent": [255 / 255, 193 / 255, 7 / 255],
#                                        "live": [128 / 255, 115 / 255, 172 / 255],
#                                        "dead": [216 / 255, 218 / 255, 235 / 255],
#                                        "diff": [0 / 255, 141 / 255, 255 / 255]}
#                            )
# with open(os.path.join(lab_project_path, "rapa","main_dict_colorschem1_basecolordark.pkl"), 'wb') as f:
#     pickle.dump(main_dict1,f,protocol=pickle.HIGHEST_PROTOCOL)

#main_dict2 = main_function(output_path="/home/eilay/PycharmProjects/lab_projects/rapa/seq/significant/figs/0602_structchral_figs/igG_vs_igM/final_red=dead_green=life/color_scheme2",to_plot=False,
# colors_map={"ambivalent":[255 / 255, 193 / 255, 7 / 255],
# #         "dead": [179 / 255, 88 / 255, 6 / 255],
# #         "live": [84 / 255, 39 / 255, 136 / 255],
# #         "diff": [228 / 255, 130 / 255, 217 / 255]
# #               })
# cmap = {"ambivalent":[255 / 255, 193 / 255, 7 / 255],
#         "dead": [179 / 255, 88 / 255, 6 / 255],
#         "live": [84 / 255, 39 / 255, 136 / 255] ,
#         "diff": [228 / 255, 130 / 255, 217 / 255]
#               }#eilay after morin
#main_dict2 = main_function(output_path="/home/eilay/PycharmProjects/lab_projects/rapa/seq/significant/figs/0602_structchral_figs/igG_vs_igM/final_red=dead_green=life/color_scheme2",to_plot=False)
# main_dict2 = main_function(output_path="/home/eilay/PycharmProjects/lab_projects/rapa/seq/significant/figs/0602_structchral_figs/igG_vs_igM/final_red=dead_green=life/color_scheme2",to_plot=False)

# if __name__ =="__main__":
#     endresulte_per_position = pd.read_csv("/home/eilay/PycharmProjects/lab_projects/rapa/endresulte_per_position.tsv",
#                                           sep="\t",index_col=0)
#     endresulte_per_peptide = pd.read_csv("/home/eilay/PycharmProjects/lab_projects/rapa/endresulte_per_peptide.tsv",
#                                          sep="\t",index_col=0)






