import os

import numpy as np
import pandas as pd
import freesasa
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from rapa_significant_peptide_color import start_end_of_peptide_by_id_number
import pdb_guide_Eilay_jj
import sys
sys.path.insert(0, "/home/eilay/PycharmProjects/color_peptides")
import coloring
import importlib
importlib.reload(pdb_guide_Eilay_jj)
importlib.reload(coloring)

"""
parcing getarea data:

parsing to consurf data:

"""

def make_aligmnent(pdb_folder_path="/home/eilay/PycharmProjects/lab_projects/rapa/pdb/pdb_fix_insertions",fasta_of_seq_from_pep="/home/eilay/PycharmProjects/lab_projects/rapa/seq/the_seq_from_pep.fasta",output_folder_path="/home/eilay/PycharmProjects/lab_projects/rapa/alighment_between_seq_to_pdb",muscle_folder_path="/home/eilay/musel"):
    """
    :param pdb_folder_path: the folder path of your rapa pdbs files
    :param fasta_of_seq_from_pep: the fasta of seq from pdb genrate by rapa_significant_peptide_color.the_seq_from_peptide
    :param output_folder_path: just the output folder
    :param muscle_folder_path: musscel path
    :return:
    """
    dict_of_seqs_db = {}
    dictonery_of_big_name_to_small_name = {"4kdm_H5N1_vietnam_H5.pdb": ("H5",["A","B"]),
                                            "2hty_H5N1_vietnam_N1.pdb": ("Vie_N1",["A"]),
                                            "2viu_H3N2_x31_H3.pdb": ("H3",["A","B"]),
                                            "3lzg_H1N1_california2009_H1.pdb": ("H1",["A","B"]),
                                            "3nss_H1N1_california2009_N1.pdb": ("Cal_N1",["A"]),
                                            "3vun_H3N2_x31_H3.pdb": ("H3",["A","B"])}
    atm_seq = {}
    for file in os.listdir(pdb_folder_path):
        the_seq = ""
        small_name,chains = dictonery_of_big_name_to_small_name[file]
        chain_col = ""
        for chain in chains:
            temp_seq = coloring.seq(pdb_file_path=os.path.join(pdb_folder_path,file),chain = chain)
            the_seq += temp_seq
            chain_col += chain*len(temp_seq)
        atm_seq[small_name] = [the_seq,chain_col]

    pep_seq = {}
    for record in SeqIO.parse(fasta_of_seq_from_pep, "fasta"):
        pep_seq[record.id.split(" ")[0]] = record.seq

    for key in atm_seq.keys():
        seq_db = coloring.generate_seq_db(pep_seq[key],atm_seq[key][0],output_folder_path,muscle_folder_path,startsswith = key+"_")
        seq_db.loc[seq_db["atm_pdb_seq"]!="-","chain"] = list(atm_seq[key][1])
        for chain in seq_db["chain"].unique():
            seq_db.loc[seq_db["chain"] == chain, "atm_count"] = range(1, len(
                seq_db.loc[seq_db["chain"] == chain, "atm_count"]) + 1)
        seq_db["atm_id"] =  seq_db["chain"] +"_" + seq_db["atm_count"].astype(int).astype(str)
        seq_db["atm_count"] = seq_db["atm_count"].astype(int)
        seq_db["pep_count"] = seq_db["pep_count"].astype(int)
        seq_db.columns =  seq_db.columns + "_" + key
        seq_db = seq_db.rename(columns = {"pep_count"+ "_" + key :"pep_count"})
        dict_of_seqs_db[key] = seq_db
    return dict_of_seqs_db


sasa_folder = "/home/eilay/PycharmProjects/lab_projects/rapa/sasa"
def make_sasa_df(sasa_folder=sasa_folder):
    """
    #http://curie.utmb.edu/getarea.html
    #<20 is i (inside)
    #>50 is o (outside)

    this is a python code for parcing http://curie.utmb.edu/getarea.html data

    OLD FUNCTION
    by the getarea site
    :param sasa_folder:
    :return:
    """
    sasa_df = pd.DataFrame()
    count=0
    for sasa_file_name in os.listdir(sasa_folder):
        protein = sasa_file_name.split("_")[-2]+"_"+sasa_file_name.split("_")[-1].rstrip(".txt")
        a = pd.read_csv(os.path.join(sasa_folder,sasa_file_name), sep="\s+",
                        engine='python')
        a = a.groupby("Position").mean()#mean of an all the chains because is homomer
        a["In/Out"] = a["Ratio(%)"].apply(lambda x: "i" if x<20 else "o" if x>50 else np.nan)
        a["Position"] = a.index
        a.index.name = ""
        #c = a.drop_duplicates("Position")
        #b["AminoAcid"] = c["AminoAcid"]
        #a = b
        if count==0:
            sasa_df["Position"] = a["Position"]
            sasa_df["Ratio(%)_"+protein] = a["Ratio(%)"]
            sasa_df["In/Out_" + protein] = a["In/Out"]
        else:
            a = a.rename(columns={"Ratio(%)":"Ratio(%)"+"_"+protein,"In/Out":"In/Out"+"_"+protein})
            sasa_df = sasa_df.merge(a.loc[:,["Ratio(%)"+"_"+protein,"In/Out"+"_"+protein,"Position"]],how="outer",left_on ="Position",right_on = "Position")
        count+=1
    sasa_df = sasa_df.sort_values(by="Position")
    return sasa_df

def get_max_sasa_for_aa(amino_acid,currect_value=""):
    """
    from https://pubmed.ncbi.nlm.nih.gov/24278298/
    :param amino_acid: in three or one letter
    :param currect_value: if gave so it will return his RSA
    :return: his max sasa
    """
    if len(amino_acid) == 3 :
        amino_acid = pdb_guide_Eilay_jj.get_3id_to_1id_or_1id_from_3id(amino_acid)
    amino_acid = amino_acid.upper()
    dict_of_max_sasa = {"A":121,"R":265,"N":187,"D":187,"C":148,"E":214,"Q":225,"G":97,
                        "H":216,"I":195,"L":191,"K":230,"M":203,"F":228,"P":154,"S":143,
                        "T":163,"W":264,"Y":255,"V":165}
    if currect_value == "":
        return dict_of_max_sasa[amino_acid]
    else:
        return currect_value/dict_of_max_sasa[amino_acid]

def sasa_df_by_dssp(sasa_folder=os.path.join(sasa_folder,"dssp")):
    """
    OLD FUNCTION
    https://www3.cmbi.umcn.nl/xssp/?tour=true
    :param sasa_folder:
    :return:
    shittrttttttttttttttttttttttttttttt
    """
    sasa_df = pd.DataFrame(columns=["protein","aa","chain","SASA","RSA"])
    for file in os.listdir(sasa_folder):
        if file.endswith(".dssp"):
            protein = file.split("_")[-1].lstrip(".dssp")
            with open(os.path.join(sasa_folder,file)) as f:
                counter = 0
                for line in f.readlines():
                    if counter < 28:# to skip the garbech
                        counter+=1
                        continue
                    SASA = int(line[35:38].replace(" ", ""))
                    aa = line[13:14].replace(" ", "")
                    try:
                        RSA = get_max_sasa_for_aa(aa, currect_value=SASA)
                    except:
                        return line
                    chain = line[11:12].replace(" ", "")
                    sasa_df = sasa_df.append([protein,aa,chain,SASA,RSA])
    return sasa_df

clean_insersion_pdb_folder_path = "/home/eilay/PycharmProjects/lab_projects/rapa/pdb/pdb_fix_insertions"
def sasa_from_free_sasa(clean_insersion_pdb_folder_path = clean_insersion_pdb_folder_path ):
    """
    :param pdb_folder_path:no insertion allowed
    :return:
    """
    # dictonery_of_big_name_to_small_name = {"4kdm_H5N1_vietnam_H5.pdb": ("H5",["A","B"]),
    #                                        "2hty_H5N1_vietnam_N1.pdb": ("Vie_N1",["A"]),
    #                                        "2viu_H3N2_x31_H3.pdb": ("H3",["A","B"]),
    #                                        "3lzg_H1N1_california2009_H1.pdb": ("H1",["A","B"]),
    #                                        "3nss_H1N1_california2009_N1.pdb": ("Cal_N1",["A"]),
    #                                        "3vun_H3N2_x31_H3.pdb": ("H3",["A","B"])}
    dictonery_of_big_name_to_small_name = {"4kdm_H5N1_vietnam_H5.pdb": ("H5",["A","B"]),
                                           "2hty_H5N1_vietnam_N1.pdb": ("Vie_N1",["A"]),
                                           "3lzg_H1N1_california2009_H1.pdb": ("H1",["A","B"]),
                                           "3nss_H1N1_california2009_N1.pdb": ("Cal_N1",["A"]),
                                           "3vun_H3N2_x31_H3.pdb": ("H3",["A","B"])}

    sasa_dict = {}
    for file in os.listdir(clean_insersion_pdb_folder_path):
        sasa_df = pd.DataFrame()
        if file.endswith(".pdb"):
            structure = freesasa.Structure(os.path.join(clean_insersion_pdb_folder_path,file))
            result = freesasa.calc(structure)
            area_res = result.residueAreas()
            for chain,resideis in area_res.items():
                if chain in dictonery_of_big_name_to_small_name[file][1]:
                    for res,sasa in resideis.items():
                        sasa_df.loc[chain+"_"+res, "sasa_"+dictonery_of_big_name_to_small_name[file][0]] = sasa.relativeTotal


            sasa_dict[dictonery_of_big_name_to_small_name[file][0]] = sasa_df
    return sasa_dict

conserved_path = "/home/eilay/PycharmProjects/lab_projects/rapa/conservation/by-peptide-seq"
def make_conserved_df(conserved_path = conserved_path):
    """
    :param conserved_path: the path that contain all the output from consurf server
    https://consurf.tau.ac.il/
    :return:
    """
    conserved_df = pd.DataFrame()
    count=0
    for folder in os.listdir(conserved_path):#name of protein
        if os.path.isdir(os.path.join(conserved_path,folder)):
            for file in os.listdir(os.path.join(conserved_path,folder)):
                if file == "msa_aa_variety_percentage.csv":
                    aa_variety_percentage = pd.read_csv(os.path.join(conserved_path,folder,file), skiprows=3)
                    #protein = folder.split("_")[2] +"_"+ folder.split("_")[3]
                    protein = folder
                    if count == 0:
                        conserved_df["Position"] = aa_variety_percentage["pos"]
                        conserved_df["conservation_" + protein] = aa_variety_percentage["ConSurf Grade"]
                    else:
                        aa_variety_percentage = aa_variety_percentage.rename(columns={"ConSurf Grade": "conservation_" + protein,"pos":"Position"})
                        conserved_df = conserved_df.merge(
                            aa_variety_percentage.loc[:, ["conservation_" + protein, "Position"]], how="outer",
                            left_on="Position", right_on="Position")
                    count += 1

                    conserved_df["conservation_" + protein] = conserved_df["conservation_" + protein].astype(str)

    for folder in os.listdir(conserved_path):#name of protein
        if os.path.isdir(os.path.join(conserved_path,folder)):
            for file in os.listdir(os.path.join(conserved_path,folder)):
                if file == "consurf.grades":
                    # if folder == "no_H3N2_x31_N2":
                    grade = pd.read_csv(
                        os.path.join(conserved_path,folder,file),
                        skiprows=13, sep="\t+", engine='python')
                    # else:
                    #     grade = pd.read_csv(
                    #         os.path.join(conserved_path,folder,file),
                    #         skiprows=12, sep="\t+", engine='python')
                    grade = grade.iloc[1:-2, :]  # delete normalized wierd row
                    #protein = folder.split("_")[2] +"_"+ folder.split("_")[3]
                    protein = folder
                    grade = grade.rename(columns={"SCORE": "conservation_nonbinery_" + protein,"POS":"Position"," SEQ":"SEQ"})
                    grade["Position"] = grade["Position"].astype(int)
                    conserved_df = conserved_df.merge(
                        grade.loc[:, ["conservation_nonbinery_" + protein, "Position"]], how="outer",
                        left_on="Position", right_on="Position")
                    count += 1
    conserved_df = conserved_df.rename(columns={"Position":"pep_count"})
    return conserved_df

pdb_folder_path = "/home/eilay/PycharmProjects/lab_projects/rapa/pdb"
def structchers_seq(pdb_folder_path=pdb_folder_path):
    """
    :param pdb_folder_path: the folder path of your rapa pdbs files
    :return: the_seq_from_pdb.fasta
    but i never use that
    """
    sequences = []
    for file_name in os.listdir(pdb_folder_path):
        if file_name.endswith(".pdb"):
            protein_name = file_name.lstrip(".pdb").split("_")[-1]#H1,H3,H5,N1
            chains = pdb_guide_Eilay_jj.get_chains(os.path.join(pdb_folder_path,file_name))
            the_seq = ""
            if "H" in protein_name:
                HA1 = pdb_guide_Eilay_jj.seq(pdb_file_path=os.path.join(pdb_folder_path,file_name),chain = chains[0])
                HA2 = pdb_guide_Eilay_jj.seq(pdb_file_path=os.path.join(pdb_folder_path,file_name),chain = chains[1])
                the_seq = HA1+HA2
            elif "N" in protein_name:
                the_seq = pdb_guide_Eilay_jj.seq(pdb_file_path=os.path.join(pdb_folder_path,file_name),chain = chains[0])
            dictonery_of_big_name_to_small_name = { "4kdm_H5N1_vietnam_H5.pdb":"H5 4kdm","2hty_H5N1_vietnam_N1.pdb":"Vie_N1 2hty",
                                                    "2viu_H3N2_x31_H3.pdb":"H3 2viu","3lzg_H1N1_california2009_H1.pdb":"H1 3lzg",
                                                    "3nss_H1N1_california2009_N1.pdb":"Cal_N1 3nss"}
            name = dictonery_of_big_name_to_small_name[file_name]
            record = SeqRecord(
                Seq(the_seq),
                id=name,
                name=name,
                description=name)
            sequences.append(record)

    with open(os.path.join(pdb_folder_path,"the_seq_from_pdb.fasta"), "w") as output_handle:
        SeqIO.write(sequences, output_handle, "fasta")

def get_endreslut_position(conserved_df,aligments_dict,sasa_dict):
    """
    :param conserved_df: generate from make_conserved_df
    :param aligments_dict: generate from make_aligmnent
    :param sasa_dict: generate from sasa_from_free_sasa
    :return: endresulte_per_position
    """
    endresulte_per_position = conserved_df
    for name in aligments_dict.keys():
        sasa_align = aligments_dict[name].loc[aligments_dict[name]["protein_seq_"+name]!="-"].merge(sasa_dict[name],how="left",right_index=True,left_on="atm_id_"+name)
        endresulte_per_position = endresulte_per_position.merge(sasa_align,on="pep_count",how="left")
    return endresulte_per_position

def get_end_result_to_pep(endresulte_per_position,first_index_of_peptides=0,len_of_peptides=20,len_of_overlap=15,last_index_of_peptides=0):
    """
    :param endresulte_per_position: endresulte_per_position
    :param first_index_of_peptides:
    :param len_of_peptides:
    :param len_of_overlap:
    :param last_index_of_peptides:
    :return: take the endresult per position and split tham to peptide
    """
    endresulte_per_peptide = pd.DataFrame()
    number_of_peptide = int(
        (endresulte_per_position.iloc[-1, 0] - (first_index_of_peptides + 1)) / (len_of_peptides - len_of_overlap) + 1)
    for peptide_num in range(1, number_of_peptide + 1):
        first_index, last_index = start_end_of_peptide_by_id_number(peptide_num, first_index_of_peptides=first_index_of_peptides,
                                                                    len_of_peptides=len_of_peptides, len_of_overlap=len_of_overlap,
                                                                    last_index_of_peptides=last_index_of_peptides)
        if last_index > endresulte_per_position.iloc[-1, 0]:
            last_index =  endresulte_per_position.iloc[-1, 0]-1

        endresulte_per_peptide[peptide_num] = endresulte_per_position.loc[first_index:last_index].mean(
            numeric_only=True)
        endresulte_per_peptide.loc["first_position", peptide_num] = first_index + 1
        endresulte_per_peptide.loc["last_position", peptide_num] = last_index + 1

    endresulte_per_peptide = endresulte_per_peptide.T
    endresulte_per_peptide = endresulte_per_peptide.drop(
        ["pep_count", "atm_count_H5", "atm_count_Cal_N1", "atm_count_H3", "atm_count_Vie_N1", "atm_count_H1"], axis=1)
    endresulte_per_peptide = endresulte_per_peptide.loc[:,
                             ['first_position', "last_position",
                              'conservation_nonbinery_H5', 'sasa_H5',
                              'conservation_nonbinery_Vie_N1', 'sasa_Vie_N1',
                              'conservation_nonbinery_H3', 'sasa_H3',
                              'conservation_nonbinery_N2',
                              'conservation_nonbinery_H1', 'sasa_H1',
                              'conservation_nonbinery_Cal_N1', 'sasa_Cal_N1', ]]
    endresulte_per_peptide['first_position'] = endresulte_per_peptide['first_position'].astype(int)
    endresulte_per_peptide['last_position'] = endresulte_per_peptide['last_position'].astype(int)
    """
    for col in endresulte_per_peptide.columns:
        if col in ['first_position', "last_position"]:
            continue
        if col.startswith("sasa"):
            continue
        the_real_last_pep_count = endresulte_per_position.loc[
            ~ pd.isna(endresulte_per_position[col]), ["pep_count"]].max().max()
        
        if col == 'conservation_nonbinery_N2':
            rows_that_are_bigger_the_last_pep_count = endresulte_per_peptide.loc[endresulte_per_peptide[
                                                                                     "last_position"] >= the_real_last_pep_count, [
                                                                                     col]].iloc[1:, :]
        else:
            rows_that_are_bigger_the_last_pep_count = endresulte_per_peptide.loc[endresulte_per_peptide[
                                                                                     "last_position"] >= the_real_last_pep_count, [
                                                                                     col, "sasa" + "_" + "_".join(
                    col.split("_")[2:])]].iloc[1:, :]
        for index in rows_that_are_bigger_the_last_pep_count.index:
            for col in rows_that_are_bigger_the_last_pep_count.columns:
                endresulte_per_peptide.loc[index, col] = np.nan
    endresulte_per_peptide.loc[
        111, 'conservation_nonbinery_H3'] = np.nan  # fix that loc because in the first input that i get from marwa H3N2_X31.xlsx should also be a peptide 111 althow it was not their
    """
    return endresulte_per_peptide

lab_project_path = os.path.dirname(os.path.abspath("change_file_type.py"))
first_index_of_peptides=0
len_of_peptides=20
len_of_overlap=15
last_index_of_peptides=0
def make_endresult_df(lab_project_path=lab_project_path,save=True,first_index_of_peptides=0,len_of_peptides=20,len_of_overlap=15,last_index_of_peptides=0):
    """
    THE MAIN FUNCTION

    :param lab_project_path: the suource of the GIT
    :param save: is to save the endresult or not defult is yes
    :param first_index_of_peptides:
    :param len_of_peptides:
    :param len_of_overlap:
    :param last_index_of_peptides:
    :return:
    endresulte_per_peptide,endresulte_per_position and can save them
    """
    conserved_df = make_conserved_df(os.path.join(lab_project_path,"rapa","conservation","by-peptide-seq"))
    aligments_dict = make_aligmnent()
    sasa_dict = sasa_from_free_sasa()
    endresulte_per_position = get_endreslut_position(conserved_df,aligments_dict,sasa_dict)
    endresulte_per_peptide = get_end_result_to_pep(endresulte_per_position,first_index_of_peptides=first_index_of_peptides,len_of_peptides=len_of_peptides,len_of_overlap=len_of_overlap,last_index_of_peptides=last_index_of_peptides)
    if save:
        endresulte_per_position.to_csv(os.path.join("rapa","endresulte_per_position.tsv"),"\t")
        endresulte_per_peptide.to_csv(os.path.join("rapa","endresulte_per_peptide.tsv"),"\t")
    return endresulte_per_peptide,endresulte_per_position

endresulte_per_peptide_cc,endresulte_per_position_cc = make_endresult_df(save=True)
#here the protein dont align one to each other!!!!!


"""
seqs_db = pd.DataFrame()
for record in SeqIO.parse("/home/eilay/PycharmProjects/lab_projects/rapa/seq/pdb_align_to_pep/H1.nexus", "nexus"):
    seqs_db[record.id] = list(str(record.seq))
seqs_db.index = seqs_db.index+1
seqs_db.columns = ["peptide_seq","atm_pdb_seq"]
for i in seqs_db.index:
    subseries_of_the_atm_pdb_count = seqs_db.loc[:i,"atm_pdb_seq"]#the subseries of the atm pdb count
    pdb_atm_count = (subseries_of_the_atm_pdb_count != "-").sum()
    seqs_db.loc[i,"atm_count"] = int(pdb_atm_count)
    subseries_of_the_atm_pdb_count = seqs_db.loc[:i,"peptide_seq"]#the subseries of the atm pdb count
    pdb_atm_count = (subseries_of_the_atm_pdb_count != "-").sum()
    seqs_db.loc[i,"pep_count"] = int(pdb_atm_count)
"""


def get_pdb_residiue_to_color_by_first_position_of_peptide(first_Position_of_peptide,seqs_db):
    """
    will get first position of peptide and return list of position in pdb atom counting.
    seqs_db is a pd DataFrame witch the rows are position, the first col is the peptide seq
    the second col is the pdb atm seq. both are align by musscel algoritem.
    ? represent disorderd residue on atm seq.
    - repersent gapping by musscel algoritem.

    :param first_Position_of_peptide:int
    :param seqs_db: pandas DataFrame
    :return: list of resiude regardent the pdb atm counting
    """
