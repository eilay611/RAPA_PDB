import os
import pickle
import rapa_structcheral_charectrize

dictonery_of_big_name_to_small_name = {"4kdm_H5N1_vietnam_H5.pdb": ("H5", ["A", "B"]),
                                       "2hty_H5N1_vietnam_N1.pdb": ("Vie_N1", ["A"]),
                                       "2viu_H3N2_x31_H3.pdb": ("H3", ["A", "B"]),
                                       "3lzg_H1N1_california2009_H1.pdb": ("H1", ["A", "B"]),
                                       "3nss_H1N1_california2009_N1.pdb": ("Cal_N1", ["A"]),
                                       "3vun_H3N2_x31_H3.pdb": ("H3", ["A", "B"])}
dictonery_of_small_name_to_big_name={'H5': ['4kdm_H5N1_vietnam_H5.pdb', ['A', 'B']],
 'Vie_N1': ['2hty_H5N1_vietnam_N1.pdb', ['A']],
 'H3': ['3vun_H3N2_x31_H3.pdb', ['A', 'B']],
 'H1': ['3lzg_H1N1_california2009_H1.pdb', ['A', 'B']],
 'Cal_N1': ['3nss_H1N1_california2009_N1.pdb', ['A']]}
rapa_folder = "/home/eilay/PycharmProjects/lab_projects/rapa"
pdb_folders = os.path.join(rapa_folder,"pdb","pdb_fix_insertions")
aligment_between_pdb_to_seq = rapa_structcheral_charectrize.make_aligmnent()
vikas_folder=os.path.join(rapa_folder,"vikas")
for file in os.listdir(vikas_folder):
    if file.endswith(".pkl"):
        igm_or_igg = file.removesuffix(".pkl")[-3:]#IgM or IgG
        a_file  = open(os.path.join(vikas_folder,file), "rb")
        picle_vikas = pickle.load(a_file)
        a_file.close()
        for protein,dict_of_bfactor in picle_vikas.items():
            if protein == 'N2':
                continue
            pdb_file_path = os.path.join(pdb_folders,dictonery_of_small_name_to_big_name[protein][0])
            OUT = open(os.path.join(vikas_folder,"b_factor","{}_{}.pdb".format(protein,igm_or_igg)), 'w')
            aligment_between_pdb_to_seq_protein = aligment_between_pdb_to_seq[protein]
            with open(pdb_file_path, "r") as IN:
                for line in IN.readlines():
                    if line.startswith('ATOM'):
                        l1 = line
                        # 0123456789012345678901234567890123456789012345678901234567890123456789
                        # ATOM      1  N   LEU A  62      55.829  64.222  83.106  1.00 62.37           N
                        res_number = int(line[22:26].replace(" ", ""))
                        chain_id = line[21:22].replace(" ", "")
                        atm_id = chain_id+"_"+str(res_number)
                        try:
                            b_factor_ceil_2 = "{:.2f}".format(dict_of_bfactor[int(aligment_between_pdb_to_seq_protein.loc[aligment_between_pdb_to_seq_protein["atm_id_"+protein]==atm_id,"pep_count"])])
                        except TypeError:
                            b_factor_ceil_2 = "{:.2f}".format(0)  # for a shorter list
                        # try:
                        #     b_factor_ceil_2 = "{:.2f}".format(list_of_new_b_factor[res_number - 1])
                        # except IndexError:
                        #     b_factor_ceil_2 = "{:.2f}".format(0)  # for a shorter list
                        #     if not flag:
                        #         flag = True
                        #         if (
                        #                 input("you enter a shorter list than the len of the chain of your input pdb. if it ok enter 1 if not enter 0") == 0):
                        #             raise ValueError("shorter list than the len of the chain")
                        #
                        REVAH = 6 - len(b_factor_ceil_2)
                        bfactor = " " * REVAH + b_factor_ceil_2
                        if not (len(bfactor) == 6):
                            raise ValueError("fuck" + res_number)
                        line = line[:60] + bfactor + line[66:]
                        OUT.write(line)
            OUT.close()

