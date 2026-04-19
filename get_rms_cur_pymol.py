import numpy as np
import pandas as pd
import os
username = os.getlogin()
if username!="eilay":
    COV2_project = "\\\\132.72.92.166\\eilay\\beackup_to_linux_server\\COV_2_copy_new" #leptop win
else:
    COV2_project = "/run/user/1000/gvfs/smb-share:server=132.72.92.166,share=eilay/beackup_to_linux_server/COV_2_copy_new" #leptop
    #COV2_project = "/run/user/1000/gvfs/smb-share:server=132.72.92.166,share=eilay/beackup_to_linux_server/COV_2_copy_new_adding_new_data" #leptop

regardent_folder_path = os.path.join(COV2_project,"regardent")
regardent_PDB_folder_path = os.path.join(regardent_folder_path ,"PDB")
regardent_contacts_folder_path = os.path.join(regardent_folder_path ,"contacts")
regardent_analyzide_data_folder_path = os.path.join(regardent_folder_path ,"analyzide_data")
regardent_checker_folder_path = os.path.join(regardent_folder_path ,"checker")
regardent_database_folder_path = os.path.join(regardent_folder_path,"database")
unregardent_folder_path = os.path.join(COV2_project,"unregardent")
unregardent_PDB_folder_path = os.path.join(unregardent_folder_path ,"PDB")
unregardent_contacts_folder_path = os.path.join(unregardent_folder_path,"contacts")
unregardent_analyzide_data_folder_path = os.path.join(unregardent_folder_path ,"analyzide_data")
unregardent_checker_folder_path = os.path.join(unregardent_folder_path ,"checker")
repair_folder_path = os.path.join(COV2_project,"unregardent_repaired")
repair_PDB_folder_path = os.path.join(repair_folder_path ,"PDB")
repair_mutate_PDB_folder_path = os.path.join(repair_folder_path ,"mutations")
repair_mutate_full_PDB_folder_path = os.path.join(repair_folder_path ,"mutaions_full")
repair_analyzide_data_folder_path =os.path.join( repair_folder_path ,"analyzide_data")
unanalyzed_folder_path = os.path.join(COV2_project ,"unanalyzed") #unanalyzide_folder_path>unanalyzed_folder_path
foldx_path = "/home/sacharen/foldxLinux64"
the_sabdab_file_path = os.path.join(regardent_database_folder_path,"sabdab_summary_all.tsv")
the_cov3d_file_path = os.path.join(regardent_database_folder_path,"spike_structures.tsv")


#merged = pd.read_excel(os.path.join(regardent_analyzide_data_folder_path, "bigger_cmpnd_db.xlsx"), index_col=0)
#tow_cluster_metrix_togther = pd.read_csv("/run/user/1000/gvfs/smb-share:server=132.72.92.166,share=eilay/beackup_to_linux_server/cov2_280822/new_data_base_on_cov2_new_data_without_adding_/temp/unregardent_analyse/analyse/fig_of_only_new_versus_old.tsv",sep="\t",index_col=0)
midi_bigger=pd.read_excel("/run/user/1000/gvfs/smb-share:server=132.72.92.166,share=eilay/beackup_to_linux_server/cov2_280822/regardent/analyzide_data/midell_bigger_cmpnd_db.xlsx",index_col=0)
a= list(range(1,1274))+["clean_Cluster"]
tow_cluster_metrix_togther = midi_bigger.loc[(midi_bigger["represnted_that_choosed_of_new_batch"]==1) | (midi_bigger["old_11_clusters_represnted_that_choosed"]==1),a].copy().dropna(subset=["clean_Cluster"])
tow_cluster_metrix_togther=tow_cluster_metrix_togther.rename(columns={"clean_Cluster":"Cluster"})
cmd.set_color( "cluster_1", [0.8941176470588236, 0.10196078431372549, 0.10980392156862745])
cmd.set_color(  "cluster_2", [0.21568627450980393, 0.49411764705882355, 0.7215686274509804])
cmd.set_color(  "cluster_3", [0.30196078431372547, 0.6862745098039216, 0.2901960784313726])
cmd.set_color(  "cluster_4", [0.596078431372549, 0.3058823529411765, 0.6392156862745098])
cmd.set_color(  "cluster_5", [1.0, 0.4980392156862745, 0.0])
cmd.set_color(  "cluster_6", [1.0, 1.0, 0.2])
cmd.set_color(  "cluster_7", [0.6509803921568628, 0.33725490196078434, 0.1568627450980392])
cmd.set_color(  "cluster_8", [0.9686274509803922, 0.5058823529411764, 0.7490196078431373])
cmd.set_color(  "cluster_9", [0.10196078431372549, 0.8196078431372549, 1.0])
cmd.set_color(  "cluster_10", [0.8, 0.6, 1.0])
cmd.set_color(  "cluster_11",[0.5215686274509804, 0.8784313725490196, 0.5215686274509804])
cmd.set_color(  "cluster_nan", [0.6, 0.23, 0.6])
cmd.set_color(  "cluster_new", [0.15, 0.23, 0.9])


groups = {}
all_obgects = cmd.get_names('objects', 0,'(all)')
cmd.color("gray70", "all", 1)
for obj in all_obgects:
    #num_cluster = merged.loc[obj.split("_")[-2],"sub_cluster"]
    #num_cluster = str(int(merged.loc[obj.split("_")[-2],"clean_Cluster"]))
    try:
        num_cluster = tow_cluster_metrix_togther.loc[obj.split("_")[-2],"Cluster"]
    except:
        cmd.delete(obj)#non in the strike metrix
        continue
    groups.setdefault(num_cluster,[])
    groups[num_cluster].append(obj)

for cluster_num,obj_group in groups.items():
    sub_cluster = cluster_num

    try:
        #cluster_num= int(cluster_num[:-1])
        cluster_num = int(cluster_num)
        if cluster_num>11:
            cluster_num = "new"
    except:
        cluster_num = "new"
    cmd.group("cluster{}".format(sub_cluster), members=" ".join(obj_group))
    cmd.color("cluster_{}".format(cluster_num),"cluster{}".format(int(sub_cluster)))
all_obgects = cmd.get_names('objects', 0,'(all)')
for obj in all_obgects:
    spike_chains = "+".join(list(obj.split("_")[-1]))
    light_chain = obj.split("_")[-2][-1]
    heavy_chain = obj.split("_")[-2][-2]
    cmd.color("palecyan","obj {} and chain {}".format(obj,spike_chains), 1)
    cmd.color("wheat","obj {} and resi 333-360+436-500 and chain {}".format(obj,spike_chains), 1)

#
# lists=[]
# for i,obj1 in enumerate(all_obgects):
#     spike_chains1 = "+".join(list(obj1.split("_")[-1]))
#     light_chain1 = obj1.split("_")[-2][-1]
#     heavy_chain1 = obj1.split("_")[-2][-2]
#     not_all = all_obgects.copy()
#     not_all.remove(obj1)
#     for j,obj2 in enumerate(not_all):
#         if i < j:
#             print(obj1,obj2)
#             spike_chains2 = "+".join(list(obj2.split("_")[-1]))
#             light_chain2 = obj2.split("_")[-2][-1]
#             heavy_chain2 = obj2.split("_")[-2][-2]
#             cmd.super("{} and chain {}".format(obj1,spike_chains1.split("+")[0]), "{} and chain {}".format(obj2,spike_chains2.split("+")[0]))
#             cmd.super("{} and chain {}+{}".format(obj1,light_chain1,heavy_chain1), "{} and chain {}+{}".format(obj2,light_chain2,heavy_chain2),transform=0,object="aln")
#             try:
#                 a = cmd.rms_cur("{} and chain {}+{} & aln".format(obj1,light_chain1,heavy_chain1),  "{} and chain {}+{} & aln".format(obj2,light_chain2,heavy_chain2), mobile_state=0, target_state=0,quiet=1, matchmaker=-1, cutoff=2.0, cycles=0,object=None)
#             except:
#                 a = np.NAN
#             lists.append([obj1,obj2,a])
#             lists.append([obj2,obj1,a])
# rmsss = pd.DataFrame(lists)
# rmsss = pd.pivot_table(rmsss, values=2, index=[0],columns=[1])
# #rmsss.to_csv("/run/user/1000/gvfs/smb-share:server=132.72.92.166,share=eilay/beackup_to_linux_server/COV_2_copy_new/regardent/analyzide_data/rms_pivot_RBD_rmss.tsv","\t")

for obj in all_obgects:
    cmd.set_name(obj,obj.replace("alighn_PyIgClassifyOutput_","").replace("alighn_alighn_ref_PyIgClassifyOutput_","").replace("alighn_ref_PyIgClassifyOutput_",""))
#run /home/eilay/PycharmProjects/lab_projects/get_rms_cur_pymol.py

