import time
import math
from Bio.PDB import *
import nglview as nv
import numpy as np
import time
import os
import concurrent.futures
import pandas as pd
import subprocess
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio import PDB
import itertools
import json
import seaborn as sns
import matplotlib.pyplot as plt
from matplotlib.patches import Patch
from scipy import stats
import matplotlib.ticker as ticker
from shutil import copyfile
from distutils.dir_util import copy_tree
import shutil
import codecs
import re
from scipy.spatial import distance
from scipy.cluster.hierarchy import dendrogram, linkage
from scipy.cluster.hierarchy import fclusterdata
from Bio.PDB.PDBExceptions import PDBConstructionWarning
import warnings
warnings.simplefilter('ignore', PDBConstructionWarning)
import multiprocessing
# from kneed import KneeLocator
from sklearn.datasets import make_blobs
from sklearn.cluster import KMeans
from sklearn.metrics import silhouette_score
from sklearn.preprocessing import StandardScaler
from matplotlib_venn import venn2, venn2_circles, venn2_unweighted
from sklearn.preprocessing import MinMaxScaler
from sklearn.manifold import TSNE

def make_cool_heatmap(data_frame,pallet="Blues",fig_size=(30,19),title = ""):
    """

    Parameters
    ----------
    data_frame : data_frame pandas
        any 2D data frame pandas.
    pallet : str
        like :Blues","mako","coolwarm".
        https://seaborn.pydata.org/tutorial/color_palettes.html
    fig_size : tupel, optional
        DESCRIPTION. The default is (30,19).
        for more squar put (30,24)

    Returns
    -------
    None.

    """
    sns.set(font_scale=3)
    fig,ax = plt.subplots(subplot_kw=dict(frameon=False),figsize=fig_size)
    figure = sns.heatmap(data_frame,linewidths=.3,ax=ax,cmap=sns.mpl_palette(pallet,100),xticklabels=True, yticklabels=True)
    figure.set_xticklabels(figure.get_xmajorticklabels(), fontsize = 15)
    figure.set_yticklabels(figure.get_ymajorticklabels(), fontsize = 15)
    figure.set_title(title)
    plt.tight_layout()
    plt.show()
    return fig

"""
########################################################################################################################
mouse data
########################################################################################################################
"""
#reading data
igm_mouse_db_with_comparing = pd.read_excel("/mnt/Disk2/marwa_eilay/IgMRapa4_Final_proteins_viruses_array_data.xlsx",index_col=0)
igg_mouse_db_with_comparing = pd.read_excel("/mnt/Disk2/marwa_eilay/IgGRapa4_Final_proteins_viruses_array_data.xlsx",index_col=0)
cols_of_real_data = igg_mouse_db_with_comparing.loc[:,:'rNA_70_H5N1_A/Hubei/1/2010'].columns

#lstrip D28 from igg index and d14 from igm
igm_mouse_db_with_comparing.index = igm_mouse_db_with_comparing.index.map(lambda x:x.lstrip("d14_"))
igg_mouse_db_with_comparing.index = igg_mouse_db_with_comparing.index.map(lambda x:x.lstrip("D28_"))

#checking for differnce cols
print("the cols that are on igg(28) but not igm(14) "+str(set(igg_mouse_db_with_comparing.columns).difference(set(igm_mouse_db_with_comparing.columns))))
print("the cols that are on igm(14) but not igg(28) "+str(set(igm_mouse_db_with_comparing.columns).difference(set(igg_mouse_db_with_comparing.columns))))
#checking for differnce rows
print("the rows that are on igg(28) but not igm(14) "+str(set(igg_mouse_db_with_comparing.index).difference(set(igm_mouse_db_with_comparing.index))))
print("the rows that are on igm(14) but not igg(28) "+str(set(igm_mouse_db_with_comparing.index).difference(set(igg_mouse_db_with_comparing.index))))

#adding rapa col True for Rapa false for PBS
igm_mouse_db_with_comparing["rapa"] = igm_mouse_db_with_comparing["group"].apply(lambda x: x.split("_")[0]=="Rapa")
igg_mouse_db_with_comparing["rapa"] = igg_mouse_db_with_comparing["group"].apply(lambda x: x.split("_")[0]=="Rapa")

#adding live col True for Live false for Dead
igm_mouse_db_with_comparing["live"] = igm_mouse_db_with_comparing["group"].apply(lambda x: x.split("_")[1]=="Live")
igg_mouse_db_with_comparing["live"] = igg_mouse_db_with_comparing["group"].apply(lambda x: x.split("_")[1]=="Live")

#adding type col True for Live false for Dead
def type_db(x):
    y=x.split("_")[-1]
    if y == "CD1d":
        return y
    elif y == "WT":
        return y
    else:
        return "rapa3"
igm_mouse_db_with_comparing["type"] = igm_mouse_db_with_comparing["group"].apply(type_db)
igg_mouse_db_with_comparing["type"] = igg_mouse_db_with_comparing["group"].apply(type_db)

#drop CD1d
igm_mouse_db_with_comparing = igm_mouse_db_with_comparing[igm_mouse_db_with_comparing["type"]!="CD1d"]
igg_mouse_db_with_comparing = igg_mouse_db_with_comparing[igg_mouse_db_with_comparing["type"]!="CD1d"]

#dont care about rapa 3/4 col
igm_mouse_db_with_comparing["group_unspesific"] = igm_mouse_db_with_comparing["group"].apply(lambda x:x.rstrip("_WT").replace("Rapa","R").replace("PBS","P").replace("Live","L").replace("Dead","D"))
igg_mouse_db_with_comparing["group_unspesific"] = igg_mouse_db_with_comparing["group"].apply(lambda x:x.rstrip("_WT").replace("Rapa","R").replace("PBS","P").replace("Live","L").replace("Dead","D"))

#tsne plot_igg_group_unspesific
tsne_em = TSNE(n_components=2,perplexity=30.0,n_iter=1000,verbose=1).fit_transform(igg_mouse_db_with_comparing.loc[:,cols_of_real_data])
from bioinfokit.visuz import cluster
color_class = igg_mouse_db_with_comparing["group_unspesific"].to_numpy()
cluster.tsneplot(score=tsne_em,colorlist=color_class,legendpos='upper right',legendanchor=(1.16,1),show=True)
#cluster.tsneplot(score=tsne_em,show=True)

#tsne plot_igm_group_unspesific
tsne_em = TSNE(n_components=2,perplexity=30.0,n_iter=1000,verbose=1).fit_transform(igm_mouse_db_with_comparing.loc[:,cols_of_real_data])
from bioinfokit.visuz import cluster
color_class = igg_mouse_db_with_comparing["group_unspesific"].to_numpy()
cluster.tsneplot(score=tsne_em,colorlist=color_class,legendpos='upper right',legendanchor=(1.16,1),show=True)
#cluster.tsneplot(score=tsne_em,show=True)


"""
########################################################################################################################
protein data
########################################################################################################################
"""

#reading data
igg_protein_db =  igg_mouse_db_with_comparing.loc[:,:'rNA_70_H5N1_A/Hubei/1/2010'].T
igm_protein_db =  igm_mouse_db_with_comparing.loc[:,:'rNA_70_H5N1_A/Hubei/1/2010'].T

#add col of subtype
def subtype(x):
    y = x.split("_")[2]
    try:
        int(int(y))
        return x.split("_")[3]
    except ValueError:
        return x.split("_")[2]
igg_protein_db["subtype"] = igg_protein_db.index.map(subtype)
igm_protein_db["subtype"] = igm_protein_db.index.map(subtype)

#add col of rHA
igg_protein_db["protein"] = igg_protein_db.index.map(lambda x:x.split("_")[0])
igm_protein_db["protein"] = igm_protein_db.index.map(lambda x:x.split("_")[0])

#sorting
igg_protein_db.sort_values(by = ["subtype","protein"],inplace=True)
igm_protein_db.sort_values(by = ["subtype","protein"],inplace=True)

#first look on data by heatmap
make_cool_heatmap(igg_protein_db.iloc[:,:-2],title="igg",fig_size=(40,19))
make_cool_heatmap(igm_protein_db.iloc[:,:-2],title="igm",fig_size=(40,19))

# #clustring by mouses
# X = igg_protein_db.iloc[:,:-2]
# Y = distance.pdist(X)
# squre_Y = distance.squareform(Y)
# figure = sns.clustermap(squre_Y, method='complete', xticklabels=True, yticklabels=True, figsize=(30, 24))
# plt.show()
#
# X = igm_protein_db.iloc[:,:-2]
# Y = distance.pdist(X)
# squre_Y = distance.squareform(Y)
# figure = sns.clustermap(squre_Y, method='complete', xticklabels=True, yticklabels=True, figsize=(30, 24))
# plt.show()

def cluster_map(protein_db,mouse_db_with_comparing):

    igg_protein_db = protein_db
    igg_mouse_db_with_comparing = mouse_db_with_comparing
    # clustering
    ##color row
    lut_row = dict(zip(igg_protein_db["subtype"].unique(), "rbgy"))
    row_colors = igg_protein_db["subtype"].map(lut_row)
    lut_row2 = dict(zip(igg_protein_db["protein"].unique(), "kwc"))
    row_colors = pd.DataFrame(row_colors).join(pd.DataFrame(igg_protein_db["protein"].map(lut_row2)))
    ##color col
    lut_col = dict(zip(igg_mouse_db_with_comparing["group_unspesific"].unique(), "kwcm"))
    col_colors = igg_mouse_db_with_comparing["group_unspesific"].map(lut_col)
    lut_col2 = dict(zip(igg_mouse_db_with_comparing["type"].unique(), "rb"))
    col_colors = pd.DataFrame(col_colors).join(pd.DataFrame(igg_mouse_db_with_comparing["type"].map(lut_col2)))
    ##making clustermap
    sns.set(font_scale=1)
    figure = sns.clustermap(igg_protein_db.iloc[:, :-2], row_colors=row_colors, col_colors=col_colors, cmap="Blues",
                            method='complete', xticklabels=True, yticklabels=True,figsize=(40,34) )
    sns.set(font_scale=3)
    figure.fig.suptitle('igm')
    # figure.gs.update(left=0.05, right=30/40,top=24/34,bottom=0.052)
    #
    # gs = figure.fig.add_gridspec(40,34)
    # ax1b = figure.fig.add_subplot(gs[30:40, 0:12])#legend1 row
    # ax1c = figure.fig.add_subplot(gs[30:40, 12:24])#legend2 row
    # ax1d = figure.fig.add_subplot(gs[0:15, 24:34])#legend1 col
    # ax1e = figure.fig.add_subplot(gs[15:30, 24:34])#legend2_col


    #legends
    #handles = [Patch(facecolor=lut_row[name]) for name in lut_row]
    #plt.legend(handles, lut_row, title='Subtype',
    #           bbox_to_anchor=(1, 1), bbox_transform=plt.gcf().transFigure, loc='center right')
    #handles = [Patch(facecolor=lut_row2[name]) for name in lut_row2]
    # plt.legend(handles, lut_row2, title='protein',
    #            bbox_to_anchor=(1, 1), bbox_transform=plt.gcf().transFigure, loc='best',ax=ax1c)
    # handles = [Patch(facecolor=lut_col[name]) for name in lut_col]
    # plt.legend(handles, lut_col, title='Subtype',
    #            bbox_to_anchor=(1, 1), bbox_transform=plt.gcf().transFigure, loc='best',ax=ax1d)
    # handles = [Patch(facecolor=lut_col2[name]) for name in lut_col2]
    # plt.legend(handles, lut_col2, title='protein',
    #            bbox_to_anchor=(1, 1), bbox_transform=plt.gcf().transFigure, loc='best',ax=ax1e)

    plt.show()

cluster_map(igg_protein_db,igg_mouse_db_with_comparing)
cluster_map(igm_protein_db,igm_mouse_db_with_comparing)
plt.close("all")