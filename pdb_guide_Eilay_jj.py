the_sabdab_file_path = "/mnt/Disk2/cov/COV_2_copy_new/regardent/database/sabdab_summary_all.tsv"
the_cov3d_file_path = "/mnt/Disk2/cov/COV_2_copy_new/regardent/database/spike_structures.tsv"
import time
import math
import pdbtools
from Bio.PDB import *
import nglview as nv
import numpy as np
import time
spike_seq = "MFVFLVLLPLVSSQCVNLTTRTQLPPAYTNSFTRGVYYPDKVFRSSVLHSTQDLFLPFFSNVTWFHAIHVSGTNGTKRFDNPVLPFNDGVYFASTEKSNIIRGWIFGTTLDSKTQSLLIVNNATNVVIKVCEFQFCNDPFLGVYYHKNNKSWMESEFRVYSSANNCTFEYVSQPFLMDLEGKQGNFKNLREFVFKNIDGYFKIYSKHTPINLVRDLPQGFSALEPLVDLPIGINITRFQTLLALHRSYLTPGDSSSGWTAGAAAYYVGYLQPRTFLLKYNENGTITDAVDCALDPLSETKCTLKSFTVEKGIYQTSNFRVQPTESIVRFPNITNLCPFGEVFNATRFASVYAWNRKRISNCVADYSVLYNSASFSTFKCYGVSPTKLNDLCFTNVYADSFVIRGDEVRQIAPGQTGKIADYNYKLPDDFTGCVIAWNSNNLDSKVGGNYNYLYRLFRKSNLKPFERDISTEIYQAGSTPCNGVEGFNCYFPLQSYGFQPTNGVGYQPYRVVVLSFELLHAPATVCGPKKSTNLVKNKCVNFNFNGLTGTGVLTESNKKFLPFQQFGRDIADTTDAVRDPQTLEILDITPCSFGGVSVITPGTNTSNQVAVLYQGVNCTEVPVAIHADQLTPTWRVYSTGSNVFQTRAGCLIGAEHVNNSYECDIPIGAGICASYQTQTNSPRRARSVASQSIIAYTMSLGAENSVAYSNNSIAIPTNFTISVTTEILPVSMTKTSVDCTMYICGDSTECSNLLLQYGSFCTQLNRALTGIAVEQDKNTQEVFAQVKQIYKTPPIKDFGGFNFSQILPDPSKPSKRSFIEDLLFNKVTLADAGFIKQYGDCLGDIAARDLICAQKFNGLTVLPPLLTDEMIAQYTSALLAGTITSGWTFGAGAALQIPFAMQMAYRFNGIGVTQNVLYENQKLIANQFNSAIGKIQDSLSSTASALGKLQDVVNQNAQALNTLVKQLSSNFGAISSVLNDILSRLDKVEAEVQIDRLITGRLQSLQTYVTQQLIRAAEIRASANLAATKMSECVLGQSKRVDFCGKGYHLMSFPQSAPHGVVFLHVTYVPAQEKNFTTAPAICHDGKAHFPREGVFVSNGTHWFVTQRNFYEPQIITTDNTFVSGNCDVVIGIVNNTVYDPLQPELDSFKEELDKYFKNHTSPDVDLGDISGINASVVNIQKEIDRLNEVAKNLNESLIDLQELGKYEQYIKWPWYIWLGFIAGLIAIVMVTIMLCCMTSCCSCLKGCCSCGSCCKFDEDDSEPVLKGVKLHYT"
spike_seq = spike_seq[:613]+"G"+spike_seq[614:]
import pickle
import os
import concurrent.futures
import pandas as pd
import subprocess
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio import PDB
from itertools import product
#import pymol
import itertools
import json
import seaborn as sns
import matplotlib.pyplot as plt
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
unique_atm_combi_dict = {'N_N': None, 'N_CA': None, 'N_C': None, 'N_O': None, 'N_CB': None, 'N_CG': None, 'N_CD': None, 'N_OE1': None, 'N_NE2': None, 'CA_N': None, 'CA_CA': None, 'CA_C': None, 'CA_O': None, 'CA_CB': None, 'CA_CG': None, 'CA_CD': None, 'CA_OE1': None, 'CA_NE2': None, 'C_N': None, 'C_CA': None, 'C_C': None, 'C_O': None, 'C_CB': None, 'C_CG': None, 'C_CD': None, 'C_OE1': None, 'C_NE2': None, 'O_N': None, 'O_CA': None, 'O_C': None, 'O_O': None, 'O_CB': None, 'O_CG': None, 'O_CD': None, 'O_OE1': None, 'O_NE2': None, 'CB_N': None, 'CB_CA': None, 'CB_C': None, 'CB_O': None, 'CB_CB': None, 'CB_CG': None, 'CB_CD': None, 'CB_OE1': None, 'CB_NE2': None, 'CG_N': None, 'CG_CA': None, 'CG_C': None, 'CG_O': None, 'CG_CB': None, 'CG_CG': None, 'CG_CD': None, 'CG_OE1': None, 'CG_NE2': None, 'CD_N': None, 'CD_CA': None, 'CD_C': None, 'CD_O': None, 'CD_CB': None, 'CD_CG': None, 'CD_CD': None, 'CD_OE1': None, 'CD_NE2': None, 'OE1_N': None, 'OE1_CA': None, 'OE1_C': None, 'OE1_O': None, 'OE1_CB': None, 'OE1_CG': None, 'OE1_CD': None, 'OE1_OE1': None, 'OE1_NE2': None, 'NE2_N': None, 'NE2_CA': None, 'NE2_C': None, 'NE2_O': None, 'NE2_CB': None, 'NE2_CG': None, 'NE2_CD': None, 'NE2_OE1': None, 'NE2_NE2': None, 'N_SG': None, 'CA_SG': None, 'C_SG': None, 'O_SG': None, 'CB_SG': None, 'CG_SG': None, 'CD_SG': None, 'OE1_SG': None, 'NE2_SG': None, 'N_CG1': None, 'N_CG2': None, 'CA_CG1': None, 'CA_CG2': None, 'C_CG1': None, 'C_CG2': None, 'O_CG1': None, 'O_CG2': None, 'CB_CG1': None, 'CB_CG2': None, 'CG_CG1': None, 'CG_CG2': None, 'CD_CG1': None, 'CD_CG2': None, 'OE1_CG1': None, 'OE1_CG2': None, 'NE2_CG1': None, 'NE2_CG2': None, 'N_OD1': None, 'CA_OD1': None, 'CA_ND2': None, 'C_OD1': None, 'C_ND2': None, 'O_OD1': None, 'O_ND2': None, 'CB_OD1': None, 'CB_ND2': None, 'CG_OD1': None, 'CA_CD2': None, 'C_CD1': None, 'C_CD2': None, 'O_CD1': None, 'O_CD2': None, 'CB_CD2': None, 'N_CD2': None, 'O_CE2': None, 'C_OD2': None, 'O_OD2': None, 'N_CD1': None, 'N_CE1': None, 'N_CE2': None, 'N_CZ': None, 'CA_CD1': None, 'CA_CE1': None, 'CA_CE2': None, 'CA_CZ': None, 'C_CE1': None, 'C_CE2': None, 'C_CZ': None, 'O_CE1': None, 'O_CZ': None, 'CB_CD1': None, 'CB_CE1': None, 'CB_CE2': None, 'CB_CZ': None, 'CG_CD1': None, 'CG_CE1': None, 'CG_CE2': None, 'CG_CZ': None, 'CD_CE1': None, 'CD_CE2': None, 'CD_CZ': None, 'OE1_CE1': None, 'OE1_CE2': None, 'OE1_CZ': None, 'N_OE2': None, 'CA_OE2': None, 'C_OE2': None, 'O_OE2': None, 'CB_OE2': None, 'CG_OE2': None, 'CD_OE2': None, 'OE1_OE2': None, 'NE2_OE2': None, 'CG_CD2': None, 'N_NE': None, 'N_NH1': None, 'N_NH2': None, 'CA_NE': None, 'CA_NH1': None, 'CA_NH2': None, 'C_NE': None, 'C_NH1': None, 'C_NH2': None, 'O_NE': None, 'O_NH1': None, 'O_NH2': None, 'CB_NE': None, 'CB_NH1': None, 'CB_NH2': None, 'CG_NE': None, 'CG_NH1': None, 'CG_NH2': None, 'CD_NE': None, 'CD_NH1': None, 'CD_NH2': None, 'OE1_NE': None, 'OE1_NH1': None, 'OE1_NH2': None, 'NE2_NE': None, 'NE2_CZ': None, 'NE2_NH1': None, 'NE2_NH2': None, 'N_OH': None, 'CA_OH': None, 'C_OH': None, 'O_OH': None, 'N_OG': None, 'CA_OG': None, 'C_OG': None, 'O_OG': None, 'CB_OG': None, 'CG_OG': None, 'CD_OG': None, 'OE1_OG': None, 'NE2_OG': None, 'SG_N': None, 'SG_CA': None, 'SG_C': None, 'SG_O': None, 'SG_CB': None, 'SG_CG': None, 'SG_CD': None, 'SG_OE1': None, 'SG_NE2': None, 'SG_SG': None, 'SG_CG1': None, 'SG_CG2': None, 'N_ND2': None, 'SG_OD1': None, 'SG_ND2': None, 'SG_CD1': None, 'SG_CD2': None, 'C_OG1': None, 'O_OG1': None, 'SG_OG': None, 'SG_CE2': None, 'SG_CE1': None, 'SG_CZ': None, 'N_OD2': None, 'CA_OD2': None, 'CB_OD2': None, 'SG_OD2': None, 'SG_OE2': None, 'SG_NE': None, 'SG_NH1': None, 'SG_NH2': None, 'SG_OH': None, 'C_NE1': None, 'CG1_N': None, 'CG1_CA': None, 'CG1_C': None, 'CG1_O': None, 'CG1_CB': None, 'CG1_CG': None, 'CG1_CD': None, 'CG1_OE1': None, 'CG1_NE2': None, 'CG2_N': None, 'CG2_CA': None, 'CG2_C': None, 'CG2_O': None, 'CG2_CB': None, 'CG2_CG': None, 'CG2_CD': None, 'CG2_OE1': None, 'CG2_NE2': None, 'CG1_SG': None, 'CG2_SG': None, 'CG1_CG1': None, 'CG1_CG2': None, 'CG2_CG1': None, 'CG2_CG2': None, 'CG1_OD1': None, 'CG1_ND2': None, 'CG2_OD1': None, 'CG2_ND2': None, 'CG1_CD1': None, 'CG1_CD2': None, 'CG2_CD1': None, 'CG2_CD2': None, 'N_OG1': None, 'CA_OG1': None, 'CB_OG1': None, 'CG1_OG1': None, 'CG2_OG1': None, 'CG1_NH2': None, 'CG2_NE': None, 'CG2_CZ': None, 'CG2_NH1': None, 'CG2_NH2': None, 'CG2_CE1': None, 'CG2_CE2': None, 'CG1_OD2': None, 'CG2_OD2': None, 'CG1_CE1': None, 'CG1_CE2': None, 'CG1_CZ': None, 'CG1_OE2': None, 'CG2_OE2': None, 'CG1_NE': None, 'CG1_NH1': None, 'CG1_OG': None, 'CG2_OG': None, 'N_NE1': None, 'N_CZ2': None, 'CA_NE1': None, 'CA_CE3': None, 'CA_CZ2': None, 'CA_CZ3': None, 'CA_CH2': None, 'C_CE3': None, 'C_CZ2': None, 'C_CZ3': None, 'C_CH2': None, 'O_NE1': None, 'O_CE3': None, 'O_CZ2': None, 'O_CZ3': None, 'O_CH2': None, 'CB_NE1': None, 'CB_CE3': None, 'CB_CZ2': None, 'CB_CZ3': None, 'CB_CH2': None, 'CG1_NE1': None, 'CG1_CE3': None, 'CG1_CZ2': None, 'CG2_NE1': None, 'CG2_CE3': None, 'CG2_CZ2': None, 'CG2_CZ3': None, 'CG2_CH2': None, 'OD1_N': None, 'OD1_CA': None, 'OD1_C': None, 'OD1_O': None, 'OD1_CB': None, 'OD1_CG': None, 'ND2_CA': None, 'ND2_C': None, 'ND2_O': None, 'ND2_CB': None, 'OD1_SG': None, 'ND2_N': None, 'ND2_SG': None, 'OD1_CG1': None, 'OD1_CG2': None, 'ND2_CG1': None, 'ND2_CG2': None, 'CG_ND2': None, 'OD1_OD1': None, 'OD1_ND2': None, 'ND2_CG': None, 'ND2_OD1': None, 'ND2_ND2': None, 'OD1_CD1': None, 'OD1_CD2': None, 'ND2_CD1': None, 'ND2_CD2': None, 'CG_OG1': None, 'OD1_OG1': None, 'ND2_OG1': None, 'OD1_CD': None, 'OD1_NE': None, 'OD1_CZ': None, 'OD1_NH1': None, 'OD1_NH2': None, 'ND2_NE': None, 'ND2_CZ': None, 'ND2_NH1': None, 'ND2_NH2': None, 'C_NZ': None, 'O_CE': None, 'O_NZ': None, 'OD1_CE2': None, 'OD1_NE2': None, 'ND2_NE2': None, 'CG_OD2': None, 'OD1_OD2': None, 'ND2_OD2': None, 'ND2_CD': None, 'OD1_CE1': None, 'ND2_CE1': None, 'ND2_CE2': None, 'OD1_OG': None, 'ND2_OG': None, 'OD1_OE1': None, 'ND2_OE1': None, 'N_CE3': None, 'N_CZ3': None, 'N_CH2': None, 'CG_NE1': None, 'CG_CZ2': None, 'CG_CH2': None, 'OD1_NE1': None, 'OD1_CZ2': None, 'OD1_CH2': None, 'ND2_NE1': None, 'CD1_C': None, 'CD1_O': None, 'CD2_CA': None, 'CD2_C': None, 'CD2_O': None, 'CD2_CB': None, 'CD1_N': None, 'CD1_CA': None, 'CD1_CB': None, 'CD1_SG': None, 'CD2_N': None, 'CD2_SG': None, 'CD1_CG1': None, 'CD1_CG2': None, 'CD2_CG1': None, 'CD2_CG2': None, 'CD1_CG': None, 'CD1_OD1': None, 'CD1_ND2': None, 'CD2_CG': None, 'CD2_OD1': None, 'CD2_ND2': None, 'CD1_CD1': None, 'CD1_CD2': None, 'CD2_CD1': None, 'CD2_CD2': None, 'CD1_OG1': None, 'CD2_OG1': None, 'CD1_CD': None, 'CD1_NE': None, 'CD1_CZ': None, 'CD1_NH1': None, 'CD1_NH2': None, 'CD2_CD': None, 'CD2_NE': None, 'CD2_CZ': None, 'CD2_NH1': None, 'CD2_NH2': None, 'CD1_ND1': None, 'CD2_ND1': None, 'N_CE': None, 'N_NZ': None, 'CA_CE': None, 'CA_NZ': None, 'C_CE': None, 'CB_CE': None, 'CB_NZ': None, 'CG_CE': None, 'CG_NZ': None, 'CD1_CE': None, 'CD1_NZ': None, 'CD2_CE': None, 'CD2_NZ': None, 'CD1_CE1': None, 'CD1_CE2': None, 'CD2_CE1': None, 'CD2_CE2': None, 'CD1_OD2': None, 'CD2_OD2': None, 'CD1_OE1': None, 'CD1_OE2': None, 'CD2_OE1': None, 'CD2_OE2': None, 'CD1_NE2': None, 'CD1_OG': None, 'CD2_OG': None, 'CG_CE3': None, 'CG_CZ3': None, 'CD1_NE1': None, 'CD1_CE3': None, 'CD1_CZ2': None, 'CD1_CZ3': None, 'CD1_CH2': None, 'CD2_NE1': None, 'CD2_CE3': None, 'CD2_CZ2': None, 'CD2_CZ3': None, 'CD2_CH2': None, 'CD1_OH': None, 'OG1_C': None, 'OG1_O': None, 'OG1_N': None, 'OG1_CA': None, 'OG1_CB': None, 'OG1_CG1': None, 'OG1_CG2': None, 'OG1_CG': None, 'OG1_OD1': None, 'OG1_ND2': None, 'OG1_CD1': None, 'OG1_CD2': None, 'OG1_OG1': None, 'OG1_CD': None, 'OG1_NE': None, 'OG1_CZ': None, 'OG1_NH1': None, 'OG1_NH2': None, 'C_ND1': None, 'O_ND1': None, 'OG1_CE': None, 'OG1_NZ': None, 'CG2_CE': None, 'CG2_NZ': None, 'OG1_CE2': None, 'OG1_OD2': None, 'OG1_CE1': None, 'OG1_OG': None, 'OG1_NE1': None, 'OG1_CE3': None, 'OG1_CZ2': None, 'OG1_CZ3': None, 'OG1_CH2': None, 'N_ND1': None, 'CA_ND1': None, 'CB_ND1': None, 'OG1_ND1': None, 'OG1_NE2': None, 'CG2_ND1': None, 'NH2_C': None, 'NE_CA': None, 'NE_C': None, 'NE_O': None, 'NE_CB': None, 'NE_CG2': None, 'CZ_N': None, 'CZ_CA': None, 'CZ_C': None, 'CZ_O': None, 'CZ_CB': None, 'CZ_CG2': None, 'NH1_N': None, 'NH1_CA': None, 'NH1_C': None, 'NH1_O': None, 'NH1_CB': None, 'NH1_CG2': None, 'NH2_N': None, 'NH2_CA': None, 'NH2_O': None, 'NH2_CB': None, 'NH2_CG1': None, 'NH2_CG2': None, 'CD_OD1': None, 'NE_N': None, 'NE_CG': None, 'NE_OD1': None, 'NE_ND2': None, 'CZ_CG': None, 'CZ_OD1': None, 'CZ_ND2': None, 'NH1_CG': None, 'NH1_OD1': None, 'NH1_ND2': None, 'NH2_CG': None, 'NH2_OD1': None, 'NH2_ND2': None, 'CD_CD1': None, 'CD_CD2': None, 'NE_CD1': None, 'NE_CD2': None, 'CZ_CD1': None, 'CZ_CD2': None, 'NH1_CD1': None, 'NH1_CD2': None, 'NH2_CD1': None, 'NH2_CD2': None, 'CD_OG1': None, 'NE_OG1': None, 'CZ_OG1': None, 'NH1_OG1': None, 'NH2_OG1': None, 'NE_CD': None, 'NE_NE': None, 'NE_CZ': None, 'NE_NH1': None, 'NE_NH2': None, 'CZ_CD': None, 'CZ_NE': None, 'CZ_CZ': None, 'CZ_NH1': None, 'CZ_NH2': None, 'NH1_CD': None, 'NH1_NE': None, 'NH1_CZ': None, 'NH1_NH1': None, 'NH1_NH2': None, 'NH2_CD': None, 'NH2_NE': None, 'NH2_CZ': None, 'NH2_NH1': None, 'NH2_NH2': None, 'NE_OE1': None, 'NE_NE2': None, 'CZ_OE1': None, 'NH2_OE1': None, 'NE_CE1': None, 'CZ_CE1': None, 'NH1_CE1': None, 'CG_ND1': None, 'CD_ND1': None, 'NE_ND1': None, 'CD_CE': None, 'CD_NZ': None, 'NE_CE': None, 'NE_NZ': None, 'CZ_CE': None, 'CZ_NZ': None, 'NH1_CE': None, 'NH1_NZ': None, 'NH2_CE': None, 'NH2_NZ': None, 'NE_CE2': None, 'CZ_CE2': None, 'NH1_CE2': None, 'NH2_CE1': None, 'NH2_CE2': None, 'CD_OD2': None, 'NE_OD2': None, 'CZ_OD2': None, 'NH1_OD2': None, 'NH2_OD2': None, 'CD_ND2': None, 'CZ_NE2': None, 'NH1_OE1': None, 'NH1_NE2': None, 'NH2_NE2': None, 'NE_OG': None, 'CZ_OG': None, 'NH1_OG': None, 'NH2_OG': None, 'CD_NE1': None, 'CD_CE3': None, 'CD_CZ2': None, 'CD_CZ3': None, 'CD_CH2': None, 'NE_NE1': None, 'NE_CE3': None, 'NE_CZ2': None, 'NE_CZ3': None, 'NE_CH2': None, 'CZ_NE1': None, 'CZ_CE3': None, 'CZ_CZ2': None, 'CZ_CZ3': None, 'CZ_CH2': None, 'NH1_NE1': None, 'NH1_CE3': None, 'NH1_CZ2': None, 'NH1_CZ3': None, 'NH1_CH2': None, 'NH2_NE1': None, 'NH2_CE3': None, 'NH2_CZ2': None, 'NH2_CZ3': None, 'NH2_CH2': None, 'CG_OH': None, 'CD_OH': None, 'NE_OH': None, 'CZ_OH': None, 'NH1_OH': None, 'NH2_OH': None, 'OG1_OE1': None, 'OE1_OG1': None, 'NE2_OG1': None, 'OE1_CD1': None, 'OE1_CD2': None, 'NE2_CD1': None, 'NE2_CD2': None, 'OE1_OD1': None, 'OE1_OD2': None, 'NE2_OD1': None, 'CD2_NE2': None, 'CB_OH': None, 'CE1_CA': None, 'CE1_C': None, 'CE1_O': None, 'CE1_CB': None, 'CE2_C': None, 'CE2_O': None, 'OH_N': None, 'OH_CA': None, 'OH_C': None, 'OH_O': None, 'OH_CB': None, 'CE1_N': None, 'CE1_CG': None, 'CE1_CD': None, 'CE2_N': None, 'CE2_CA': None, 'CE2_CB': None, 'CE2_CG': None, 'CE2_CD': None, 'OH_CG': None, 'OH_CD': None, 'CD2_OH': None, 'CE1_CD1': None, 'CE1_CD2': None, 'CE1_CE1': None, 'CE1_CE2': None, 'CE1_CZ': None, 'CE1_OH': None, 'CE2_CD1': None, 'CE2_CD2': None, 'CE2_CE1': None, 'CE2_CE2': None, 'CE2_CZ': None, 'CE2_OH': None, 'OH_CD1': None, 'OH_CD2': None, 'OH_CE1': None, 'OH_CE2': None, 'OH_CZ': None, 'OH_OH': None, 'CE1_OG1': None, 'CE1_CG2': None, 'CE2_OG1': None, 'CE2_CG2': None, 'OH_CG2': None, 'CE2_OD1': None, 'CE1_OD1': None, 'CE1_ND2': None, 'CE2_ND2': None, 'OH_OD1': None, 'OH_ND2': None, 'CE1_CG1': None, 'CE2_CG1': None, 'CZ_CG1': None, 'OH_CG1': None, 'OH_OG1': None, 'CE1_NE1': None, 'CE1_CE3': None, 'CE1_CZ2': None, 'CE1_CZ3': None, 'CE1_CH2': None, 'CE2_NE1': None, 'CE2_CE3': None, 'CE2_CZ2': None, 'CE2_CZ3': None, 'CE2_CH2': None, 'OH_CE3': None, 'OH_CZ3': None, 'CE1_NE2': None, 'CE2_NE2': None, 'CE1_NH2': None, 'CE2_NH2': None, 'OH_NE': None, 'OH_NH2': None, 'CE1_OD2': None, 'CE2_OD2': None, 'OH_OD2': None, 'CG2_OH': None, 'OG1_OH': None, 'OG1_OE2': None, 'OD1_CE3': None, 'OD1_CZ3': None, 'ND2_CE3': None, 'ND2_CZ3': None, 'ND2_OH': None, 'OG_CA': None, 'OG_C': None, 'OG_O': None, 'OG_CB': None, 'OG_CG': None, 'OG_CD2': None, 'OG_N': None, 'OG_OG1': None, 'OG_CG2': None, 'OG_OD1': None, 'OG_ND2': None, 'OG_OG': None, 'OG_CD1': None, 'OG_CE1': None, 'OG_CE2': None, 'OG_CZ': None, 'OG_CD': None, 'OG_NE': None, 'OG_NH1': None, 'OG_NH2': None, 'OG_CG1': None, 'OG_OH': None, 'OG_OE1': None, 'OG_OE2': None, 'OG_OD2': None, 'OG_NE2': None, 'CE1_OG': None, 'CE2_OG': None, 'CE1_NE': None, 'CE1_NH1': None, 'CE2_NE': None, 'CE2_NH1': None, 'CE2_OE2': None, 'CE1_OE1': None, 'CE2_OE1': None, 'NE_CG1': None, 'NH1_CG1': None, 'NE_OE2': None, 'CZ_OE2': None, 'NH1_OE2': None, 'NH2_OE2': None, 'NH1_ND1': None, 'CG1_OH': None, 'CG1_ND1': None, 'CG1_NZ': None, 'CG1_CE': None, 'CE2_CE': None, 'CE2_NZ': None, 'CE2_ND1': None, 'CZ_ND1': None, 'OH_OE1': None, 'OH_NE2': None, 'CE1_CE': None, 'CE1_NZ': None, 'OH_CE': None, 'OH_NZ': None, 'OH_OG': None, 'CE1_ND1': None, 'OH_ND1': None, 'CE1_OE2': None, 'OH_OE2': None, 'OD2_CG1': None, 'OD1_OH': None, 'OD2_N': None, 'OD2_CA': None, 'OD2_C': None, 'OD2_O': None, 'OD2_CB': None, 'OD2_CG': None, 'OD2_CD1': None, 'OD2_CD2': None, 'OD2_CE1': None, 'OD2_CE2': None, 'OD2_CZ': None, 'OD2_OH': None, 'OD2_CD': None, 'OD2_OD1': None, 'OD2_OD2': None, 'OD1_CE': None, 'OD1_NZ': None, 'OD2_CE': None, 'OD2_NZ': None, 'OD2_CG2': None, 'OD2_NE': None, 'OD2_NH1': None, 'OD2_NH2': None, 'OD1_ND1': None, 'OD2_ND1': None, 'OD2_NE2': None, 'OD2_OG': None, 'OD2_OG1': None, 'CE_CA': None, 'CE_C': None, 'CE_O': None, 'CE_CB': None, 'CE_CG': None, 'CE_CD2': None, 'CE_CE2': None, 'CE_CZ': None, 'NZ_C': None, 'NZ_CB': None, 'NZ_CG': None, 'NZ_CD2': None, 'NZ_CE2': None, 'CE_N': None, 'CE_CD1': None, 'CE_CE1': None, 'NZ_N': None, 'NZ_CA': None, 'NZ_O': None, 'NZ_CD1': None, 'NZ_CE1': None, 'CE_CD': None, 'NZ_CD': None, 'CE_OD1': None, 'CE_OD2': None, 'NZ_OD1': None, 'NZ_OD2': None, 'CE_CE': None, 'CE_NZ': None, 'NZ_CE': None, 'NZ_NZ': None, 'CE_CG1': None, 'CE_CG2': None, 'NZ_CG2': None, 'NZ_CZ': None, 'NZ_OH': None, 'NZ_CG1': None, 'CE_OH': None, 'CE_OE1': None, 'CE_OE2': None, 'NZ_OE1': None, 'NZ_OE2': None, 'NH2_ND1': None, 'OG_CE': None, 'OG_NZ': None, 'ND1_CG1': None, 'ND1_N': None, 'ND1_CA': None, 'ND1_C': None, 'ND1_O': None, 'ND1_CD2': None, 'ND1_CE2': None, 'ND1_CZ': None, 'ND1_CB': None, 'ND1_CG': None, 'ND1_CD1': None, 'ND1_CE1': None, 'ND1_OH': None, 'NE2_CE1': None, 'ND1_CD': None, 'ND1_OD1': None, 'ND1_OD2': None, 'NE2_OD2': None, 'ND1_CG2': None, 'ND1_NE': None, 'ND1_NH1': None, 'ND1_NH2': None, 'ND1_ND1': None, 'ND1_NE2': None, 'NE2_ND1': None, 'ND1_OG': None, 'ND1_OG1': None, 'NE2_CE2': None, 'NE2_OH': None, 'ND1_SG': None, 'CE1_SG': None, 'ND1_CE': None, 'ND1_NZ': None, 'NE2_CE': None, 'NE2_NZ': None, 'NE2_ND2': None, 'OG_ND1': None, 'OG_SG': None, 'OG1_SG': None, 'OE1_OH': None, 'OE1_CE': None, 'OE1_NZ': None, 'OD2_OE1': None, 'CE2_SG': None, 'CZ_SG': None, 'CG1_CZ3': None, 'CG1_CH2': None, 'NE1_O': None, 'NE1_CB': None, 'NE1_CG': None, 'NE1_CD1': None, 'CE3_C': None, 'CE3_O': None, 'CE3_CB': None, 'CE3_CD1': None, 'CZ2_CB': None, 'CZ2_CD1': None, 'CZ3_CB': None, 'CZ3_CD1': None, 'CH2_CB': None, 'CH2_CD1': None, 'NE1_C': None, 'CE3_N': None, 'CE3_CA': None, 'CE3_CG': None, 'CZ2_C': None, 'CZ2_O': None, 'CZ3_CA': None, 'CZ3_C': None, 'CZ3_O': None, 'CZ3_CG': None, 'CH2_C': None, 'CH2_O': None, 'NE1_N': None, 'NE1_CA': None, 'CE3_CD': None, 'CZ2_N': None, 'CZ2_CA': None, 'CZ3_N': None, 'CZ3_CD': None, 'CH2_N': None, 'CH2_CA': None, 'NE1_CD2': None, 'NE1_CE1': None, 'NE1_CE2': None, 'NE1_CZ': None, 'CE3_CD2': None, 'CE3_CE1': None, 'CE3_CE2': None, 'CE3_CZ': None, 'CE3_OH': None, 'CZ2_CG': None, 'CZ2_CD2': None, 'CZ2_CE1': None, 'CZ2_CE2': None, 'CZ3_CD2': None, 'CZ3_CE1': None, 'CZ3_CE2': None, 'CZ3_CZ': None, 'CZ3_OH': None, 'CH2_CG': None, 'CH2_CD2': None, 'CH2_CE1': None, 'CH2_CE2': None, 'CH2_CZ': None, 'NE1_OG1': None, 'NE1_CG2': None, 'CE3_OG1': None, 'CE3_CG2': None, 'CZ2_OG1': None, 'CZ2_CG2': None, 'CZ3_OG1': None, 'CZ3_CG2': None, 'CH2_OG1': None, 'CH2_CG2': None, 'CE3_OD1': None, 'CE3_ND2': None, 'CZ3_OD1': None, 'CZ3_ND2': None, 'CZ2_CZ': None, 'NE1_OD1': None, 'CZ2_OD1': None, 'CH2_OD1': None, 'NE1_CG1': None, 'CE3_CG1': None, 'CZ2_CG1': None, 'CZ3_CG1': None, 'CH2_CG1': None, 'NE1_NE1': None, 'NE1_CE3': None, 'NE1_CZ2': None, 'NE1_CZ3': None, 'NE1_CH2': None, 'CE3_NE1': None, 'CE3_CE3': None, 'CE3_CZ2': None, 'CE3_CZ3': None, 'CE3_CH2': None, 'CZ2_NE1': None, 'CZ2_CE3': None, 'CZ2_CZ2': None, 'CZ2_CZ3': None, 'CZ2_CH2': None, 'CZ3_NE1': None, 'CZ3_CE3': None, 'CZ3_CZ2': None, 'CZ3_CZ3': None, 'CZ3_CH2': None, 'CH2_NE1': None, 'CH2_CE3': None, 'CH2_CZ2': None, 'CH2_CZ3': None, 'CH2_CH2': None, 'NE1_ND1': None, 'NE1_NE2': None, 'CE3_ND1': None, 'CE3_NE2': None, 'CZ2_ND1': None, 'CZ2_NE2': None, 'CZ3_ND1': None, 'CZ3_NE2': None, 'CH2_ND1': None, 'CH2_NE2': None, 'NE1_CD': None, 'NE1_NE': None, 'NE1_NH2': None, 'CE3_NE': None, 'CE3_NH2': None, 'CZ2_CD': None, 'CZ2_NE': None, 'CZ3_NE': None, 'CH2_NE': None, 'NE1_OD2': None, 'CE3_OD2': None, 'CZ2_OD2': None, 'CZ3_OD2': None, 'CH2_OD2': None, 'NE1_OG': None, 'NE1_NH1': None, 'CE3_NH1': None, 'CZ2_NH1': None, 'CZ2_NH2': None, 'CZ3_NH1': None, 'CZ3_NH2': None, 'CH2_CD': None, 'CH2_NH1': None, 'CH2_NH2': None, 'NE1_OH': None, 'CZ2_OH': None, 'CH2_OH': None, 'ND1_NE1': None, 'ND1_CE3': None, 'ND1_CZ2': None, 'ND1_CZ3': None, 'ND1_CH2': None, 'NE2_NE1': None, 'NE2_CE3': None, 'NE2_CZ2': None, 'NE2_CZ3': None, 'NE2_CH2': None, 'ND1_OE1': None, 'CE_OG1': None, 'NZ_OG1': None, 'CE_NE': None, 'CE_NH1': None, 'CE_NH2': None, 'NZ_NE': None, 'NZ_NH1': None, 'NZ_NH2': None, 'CE_ND1': None, 'CE_NE2': None, 'NZ_ND1': None, 'NZ_NE2': None, 'CE_OG': None, 'NZ_OG': None, 'CE_NE1': None, 'CE_CE3': None, 'CE_CZ2': None, 'CE_CZ3': None, 'CE_CH2': None, 'NZ_NE1': None, 'NZ_CE3': None, 'NZ_CZ2': None, 'NZ_CZ3': None, 'NZ_CH2': None, 'OD2_NE1': None, 'OD2_CE3': None, 'OD2_CZ2': None, 'OD2_CZ3': None, 'OD2_CH2': None, 'OD2_ND2': None, 'ND2_CZ2': None, 'ND2_CH2': None, 'OH_NH1': None, 'OH_CZ2': None, 'OG_NE1': None, 'OG_CE3': None, 'OG_CZ2': None, 'OG_CZ3': None, 'OG_CH2': None, 'OG_SD': None, 'OE2_CA': None, 'OE2_C': None, 'OE2_O': None, 'OE2_CB': None, 'OE2_N': None, 'OE2_CG': None, 'OE2_CD1': None, 'OE2_CD2': None, 'OE2_CE1': None, 'OE2_CE2': None, 'OE2_CZ': None, 'OE2_OG': None, 'OE2_OG1': None, 'OE2_CG2': None, 'OE2_CD': None, 'OE2_OE1': None, 'OE2_OE2': None, 'OE2_CE': None, 'OE2_NZ': None, 'OE1_ND2': None, 'OE2_OD1': None, 'OE2_ND2': None, 'OE2_CG1': None, 'OE2_NE': None, 'OE2_NH1': None, 'OE2_NH2': None, 'OE1_NE1': None, 'OE1_CE3': None, 'OE1_CZ2': None, 'OE1_CZ3': None, 'OE1_CH2': None, 'OE2_NE1': None, 'OE2_CE3': None, 'OE2_CZ2': None, 'OE2_CZ3': None, 'OE2_CH2': None, 'CA_SD': None, 'C_SD': None, 'O_SD': None, 'CB_SD': None, 'CG_SD': None, 'CD_SD': None, 'OE1_SD': None, 'OE2_SD': None, 'OE1_ND1': None, 'OE2_ND1': None, 'OE2_NE2': None, 'OE2_OH': None, 'CE_ND2': None, 'NZ_ND2': None, 'N_SD': None, 'OD1_OE2': None, 'ND2_OE2': None, 'ND2_CE': None, 'ND2_NZ': None, 'OD1_SD': None, 'ND2_SD': None, 'ND2_ND1': None, 'CG1_SD': None, 'CG2_SD': None, 'CD1_SD': None, 'NH2_SD': None, 'NE_SD': None, 'CZ_SD': None, 'NH1_SD': None, 'NE1_ND2': None, 'CE3_OG': None, 'CZ2_OG': None, 'CZ3_OG': None, 'CH2_OG': None, 'NE1_OE1': None, 'NE1_OE2': None, 'CE3_OE1': None, 'CE3_OE2': None, 'CZ2_OE1': None, 'CZ2_OE2': None, 'CZ3_OE1': None, 'CZ3_OE2': None, 'CH2_OE1': None, 'CH2_OE2': None, 'CZ2_ND2': None, 'CH2_ND2': None, 'CZ3_CE': None, 'CH2_CE': None, 'OD2_OE2': None, 'OG1_SD': None, 'CE_SG': None, 'NZ_SG': None, 'SG_CE': None, 'SG_NZ': None, 'SG_OG1': None, 'OE2_OD2': None, 'OE2_SG': None, 'OD2_SG': None, 'CD2_SD': None, 'OH_NE1': None, 'OH_CH2': None, 'CE1_SD': None, 'CE2_SD': None, 'OH_SD': None, 'ND1_ND2': None, 'ND1_SD': None, 'NE2_SD': None, 'NE1_NZ': None, 'CE3_NZ': None, 'CZ2_CE': None, 'CZ2_NZ': None, 'CZ3_NZ': None, 'CH2_NZ': None, 'NE1_CE': None, 'CE3_CE': None, 'NE1_SD': None, 'CE3_SD': None, 'CZ2_SD': None, 'CZ3_SD': None, 'SD_NH2': None, 'SD_N': None, 'SD_CA': None, 'SD_C': None, 'SD_CB': None, 'SD_CG2': None, 'SD_O': None, 'SD_CG1': None, 'SD_CG': None, 'SD_CD1': None, 'SD_CD2': None, 'SD_CE1': None, 'SD_CE2': None, 'SD_CZ': None, 'SD_OH': None, 'SD_ND1': None, 'SD_NE2': None, 'SD_OD1': None, 'SD_ND2': None, 'SD_OG': None, 'SD_NE1': None, 'SD_CE3': None, 'SD_CZ2': None, 'SD_CZ3': None, 'SD_SD': None, 'SD_CE': None, 'CE_SD': None, 'SD_CD': None, 'SD_OE1': None, 'SD_OE2': None, 'NE_SG': None, 'NH1_SG': None, 'NH2_SG': None, 'OH_SG': None, 'SD_NE': None, 'SD_NH1': None, 'SD_OG1': None, 'SD_OD2': None, 'SD_NZ': None, 'OD2_SD': None, 'NZ_SD': None, 'ND1_OE2': None, 'SG_ND1': None, 'SG_NE1': None, 'SG_CE3': None, 'SG_CZ2': None, 'SG_CZ3': None, 'SG_CH2': None, 'NE1_SG': None, 'CE3_SG': None, 'CZ2_SG': None, 'CZ3_SG': None, 'CH2_SG': None, 'SG_SD': None, 'SD_SG': None, 'CH2_SD': None, 'SD_CH2': None}
unique_atm_combi_equal=[]
for a in unique_atm_combi_dict.keys():
    if a.split("_")[0] == a.split("_")[1]:
        unique_atm_combi_equal.append(a)
unique_atm_combi_equal_dict = {'N_N': None, 'CA_CA': None, 'C_C': None, 'O_O': None, 'CB_CB': None, 'CG_CG': None, 'CD_CD': None, 'OE1_OE1': None, 'NE2_NE2': None, 'SG_SG': None, 'CG1_CG1': None, 'CG2_CG2': None, 'OD1_OD1': None, 'ND2_ND2': None, 'CD1_CD1': None, 'CD2_CD2': None, 'OG1_OG1': None, 'NE_NE': None, 'CZ_CZ': None, 'NH1_NH1': None, 'NH2_NH2': None, 'CE1_CE1': None, 'CE2_CE2': None, 'OH_OH': None, 'OG_OG': None, 'OD2_OD2': None, 'CE_CE': None, 'NZ_NZ': None, 'ND1_ND1': None, 'NE1_NE1': None, 'CE3_CE3': None, 'CZ2_CZ2': None, 'CZ3_CZ3': None, 'CH2_CH2': None, 'OE2_OE2': None, 'SD_SD': None}
"""
#####################################################################################################################
assist functions
#####################################################################################################################
"""

"""
input: pdb_file_path

output: the atom features of a pbd file
lists of list that represents various properties of each of the atoms 
in the protein structure file. For each atom eight features are stored as follows:
0 ATOM_serial Int The serial number of the atom in the pdb
1 ATOM_symbol Str The name of the atom in abbreviated manner
2 Res_number Int Residue sequence number
3 RES_name Str Abbreviated residue name, ARG for Argenine,GLY for Glycine ..
4 Chain_id Str The chain identifier. Some protein complexes composed of several chains
5 X Float Orthogonal coordinates for X in Angstroms
6 Y Float Orthogonal coordinates for Y in Angstroms
7 Z Float Orthogonal coordinates for Z in Angstroms
"""
def get_atom_features(pdb_file_path):
    if not pdb_file_path.endswith(".pdb"):
        pdb_file_path = pdb_file_path + ".pdb"
    features_matrix=[]#define empty tuple
    IN=open(pdb_file_path,'r')
    for line in IN.readlines():
        one_line_vector=()
        if line.startswith('ATOM'):
            #0123456789012345678901234567890123456789012345678901234567890123456789
            #ATOM      1  N   LEU A  62      55.829  64.222  83.106  1.00 62.37           N
            atom_serial=int(line[6:11].replace(" ",""))
            atom_simbol=line[12:16].replace(" ","")

            res_number = int(line[22:26].replace(" ", ""))
            res_name=line[17:20].replace(" ","")

            chain_id = line[21:22].replace(" ", "")

            x = float(line[30:38].replace(" ", ""))
            y = float(line[38:46].replace(" ", ""))
            z = float(line[46:54].replace(" ", ""))
            one_line_vector = (atom_serial,atom_simbol,res_number,res_name,chain_id,x,y,z)  # define atom features tuple
            #print one_line_vector
            features_matrix.append(one_line_vector)
    features_matrix=tuple(features_matrix)
    IN.close()
    return features_matrix

def get_atom_features_db(pdb_file_path):
    """
    input:  the pdb file path
    output: the atom features of a pbd file
    lists of list that represents various properties of each of the atoms
    in the protein structure file. For each atom eight features are stored as follows:
    0 ATOM_serial Int The serial number of the atom in the pdb
    1 ATOM_symbol Str The name of the atom in abbreviated manner
    2 Res_number Int Residue sequence number
    3 RES_name Str Abbreviated residue name, ARG for Argenine,GLY for Glycine ..
    4 Chain_id Str The chain identifier. Some protein complexes composed of several chains
    5 X Float Orthogonal coordinates for X in Angstroms
    6 Y Float Orthogonal coordinates for Y in Angstroms
    7 Z Float Orthogonal coordinates for Z in Angstroms
    """
    features_matrix=[]#define empty tuple
    IN=open(pdb_file_path,'r')
    for line in IN.readlines():
        one_line_vector=()
        if line.startswith('ATOM'):
            #0123456789012345678901234567890123456789012345678901234567890123456789
            #ATOM      1  N   LEU A  62      55.829  64.222  83.106  1.00 62.37           N
            atom_serial=int(line[6:11].replace(" ",""))
            atom_simbol=line[12:16].replace(" ","")

            res_name=line[17:20].replace(" ","")
            chain_id = line[21:22].replace(" ", "")

            res_number = int(line[22:26].replace(" ", ""))

            x = float(line[30:38].replace(" ", ""))
            y = float(line[38:46].replace(" ", ""))
            z = float(line[46:54].replace(" ", ""))

            condormation_rate = float(line[55:61].replace(" ", ""))
            try:
                b_factor =  "{:.2f}".format(float(line[61:67].replace(" ", "")))
            except ValueError:
                b_factor = None
            one_line_vector = (atom_serial,atom_simbol,res_name,res_number,chain_id,x,y,z,condormation_rate,b_factor)  # define atom features tuple
            #print one_line_vector
            features_matrix.append(one_line_vector)
    features_matrix=tuple(features_matrix)
    IN.close()
    return pd.DataFrame(features_matrix,columns = ["ATOM_serial","ATOM_symbol","RES_name","Res_number","Chain_id","X","Y","Z","condormation_rate","b_factor"])

"""
input:
1. all_atom_features - atom features in tuple of tuples (e.g.:atom_features_1jyf that was defined in 7.1.a)
2. Chain_id – a string
3. Res_number – an int
4. ATOM_symbol – a string 

output:
tuple which its items are the values of x,y and z coordinates of the specified 
atom. 
"""
def get_xyz_from_atom_features(all_atom_features,Chain_id,Res_number,ATOM_symbol):
    for tup in all_atom_features:
        if tup[4] == Chain_id and tup[2] == Res_number and tup[1] == ATOM_symbol:
            return (tup[5],tup[6],tup[7])



def calculate_distance(atom1_coor,atom2_coor,Distance_cutoff):
    """
    input:
    1- atom1_coor- Tuple with the X,Y,Z coordinates for atom1
    2- atom2_coor-Tuple with the X,Y,Z coordinates for atom2
    3- Distance_cutoff- The distance cutoff in Angstrom(A0), a float value, to determine whether the two atoms interact with each other.

    output:
    two values (what it really returns is a tuple with two values) :
    1- The calculated distance (a float value)
    2- A Boolean value: True if the distance between the two atoms is within the accepted distance cutoff. Else returns False
    """
    dis = ((atom1_coor[0]-atom2_coor[0])**2+(atom1_coor[1]-atom2_coor[1])**2+(atom1_coor[2]-atom2_coor[2])**2)**0.5
    flag =True
    if dis < Distance_cutoff:
        flag = True
    else:
        flag = False
    return (dis,flag)

def get_res_dist(pdb_file_path1,chain1,pdb_file_path2,chain2,how="all"):
    atm_db1 = get_atom_features_db(pdb_file_path1)
    atm_db1 = atm_db1.loc[atm_db1["Chain_id"] == chain1]
    atm_db2 = get_atom_features_db(pdb_file_path2)
    atm_db2 = atm_db2.loc[atm_db2["Chain_id"].isin(chain2)]
    ress1 = list(atm_db1["Res_number"].unique())
    ress2 = list(atm_db2["Res_number"].unique())
    output = list(product(ress1, ress2))
    res_dist = pd.DataFrame(output,columns=["res1","res2"])
    counter=1
    index_of_combi = 0
    percent=0
    list_of_end = []
    last_res = atm_db2.at[0,"Res_number"]
    unique_atm_combi_equal_dict = {'N_N': None, 'CA_CA': None, 'C_C': None, 'O_O': None, 'CB_CB': None, 'CG_CG': None, 'CD_CD': None, 'OE1_OE1': None, 'NE2_NE2': None, 'SG_SG': None, 'CG1_CG1': None, 'CG2_CG2': None, 'OD1_OD1': None, 'ND2_ND2': None, 'CD1_CD1': None, 'CD2_CD2': None, 'OG1_OG1': None, 'NE_NE': None, 'CZ_CZ': None, 'NH1_NH1': None, 'NH2_NH2': None, 'CE1_CE1': None, 'CE2_CE2': None, 'OH_OH': None, 'OG_OG': None, 'OD2_OD2': None, 'CE_CE': None, 'NZ_NZ': None, 'ND1_ND1': None, 'NE1_NE1': None, 'CE3_CE3': None, 'CZ2_CZ2': None, 'CZ3_CZ3': None, 'CH2_CH2': None, 'OE2_OE2': None, 'SD_SD': None}
    for atm_db1_index in atm_db1.index:
        for atm_db2_index in atm_db1.index:
            rate_progres = counter/res_dist.shape[0]
            if int(rate_progres)>percent:
                print(str(int(rate_progres))+"%")
                percent = int(rate_progres)
            counter += 1
            atm1 = atm_db1.at[atm_db1_index,"ATOM_symbol"]
            atm2 = atm_db2.at[atm_db2_index,"ATOM_symbol"]
            if atm2!=atm1:
                continue
            res1 = atm_db1.at[atm_db2_index,"Res_number"]
            res2 = atm_db2.at[atm_db2_index,"Res_number"]
            atom1_coor = tuple(atm_db1.loc[atm_db1_index,"X":"Z"])
            atom2_coor = tuple(atm_db2.loc[atm_db2_index,"X":"Z"])
            #atom2_coor
            dist,flag = calculate_distance(atom1_coor, atom2_coor, 12)

            if last_res != res2:
                list_of_end.append([res1,last_res]+list(unique_atm_combi_equal_dict.values()))
                index_of_combi+=1
                last_res = res2
                unique_atm_combi_equal_dict = {'N_N': None, 'CA_CA': None, 'C_C': None, 'O_O': None, 'CB_CB': None, 'CG_CG': None, 'CD_CD': None, 'OE1_OE1': None, 'NE2_NE2': None, 'SG_SG': None, 'CG1_CG1': None, 'CG2_CG2': None, 'OD1_OD1': None, 'ND2_ND2': None, 'CD1_CD1': None, 'CD2_CD2': None, 'OG1_OG1': None, 'NE_NE': None, 'CZ_CZ': None, 'NH1_NH1': None, 'NH2_NH2': None, 'CE1_CE1': None, 'CE2_CE2': None, 'OH_OH': None, 'OG_OG': None, 'OD2_OD2': None, 'CE_CE': None, 'NZ_NZ': None, 'ND1_ND1': None, 'NE1_NE1': None, 'CE3_CE3': None, 'CZ2_CZ2': None, 'CZ3_CZ3': None, 'CH2_CH2': None, 'OE2_OE2': None, 'SD_SD': None}
            if flag == True:
                unique_atm_combi_equal_dict[atm1+"_"+atm2] = dist
                #res_dist.loc[index_of_combi,atm1+"_"+atm2] = dist

    res_dist = pd.DataFrame(list_of_end,columns=["res1","last_res"]+list(unique_atm_combi_equal_dict.keys()))
    return res_dist

def get_coordinates_of_chain(atom_features,chain_ID):
    """
    input:
     atom_features tuple and a chain ID

    output:
    dictionary of dictionaries (D_O_D) where the outer keys are the tuples (Res_number,
    amino_acid_type) and the values are dictionaries where the keys are the atoms names
    and the values are tuples with the x,y,z coordinates of the atom from the given chain.
    """
    end_dict_of_dict = {}
    for tup in atom_features:
        if tup[4] != chain_ID:
            continue
        outer_key = (tup[2],tup[3])
        inner_key = tup[1]
        if(outer_key in end_dict_of_dict):
            end_dict_of_dict[outer_key][inner_key] = (tup[5],tup[6],tup[7])
        else:
            end_dict_of_dict[outer_key]={inner_key:(tup[5],tup[6],tup[7])}
    return end_dict_of_dict


"""
input:
1. L_res1- list of tuples [(xatom1,yatom1,zatom1),( xatom2,yatom2,zatom2)…] holding all of the
atoms coordinates of a given residue.
2. L_res2- list of tuples holding the atoms coordinate of additional residue
3. distance- distance parameter

output:
True if it finds a pair of atoms (each different residue) that are in contact or return
 False otherwise.
"""
def check_if_two_res_contact (L_res1,L_res2,distance):
    flag = False
    for dot1 in L_res1:
        for dot2 in L_res2:
            flag = flag or calculate_distance(dot1,dot2,distance)[1]
    return flag


"""
input:
1. DoD_chain1- a dictionary of dictionaries holding the atoms coordinate for each residue from chain1
2. DoD_chain2- a dictionary of dictionaries holding the atoms coordinate for each residue from chain2
3. def_dist - Interaction distance parameter

output:
dictionary with the contacting residues, where:
- The keys are the contact residues from chain 1 (residues not in contact should not be in the keys). The contact residues are represented by tuples (RES_NUMBER, RES_TYPE)
- The value is list of tuples (RES_NUMBER, RES_TYPE) representing all the contact residues pairs from chain2.
For example: if TYR number 122 from chain 1 interacts with ARG13,HIS27 and TRP20 from chain 2, it will
be stored as:
contact_dic[(122,’TYR’)]=[(13,’ARG’),(‘27,HIS’), (‘20,TRP’)]

"""
def get_contacts_between_2_chains(D_of_D_chain1,D_of_D_chain2,def_dist):
    contact_dict = {}
    for res1,coor_of_res1 in D_of_D_chain1.items():
        temp_contact_list = []
        for res2,coor_of_res2 in D_of_D_chain2.items():
            if check_if_two_res_contact(list(coor_of_res1.values()),list(coor_of_res2.values()),def_dist):
                temp_contact_list.append(res2)

        if not temp_contact_list == []:
                contact_dict[res1] = temp_contact_list
    return contact_dict


def get_contacts_distance_between_2_chains(D_of_D_chain1,D_of_D_chain2,def_dist):
    contact_dict = {}
    temp_contact_list = []
    for res1,coor_of_res1 in D_of_D_chain1.items():
        # temp_contact_list = []
        for res2,coor_of_res2 in D_of_D_chain2.items():
            for atm_sym1,dot1 in coor_of_res1.items():
                for atm_sym2,dot2 in coor_of_res2.items():
                    if atm_sym1 == atm_sym2 == "CA":
                    # if flag==True:
                    #     print(dist)
                    #     residue_db.append([res1[0],res1[1],res2[0],res2[1],atm_sym1+"_"+atm_sym2,dist])
                    #temp_contact_list.append((res2,atm_sym1+"_"+atm_sym2,calculate_distance(dot1, dot2, def_dist)[0]))
                        dist,flag=calculate_distance(dot1, dot2, def_dist)
                        temp_contact_list.append([res1[0],res1[1],res2[0],res2[1],atm_sym1+"_"+atm_sym2,dist])
        # if not temp_contact_list == []:
        #         contact_dict[res1] = temp_contact_list
    residue_db = pd.DataFrame(temp_contact_list,columns=["res1_position","res1_symbol","res2_position","res2_symbol","atm_sym1_atm_sym2","distance"])
    return contact_dict,residue_db,temp_contact_list

# starttime = time.time()
# atm = get_atom_features("/home/eilay/Downloads/7a98_open.pdb")
# dod_A = get_coordinates_of_chain(atm,"A")
# dod_B = get_coordinates_of_chain(atm,"B")
# dod_C = get_coordinates_of_chain(atm,"C")
#
# contact_dict_A_A,residue_db_A_A,temp_contact_list_A_A = get_contacts_distance_between_2_chains(dod_A,dod_A,12)
# temp_contact_db_A_A = pd.DataFrame(temp_contact_list_A_A)
# #temp_contact_db_A_A.to_csv("/home/eilay/Downloads/A_A.tsv","\t")
#
# contact_dict_A_B,residue_db_A_B,temp_contact_list_A_B  = get_contacts_distance_between_2_chains(dod_A,dod_B,12)
# temp_contact_db_A_B = pd.DataFrame(temp_contact_list_A_B)
# #temp_contact_db_A_A.to_csv("/home/eilay/Downloads/A_B.tsv","\t")
#
# contact_dict_A_C,residue_db_A_C,temp_contact_list_A_C = get_contacts_distance_between_2_chains(dod_A,dod_C,12)
# temp_contact_db_A_C = pd.DataFrame(temp_contact_list_A_C)
# #temp_contact_db_A_A.to_csv("/home/eilay/Downloads/A_C.tsv","\t")
# temp_contact_db_A_A["chain"] ="A_A"
# temp_contact_db_A_B["chain"] ="A_B"
# temp_contact_db_A_C["chain"] ="A_C"
# final = pd.concat([temp_contact_db_A_A,temp_contact_db_A_B,temp_contact_db_A_C])
# final.columns =["res1_position","res1_symbol","res2_position","res2_symbol","atm_sym1_atm_sym2","distance","chain"]
# final.to_csv("/home/eilay/Downloads/open_confermation.tsv","\t")
#

"""
The function gets a string of amino acid name in one(‘V’) or in three letters(‘VAL’) and returns a string of an
id of amino acids in three letters(‘VAL’) or in one letter(‘V’), respectively
"""
def get_3id_to_1id_or_1id_from_3id(amino_str):
    dict_amino_3to1 ={'CYS': 'C', 'ASP': 'D', 'SER': 'S', 'GLN': 'Q', 'LYS': 'K',
     'ILE': 'I', 'PRO': 'P', 'THR': 'T', 'PHE': 'F', 'ASN': 'N',
     'GLY': 'G', 'HIS': 'H', 'LEU': 'L', 'ARG': 'R', 'TRP': 'W',
     'ALA': 'A', 'VAL':'V', 'GLU': 'E', 'TYR': 'Y', 'MET': 'M'}

    if len(amino_str) == 3:
         return dict_amino_3to1[amino_str]
    elif len(amino_str) == 1:
        for key, value in dict_amino_3to1.items():
            if value == amino_str:
                return key
    else:
        raise ValueError("not 3 or 1 len")

"""
input:
DoD_chain1- a dictionary of dictionaries holding the atoms coordinate for each residue from chain1
output:
tupel:
1 the seq of the DOD of the ordered_residues
2. the first index of the ordered_residues
3. the last place of the ordered_residues
"""
def get_seq_and_first_index_and_last_place(D_of_D_chain):

    keys = list(D_of_D_chain.keys())
    resideues = ""
    first_index = -8
    last_place = -8
    for place , res in keys:
        last_place = place
        if first_index == -8:
            first_index = place-1
        if place < 0:
            return ("neg_place_in_struc",-1,-1)
        resideues += get_3id_to_1id_or_1id_from_3id(res)
    return (resideues,first_index,last_place)



def seq(D_of_D_chain="",pdb_file_path="",chain =""):
    """
    input:
    DoD_chain1- a dictionary of dictionaries holding the atoms coordinate for each residue from chain1
    or:
    pdb_file_path - path of pdb
    chain - chain of pdb
    output:
    tupel:
    1 the seq(of orderd residued) of the DOD of the ordered_residues
    """
    if D_of_D_chain != "":
        assert chain =="" and pdb_file_path ==""
        keys = list(D_of_D_chain.keys())
        resideues = ""
        privius_place = keys[0][0]
        counter = 0
        for place , res  in keys :
            if counter==0:
                try:
                    resideues += "X"* (place-1)# fixing aligmnent on the beggining
                except TypeError:
                    resideues += "X"* (int(place[:-1])-1)# fixing 1A
                    print("i use "+ str(place[:-1])+ "insted" + str(place)  )
            jump = place - privius_place
            resideues += "X"*(jump-1) + get_3id_to_1id_or_1id_from_3id(res)
            privius_place = place
            counter+=1
    else:
        assert  chain != "" and pdb_file_path != ""
        atm = get_atom_features(pdb_file_path.rstrip(".pdb"))
        D_of_D_chain = get_coordinates_of_chain(atm,chain)
        keys = list(D_of_D_chain.keys())
        resideues = ""
        privius_place = keys[0][0]
        counter = 0
        for place , res  in keys :
            if counter==0:
                try:
                    resideues += "X"* (place-1)# fixing aligmnent on the beggining
                except TypeError:
                    resideues += "X"* (int(place[:-1])-1)# fixing aligmnent on the beggining
                    print("i use "+ str(place[0])+ "insted" + str(place)  )
            jump = place - privius_place
            resideues += "X"*(jump-1) + get_3id_to_1id_or_1id_from_3id(res)
            privius_place = place
            counter+=1

    return resideues


def get_chains(pdb_file_path):
    """

    :param pdb_file_path: a file path of a pdb
    :return: list of the chain_id (exp. ["A","L","H","K"]
    """
    parser = PDBParser()
    structure = parser.get_structure("", pdb_file_path)
    list_of_chains = []
    for i in range(len(list(structure.get_chains()))):
        chain_id = list(structure.get_chains())[i].get_id()
        list_of_chains.append(chain_id)
    return list_of_chains


"""
input:
DoD_chain1- a dictionary of dictionaries holding the atoms coordinate for each residue from chain1
chain_id - the chain you want the check will be
pdb_name - exp. 6XDG the four letter from the pdb data

output:
1 -1 if its not a protein at all(like chemistry molecoles arond)
2 "neg_place_in_struc" if their is any negetive place in the sturctcher(bad data)
3 "antigen" if it is antigen
4 "heavy" if it heavy chain of antibody
5 "light" if it light chain of antibody
6 "unknown" if it is not in sabdab
"""
def what_is_the_chain(D_of_D_chain,chain_id,pdb_name,sabdab_data_file_path=the_sabdab_file_path):
    resideues,first_index,last_place = get_seq_and_first_index_and_last_place(D_of_D_chain)
    if D_of_D_chain == {}:
        return -1
    if resideues =="neg_place_in_struc":
        return "neg_place_in_struc"
    try:
        sabdab_db = pd.read_csv(sabdab_data_file_path,"\t")
    except:
        raise ("cheak your sabdab_database in line 7,8")

    is_this_antigen = chain_id in str(list(sabdab_db[(sabdab_db.pdb == pdb_name)].antigen_chain))
    is_this_Hchain = chain_id in str(list(sabdab_db[(sabdab_db.pdb == pdb_name)].Hchain))
    is_this_Lchain = chain_id in str(list(sabdab_db[(sabdab_db.pdb == pdb_name)].Lchain))

    if is_this_antigen :
        return "antigen"
    elif is_this_Hchain:
        return "heavy"
    elif is_this_Lchain:
        return "light"

    return "unknown"



def delete_chain(pdb_file_path,chain_to_delete,clean):
    """
    :param pdb_file_path: file path of pdb
    :param chain_to_delete: the chain to delete exp: "A"
    :param clean: True- write over the original
                  False - make anthor file with the smae name + .myout
    :return:
    """
    pdb_io = PDB.PDBIO()
    pdb_parser = PDB.PDBParser()
    structure = pdb_parser.get_structure(" ", pdb_file_path)

    for model in structure:
        model.detach_child(chain_to_delete)
    pdb_io.set_structure(structure)
    if clean:
        pdb_io.save(pdb_file_path)
    else:
        pdb_io.save(pdb_file_path+".myout")




def make_chain(pdb_file_path,chains):
    """

    :param pdb_file_path:
    :param chains: ???????????????????????
    :return:
    """
    pdb_io = PDB.PDBIO()
    pdb_parser = PDB.PDBParser()
    structure = pdb_parser.get_structure(" ", pdb_file_path)
    for model in structure:
        for chain in chains:
            model.add(chain)
    pdb_io.set_structure(structure)
    pdb_io.save(pdb_file_path+".myout")


def change_shifting_of_pdb(shift_right,pdb_file_path,chain_to_shift,clean=False):
    """
    :param shift_right: how many index to shift the indexing of chain
    :param pdb_file_path:
    :param chain_to_shift: which chain to shift
    :param clean: True- write over the original
                  False - make anthor file with the smae name + .myout
    :return:
    """
    pdb_io = PDB.PDBIO()
    pdb_parser = PDB.PDBParser()
    structure = pdb_parser.get_structure(" ", pdb_file_path)

    for model in structure:
        for chain in model:
            if chain.id == chain_to_shift:
                for i, residue in enumerate(chain.get_residues()):
                    res_id = list(residue.id)
                    res_id[1] = int(res_id[1]) + shift_right
                    residue.id = tuple(res_id)
    pdb_io.set_structure(structure)
    if clean:
        pdb_io.save(pdb_file_path)
    else:
        pdb_io.save(pdb_file_path + ".myout")

def remain_only_those_resi(pdb_file_path,resi_to_keep,clean = False,chains_to_manipulate="all",chains_to_skip=[]):
    """
    :param pdb_file_path:
    :param resi_to_keep: list of position of resi to keep eg [1,2,3,4,5,6,7,11]
    :param clean: is to delete the original file
    :param chains_to_manipulate:
    :param chains_to_skip:
    :return:
    """
    OUT = open(pdb_file_path+".myout",'w')
    with open(pdb_file_path,"r")as IN:
        for line in IN.readlines():
            if line.startswith('ATOM'):
                l1 = line
                #0123456789012345678901234567890123456789012345678901234567890123456789
                #ATOM      1  N   LEU A  62      55.829  64.222  83.106  1.00 62.37           N
                res_number = int(line[22:26].replace(" ", ""))
                chain = line[21: 22].replace(" ", "")
                if (chain in chains_to_skip) or ((chains_to_manipulate == "all" or chain in chains_to_manipulate) and res_number in resi_to_keep):
                    OUT.write(line)
    OUT.close()
    if clean:
        os.remove(pdb_file_path)
        os.rename(pdb_file_path+".myout",pdb_file_path)


def change_first_zero_of_pdb_and_get_rid_anything_is_not_atom(pdb_file_path):
    """
    :param pdb_file_path:
    :return:  sometimes the first residiue is in index 0 and it make a lot of mess after it so you change it to 1A and
     it get rid from anithing is not ATOM symbool to save space note you will lost data like resulotion and more
    """
    IN=open(pdb_file_path+'.myout','w')
    with open(pdb_file_path,"r") as f:
        for line in f.readlines():
            if line.startswith('ATOM'):
                res_number = int(line[22:26].replace(" ", ""))
                if res_number == 0:
                    line = line[:25]+"1A"+line[27:]
                    IN.write(line)
                else:
                    IN.write(line)

def change_first_zero_of_pdb(pdb_file_path):
    """
    :param pdb_file_path:
    :return:  sometimes the first residiue is in index 0 and it make a lot of mess after it so you change it to 1A
    """
    IN=open(pdb_file_path+'.myout','w')
    with open(pdb_file_path,"r") as f:
        for line in f.readlines():
            if line.startswith('ATOM'):
                res_number = int(line[22:26].replace(" ", ""))
                if res_number == 0:
                    line = line[:25]+"1A"+line[27:]
                    IN.write(line)
                else:
                    IN.write(line)
            else:
                IN.write(line)
    IN.close()


def change_first_zero_of_pdb_and_get_rid_anything_is_not_atom_by_folder(pdb_path_folder):
    """
    :param pdb_path_folder:
    :return: same as above but to folder
    """
    for file in os.listdir(pdb_path_folder):
        filename = os.fsdecode(file)
        if filename.endswith(".pdb"):
            change_first_zero_of_pdb_and_get_rid_anything_is_not_atom(pdb_path_folder+"/"+filename)

def change_first_zero_of_pdb_by_folder(pdb_path_folder):
    """
    :param pdb_path_folder:
    :return: same as above but to folder
    """
    for file in os.listdir(pdb_path_folder):
        filename = os.fsdecode(file)
        if filename.endswith(".pdb"):
            change_first_zero_of_pdb(pdb_path_folder+"/"+filename)


def b_factor_text(pdb_file_path,list_of_new_b_factor=[],list_of_tupels_new_b_factor=[],a="",output_folder_path=""):
    """
    :param pdb_file_path:
    :param list_of_new_b_factor:  list of numbers that represnt each the b-factor to enter base on the place in the list
    :param list_of_tupels_new_b_factor:  list of tupels e.g [(18,100),...](witch 18 is the residue and 100 is the bfactor num)
     that represnt each the b-factor to enter base on the place in the list
    :parm a: is a str to add before .myout
    :parm output_folder_path: where to put the output
    :return: used to color the pdb on discovery stodiu >>> ctrl+d >>> color by isotropic displacement , display style CPK
    """
    assert list_of_new_b_factor == [] or list_of_tupels_new_b_factor == []

    flag = False
    if output_folder_path == "":
        OUT=open(pdb_file_path+a+".myout",'w')
    else:
        OUT = open(output_folder_path+"/" +a+pdb_file_path.split("/")[-1] + ".myout", 'w')
    with open(pdb_file_path,"r")as IN:
        for line in IN.readlines():
            if line.startswith('ATOM'):
                l1 = line
                #0123456789012345678901234567890123456789012345678901234567890123456789
                #ATOM      1  N   LEU A  62      55.829  64.222  83.106  1.00 62.37           N
                res_number = int(line[22:26].replace(" ", ""))
                if list_of_new_b_factor != []:
                    try:
                        b_factor_ceil_2 = "{:.2f}".format(list_of_new_b_factor[res_number-1])
                    except IndexError:
                        b_factor_ceil_2 = "{:.2f}".format(0)#for a shorter list
                        if not flag:
                            flag = True
                            if (input("you enter a shorter list than the len of the chain of your input pdb. if it ok enter 1 if not enter 0") == 0):
                                raise ValueError("shorter list than the len of the chain")
                    REVAH = 6 - len(b_factor_ceil_2)
                    bfactor = " " * REVAH + b_factor_ceil_2
                    if not(len (bfactor) == 6 ):
                        raise ValueError("fuck"+ str(res_number))
                    line = line[:60]+bfactor+line[66:]
                else:
                    found_it = False
                    for (res_to_set,b_factor_to_set) in list_of_tupels_new_b_factor:
                        if str(res_to_set) == str(res_number):
                            b_factor_ceil_2 = "{:.2f}".format(b_factor_to_set)
                            REVAH = 6 - len(b_factor_ceil_2)
                            bfactor = " " * REVAH + b_factor_ceil_2
                            if not(len (bfactor) == 6 ):
                                raise ValueError("fuck"+ str(res_number))
                            line = line[:60]+bfactor+line[66:]
                            found_it = True
                            break
                    if found_it == False:
                        b_factor_ceil_2 = "{:.2f}".format(0)
                        REVAH = 6 - len(b_factor_ceil_2)
                        bfactor = " " * REVAH + b_factor_ceil_2
                        if not (len(bfactor) == 6):
                            raise ValueError("fuck" + str(res_number))
                        line = line[:60] + bfactor + line[66:]
                        found_it = True

                OUT.write(line)
            else:
                OUT.write(line)
    OUT.close()

# def b_factor_text(pdb_file_path,list_of_new_b_factor=[],list_of_tupels_new_b_factor=[]):
#     """
#     :param pdb_file_path:
#     :param list_of_new_b_factor:  list of numbers that represnt each the b-factor to enter base on the place in the list
#     :param list_of_tupels_new_b_factor:  list of tupels e.g [(18,100),...](witch 18 is the residue and 100 is the bfactor num)
#      that represnt each the b-factor to enter base on the place in the list
#     :return: used to color the pdb on discovery stodiu >>> ctrl+d >>> color by isotropic displacement , display style CPK
#     """
#     assert list_of_new_b_factor == [] or list_of_tupels_new_b_factor == []
#     OUT=open(pdb_file_path+".myout",'w')
#     with open(pdb_file_path,"r")as IN:
#         for line in IN.readlines():
#             if line.startswith('ATOM'):
#                 l1 = line
#                 #0123456789012345678901234567890123456789012345678901234567890123456789
#                 #ATOM      1  N   LEU A  62      55.829  64.222  83.106  1.00 62.37           N
#                 res_number = int(line[22:26].replace(" ", ""))
#                 if list_of_new_b_factor != []:
#                     b_factor_ceil_2 = "{:.2f}".format(list_of_new_b_factor[res_number-1])
#                     REVAH = 6 - len(b_factor_ceil_2)
#                     bfactor = " " * REVAH + b_factor_ceil_2
#                     if not(len (bfactor) == 6 ):
#                         raise ValueError("fuck"+ res_number)
#                     line = line[:60]+bfactor+line[66:]
#                 else:
#                     for (res_to_set,b_factor_to_set) in list_of_tupels_new_b_factor:
#                         if res_to_set ==
#                     b_factor_ceil_2 = "{:.2f}".format(list_of_new_b_factor[res_number-1])
#                     REVAH = 6 - len(b_factor_ceil_2)
#                     bfactor = " " * REVAH + b_factor_ceil_2
#                     if not(len (bfactor) == 6 ):
#                         raise ValueError("fuck"+ res_number)
#                     line = line[:60]+bfactor+line[66:]
#                 OUT.write(line)
#     OUT.close()

def is_atom_lines_same(small_pdb_file_path,big_pdb_file_path):
    """
    :param small_pdb_file_path:
    :param big_pdb_file_path:
    :return: chack if tow pdbs have the same ATOM lines
    """
    small = open(small_pdb_file_path,"r")
    big = open(big_pdb_file_path,"r")
    counter = 0
    for line_big in big.readlines():
        if line_big.startswith("ATOM"):
            counter+=1
            line_small = small.readline()
            if not( line_small == line_big ):
                return False,counter
    small.close()
    big.close()
    return True


def resolotion(pdb_file_path):
    """
    :param pdb_file_path:
    :return: get resulotion
    """
    pdb_io = PDB.PDBIO()
    pdb_parser = PDB.PDBParser()
    structure = pdb_parser.get_structure(" ", pdb_file_path)
    return structure.header["resolution"]


def is_it_full(pdb_file_path,the_max_begining=90,the_min_ending=1000,the_min_len=700,the_len_of_the_protein=1273):
    """
    :param pdb_file_path:
    :param the_max_begining:
    :param the_min_ending:
    :param the_min_len:
    :param the_len_of_the_protein:
    :return:
    """
    pdb_parser = PDB.PDBParser()
    structure = pdb_parser.get_structure(" ", pdb_file_path)
    atm = get_atom_features(pdb_file_path.rstrip(".pdb"))
    filename=pdb_file_path.split("/")[-1]
    for i in range(len(list(structure.get_chains()))):
        chain_id = list(structure.get_chains())[i].get_id()
        d_o_d = get_coordinates_of_chain(atm,chain_id)
        is_it = what_is_the_chain(d_o_d,chain_id,filename.lstrip("PyIgClassifyOutput_").rstrip(".pdb"))
        if is_it == "antigen":
            seq,first,last = get_seq_and_first_index_and_last_place(d_o_d)
            if (first in range(0,the_max_begining) and last in range (the_min_ending,the_len_of_the_protein) and len(seq)>the_min_len):
                return True,first,last,len(seq)
            else :
                return False,first,last,len(seq)

"""
#####################################################################################################################
contacts makerrrrrrrrr
#####################################################################################################################
note:
thier is an assumption that every protein that is not spike is an antibodie i.e
the PDB should contain only small chimical molecols,antibosies, spike

0.85 is the trashhold i recomand of matching between spike WT to the spesific DOD that was given. no longer relevant since the last update.

input:
1. file_name like 'PyIgClassifyOutput_7K8X.pdb'
2. def_dist - Interaction distance parameter

output:
will write a txt file with the name: name of pdb_spike chain_antibody chain.txt 
for example: '6W41_C_H.txt' 6W41 is the pdb id ,C is the chain id of the spike,H is the chain ID of the antibodie
"""
def contacts_between_spike_and_antibodies(file_name,def_dist,input_directory,output_directory,sabdab_data_file_path=the_sabdab_file_path):
    atom_features = get_atom_features(input_directory + "/" + file_name.rstrip(".pdb"))

    parser = PDBParser()
    structure = parser.get_structure(file_name.rstrip(".pdb"), input_directory + "/" + file_name)

    spike_chain_dods = []
    antibody_chain_dods= []
    spike_chain_ids = [] #self check
    antibody_chain_ids= [] #self check
    is_brake = False
    for i in range(len(list(structure.get_chains()))):
        chain_id = list(structure.get_chains())[i].get_id()
        d_o_d = get_coordinates_of_chain(atom_features,chain_id)
        try:
            is_it = what_is_the_chain(d_o_d,chain_id,file_name.lstrip("PyIgClassifyOutput_").rstrip(".pdb"),sabdab_data_file_path)
        except UnboundLocalError:
            raise UnboundLocalError(file_name+chain_id)
        if(is_it=="antigen"):
            spike_chain_ids.append(chain_id)
            spike_chain_dods.append(d_o_d)
        elif(is_it=="heavy" or is_it=="light" ):
            antibody_chain_ids.append(chain_id)
            antibody_chain_dods.append(d_o_d)
        elif(is_it=="neg_place_in_struc"):
            continue

    for spike_chain_id,spike_chain_dod in zip(spike_chain_ids,spike_chain_dods):
        for antibody_chain_id,antibody_chain_dod in zip(antibody_chain_ids,antibody_chain_dods):
            contacts = get_contacts_between_2_chains(spike_chain_dod,antibody_chain_dod,def_dist)
            if not contacts == {}:
                pd.Series(contacts,name="antibody res").to_csv(output_directory+"//"+file_name.lstrip("PyIgClassifyOutput_").rstrip(".pdb")+"_"+spike_chain_id+"_"+antibody_chain_id+".txt","\t")




"""
''' this is the directory(of the folder) of your input pdb spike+antibodies structchers!!!!!!!!!!!!!!!!!!!!!!'''
inputt_directory = "E:/dowloads/python/tomer_lab/structural_bio/COMLEX/new_29.12_spike_ab_complex/spike_ab_complex/COV_2/PDB"       
''' this is the directory(of the folder choose) of your input pdb spike+antibodies structchers!!!!!!!!!!!!!!!!!!!!!!'''

outputt_directory = "E:/dowloads/python/tomer_lab/structural_bio/COMLEX/new_29.12_spike_ab_complex/spike_ab_complex/COV_2/contacts"
"""
#input PDB FOLDER
#output CONTACTS FOLDER
def make_contacts_folder(PDB_folder_path,contacts_folder_path):
    input_directory = PDB_folder_path
    output_directory = contacts_folder_path
    #example PyIgClassifyOutput_6XKP.pdb
    antibodies_spike_complex_names = []
    for file in os.listdir(input_directory):
         filename = os.fsdecode(file)
         if filename.endswith(".pdb"):
             antibodies_spike_complex_names.append(filename)
    #example 6XKP
    pdb_that_already_compiled = []
    for file in os.listdir(output_directory):
         filename = os.fsdecode(file)
         if len(filename) == 12:
             pdb_that_already_compiled.append(filename[0:4])

    for file_name in antibodies_spike_complex_names:
        binarry_already_compiled = False

        for already_compiled in pdb_that_already_compiled:
            if (file_name.endswith(already_compiled+".pdb")):
                binarry_already_compiled = True
                break

        if(binarry_already_compiled == False):
            contacts_between_spike_and_antibodies(file_name,5,input_directory,output_directory)



"""
#####################################################################################################################
analyzer
#####################################################################################################################
"""
"""
input:
1. file_name like 7K90_C_G(before).txt
output:
number of residue contacts
"""
def num_of_contacts(contacts_files_mame):
    num_of_contacts = 0
    with open(contacts_files_mame,"r") as f:
        for line in f.readlines():
            line.rstrip("\n")
            line_sep = line.split("\t")
            num_of_contacts += line_sep[2].count("(")
    return num_of_contacts

"""
input:
1. file_name like 7K90_C_G(before).txt
2. file name like 7K90_C_G(after).txt
output:
the delta between the number of residue contacts
"""
def difference_between_contacts_after_mutaion(before_mutaion_contacts,after_mutaion_contacts):
    return num_of_contacts(after_mutaion_contacts) - num_of_contacts(before_mutaion_contacts)


def make_pfm_from_pdbs(pdbs_folder_path,sabdab_data_file_path=the_sabdab_file_path):
    """
    :param pdbs_folder_path: the path of the pdb folder
    :return: pfm : https://en.wikipedia.org/wiki/Position_weight_matrix
    """

    PFM = pd.DataFrame(np.zeros([20,1273],dtype=int), index = ["A","C","D","E","F","G","H","I","K","L","M","N","P","Q","R","S","T","V","W","Y"])
    PFM_info = pd.DataFrame(np.zeros([20,1273],dtype=str), index = ["A","C","D","E","F","G","H","I","K","L","M","N","P","Q","R","S","T","V","W","Y"])
    pdbs_files_mames = []
    for file in os.listdir(pdbs_folder_path):
        filename = os.fsdecode(file)
        if filename.endswith(".pdb"):
            pdbs_files_mames.append(filename.rstrip(".pdb"))

    for filename in pdbs_files_mames:
        atom_features = get_atom_features(pdbs_folder_path + "/" + filename)

        parser = PDBParser()
        structure = parser.get_structure(filename, pdbs_folder_path + "/"+ filename+".pdb")

        spike_chain_dods = []
        spike_chain_ids = [] #self check
        is_brake = False
        for i in range(len(list(structure.get_chains()))):
            chain_id = list(structure.get_chains())[i].get_id()
            d_o_d = get_coordinates_of_chain(atom_features,chain_id)
            is_it = what_is_the_chain(d_o_d,chain_id,filename.lstrip("PyIgClassifyOutput_").rstrip(".pdb"),sabdab_data_file_path)
            if(is_it=="antigen"):
                spike_chain_ids.append(chain_id)
                spike_chain_dods.append(d_o_d)
            else:
                continue
        for d_o_d in spike_chain_dods:
            for (Res_number,amino_acid_type) in d_o_d.keys():
                PFM_info.loc[get_3id_to_1id_or_1id_from_3id(amino_acid_type),int(Res_number)-1] += (filename[19:23] +",")
                PFM.loc[get_3id_to_1id_or_1id_from_3id(amino_acid_type),int(Res_number)-1] += 1
    return (PFM.T,PFM_info.T)


def contact_metrix(contacts_folder_path):#work good only on non regardent data
    """
    :param contacts_folder_path: only non regardent contact folder
    :return: it will make contact metrix
    """
    contact_metrix = pd.DataFrame(np.zeros([20,1273],dtype=int), index = ["A","C","D","E","F","G","H","I","K","L","M","N","P","Q","R","S","T","V","W","Y"])
    contact_metrix_info = pd.DataFrame(np.zeros([20,1273],dtype=str), index = ["A","C","D","E","F","G","H","I","K","L","M","N","P","Q","R","S","T","V","W","Y"])
    contacts_files_mames = []
    for file in os.listdir(contacts_folder_path):
        filename = os.fsdecode(file)
        if filename.endswith(".txt") and filename != "contact_metrix.txt" and filename != "contact_metrix_info.txt":
            contacts_files_mames.append(filename)
    for filename in contacts_files_mames:
        with open(contacts_folder_path + "/" + filename  ,"r") as f:
            counter = 0
            for line in f.readlines():
                if counter == 0 :
                    counter += 1
                    continue
                line.rstrip("\n")
                line_sep = line.split("\t")
                regardent_contact_key = ""
                if len(filename.split("_")[0]) == 4:
                    regardent_contact_key = filename.split("_")[0]+filename.split("_")[1] #6W41_C_H.txt  >>> 6W41C
                elif len(filename.split("_")[0]) == 6:
                    regardent_contact_key = filename.split("_")[0]# 7NX8BA_E_A.txt >>>  7NX8BA
                else:
                    raise ValueError("why do you have a contact first split by _ in len differnt than 4,6")
                if not(regardent_contact_key in contact_metrix_info.loc[get_3id_to_1id_or_1id_from_3id(line_sep[1]),int(line_sep[0]) -1]):
                    contact_metrix_info.loc[get_3id_to_1id_or_1id_from_3id(line_sep[1]),int(line_sep[0]) -1] += (regardent_contact_key +",")
                    contact_metrix.loc[get_3id_to_1id_or_1id_from_3id(line_sep[1]),int(line_sep[0]) -1] += 1
    return (contact_metrix.T,contact_metrix_info.T)


def analyzer(PDB_folder_path,contacts_folder_path,analyzide_data_folder_path):
    (contact_metrixs,contact_metrix_info) = contact_metrix(contacts_folder_path)
    (PFM,PFM_info) = make_pfm_from_pdbs(PDB_folder_path)
    contacts_over_PFM = contact_metrixs/PFM

    contact_metrixs.to_csv(analyzide_data_folder_path + "/contact_metrix.txt","\t")
    contact_metrix_info.to_csv(analyzide_data_folder_path + "/contact_metrix_info.txt","\t")
    PFM.to_csv(analyzide_data_folder_path + "/PFM.txt","\t")
    PFM_info.to_csv(analyzide_data_folder_path + "/PFM_info.txt","\t")
    contacts_over_PFM.to_csv(analyzide_data_folder_path + "/contacts_over_PFM.txt","\t")


def the_big_picture(PDB_folder_path,contacts_folder_path,analyzide_data_folder_path):
    make_contacts_folder(PDB_folder_path,contacts_folder_path)
    analyzer(PDB_folder_path,contacts_folder_path,analyzide_data_folder_path)

"""
#################################################################################
clstringggggggggggggggggggg
##################################################################################
"""
"""
(h,h)/(l,h)
"""
def make_pair_of_antibodies(chain_id_list,pdb_name,sabdab_data_base_path):
    sabdab_data_base = pd.read_csv(sabdab_data_base_path,"\t")
    pairs_index = []
    pairs_chine_id = []
    sabdab_data_base_true_false = sabdab_data_base[sabdab_data_base.pdb == pdb_name].copy().isnull()

    Lchains_true_false = list(sabdab_data_base_true_false.Lchain)
    Lchains = list(sabdab_data_base[sabdab_data_base.pdb == pdb_name].Lchain)
    Hcahins = list(sabdab_data_base[sabdab_data_base.pdb == pdb_name].Hchain)

    for index in range(len(Hcahins)):
        if (Lchains_true_false[index]):
            pairs_chine_id.append((Hcahins[index],Hcahins[index]))
        else:
            pairs_chine_id.append((Lchains[index],Hcahins[index]))

    """
    for row in sabdab_data_base[sabdab_data_base.pdb == pdb_name].iterrows():
        if ( list(row.Lchain.isnull())[0] ):
            pairs_chine_id.append((list(sabdab_data_base[sabdab_data_base.pdb == pdb_name].Hchain)[0],(list(sabdab_data_base[sabdab_data_base.pdb == pdb_name].Hchain)[0])))          
        else:
            pairs_chine_id.append((list(sabdab_data_base[sabdab_data_base.pdb == pdb_name].Lchain)[0],list(sabdab_data_base[sabdab_data_base.pdb == pdb_name].Hchain)[0]))
    """
    for (chainA,chainB) in pairs_chine_id:
        try:
            pairs_index.append((chain_id_list.index(chainA),chain_id_list.index(chainB)))
        except ValueError:
            print(pdb_name,pairs_chine_id,chainA,chainB,chain_id_list)

    return pairs_index

def fasta_of_antibodies_pairs(pdbs_folder_path,analyzide_data_folder_path,sabdab_data_base_path):
    """
    :param pdbs_folder_path:
    :param analyzide_data_folder_path:
    :param sabdab_data_base_path:
    :return: will write fasta that will combine LH antibodie together
    """
    pdbs_files_mames = []
    for file in os.listdir(pdbs_folder_path):
        filename = os.fsdecode(file)
        if filename.endswith(".pdb"):
            pdbs_files_mames.append(filename.rstrip(".pdb"))
    pairess = {}
    sequences = []
    for filename in pdbs_files_mames:
        atom_features = get_atom_features(pdbs_folder_path + "/" + filename)

        parser = PDBParser()
        structure = parser.get_structure(filename, pdbs_folder_path + "/"+ filename+".pdb")

        antibody_list = []
        chain_id_list = []
        for i in range(len(list(structure.get_chains()))):
            chain_id = list(structure.get_chains())[i].get_id()
            d_o_d = get_coordinates_of_chain(atom_features,chain_id)
            try:
                is_it = what_is_the_chain(d_o_d,chain_id,filename.lstrip("PyIgClassifyOutput_").rstrip(".pdb"))
            except UnboundLocalError:
                return filename,chain_id

            if(is_it=="heavy" or is_it=="light"):
                antibody_list.append(d_o_d)
                chain_id_list.append(chain_id)

        couples = make_pair_of_antibodies(chain_id_list,filename.lstrip("PyIgClassifyOutput_").lower(),sabdab_data_base_path)
        pairs = couples ,chain_id_list
        pairess[filename] = pairs

        for couple in couples:

            record = SeqRecord(
            Seq(seq(antibody_list[couple[0]]) + seq(antibody_list[couple[1]])),
            id = filename.lstrip("PyIgClassifyOutput_") + chain_id_list[couple[0]] + chain_id_list[couple[1]],
            name = filename.lstrip("PyIgClassifyOutput_") + chain_id_list[couple[0]] + chain_id_list[couple[1]],
            description = filename.lstrip("PyIgClassifyOutput_") + chain_id_list[couple[0]] + chain_id_list[couple[1]])

            sequences.append(record)
    with open(analyzide_data_folder_path + "/fasta_of_antibodies_pairs.fasta", "w") as output_handle:
            SeqIO.write(sequences, output_handle, "fasta")


"""
make clustring by cd-hit
"""
def read_cluster_of_pairs(clsr_pairs_file_path):
    dictionary = {}
    with open(clsr_pairs_file_path, "r") as f:
        is_first_time = True
        outer_key = ""
        inner_key = ""
        pdb_file = ""
        chain = ""
        description = ""
        for line in f.readlines():
            line = line.rstrip("\n")
            if line.startswith(">Cluster"):
                if(not is_first_time):
                    dictionary[outer_key] = pd.DataFrame(dictionary[outer_key]).T
                is_first_time = False
                outer_key = line[9:]
            elif(not line == ""):
                inner_key = line[0]
                split_line = line.split(">")
                pdb_file = split_line[1][0:4]
                chain = split_line[1][4:6]
                split_split_line = split_line[1].split(" ")
                description = split_split_line[-1]
                inner_inner_key = {"pdb_file":pdb_file, "chain":chain, "description": description }
                if(outer_key in dictionary):
                    dictionary[outer_key][inner_key] = inner_inner_key
                else:
                    dictionary[outer_key] = {inner_key: inner_inner_key}
        dictionary[outer_key] = pd.DataFrame(dictionary[outer_key]).T
        return dictionary


def choose_represent(list_of_pdb_name,PDB_folder_path):
    """
    :param list_of_pdb_name:
    :param PDB_folder_path:
    :return: assist to list_of_represnters
    """
    list_of_pdb_file_path = []

    for i in list_of_pdb_name:
        list_of_pdb_file_path.append(PDB_folder_path+"/PyIgClassifyOutput_"+i+".pdb")
    full = []
    not_full = []
    for i in list_of_pdb_file_path:
        #print("the is " + str(i) + "and"+str(list_of_pdb_file_path))
        if is_it_full(i)[0]:
            full.append(i)
        else:
            not_full.append(i)

    min_res = 10000000000000000000000000000000000000000000000000000000000000000000000000000000
    min_path = ""
    for i in full:
        r =  resolotion(i)
        if min_res > r:
            min_res = r
            min_path = i

    if min_path == "":
        min_res = 10000000000000000000000000000000000000000000000000000000000000000000000000000000
        min_path = ""
        for i in not_full:
            r =  resolotion(i)
            if min_res > r:
                min_res = r
                min_path = i

    return min_path,min_path[-8:-4]

def list_of_represnters(clstr_file_path,PDB_folder_path):
    """
    :param clstr_file_path: from CD-HIT
    :param PDB_folder_path:
    :return: list of selected pdb with the best resulotion from the full group if exist
    """
    list_of_represnters_to_return = []
    clst = read_cluster_of_pairs(clstr_file_path)
    for i in clst.values():
        pdb_names = list(i.pdb_file)
        path_rep,name_rep = choose_represent(pdb_names,PDB_folder_path)
        list_of_represnters_to_return.append(name_rep+list(i[i.pdb_file == name_rep].chain)[0])
    return list_of_represnters_to_return


"""
input:
    tow paths
    list of selected exp [6XPDLH,6XPRHH] meaning first 4 letters is pdbID next 2 letters represent tow Light + Heavy chain or Heavy chain alone
output:
    fill the unregardent folder in the relevant selected from the list
"""
def regardner(regardent_PDB_folder_path,unregardent_PDB_folder_path,list_of_selected,sabdab_data_base_path):
    sabdab = pd.read_csv(sabdab_data_base_path,"\t")
    for selected in list_of_selected:
        antigen_chain2 = "flag"
        first_chain = selected[4]
        second_chain = selected[5]
        #NO ALONE ANTIBODI
        if first_chain == second_chain:
            continue
        pdb_id = selected[0:4]
        antigene_list = list(sabdab[(sabdab.pdb == pdb_id.lower()) & (sabdab.Hchain == second_chain)].antigen_chain)
        if len(antigene_list)>1 :
            raise ValueError("thier is more than one pdbs with the same name and the same hevy chain"+ pdb_id.lower())
        if len(antigene_list) == 0:
            raise ValueError("thier is no antigen in sabdab for that " + pdb_id.lower() + " the heavy cahin " + second_chain)

        if len(antigene_list[0])>1:
            split = antigene_list[0].replace(" ","").split("|")
            if (len(split) == 2) and (len(split[0]) == 1) and (len(split[1]) == 1):
                antigene_chain = split[0]
                antigen_chain2 =  split[1]
            else:
                raise ValueError("weird antigene please chack the data "+ pdb_id.lower() + " the heavy cahin " + second_chain )
        else:
            antigene_chain = antigene_list[0]
        file_name = "PyIgClassifyOutput_" + pdb_id.upper() + ".pdb"
        IN = open(regardent_PDB_folder_path + "/" + file_name,'r')
        OUT = open(unregardent_PDB_folder_path + "/PyIgClassifyOutput_" + selected + ".pdb", "w")
        #copy and peste
        for line in IN.readlines():
            if line.startswith('ATOM'):
                chain_id = line[21:22].replace(" ", "")
                if (chain_id == first_chain) or (chain_id == second_chain) or (chain_id == antigene_chain) or (chain_id == antigen_chain2):
                    OUT.write(line)
        IN.close()
        OUT.close()

def cif_to_pdb(input_cif_file):
    parser = MMCIFParser()
    structure = parser.get_structure("c", input_cif_file)
    io = PDB.PDBIO()
    io.set_structure(structure[0])
    io.save( input_cif_file.replace(".cif",".myout.pdb"))


########################################################################################################################
#########################################################pymol##########################################################
########################################################################################################################


def pymol_align_native(native_path,receptor_path):
    """
    from nirit need some extra work
    :param native_path:
    :param receptor_path:
    :return:
    """
  #  import __main__
  #  __main__.pymol_argv=['pymol','-gc']
    pymol.finish_launching()
    #load Native and recptor
    path="/".join(native_path.split("/")[0:-1])
    os.chdir(path)
    n=native_path.split("/")[-1].split(".")[0]
    r=receptor_path.split("/")[-1].split(".")[0]
    pymol.cmd.load(native_path,n)
    pymol.cmd.load(receptor_path,r)
    #align 1G9I_EI, 1BTP_A
    pymol.cmd.align(r,n,transform=0,object='aln')
    #select native, 1BTP_A or (1G9I_EI and chain I)
    pymol.cmd.select("native","%s or ( %s and chain B)"%(r,n))
    #save 1BTP_A_native.pdb, native
    pymol.cmd.save("n.pdb","native")
   # pymol.cmd.save("native.pdb")



