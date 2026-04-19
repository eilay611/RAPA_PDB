# TODO replace 7DK7fe
# TODO check CDR3
import sys
sys.path.insert(0,"/home/eilay/PycharmProjects/lab_projects")
import anarci
from Bio import SeqIO
import pickle
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import numpy as np
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq


##### -------------- Functions ------------------ ##### <editor-fold>

def dump_to_pickle(path, my_info):
    with open(path, 'wb') as f:
        pickle.dump(my_info, f)


def load_from_pickle(path):
    def open_file():
        return open(path, 'rb')

    with open_file() as f:
        return pickle.load(f)


def anarci_results_report_orig(sequences, numbering, alignment_details, hit_tables):
    # Iterate over the sequences
    for i in range(len(sequences)):
        if numbering[i] is None:
            print('ANARCI did not number', sequences[i][0])
        else:
            print('ANARCI numbered', sequences[i][0])
            print('It identified %d domain(s)' % len(numbering[i]))

            # Iterate over the domains
            for j in range(len(numbering[i])):
                domain_numbering, start_index, end_index = numbering[i][j]
                print('This is the IMGT numbering for the %d\th domain:' % j, domain_numbering)
                print('This is the bit of the sequence it corresponds to:',
                      sequences[i][1][start_index:end_index + 1])
                print('These are the details of the alignment:')
                for (key, value) in alignment_details[i][j].items():
                    print(key, ':', value)
                print('This is the summary of the hits that HMMER found')
                for line in hit_tables[i]:
                    print(line)
        print('\n', '_' * 40)

    print('Do with this infomation as you wish')

    print('\n', '*' * 40)


def anarci_results_report_mine(numbering, alignment_details, hit_tables, sequences):
    # Each should have the same number of elements as the number of sequences submitted
    assert len(numbering) == len(alignment_details) == len(hit_tables) == len(sequences)

    #### Print anarci results report
    # anarci_results_report(sequences, numbering, alignment_details, hit_tables)

    bad_chains_list = []

    species_count = {}
    chains_count = {}
    chains_type_count = {}
    chains_type_length_count = {}

    for i in range(len(sequences)):
        seq_name = sequences[i][0]

        if numbering[i] is None:
            print('@@@---------- Sequence', seq_name)
            print('###!!! ANARCI did not number this sequence')
            bad_chains_list.append(seq_name)
        else:
            chains_num = len(numbering[i])
            chains_count[chains_num] = chains_count.get(chains_num, 0) + 1

            if chains_num != 1:
                print('@@@---------- Sequence', seq_name)
                print('###!!! ANARCI identified %d domain(s)' % len(numbering[i]))
                bad_chains_list.append(seq_name)

            # Iterate over the domains
            for j in range(len(numbering[i])):
                assert alignment_details[i][j]['query_name'] == seq_name

                assert alignment_details[i][j]['scheme'] == 'imgt'

                species = alignment_details[i][j]['species']
                species_count[species] = species_count.get(species, 0) + 1

                chain_type = alignment_details[i][j]['chain_type']

                chains_type_count[chain_type] = chains_type_count.get(chain_type, 0) + 1

                chains_type_length = len(numbering[i][j][0])
                if chain_type not in chains_type_length_count:
                    chains_type_length_count[chain_type] = {}
                chains_type_length_count[chain_type][chains_type_length] = chains_type_length_count[chain_type].get(
                    chains_type_length, 0) + 1

    print('\nchains_count:', chains_count)
    print('\nless_chains_list:', bad_chains_list)
    print('\nchains_type_count:', chains_type_count)
    print('\nchains_type_length_count:', chains_type_length_count)
    print('\nspecies_count:', species_count)


def compare_anarci_df_seqs_to_orig(anarci_df, seqs):
    """
    Gets the anarci results df containing IMGT numbered seqs,
    and the original sequences.
    Concatenates the IMGT numbered seqs from the anarci_df,
    deletes the added '-' symbols, and compares with the original sequence.
    Since in some cases the first letters from the original sequence are
    dropped, looks for the beginning of the IMGT numbered seq.
    :param anarci_df:
    :param seqs:
    :return:
    """
    correctness_list = [] # correct = True, error = False
    errors_list = []
    for i, i_name in enumerate(anarci_df.index):
        print(i_name)
        id = anarci_df.loc[i_name, 'Id']
        seq = ''.join(anarci_df.loc[i_name, '1':])
        seq = seq.replace('-', '')

        if seqs[id][:len(seq)] == seq:
            correctness_list.append(True)
        else:
            correct = False

            print()
            print('@@@', id, ' - oops, sequences are different...')
            print('original seq:           ', seqs[id])
            print('concatenated anarci seq:', seq)
            print("eilay")
            for start in range(len(seqs[id])):
                if seqs[id][start:start + len(seq)] == seq:
                    print('# But - Yay! They are identical if starting at position', start + 1)
                    correct = True

            correctness_list.append(correct)
            if not correct:
                errors_list.append(id)

    # Check if all correctness_list items are true.
    if sum(correctness_list) == anarci_df.shape[0]:
        print('\nFinal result: all anarci seqs are identical to orig\n')
    else: # If not, there's an error in at least 1 sequence! print it.
        raise Exception('\nFinal result: not all anarci seqs are identical to orig!!!:\n', errors_list)


def hamming_distance(str1, str2):
    """
    Gets 2 strings with the same length.
    Returns the number of mismatched letters between them - hamming distance.
    :param str1:
    :param str2:
    :return: int
    """
    if len(str1) != len(str2):
        print('length mismatch!')
        return 'length mismatch!'

    hamming = 0
    for i in range(len(str1)):
        if str1[i] != str2[i]:
            hamming += 1

    return hamming


def get_seqs_distances_df(anarci_df, all_ab_names, first_aa_col='1', last_aa_col=None,
                          write_csv=None, write_length_to_text=None):
    """
    Gets the anarci results df containing IMGT numbered seqs.
    Returns a pairwise distances matrix (pandas df) and the
    sequences length (int).
    :param anarci_df:
    :param all_ab_names:
    :return: dist_df, seq_len
    """
    dist_df = pd.DataFrame(index=all_ab_names, columns=all_ab_names, dtype=np.float64)

    seqs = {}
    for i in range(len(all_ab_names)):
        for j in range(len(all_ab_names)):
            if j < i:
                ab_1 = dist_df.index[i]
                ab_2 = dist_df.columns[j]

                if (ab_1 != '7DK7fe') and (ab_2 != '7DK7fe'): # dont analyze '7DK7fe' # TODO add it manually

                    seq_1 = ''.join(anarci_df.loc[ab_1, first_aa_col:last_aa_col])
                    seq_2 = ''.join(anarci_df.loc[ab_2, first_aa_col:last_aa_col])

                    hamming = np.float64(hamming_distance(seq_1, seq_2))
                    dist_df.loc[ab_1, ab_2] = hamming
                    dist_df.loc[ab_2, ab_1] = hamming

                    seqs[ab_1] = seq_1
                    seqs[ab_2] = seq_2

            elif j == i:
                dist_df.iloc[i, i] = 0


    seq_len = len(anarci_df.loc[all_ab_names[0], first_aa_col:last_aa_col])

    # drop '7DK7fe' # TODO add it manually
    if '7DK7fe' in dist_df.index:
        dist_df.drop(index='7DK7fe', inplace=True)
    if '7DK7fe' in dist_df.columns:
        dist_df.drop(columns='7DK7fe', inplace=True)

    if write_csv is not None:
        dist_df.to_csv(write_csv)
    if write_length_to_text is not None:
        write_list_to_txt([seq_len], write_length_to_text)

    return dist_df, seq_len, seqs


def get_lower_matrix(mat_df):
    """
    Gets n*n matrix (pandas df format). Returns
    the same matrix only with the lower triangle.
    Other values are turned to None

    :param mat_df: pandas df
    :return: pandas df
    """
    mat_lower_df = mat_df.copy()
    for i in range(mat_lower_df.shape[0]):
        for j in range(mat_lower_df.shape[0]):
            if j >= i:
                mat_lower_df.iloc[i, j] = None

    return mat_lower_df


def plot_histogram_from_pairwise_df(dist_lower_df, title='', line_x=None,
                                    figsize=(18, 8),output_to_save=""):
    """
    Gets a pairwise df. Plots a histogram of the values
    with value-counts annotated.

    :param dist_lower_df: pandas df
    :param title: plot title
    :return:
    """
    pair_counts = dist_lower_df.melt()['value'].value_counts()

    plt.figure(figsize=figsize)
    plt.bar(pair_counts.index, pair_counts)

    x_vals = list(pair_counts.index)
    y_vals = list(pair_counts.values)

    for i in range(len(y_vals)):
        plt.annotate(str(y_vals[i]), xy=(x_vals[i], y_vals[i]),
                     ha='center', va='bottom', fontsize=8)

    if line_x is not None:
        plt.axvline(x=line_x, color='black')
        plt.annotate(str(line_x), xy=(line_x+line_x/10, max(y_vals)),
                     ha='left', va='bottom', fontsize=13)

    plt.ylabel('Pairs Count')
    plt.title(title)
    if output_to_save != "":
        pair_counts.to_csv(output_to_save+title+".csv")
        plt.savefig(output_to_save+title+".png")
    plt.show()

def plot_dist_info(dist_df, seq_len, title_add='',output_to_save=""):

    # lower triangle matrix
    dist_lower_df = get_lower_matrix(dist_df)
    # lower triangle matrix - percentage of difference
    dist_lower_perc_df = 100 * dist_lower_df / seq_len

    ### plot clustermap
    plt.figure(figsize=(18, 15))
    g = sns.clustermap(dist_df, mask=(dist_df.isna()))
    g.fig.suptitle('Distances - {}'.format(title_add))
    if output_to_save != "":
        dist_df.to_csv(output_to_save+'Distances - {}'.format(title_add)+".csv")
        plt.savefig(output_to_save+'Distances - {}'.format(title_add)+".png")
    plt.show()

    ### plot histograms
    plot_histogram_from_pairwise_df(dist_lower_df, line_x=(seq_len*0.05),
                                    title='Distances Histogram - {}'.format(title_add),output_to_save = output_to_save)

    # plot_histogram_from_pairwise_df(dist_lower_perc_df, line_x=5,
    #                                 title='Percentage Distances Histogram - {}'.format(title_add))


def print_less_than_thresh_report(dist_df, seq_len,
                                  info_heavy, info_light, thresh=0.05, print_chains=True):
    suspected_abs = set()
    counter = 0

    print('\n#########------------- Similar pairs report----------#########\n'.format(counter))
    list_of_tupel_to_return = []
    for ab1 in dist_df.index:
        for ab2 in dist_df.columns:
            if dist_df.loc[ab1, ab2] <= thresh * seq_len:
                counter += 1
                print('\n###------------------{}----------------###\n'.format(counter))

                print('Abs pair:', ab1, ab2)
                print('Different amino acids:', int(dist_df.loc[ab1, ab2]),
                      'out of', seq_len)
                list_of_tupel_to_return.append((ab1,ab2))

                if print_chains:
                    print('\nHeavy chains:')
                    print('Ab', ab1, info_heavy['seqs'][ab1])
                    print('Ab', ab2, info_heavy['seqs'][ab2])

                    print('\nLight chains:')
                    print('Ab', ab1, info_light['seqs'][ab1])
                    print('Ab', ab2, info_light['seqs'][ab2])

                suspected_abs.update({ab1, ab2})

    print('\nFound a total of {} pairs'.format(counter))
    print('containing {} unique abs:'.format(len(suspected_abs)))
    print(suspected_abs)
    return list_of_tupel_to_return

def write_list_to_txt(listToWrite, path, encoding='utf-8'):
    with open(path, 'w', encoding=encoding) as thefile:
        for item in listToWrite:
            thefile.write("%s\n" % str(item))

def write_fasta_of_imgt_seq(info_heavy,info_light,root_path):
    sequences = []
    for id_hertz, seq_light_imgt in info_light["light_seqs"].items():
        record = SeqRecord(
            Seq(seq_light_imgt.replace("-","X")+info_heavy["heavy_seqs"][id_hertz].replace("-","X")),
            id=id_hertz,
            name=id_hertz,
            description=id_hertz)
        sequences.append(record)
    with open(root_path + "fasta_of_antibodies_pairs_imgt.fasta", "w") as output_handle:
            SeqIO.write(sequences, output_handle, "fasta")
def write_fasta_of_imgt_seq_seprately(info_heavy,info_light,root_path):
    sequences = []
    for id_hertz, seq_light_imgt in info_light["light_seqs"].items():
        record = SeqRecord(
            Seq(seq_light_imgt.replace("-","X")),
            id=id_hertz,
            name=id_hertz,
            description=id_hertz)
        sequences.append(record)
    with open(root_path + "fasta_of_antibodies_light_imgt.fasta", "w") as output_handle:
            SeqIO.write(sequences, output_handle, "fasta")

    sequences = []
    for id_hertz, seq_light_imgt in info_light["light_seqs"].items():
        record = SeqRecord(
            Seq(info_heavy["heavy_seqs"][id_hertz].replace("-","X")),
            id=id_hertz,
            name=id_hertz,
            description=id_hertz)
        sequences.append(record)
    with open(root_path + "fasta_of_antibodies_heavy_imgt.fasta", "w") as output_handle:
            SeqIO.write(sequences, output_handle, "fasta")


# def get_imgt_num_from_numbering(num_aa_tuple):
#     imgt_num_tup = num_aa_tuple[0]
#
#     letter = '.' + imgt_num_tup[1]
#     if letter == '. ':
#         letter = ''
#
#     return str(imgt_num_tup[0]) + letter


# def get_all_imgt_nums_list():
#     letters = ['A', 'B', 'C', 'D', 'E', 'F', 'G', 'H', 'I']
#
#     all_imgt_nums = []
#     for num in range(1, 128 + 1):
#
#         if num not in [111, 112]:
#             all_imgt_nums.append(str(num))
#
#         elif num == 111:
#             all_imgt_nums.append(str(111))
#             for letter in letters:
#                 all_imgt_nums.append(str(111) + '.' + letter)
#
#         if num == 112:
#             for letter in letters[::-1]:
#                 all_imgt_nums.append(str(112) + '.' + letter)
#             all_imgt_nums.append(str(112))

##### -------------------------------- ##### </editor-fold>

##### ------------- Chain analysis function ------------- ##### <editor-fold>

def get_chain_analysis(root_path, chain): # chain = 'heavy' / 'light'
    """
    :param root_path: path +/ at the end
    :param chain:
    :return:
    """

    print('\n@@@@@@ Starting to analyze {} chain\n'.format(chain))

    ### Calculated params
    if chain == 'light':
        fasta_name = 'fasta_of_antibodies_light_chains'
        chains = ['K', 'L']
        type_letters = 'KL'
    elif chain == 'heavy':
        fasta_name = 'fasta_of_antibodies_heavy_chains'
        chains = ['H']
        type_letters = 'H'
    else:
        raise ValueError("chain have to be  'heavy' / 'light'")

    ##### -------------- Get IMGT numbering with anarci ------------------ ##### <editor-fold>

    ###--- Get sequences
    seqs_fasta = SeqIO.to_dict(SeqIO.parse(open(root_path + fasta_name + '.fasta'), 'fasta'))
    seqs = {key: str(seqs_fasta[key].seq) for key in seqs_fasta}

    sequences = list(seqs.items())

    all_ab_names = list(seqs_fasta.keys()) # all names. Including these that anarci fails on.

    ###--- Run anarci

    results = anarci.anarci(sequences, scheme="imgt", allowed_species=['human'],
                     allow=set(chains), ncpu=1, output=False, assign_germline=True)

    numbering, alignment_details, hit_tables = results

    # dump_to_pickle(root_path + fasta_name + '__anarci_output.pkl', results)

    anarci_results_report_mine(numbering, alignment_details, hit_tables, sequences)

    ###--- get table (create csv and read it)

    anarci.csv_output(sequences, numbering, alignment_details, root_path + '120_anarci_output__')

    anarci_df = pd.read_csv(root_path + '120_anarci_output___{}.csv'.format(type_letters))
    anarci_df.index = anarci_df['Id']

    ###--- Make sure original sequence and anarci sequence (without '-') are identical
    compare_anarci_df_seqs_to_orig(anarci_df, seqs)

    ##### -------------------------------- ##### </editor-fold>

    ##### -------------- Compare sequences ------------------ ##### <editor-fold>

    # num_comparisons = (len(ab_names) * len(ab_names) - len(ab_names)) / 2

    dist_df, seq_len, seqs_imgt = get_seqs_distances_df(anarci_df, all_ab_names,
                                   first_aa_col='1', last_aa_col=None,
                                   write_csv=root_path + 'comparisons_df__{}.csv'.format(chain),
                                   write_length_to_text=root_path + 'length__{}.csv'.format(chain))

    dist_cdr3_df, seq_cdr3_len, cdr3_seqs = get_seqs_distances_df(anarci_df, all_ab_names,
                                             first_aa_col='105', last_aa_col='117',
                                             write_csv=root_path + 'comparisons_cdr3_df__{}.csv'.format(chain),
                                             write_length_to_text=root_path + 'length__cdr3_{}.csv'.format(chain))

    plot_dist_info(dist_df, seq_len, title_add=chain + ' chain',output_to_save=root_path)

    plot_dist_info(dist_cdr3_df, seq_cdr3_len, title_add='CDR3 ' + chain + ' chain',output_to_save=root_path)

    ##### -------------------------------- ##### </editor-fold>

    return {'seqs': seqs,
            'all_ab_names': all_ab_names,
            'anarci_df': anarci_df,
            'dist_df': dist_df,
            'seq_len': seq_len,
            'dist_cdr3_df': dist_cdr3_df,
            'seq_cdr3_len': seq_cdr3_len,
            'cdr3_seqs': cdr3_seqs,
            chain+'_seqs':seqs_imgt
            }
def list_of_tupels_of_similarty_abs(root_path):
    info_heavy = get_chain_analysis(root_path, 'heavy')
    info_light = get_chain_analysis(root_path, 'light')
    write_fasta_of_imgt_seq(info_heavy, info_light, root_path)
    write_fasta_of_imgt_seq_seprately(info_heavy, info_light, root_path)
    dist_df_h = info_heavy['dist_df']
    dist_df_l = info_light['dist_df']
    assert (dist_df_h.index == dist_df_l.index).sum() == dist_df_h.shape[0]
    assert (dist_df_h.columns == dist_df_l.columns).sum() == dist_df_h.shape[1]
    dist_df__hl = dist_df_h + dist_df_l
    seq_len_hl = info_heavy['seq_len'] + info_light['seq_len']
    plot_dist_info(dist_df__hl, seq_len_hl, title_add='both chains',output_to_save = root_path)
    list_of_tupels = print_less_than_thresh_report(get_lower_matrix(dist_df__hl), seq_len_hl, info_heavy, info_light)
    return info_light,info_heavy,list_of_tupels

def check_dis():
    print("woel")
##### -------------------------------- ##### </editor-fold>

"""
##### -------------- Run Heavy + Light ------------------ ##### <editor-fold>

#root_path = '/home/liel/PycharmProjects/covid_ab_seqs/'
root_path = '/home/eilay/Downloads/'


# Get heavy + light chain analysis
info_heavy = get_chain_analysis(root_path, 'heavy')
info_light = get_chain_analysis(root_path, 'light')

# Get heavy + light chain analysis
dist_df_h = info_heavy['dist_df']
dist_df_l = info_light['dist_df']

assert (dist_df_h.index == dist_df_l.index).sum() == dist_df_h.shape[0]
assert (dist_df_h.columns == dist_df_l.columns).sum() == dist_df_h.shape[1]

#### Plot and get report - entire chain

# Get united distances matrix and seq_len
dist_df__hl = dist_df_h + dist_df_l
seq_len_hl = info_heavy['seq_len'] + info_light['seq_len']

plot_dist_info(dist_df__hl, seq_len_hl, title_add='both chains')

print_less_than_thresh_report(get_lower_matrix(dist_df__hl), seq_len_hl, info_heavy, info_light)

dist_df__hl.to_csv(root_path + 'comparisons_df__both_chains.csv')
write_list_to_txt([seq_len_hl], root_path + 'length__both_chains.csv')

# Investigate similar pairs
pairs_hl = get_lower_matrix(dist_df__hl).stack().reset_index().rename(columns={'level_0': 'ab1', 'level_1': 'ab2', 0: 'distance'})

pairs_similar_hl = pairs_hl.loc[(pairs_hl['distance'] <= 20) &
                                (pairs_hl['distance'] >= seq_len_hl * 0.05)].sort_values('distance')

# Plot and get report - CDR3
dist_cdr3_df__hl = info_heavy['dist_cdr3_df'] + info_light['dist_cdr3_df']
seq_cdr3_len__hl = info_heavy['seq_cdr3_len'] + info_light['seq_cdr3_len']

plot_dist_info(dist_cdr3_df__hl, seq_cdr3_len__hl,
               title_add='CDR3 both chains')

print_less_than_thresh_report(get_lower_matrix(dist_cdr3_df__hl), seq_cdr3_len__hl,
                              info_heavy, info_light, print_chains=False)

dist_cdr3_df__hl.to_csv(root_path + 'comparisons_df__cdr3__both_chains.csv')
write_list_to_txt([seq_cdr3_len__hl], root_path + 'length__cdr3__both_chains.csv')

# Investigate similar pairs CD3
pairs_cdr3 = get_lower_matrix(dist_cdr3_df__hl).stack().reset_index().rename(columns={'level_0': 'ab1', 'level_1': 'ab2', 0: 'distance'})

pairs_similar_cdr3 = pairs_cdr3.loc[(pairs_cdr3['distance'] <= 5) &
                                    (pairs_cdr3['distance'] >= seq_cdr3_len__hl * 0.05)].sort_values('distance')


##### Investigate similar pairs

## similar hl - check CDR3 similarity
for i in pairs_similar_hl.index:
    print('\n## ------------------------------------------------------------- ##')
    print(pairs_similar_hl.loc[i])
    ab1 = pairs_similar_hl.loc[i, 'ab1']
    ab2 = pairs_similar_hl.loc[i, 'ab2']

    print('\n#------ CDR3 dist: \n')
    print(pairs_cdr3.loc[(pairs_cdr3['ab1'] == ab1) & (pairs_cdr3['ab2'] == ab2)])
    dist_cdr3 = pairs_cdr3.loc[(pairs_cdr3['ab1'] == ab1) & (pairs_cdr3['ab2'] == ab2), 'distance']
    print(float(dist_cdr3) / seq_cdr3_len__hl)

    print('\n#------ CDR3s: heavy\n')
    print('Ab', ab1, info_heavy['cdr3_seqs'][ab1])
    print('Ab', ab2, info_heavy['cdr3_seqs'][ab2])

    print('\n#------ CDR3s: light\n')
    print('Ab', ab1, info_light['cdr3_seqs'][ab1])
    print('Ab', ab2, info_light['cdr3_seqs'][ab2])

    print('\n#------ heavy\n')
    print('Ab', ab1, info_heavy['anarci_df'].loc[ab1, 'v_gene'])
    print('Ab', ab1, info_heavy['anarci_df'].loc[ab1, 'j_gene'], '\n')

    print('Ab', ab2, info_heavy['anarci_df'].loc[ab2, 'v_gene'])
    print('Ab', ab2, info_heavy['anarci_df'].loc[ab2, 'j_gene'], '\n')
    print('Ab', ab2, info_heavy['seqs'][ab2], '\n')

    print('\n#------ light\n')
    print('Ab', ab1, info_light['anarci_df'].loc[ab1, 'v_gene'])
    print('Ab', ab1, info_light['anarci_df'].loc[ab1, 'j_gene'], '\n')

    print('Ab', ab2, info_light['anarci_df'].loc[ab2, 'v_gene'])
    print('Ab', ab2, info_light['anarci_df'].loc[ab2, 'j_gene'], '\n')

## similar CDR3 - check HL similarity
for i in pairs_similar_cdr3.index:
    print('\n## ------------------------------------------------------------- ##')
    print(pairs_similar_cdr3.loc[i])
    ab1 = pairs_similar_cdr3.loc[i, 'ab1']
    ab2 = pairs_similar_cdr3.loc[i, 'ab2']

    print('\n#------ HL dist: \n')
    print(pairs_hl.loc[(pairs_hl['ab1'] == ab1) & (pairs_hl['ab2'] == ab2)])
    dist_hl = pairs_hl.loc[(pairs_hl['ab1'] == ab1) & (pairs_hl['ab2'] == ab2), 'distance']
    print(float(dist_hl) / seq_len_hl)

    print('\n#------ CDR3s: heavy\n')
    print('Ab', ab1, info_heavy['cdr3_seqs'][ab1])
    print('Ab', ab2, info_heavy['cdr3_seqs'][ab2])

    print('\n#------ CDR3s: light\n')
    print('Ab', ab1, info_light['cdr3_seqs'][ab1])
    print('Ab', ab2, info_light['cdr3_seqs'][ab2])

    print('\n#------ heavy\n')
    print('Ab', ab1, info_heavy['anarci_df'].loc[ab1, 'v_gene'])
    print('Ab', ab1, info_heavy['anarci_df'].loc[ab1, 'j_gene'])
    print('Ab', ab1, info_heavy['seqs'][ab1], '\n')

    print('Ab', ab2, info_heavy['anarci_df'].loc[ab2, 'v_gene'])
    print('Ab', ab2, info_heavy['anarci_df'].loc[ab2, 'j_gene'])
    print('Ab', ab2, info_heavy['seqs'][ab2], '\n')

    print('\n#------ light\n')
    print('Ab', ab1, info_light['anarci_df'].loc[ab1, 'v_gene'])
    print('Ab', ab1, info_light['anarci_df'].loc[ab1, 'j_gene'])
    print('Ab', ab1, info_light['seqs'][ab1], '\n')

    print('Ab', ab2, info_light['anarci_df'].loc[ab2, 'v_gene'])
    print('Ab', ab2, info_light['anarci_df'].loc[ab2, 'j_gene'])
    print('Ab', ab2, info_light['seqs'][ab2], '\n')


# info_light['seqs']['7DK7fe']
# info_heavy['seqs']['7DK7fe']


# numbers_list = numbering[max_len_seq_index][0][0]
#
# for tup in numbers_list:
#     print(get_imgt_num_from_numbering(tup))
#
#
#
# all_imgt_nums_list = get_all_imgt_nums_list()
#
#
# seqs_df = pd.DataFrame(index=list(seqs_fasta.keys()),
#                        columns=['species', 'chain_type', 'v_gene', 'j_gene'] + all_imgt_nums_list)


##### -------------------------------- ##### </editor-fold>
"""
