import pandas as pd
from Bio import SeqIO
import Bio.Data.IUPACData
import re
import collections
import numpy as np
import itertools
import argparse

# file = "/Users/sarapour/Downloads/orf_trans.fasta"
# NOTE: all stop codons (*) were removed from the string, so A*B becomes AB
aa_list = list(Bio.Data.IUPACData.protein_letters)


def read_sequences(file):
    """
    reads in fasta text file into a fasta generator object with biopython seqio
    :param file: fasta format text file
    :return: generator object of fasta
    """
    fasta_sequences = SeqIO.parse(open(file), 'fasta')
    return fasta_sequences


# q1: empirically count the fraction of proteins which match the pattern R-x(2)-[ST]-x-[ST], the same kinase binding
# motif investigated in Part I

def sliding(a, n):
    """
    Divides string a into smaller subsequences of size n.
    :param a: string
    :param n: integer
    :return: returns object representing all overlapping string windows in a of size n.
    """
    return (a[i:i + n] for i in range(len(a) - n + 1))


def substring_count(a, b):
    """
    Counts how many matches to b there are in string a.
    :param a: string
    :param b: string of shorter length you want to find and count in a
    :return: count of how many times b appears in a
    """
    return sum(s == b for s in sliding(a, len(b)))


def count_overlapping(file, seq_to_find):
    """
    HELPER: counts overlapping matches of seq_to_find to file
    :param file: fasta text
    :param seq_to_find: string of aa
    :return: integer of counts
    """
    fasta_seqs = read_sequences(file)
    total = 0
    for fasta in fasta_seqs:
        name, sequence = fasta.id, str(fasta.seq)
        sequence = sequence.replace("*", "")
        counts = substring_count(sequence, seq_to_find)
        total += counts
    return total


# def count_indiv_matches(pattern, file):
#     """
#     counts all non overlapping matches in a sequence
#     """
#     fasta_seqs = read_sequences(file)
#     regex = re.compile(pattern)
#     total = 0
#     for fasta in fasta_seqs:
#         name, sequence = fasta.id, str(fasta.seq)
#         counts = len(re.findall(regex, sequence))
#         total += counts
#     return total


def count_all_matches(file, pattern='R..[ST].[ST]'):
    """
    counts regex-like pattern of amino acid sequences in a fasta file. it is counting number of proteins that have the
     match, not the total count. q1.
    :param pattern: regex exp you want to count
    :param file: text file of fasta sequences
    :return: num proteins that have regex/num total inputted, frequency
    """
    pattern = re.compile(pattern)
    fasta_seqs = read_sequences(file)
    has_match = 0
    total = 0
    for fasta in fasta_seqs:
        name, sequence = fasta.id, str(fasta.seq)
        sequence = sequence.replace("*", "")
        total += 1
        if bool(pattern.search(sequence)):
            has_match += 1
    fraction_match = has_match / total
    return fraction_match


def measure_marginal_freq_aa(file):
    """
    calculates the empirical marginal frequencies of amino acids in the fasta file
    use like measure_marginal_freq_aa(file), q2
    :param file: text fasta file
    :return: vector of the amino acid freqs in alph order
    """
    fasta_seqs = read_sequences(file)
    aa_dict = {key: 0 for key in aa_list}
    aa_dict = collections.Counter(aa_dict)
    total_aa = 0
    for fasta in fasta_seqs:
        name, sequence = fasta.id, str(fasta.seq)
        sequence = sequence.replace("*", "")
        total_aa += len(sequence)
        aa_counts = {e: sequence.count(e) for e in set(sequence)}
        aa_counts = collections.Counter(aa_counts)
        aa_dict = aa_dict + aa_counts
    aa_dict = {k: v / total for total in (sum(aa_dict.values()),) for k, v in aa_dict.items()}
    aa_array = pd.DataFrame([collections.OrderedDict(sorted(aa_dict.items()))]).T
    return aa_array


# m = measure_marginal_freq_aa(file1)


# Q4
def expected_freq_diaa(file):
    """
    Using your marginal amino acid frequencies from Step 2, calculate the expected frequencies of di-amino acid ‘words’,
    e.g., Alanine followed by Glycine forms the di-amino acid word AG. Here you should assume independence between
     amino acids. q4
    :param file: text fasta file
    :return: 20-by-20 frequency matrix.
    """
    single_freq = measure_marginal_freq_aa(file)
    single_freq_array = np.array(list(single_freq.values))
    difreq_freq_array = np.outer(single_freq_array, single_freq_array.transpose())
    difreq_freq_array = pd.DataFrame(difreq_freq_array, aa_list, aa_list)
    return difreq_freq_array


def chunk_string(s, n):
    """
    divides string s into list of smaller substrings of size n
    :param s:string
    :param n: integer
    :return: list of overlapping substrings of size n.
    """
    return [s[i:i + n] for i in range(len(s) - n + 1)]


# q5
def count_diaa(file):
    """
    calculate total number of diamino acids in fasta
    :param file: fasta file
    :return: integer
    """
    fasta_seqs = read_sequences(file)
    total_di = 0
    for fasta in fasta_seqs:
        name, sequence = fasta.id, str(fasta.seq)
        sequence = sequence.replace("*", "")
        total_di += len(sequence) - 2 + 1
    return total_di


def empiric_freq_diaa(file):
    """
    Empirically measure di-amino acid frequencies, q5
    :param file: the FASTA file of aas
    :return:dictionary of diamino acid counts in the whole file.
    """
    fasta_seqs = read_sequences(file)
    aa_list = Bio.Data.IUPACData.protein_letters
    aa_list = sorted(list(itertools.product(aa_list, repeat=2)))
    aa_list = ["".join(a) for a in aa_list]
    aa_dict = {key: 0 for key in aa_list}
    aa_dict = collections.Counter(aa_dict)
    total_aa = 0
    for fasta in fasta_seqs:
        name, sequence = fasta.id, str(fasta.seq)
        sequence = sequence.replace("*", "")
        total_aa += len(sequence)
        aa_counts = collections.Counter(chunk_string(sequence, 2))
        aa_dict = aa_dict + aa_counts
    return aa_dict


def order_diaa_to_array(file, normalize=True):
    """
    use: matrix = order_dinuc_to_array(file)
    :param file: fasta text file
    :param normalize: normalize makes it marginal frequencies by dividing by total number of diaa, normalize is just
    counts from empiric_freq_diaa above, but in df format
    :return: 20-20 frequency matrix in alphabetical order, rows are first aa, columns are second aa. outputs frequency matrix
    divided by the total number of diaa in the file. q5 contd,
    """
    aa_list = list(Bio.Data.IUPACData.protein_letters)
    dict_diaa = empiric_freq_diaa(file)
    od = collections.OrderedDict(sorted(dict_diaa.items()))
    total_di = count_diaa(file)
    arr = np.array(list(od.values()))
    if normalize:
        arr = arr / total_di
        arr2 = np.reshape(arr, (-1, 20))
    else:
        arr2 = np.reshape(arr, (-1, 20))
    arr2 = pd.DataFrame(arr2)
    arr2.columns = aa_list
    arr2.rename(dict(enumerate(arr2.columns)), inplace=True)
    return arr2


# q6
def conditional_freq(file, aa=None):
    """
     Using your empirically measured di-amino acid frequencies, calculate conditional amino acid frequencies where,
    for each vector, the conditioning is on a different one of the 20 possible preceding amino acids.q6
    :param file: fasta text
    :param aa: given aa that you want to see the conditional freq of the 20 aa after. if aa is none then it returns df
    of all conditional probabilities.
    :return: vector size 20
    """
    if aa is None:
        empiric_freq = order_diaa_to_array(file, normalize=False)
        empiric_freq = empiric_freq.div(empiric_freq.sum(axis=1), axis=0)
        return empiric_freq
    else:
        aa = aa.upper()
        empiric_freq = order_diaa_to_array(file, normalize=False)
        aa_select = empiric_freq.loc[aa]
        aa_sum = sum(aa_select)
        aa_norm = aa_select / aa_sum
        return aa_norm


def conditional_I_Q(file):
    """
    specific instance of conditional_freq but tells you conditional frequenceis for all aa if they start with I or Q
    :param file: fasta text
    :return: array of 2x20, first correspond to conditional prob of each aa (alphabetical) if start with I, second is
    same but if it start with Q. q6
    """
    I = np.array(conditional_freq(file, 'I'))
    Q = np.array(conditional_freq(file, 'Q'))
    merge = np.vstack((I, Q))
    merge = pd.DataFrame(merge, index=['I', 'Q'], columns=aa_list)
    return merge


def main():
    parser = argparse.ArgumentParser(description='Process and analyze amino acid compositions of a proteome fasta file')
    parser.add_argument("-f", "--file", help="path to fasta file", required=True)

    parser.add_argument("-frac", "--fraction", nargs='?', type=str, const='R..[ST].[ST]', help="get count of proteins "
                                                                                               "with regex match. if no "
                                                                                               "argument passed it will "
                                                                                               "find R X2 [ST] X [ST]. ")
    parser.add_argument("-marg", "--marginal", action='store_true', help="empirically measuring the marginal frequency "
                                                                         "of each of the 20 amino acids in the yeast "
                                                                         "proteome ")
    parser.add_argument("-exp_di", "--expected_difrequencies", action='store_true', help="the expected frequencies of "
                                                                                         "di-amino acid ‘words’, e.g., "
                                                                                         "Alanine followed by Glycine "
                                                                                         "forms the di-amino acid word "
                                                                                         "AG.Here you should assume "
                                                                                         "independence between amino "
                                                                                         "acids.")
    parser.add_argument("-emp_di", "--empirical_difrequencies", action='store_true',
                        help="Empirically measure di-amino "
                             "acid frequencies.")
    parser.add_argument("-cond_iq", "--conditiona_iso_gln", action='store_true', help="the conditioning is on a "
                                                                                      "I or Q.")

    args = parser.parse_args()
    if args.fraction:
        print(count_all_matches(args.file, args.fraction))
    if args.marginal:
        print(measure_marginal_freq_aa(args.file))
    if args.expected_difrequencies:
        print(expected_freq_diaa(args.file))
    if args.empirical_difrequencies:
        print(order_diaa_to_array(args.file))
    if args.conditiona_iso_gln:
        print(conditional_I_Q(args.file))


if __name__ == "__main__":
    main()
