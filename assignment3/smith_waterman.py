import pandas as pd
import numpy as np
import itertools
import builtins
from Bio import SeqIO

fa = '/path/to/fasta'

bases = ['A', 'C', 'G', 'U']
sub_mat = pd.DataFrame(-1, index=bases, columns=bases)
sub_mat.at['G', 'C'] = 3
sub_mat.at['C', 'G'] = 3
sub_mat.at['A', 'U'] = 2
sub_mat.at['U', 'A'] = 2
sub_mat.at['G', 'U'] = 1
sub_mat.at['U', 'G'] = 1

gap_cost = 1

align_hash = {1: ':', 2: '|', 3: '|', -1: ' '}

sub_mat


def smith_waterman(str_a, str_b):
    """Smith–Waterman algorithm. Local alignment + linear gap cost.

    INPUT: 2x strings
    OUTPUT: None
    """
    if not str_a or not str_b:  # base case
        return
    # print(str_a, str_b)
    # print('WS started')
    str_a = str_a[::-1]
    rows = len(str_a) + 1
    cols = len(str_b) + 1
    results = np.zeros((rows, cols))  # result matrix
    arr = np.full([rows, cols], fill_value='-', dtype=str)  # Traceback matrix [[" "] * cols for _ in range(rows)]
    track = arr.tolist()
    align = np.full([rows, cols], fill_value=' ', dtype=str)

    alignments = []

    def back_iter(i, j, current_result_a, current_result_b, current_align):
        """Reconstruct alignment for Smith–Waterman algo."""
        # print(i, j, track[i][j], results[i][j])

        if i == 0 or j == 0 or track[i][j] == "-":
            # return reversed strings
            string_to_print = str(max_i) + "\t3' " + current_result_a[::-1] + " 5'\t" + str(i + 1)
            string_to_print += '\n\t   ' + current_align[::-1]
            string_to_print += '\n' + str(j + 1) + "\t5' " + current_result_b[::-1] + " 3'\t" + str(max_j)

            alignments.append(string_to_print)

        else:
            if "1" in track[i][j]:
                # match - go recursive by diagonale
                back_iter(i - 1, j - 1, current_result_a + str_a[i - 1], current_result_b + str_b[j - 1],
                          current_align + align[i][j])

            if "2" in track[i][j] and current_result_b[-1] != '-':
                # delete
                back_iter(i - 1, j, current_result_a + str_a[i - 1], current_result_b + '-', current_align + ' ')

            if "3" in track[i][j] and current_result_a[-1] != '-':
                # insertion
                back_iter(i, j - 1, current_result_a + '-', current_result_b + str_b[j - 1], current_align + ' ')

    for i, j in itertools.product(range(1, rows), range(1, cols)):
        match_score = sub_mat.loc[str_a[i - 1]][str_b[j - 1]]

        align[i][j] = align_hash[match_score]

        match = results[i - 1][j - 1] + match_score
        delete = results[i - 1][j] - gap_cost
        insert = results[i][j - 1] - gap_cost

        # update results with the score;
        if match == delete and match > insert and match > 0:
            results[i][j] = match
            track[i][j] = "12"

        elif match > delete and match == insert and match > 0:
            results[i][j] = match
            track[i][j] = "13"

        elif match == delete and match == insert and match > 0:
            results[i][j] = match
            track[i][j] = "123"

        elif insert == delete and insert > match and insert > 0:
            results[i][j] = insert
            track[i][j] = "23"

        else:
            results[i][j], track[i][j] = max((match, "1"),  # ↖
                                             (delete, "2"),  # ↑
                                             (insert, "3"),  # ←
                                             (0, "-"))
    optimum_locations = [tuple(coords) for coords in np.argwhere(results == np.max(results))]

    for x in range(len(optimum_locations)):
        max_i, max_j = optimum_locations[x]
        # print(max_i, max_j)
        back_iter(max_i, max_j, '', '', '')

    print('Score ', results[optimum_locations[0][0]][optimum_locations[0][1]])
    for k in range(len(alignments)):
        print('\nAlignment', k + 1, ':')
        print(alignments[k])

def read_do_sw(file):
    """
    RUN
    :param file: fasta file like test_input.fa
    :return: None, prints alignment
    """
    list_seq =[]
    fasta_sequences = SeqIO.parse(open(file), 'fasta')
    try:
        for fasta in fasta_sequences:
            name, sequence = fasta.id, str(fasta.seq).split('#', 1)[0]
            list_seq.append(sequence)
        list_seq_1 = list_seq[::2]
        list_seq_2 = list_seq[1::2]
    except ValueError:
        print('fasta names must all be either seq1 or seq2, paired, and should be of equal length')
    else:
        get_sw(list_seq_1, list_seq_2)


def get_sw(list_seq_1, list_seq_2):
    assert len(list_seq_1) == len(list_seq_2), "seq1 and seq2 have different lengths"
    for i in range(len(list_seq_1)):
        print(i+1)
        str_a = list_seq_1[i]
        str_b = list_seq_2[i]
        smith_waterman(str_a, str_b)

smith_waterman('GGGGGAAAAAUUUUUCCCCCGGGG', 'AAAAGGGGGAAAAUUUUUUUAAA')
smith_waterman('CACGGGGGAAAAAUUUUUCCCCCCCGGGG', 'UGUAAAAGGGGGAAAUUUUUUUAAA')

smith_waterman('CAAGUGACGGUUGAAA', 'UUUCAACCGCACUUG')

smith_waterman('UGAGGUAGUAGGUUGUAUA', 'GCAAUGAUGCCUACCAAACAUUUCCAGACUUAACAUUUUGGUCUCUG')

smith_waterman('GUGAAAUGUU', 'AUUUCCAGGAAUUUAUUCCCCUUCAUAAUUUGUCUCAUUUCAUUUUAUUUCAUCCACUUGGUAGAUGAAGUCACG')
smith_waterman('AAAGAAUUC','CAUGAAUGAAGAUAGGUUGUAAACUGAAUGCUGUGAUAAUACUCUGUAUUCUUUAUGGAAAAUGUUGUCCUGU' )
smith_waterman('AUAUGUUGGAUGAUGGAGU', 'GCAAUGAUGCCUACCAAACAUUUCCAGACUUAACAUUUUGGUCUCUG')
smith_waterman('UUGUAAAGUG', 'AUUUCCAGGAAUUUAUUCCCCUUCAUAAUUUGUCUCAUUUCAUUUUAUUUCAUCCACUUGGUAGAUGAAGUCACG')

# str_a, str_b = ('GGGGGAAAAAUUUUUCCCCCGGGG', 'AAAAGGGGGAAAAUUUUUUUAAA']
