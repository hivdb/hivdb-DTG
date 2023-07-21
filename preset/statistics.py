import numpy as np
from sklearn.metrics import jaccard_score
from scipy.stats import percentileofscore
from collections import defaultdict
from itertools import combinations
from copy import deepcopy
from scipy import stats
from operator import itemgetter


# TODO, accept both contigent table and two list

def calc_jaccard(listA, listB):

    jaccard = jaccard_score(listA, listB)

    permuted_jaccards = []

    num_repeat = 2000
    for _ in range(num_repeat):
        random_list = np.random.permutation(listA)

        permuted_jaccards.append(jaccard_score(listB, random_list))

    p_value = 1.0 - (
        percentileofscore(permuted_jaccards, jaccard) / 100.0)

    return {
        'jaccard': round(jaccard, 2),
        'jaccard_p_value': p_value,
    }


def get_cross_table(records):

    items = set([
        j
        for i in records
        for j in i
    ])

    cross_table = defaultdict(dict)

    # TODO: a table for string index, contingency table
    for i in items:
        cross_table[i] = {
            i: 0
            for i in items
        }

    for i in records:
        for j in i:
            cross_table[j][j] += 1

        for j, k in combinations(i, 2):
            cross_table[j][k] += 1
            cross_table[k][j] += 1

    return cross_table


def get_binary_mark_list(records):

    items = set([
        j
        for i in records
        for j in i
    ])

    result = defaultdict(list)

    for i in items:
        for j in records:
            result[i].append(i in j)

    return result


def calc_contigency_table(listA, listB):

    match_list = list(zip(listA, listB))

    both = [
        1
        for i, j in match_list
        if i and j
    ]

    i_only = [
        1
        for i, j in match_list
        if i and not j
    ]

    j_only = [
        1
        for i, j in match_list
        if not i and j
    ]

    none = [
        1
        for i, j in match_list
        if not i and not j
    ]

    return {
        'both': len(both),
        'i_only': len(i_only),
        'j_only': len(j_only),
        'none': len(none)
    }


def calc_spearman(listA, listB):
    sp = stats.spearmanr(listA, listB)

    return {
        'spearman_rho': round(sp.statistic, 3),
        'spearman_p_value': sp.pvalue,
    }


def calc_odds_ratio(contigency_table):
    both = contigency_table['both']
    i_only = contigency_table['i_only']
    j_only = contigency_table['j_only']
    none = contigency_table['none']

    i_total = both + i_only
    none_i_total = j_only + none
    # j_total = both + j_only
    # none_j_total = i_only + none

    try:
        r_ratio = (both / i_total) / (j_only / none_i_total)
    except ZeroDivisionError:
        r_ratio = 'inf'

    try:
        odds_ratio = (both / i_only) / (j_only / none)
    except ZeroDivisionError:
        odds_ratio = 'inf'

    return {
        'r_ratio': r_ratio,
        'odds_ratio': odds_ratio,
    }


def calc_holm_Bonferroni(records, column, alpha):

    records.sort(key=itemgetter(column))

    num_total = len(records)

    for idx, r in enumerate(records):
        p = r[column]
        correct_p = alpha / (num_total - idx)
        if p >= correct_p:
            break
        r[f'corrected_{column}'] = correct_p

    return records
