#!/usr/bin/env python
# coding: utf-8


import pandas as pd
import rapidfuzz
import numpy as np
from scipy.stats import binom
import sys
import logging

class TrieNode:
    def __init__(self):
        self.children = {}
        self.is_end_of_word = False


class Trie:
    def __init__(self):
        self.root = TrieNode()

    def insert(self, word):
        current = self.root
        for letter in word:
            if letter not in current.children:
                current.children[letter] = TrieNode()
            current = current.children[letter]
        current.is_end_of_word = True

    def build_trie(self, words):
        for word in words:
            self.insert(word)


def trie_search(trie1, trie2, max_distance):
    """Given two tries, returns all pairs of words that are within a given edit distance of each other. Words must all be of the same length."""
    return recursive_trie_search(trie1.root, trie2.root, '', '', max_distance)


def recursive_trie_search(node1, node2, prefix1, prefix2, max_distance):
    matches = []
    d = rapidfuzz.distance.Levenshtein.distance(prefix1, prefix2) # should be able to make more efficient by caching a Levensthein computation matrix
    
    if d > max_distance:
        return matches

    if node1.is_end_of_word and node2.is_end_of_word:
        if d <= max_distance:
            matches.append((prefix1, prefix2, d))
    for letter1 in node1.children:
        for letter2 in node2.children:
            matches += recursive_trie_search(node1.children[letter1], node2.children[letter2], prefix1 + letter1, prefix2 + letter2, max_distance)
    
    return matches


def find_matches_rapidfuzz(whitelist, barcodes, max_distance, batch_size=1000):
    if not isinstance(whitelist, list):
        raise ValueError('whitelist must be a list')
    if not isinstance(barcodes, list):
        raise ValueError('barcodes must be a list')
    if not isinstance(max_distance, int):
        raise ValueError('max_distance must be an integer')
    
    to_correct = barcodes.copy()

    matches = []

    while len(to_correct) > 0:
        next_batch = to_correct[:min([batch_size, len(to_correct)])]
        to_correct = [i for i in to_correct if i not in next_batch]

        distances = rapidfuzz.process.cdist(next_batch, whitelist, scorer=rapidfuzz.distance.Levenshtein.distance, score_cutoff=max_distance)
        x = pd.DataFrame(np.argwhere(distances<=max_distance), columns=['CR_index', 'WL_index'])
        x['ED'] = [distances[cr_index,wl_index] for cr_index, wl_index in zip(x.CR_index, x.WL_index)]
        x['CR'] = np.array(next_batch)[x.CR_index.to_list()]
        x['WL_match'] = np.array(whitelist)[x.WL_index.to_list()]
        matches += [(wl_match, cr, ed) for wl_match, cr, ed in zip(x.WL_match, x.CR, x.ED)]

    return matches



def find_matches(whitelist, barcodes, max_distance):
    if not isinstance(whitelist, list):
        raise ValueError('whitelist must be a list')
    if not isinstance(barcodes, list):
        raise ValueError('barcodes must be a list')
    if not isinstance(max_distance, int):
        raise ValueError('max_distance must be an integer')
    
    # infer expected barcode length
    expected_barcode_length = len(whitelist[0])
    for w in whitelist:
        if len(w) != expected_barcode_length:
            raise ValueError('All barcodes in the whitelist must be the same length')

    # the trie search is fast if have many items in the whitelist and many barcodes to correct
    # but it only works when the barcodes to correct are all the same length as the expected barcode length
    # so we split the barcodes into two groups: those that are the correct length and those that are not
    # we use the trie search to find matches for the barcodes that are the correct length
    # and use rapidfuzz to find matches for the barcodes that are not the correct
    all_to_correct_full_length = [x for x in barcodes if len(x) == expected_barcode_length]
    all_to_correct_incorrect_length = [x for x in barcodes if len(x) != expected_barcode_length]
    
    whitelist_trie = Trie()
    whitelist_trie.build_trie(whitelist)

    barcodes_trie = Trie()
    barcodes_trie.build_trie(all_to_correct_full_length)
    
    matches_full_length = trie_search(whitelist_trie, barcodes_trie, max_distance)
    matches_incorrect_length = find_matches_rapidfuzz(whitelist, all_to_correct_incorrect_length, max_distance)
    matches_all = {i: [] for i in barcodes} # barcode --> [(potential_correction_1, edit_distance_1), (potential_correction_2, edit_distance_2), ...]

    for w, c, d in (matches_full_length + matches_incorrect_length):
        matches_all[c].append((w, d))
    
    for i in barcodes:
        matches_all[i] = sorted(matches_all[i], key=lambda x: x[1])
    
    return matches_all


def correct_barcode_edit_distance_only(barcode, potential_corrections, max_ed=2, min_ed_diff=2):
    for i in potential_corrections:
        if len(i) != 2:
            raise ValueError('Each potential correction must be a tuple of length 2 (correction, edit_distance)')
    if len(potential_corrections) == 0:
        return None
    elif len(potential_corrections) == 1:
        return potential_corrections[0][0] if potential_corrections[0][1] <= max_ed else None
    elif len(potential_corrections) > 1:
        return potential_corrections[0][0] if ((potential_corrections[0][1] <= max_ed) and (potential_corrections[0][1] - potential_corrections[1][1] >= min_ed_diff)) else None


def prob_of_error(uncorrected, potential_correction, phred, indel_error_rate=0.05):
    ops = rapidfuzz.distance.Levenshtein.editops(uncorrected, potential_correction)
    n_insertions = sum([i.tag == 'insert' for i in ops])
    n_deletions = sum([i.tag == 'delete' for i in ops])
    substition_positions = [i.src_pos for i in ops if i.tag == 'replace']

    total_p = 1
    for i in range(n_insertions + n_deletions):
        total_p *= indel_error_rate
    for i in substition_positions:
        q = ord(phred[i]) - 33
        p = 10**(-q/10)
        total_p *= p
    return total_p if total_p != 1 else None


def correct_barcode_cellranger(barcode, potential_corrections, umis_per_cb, phred=None, max_ed=2, error_rate=0.05, prob_threshold=0.95):
    for i in potential_corrections:
        if len(i) != 2:
            raise ValueError('Each potential correction must be a tuple of length 2 (correction, edit_distance)')
        if i[0] not in umis_per_cb:
            raise ValueError(f'UMI counts must be provided for each potential correction. barcode={barcode}, corrections={potential_corrections}')
    if len(potential_corrections) == 0:
        return None
    elif len(potential_corrections) == 1:
        return potential_corrections[0][0] if potential_corrections[0][1] <= max_ed else None
    elif len(potential_corrections) > 1:
        cb_length = len(potential_corrections[0][0])
        umis = np.array([umis_per_cb[i[0]] for i in potential_corrections])
        ed = np.array([i[1] for i in potential_corrections])
        error_probability = binom.pmf(ed, cb_length, error_rate) if phred is None else [prob_of_error(barcode, i[0], phred, indel_error_rate=error_rate) for i in potential_corrections]
        prior = umis / sum(umis)
        prob_unnormalized = error_probability * prior
        prob_normalized = prob_unnormalized / sum(prob_unnormalized)
        if max(prob_normalized) >= prob_threshold:
            return potential_corrections[np.argmax(prob_normalized)][0]
        else:
            return None

