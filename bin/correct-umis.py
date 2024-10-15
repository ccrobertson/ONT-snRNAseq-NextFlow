#!/usr/bin/env python
# coding: utf-8

import sys
import collections
import itertools
from pathlib import Path
import argparse

from editdistance import eval as edit_distance
import numpy as np
import pandas as pd
from umi_tools import UMIClusterer


parser = argparse.ArgumentParser()
parser.add_argument('--corrected-barcodes', required=True)
parser.add_argument('--read-to-gene-assignments', required=True)
args = parser.parse_args()


def get_adj_list_directional_lev(self, umis, counts, threshold=2):
    """Use Levenshtein distance for UMIclustering instead of hamming.

    This function is to monkey-patch UMIClusterer._get_adj_list_directional
    """
    adj_list = {umi: [] for umi in umis}
    iter_umi_pairs = itertools.combinations(umis, 2)
    for umi1, umi2 in iter_umi_pairs:
        if edit_distance(umi1, umi2) <= threshold:
            if counts[umi1] >= (counts[umi2] * 2) - 1:
                adj_list[umi1].append(umi2)
            if counts[umi2] >= (counts[umi1] * 2) - 1:
                adj_list[umi2].append(umi1)

    return adj_list


def umi_clusterer_call(self, umis, threshold):
    """To replace UMIClusterer.__call__.

    Use this method to monkey-patch the UMICluterer.__call__ in order to remove the
    nessesity for all UMIs to be the same length, allowing for deletions in the UMIs.

    https://github.com/CGATOxford/UMI-tools/blob/c3ead0792ad590822ca72239ef01b8e559802d
    a9/umi_tools/network.py#L350
    """
    counts = umis
    umis = list(umis.keys())

    self.positions += 1
    number_of_umis = len(umis)
    self.total_umis_per_position += number_of_umis

    if number_of_umis > self.max_umis_per_position:
        self.max_umis_per_position = number_of_umis

    adj_list = self.get_adj_list(umis, counts, threshold)
    clusters = self.get_connected_components(umis, adj_list, counts)
    final_umis = [list(x) for x in self.get_groups(clusters, adj_list, counts)]

    return final_umis


# Monkey-patch the umi-tools clusterer with a modified method
# using Levenshtein instead of Hamming distance.
UMIClusterer._get_adj_list_directional = get_adj_list_directional_lev
UMIClusterer.__call__ = umi_clusterer_call


def cluster(umis):
    """Cluster UMIs.

    Search for UMI clusters within subsets of reads sharing the same corrected barcode
    and gene. In this way the search space is dramatically reduced.

    We are using the UMI-tools directional deduplication method (modified to use
    Levenshtein distance). Connections between nodes within a cluster are generated
    based on edit distance threshold and whether node A counts >= (2* node B counts).
    https://umi-tools.readthedocs.io/en/latest/the_methods.html

    """
    if len(umis) == 1:  # early return
        return umis
    clusterer = UMIClusterer(cluster_method="directional")
    umi_counts = collections.Counter(umis)
    clusters = clusterer(umi_counts, threshold=2)

    if len(clusters) == len(umis):  # no corrections
        return umis

    # create list of corrections
    umi_map = dict()
    for clust in clusters:
        if len(clust) > 1:
            for umi in clust[1:]:
                umi_map[umi] = clust[0]
    if len(umi_map) > 0:  # pd.Series.replace is weird slow
        umis = umis.replace(umi_map)
    return umis


def create_region_name(row, ref_interval):
    """Create a fake gene name from alignment coordinates."""
    # The idea here is to slice the reference into a grid and label reads
    # with the chunk that they overlap. Reads intersecting the same chunk
    # are then grouped together.
    midpoint = int((row.start + row.end) / 2)
    interval_start = int(np.floor(midpoint / ref_interval) * ref_interval)
    interval_end = int(np.ceil(midpoint / ref_interval) * ref_interval)
    gene = f"{row.chr}_{interval_start}_{interval_end}"
    return gene


def cluster_dataframe(df, ref_interval):
    """Process records from tags file."""
    # Create column to keep track of non-assigned genes
    df['no_gene'] = False
    df_no_gene = df.loc[df.gene == '-']
    if len(df_no_gene) > 0:
        # Create a temporary gene name based on chr and location.
        regions = df_no_gene.apply(
            create_region_name, axis=1, args=(ref_interval,))
        df.loc[regions.index, 'gene'] = regions
        df.loc[df.index.isin(regions.index), 'no_gene'] = True
    # Create gene/cell index for subsetting reads prior to clustering.
    df["gene_cell"] = df["gene"] + ":" + df["CB"]
    df['read_id'] = df.index
    df.set_index('gene_cell', inplace=True, drop=True)
    # UB: corrected UMI tag
    groups = df.groupby("gene_cell")["UR"]
    df["UB"] = groups.transform(cluster)
    df.set_index('read_id', drop=True, inplace=True)
    # Reset unassigned genes to '-'
    df.loc[df.no_gene, 'gene'] = '-'
    df.drop(columns='no_gene', inplace=True)
    return df



REF_INTERVAL = 1000


barcodes = pd.read_csv(args.corrected_barcodes, sep='\t', index_col=None)
barcodes.CB = barcodes.CB.fillna('-')
barcodes.UR = barcodes.UR.fillna('-')


gene_assignments = pd.read_csv(args.read_to_gene_assignments, sep='\t', header=None, names=['read_id', 'gene', 'overlaps_exon'])


assert(barcodes.read_id.value_counts().max() == 1)
assert(gene_assignments.read_id.value_counts().max() == 1)

assert(barcodes.read_id.isin(gene_assignments.read_id).all())


# TODO: this should be true, after fixes
# assert(len(barcodes) == len(gene_assignments))
# some are missing from the barcodes, not sure why. Need to check.
barcodes = barcodes.merge(gene_assignments, on='read_id', how='left')


for_umi_correction = barcodes[(barcodes.CB != '-') & (barcodes.UR != '-')]

for_umi_correction = cluster_dataframe(for_umi_correction.set_index('read_id'), REF_INTERVAL).reset_index()

output = pd.concat([for_umi_correction, barcodes[(barcodes.CB == '-') | (barcodes.UR == '-')]])
assert(len(output) == len(barcodes))

output.to_csv(sys.stdout, sep='\t')