#!/usr/bin/env python
# coding: utf-8


import pandas as pd
import time
import correct_barcodes_levenshtein
import sys
import logging
import argparse

parser = argparse.ArgumentParser()
parser.add_argument('--tag-file', required=True)
parser.add_argument('--umi-counts', nargs='+', required=True, help='UMI counts file. Must be a tab-separated file with two columns: CR and n_umis (include header)')
parser.add_argument('--method', required=True)
args = parser.parse_args()

logging.basicConfig(level=logging.DEBUG, format='%(asctime)s - %(levelname)s: %(message)s')

TAG_FILE = args.tag_file
UMI_COUNTS = args.umi_counts
METHOD = args.method

METHOD_CHOICES = ['phred', 'error_rate', 'edit_distance']
if not METHOD in METHOD_CHOICES:
    raise ValueError('--method must be one of: '.format(', '.join(METHOD_CHOICES)))


umi_counts = pd.concat([pd.read_csv(f, sep='\t', header=0) for f in UMI_COUNTS])
umi_counts.columns = ['CR', 'umis']
umi_counts = umi_counts.groupby('CR').sum().reset_index()
umi_counts = umi_counts.set_index('CR').umis.to_dict()

whitelist = list(sorted(umi_counts.keys()))
whitelist_set = set(whitelist)

tags = pd.read_csv(TAG_FILE, sep='\t')

whitelist_match = tags[tags.CR.isin(whitelist_set)]
no_whitelist_match = tags[~tags.CR.isin(whitelist_set)]


logging.info('Finding close whitelist barcodes for each non-whitelisted barcode')
start = time.time()
matches = correct_barcodes_levenshtein.find_matches(whitelist, no_whitelist_match.CR.unique().tolist(), 3)
end = time.time()
logging.info('Elapsed time: {:,} seconds'.format(end-start))



logging.info('Correcting barcodes...')
start = time.time()

if METHOD == 'edit_distance':
    no_whitelist_match['CB'] = no_whitelist_match.CR.map(lambda x: correct_barcodes_levenshtein.correct_barcode_edit_distance_only(x, matches[x], max_ed=2, min_ed_diff=2))
elif METHOD == 'error_rate':
    no_whitelist_match['CB'] = no_whitelist_match.CR.map(lambda x: correct_barcodes_levenshtein.correct_barcode_cellranger(x, matches[x], {i[0]: umi_counts[i[0]] for i in matches[x]}, max_ed=2, error_rate=0.05, prob_threshold=0.95))
elif METHOD == 'phred':
    no_whitelist_match['CB'] = [correct_barcodes_levenshtein.correct_barcode_cellranger(cr, matches[cr], {i[0]: umi_counts[i[0]] for i in matches[cr]}, phred=cy, max_ed=2, error_rate=0.05, prob_threshold=0.95) for cr, cy in zip(no_whitelist_match.CR, no_whitelist_match.CY)] 

end = time.time()
logging.info('Elapsed time: {:,} seconds'.format(end-start))



whitelist_match['CB'] = whitelist_match.CR
corrected = pd.concat([whitelist_match, no_whitelist_match])

corrected.to_csv(sys.stdout, sep='\t', index=False)