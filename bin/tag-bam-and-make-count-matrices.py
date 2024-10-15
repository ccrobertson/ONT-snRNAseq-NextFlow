#!/usr/bin/env python
# coding: utf-8

import argparse
import logging

import pysam
import pandas as pd

logging.basicConfig(level=logging.DEBUG, format='%(asctime)s - %(levelname)s: %(message)s')


def write_mm(prefix, df, features_col, barcodes_col, count_col):
    features = df[features_col].unique()
    barcodes = df[barcodes_col].unique()
    feature_to_feature_index = {feature: i for i, feature in enumerate(features, 1)}
    barcode_to_barcode_index = {barcode: i for i, barcode in enumerate(barcodes, 1)}

    df['feature_index'] = df[features_col].map(feature_to_feature_index)
    df['barcode_index'] = df[barcodes_col].map(barcode_to_barcode_index)

    with open(f'{prefix}barcodes.tsv', 'w') as fh:
        for b in barcodes:
            fh.write(b + '\n')

    with open(f'{prefix}features.tsv', 'w') as fh:
        for f in features:
            fh.write(f + '\n')

    with open(f'{prefix}matrix.mtx', 'w') as fh:
        # generate the header
        header = ['%%MatrixMarket matrix coordinate integer general', '%', ' '.join([str(len(features)), str(len(barcodes)), str(len(df))])]
        fh.write('\n'.join(header) + '\n')

        # write the final file
        for feature_index, barcode_index, count in zip(df['feature_index'], df['barcode_index'], df[count_col]):
            fh.write(f'{feature_index} {barcode_index} {count}\n')

    return True


parser = argparse.ArgumentParser()
parser.add_argument('--gene-assignments', required=True)
parser.add_argument('--transcript-assignments', required=True)
parser.add_argument('--barcodes', required=True)
parser.add_argument('--bam-in', required=True)
parser.add_argument('--bam-out', required=True)
parser.add_argument('--count-matrix-prefix', required=True)
# a = ['--gene-assignments', '/gpfs/accounts/scjp_root/scjp0/porchard/2024-09-ONT-pipeline-updates/work/small-test-pipeline/output/9266-VD-1/9266-VD-1.gene-assignments.txt']
# a += ['--transcript-assignments', '/gpfs/accounts/scjp_root/scjp0/porchard/2024-09-ONT-pipeline-updates/work/small-test-pipeline/output/9266-VD-1/9266-VD-1.transcript-assignments.txt']
# a += ['--barcodes', '/gpfs/accounts/scjp_root/scjp0/porchard/2024-09-ONT-pipeline-updates/work/small-test-pipeline/output/9266-VD-1/9266-VD-1.corrected-umis.txt']
# a += ['--bam-in', '/gpfs/accounts/scjp_root/scjp0/porchard/2024-09-ONT-pipeline-updates/work/small-test-pipeline/work/24/6d99c52098d376851f93fe6ce4084b/merged.sorted.bam']
# a += ['--bam-out', 'test.bam']
# a += ['--count-matrix-prefix', 'test.']
#args = parser.parse_args(a)
args = parser.parse_args()



logging.info('Reading corrected barcodes and UMIs')
barcode_info = pd.read_csv(args.barcodes, sep='\t', usecols=['read_id', 'CR', 'CB', 'CY', 'UR', 'UB', 'UY'])

logging.info('Reading gene assignments')
gene_assignments = pd.read_csv(args.gene_assignments, sep='\t', header=None, names=['read_id', 'gene_id', 'overlaps_exon']).drop(columns=['overlaps_exon'])

logging.info('Reading transcript assignments')
transcript_assignments = pd.read_csv(args.transcript_assignments, sep='\t')


assert(all(transcript_assignments.read_id.isin(gene_assignments.read_id)))
assert(all(barcode_info.read_id.isin(gene_assignments.read_id)))


tags = gene_assignments.merge(transcript_assignments, how='left').merge(barcode_info, how='left').fillna('-').set_index('read_id').rename(columns={'gene_id': 'GX', 'isoform_id': 'TR'})

gene_count_matrix = tags[(tags.GX!='-') & (tags.CB!='-')].groupby(['GX', 'CB']).UB.nunique().rename('umis').reset_index()
write_mm(args.count_matrix_prefix + 'genes.', gene_count_matrix, 'GX', 'CB', 'umis')
transcript_count_matrix = tags[(tags.TR!='-') & (tags.CB!='-')].groupby(['TR', 'CB']).UB.nunique().rename('umis').reset_index()
write_mm(args.count_matrix_prefix + 'transcripts.', transcript_count_matrix, 'TR', 'CB', 'umis')


# tag bam
tags = tags.to_dict()

with pysam.AlignmentFile(args.bam_in, 'rb') as bam_in:
    with pysam.AlignmentFile(args.bam_out, 'wb', template=bam_in) as bam_out:
        count = 0
        for read in bam_in.fetch(until_eof=True):
            count += 1
            if count % 1000000 == 0:
                logging.info('Processed {:,} reads'.format(count))
            if read.is_secondary or read.is_supplementary:
                pass
            else:
                for t in ['GX', 'TR', 'CR', 'CB', 'CY', 'UR', 'UB', 'UY']:
                    read.set_tag(t, tags[t][read.query_name])
            bam_out.write(read)