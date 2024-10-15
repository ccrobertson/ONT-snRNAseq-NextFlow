#!/usr/bin/env python
# coding: utf-8


import argparse
import logging

import pysam
import pandas as pd
from scipy.io import mmread

logging.basicConfig(level=logging.DEBUG, format='%(asctime)s - %(levelname)s: %(message)s')


class Cell:

    def __init__(self, barcode):
        self.barcode = barcode
        self.supplementary_alignments = 0
        self.secondary_alignments = 0
        self.total_reads = 0 # excludes secondary/supplementary alignments
        self.uniquely_mapped_reads = 0
        self.chromosome_read_counts = dict() # chrom -> count


    def record_alignment(self, read):
        if read.is_secondary or read.is_supplementary:
            if read.is_secondary:
                self.secondary_alignments += 1
            if read.is_supplementary:
                self.supplementary_alignments += 1
            return 0
        self.total_reads += 1
        if read.mapping_quality > 0:
            self.uniquely_mapped_reads += 1
            chrom = read.reference_name
            if chrom not in self.chromosome_read_counts:
                self.chromosome_read_counts[chrom] = 0
            self.chromosome_read_counts[chrom] += 1
        return 0


    def gather_metrics(self):
        metrics = dict()
        metrics['barcode'] = self.barcode
        metrics['total_reads'] = self.total_reads
        metrics['uniquely_mapped_reads'] = self.uniquely_mapped_reads
        metrics['secondary_alignments'] = self.secondary_alignments
        metrics['supplementary_alignments'] = self.supplementary_alignments
        metrics['fraction_mitochondrial'] = self.chromosome_read_counts['chrM'] / self.uniquely_mapped_reads if 'chrM' in self.chromosome_read_counts else 0
        return metrics

    
def write_mm(prefix, df, features_col, barcodes_col, count_col, all_features=None, all_barcodes=None):
    features = df[features_col].unique() if all_features is None else all_features
    barcodes = df[barcodes_col].unique() if all_barcodes is None else all_barcodes
    if all_features is not None:
        assert(df[features_col].isin(all_features).all())
    if all_barcodes is not None:
        assert(df[barcodes_col].isin(all_barcodes).all())
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


def parse_attribute(attribute_series: pd.Series, attribute_name: str) -> pd.Series:
    """
    Parse the attributes column of a (GENCODE/RefSeq) GTF file.

    Input:
    * a [str]: the attributes element (column 9 of the GTF file)
    * regex [str]: a regular expression that will be iteratively applied to the attribute string to capture attribute key, val pairs. Default should work for GENCODE/RefSeq
    """
    if not isinstance(attribute_series, pd.Series):
        raise TypeError('attribute_series must be a pandas Series')
    if not isinstance(attribute_name, str):
        raise TypeError('attribute_name must be a string')

    return attribute_series.str.extract(f'{attribute_name} "(.*?)"')


def gtf_to_df(gtf: str, parse_attributes: list=None) -> pd.DataFrame:
    df = pd.read_csv(gtf, sep='\t', low_memory=False, header=None, names=['chrom', 'source', 'feature', 'start', 'end', 'score', 'strand', 'frame', 'attributes'], comment='#')
    if parse_attributes is not None:
        for a in parse_attributes:
            df[a] = parse_attribute(df.attributes, a)
    return df


parser = argparse.ArgumentParser()
parser.add_argument('--gene-assignments', required=True)
parser.add_argument('--transcript-assignments', required=True)
parser.add_argument('--barcodes', required=True)
parser.add_argument('--gtf', required=True)
parser.add_argument('--bam-in', required=True)
parser.add_argument('--bam-out', required=True)
parser.add_argument('--count-matrix-prefix', required=True)
parser.add_argument('--qc-metrics-out', required=True)
# a = ['--gene-assignments', '/gpfs/accounts/scjp_root/scjp0/porchard/2024-09-ONT-pipeline-updates/work/small-test-pipeline/output/9266-VD-1/9266-VD-1.gene-assignments.txt']
# a += ['--transcript-assignments', '/gpfs/accounts/scjp_root/scjp0/porchard/2024-09-ONT-pipeline-updates/work/small-test-pipeline/output/9266-VD-1/9266-VD-1.transcript-assignments.txt']
# a += ['--barcodes', '/gpfs/accounts/scjp_root/scjp0/porchard/2024-09-ONT-pipeline-updates/work/small-test-pipeline/output/9266-VD-1/9266-VD-1.corrected-umis.txt']
# a += ['--bam-in', '/gpfs/accounts/scjp_root/scjp0/porchard/2024-09-ONT-pipeline-updates/work/small-test-pipeline/work/24/6d99c52098d376851f93fe6ce4084b/merged.sorted.bam']
# a += ['--bam-out', 'test.bam']
# a += ['--count-matrix-prefix', 'test.']
# a += ['--qc-metrics-out', 'test.metrics.txt']
# args = parser.parse_args(a)
args = parser.parse_args()




logging.info('Reading corrected barcodes and UMIs')
barcode_info = pd.read_csv(args.barcodes, sep='\t', usecols=['read_id', 'CR', 'CB', 'CY', 'UR', 'UB', 'UY'])

logging.info('Reading gene assignments')
gene_assignments = pd.read_csv(args.gene_assignments, sep='\t', header=None, names=['read_id', 'gene_id', 'overlaps_exon'])

logging.info('Reading transcript assignments')
transcript_assignments = pd.read_csv(args.transcript_assignments, sep='\t')

assert(all(transcript_assignments.read_id.isin(gene_assignments.read_id)))
assert(all(barcode_info.read_id.isin(gene_assignments.read_id)))


gtf_df = gtf_to_df(args.gtf)
gtf_df = gtf_df[gtf_df.feature=='gene']
gtf_df['gene_id'] = parse_attribute(gtf_df.attributes, 'gene_id')
gtf_df['gene_name'] = parse_attribute(gtf_df.attributes, 'gene_name')

gene_id_to_gene_name = dict(zip(gtf_df.gene_id, gtf_df.gene_name))
gene_id_to_feature = {gene_id: '{}\t{}\tGene Expression'.format(gene_id, gene_id_to_gene_name[gene_id]) for gene_id in gene_id_to_gene_name}


tags = gene_assignments.merge(transcript_assignments, how='left').merge(barcode_info, how='left').fillna('-').set_index('read_id').rename(columns={'gene_id': 'GX', 'isoform_id': 'TR'})

gene_count_matrix = tags[(tags.GX!='-') & (tags.CB!='-')].groupby(['GX', 'CB']).UB.nunique().rename('umis').reset_index()
gene_count_matrix['feature'] = gene_count_matrix.GX.map(gene_id_to_feature)
write_mm(args.count_matrix_prefix + 'genes.', gene_count_matrix, 'feature', 'CB', 'umis', all_features=list(gene_id_to_feature.values()))

transcript_count_matrix = tags[(tags.TR!='-') & (tags.CB!='-')].groupby(['TR', 'CB']).UB.nunique().rename('umis').reset_index()
write_mm(args.count_matrix_prefix + 'transcripts.', transcript_count_matrix, 'TR', 'CB', 'umis')


umis_per_barcode = gene_count_matrix.groupby('CB').umis.sum().to_dict()
exon_gene_body_ratio = tags.loc[(tags.GX!='-') & (tags.CB!='-'),['GX', 'CB', 'UB', 'overlaps_exon']].groupby(['GX', 'CB', 'UB']).overlaps_exon.any().astype(int).reset_index().groupby('CB').overlaps_exon.mean().rename('exon_gene_body_ratio').to_dict()


# tag bam, and calculate QC metrics
tags = tags.to_dict()

# QC metrics
cells = dict()
no_cell_tag = 0
total_reads = 0

with pysam.AlignmentFile(args.bam_in, 'rb') as bam_in:
    with pysam.AlignmentFile(args.bam_out, 'wb', template=bam_in) as bam_out:
        for read in bam_in.fetch(until_eof=True):
            total_reads += 1
            if total_reads % 1000000 == 0:
                logging.info('Processed {:,} reads'.format(total_reads))                
            if read.is_secondary or read.is_supplementary:
                # secondary reads are filtered out early in the pipeline.
                # supplementary alignments are ignored in e.g. gene assignment, so won't have a GX tag.
                pass
            else:
                for t in ['GX', 'TR', 'CR', 'CB', 'CY', 'UR', 'UB', 'UY']:
                    read.set_tag(t, tags[t][read.query_name])
                barcode = read.get_tag('CB')
                if barcode not in cells:
                    cells[barcode] = Cell(barcode)
                cells[barcode].record_alignment(read)
            bam_out.write(read)


logging.info('Processed {:,} reads and finished writing {}'.format(total_reads, args.bam_out))

logging.info('Writing QC metrics')

with open(args.qc_metrics_out, 'w') as fh:
    print_metrics = ['barcode', 'total_reads', 'uniquely_mapped_reads', 'umis', 'fraction_mitochondrial', 'exon_gene_body_ratio']
    fh.write('\t'.join(print_metrics) + '\n')
    for cell in cells.values():
        metrics = cell.gather_metrics()
        metrics['exon_gene_body_ratio'] = exon_gene_body_ratio[metrics['barcode']] if metrics['barcode'] in exon_gene_body_ratio else 'NA'
        metrics['umis'] = umis_per_barcode[metrics['barcode']] if metrics['barcode'] in umis_per_barcode else 'NA'
        to_print = [str(metrics[i]) for i in print_metrics]
        fh.write('\t'.join(to_print) + '\n')

logging.info('Done')
