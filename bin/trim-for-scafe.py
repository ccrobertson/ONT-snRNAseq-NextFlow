#!/usr/bin/env python
# coding: utf-8


import pysam
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import logging
import argparse

logging.basicConfig(level=logging.DEBUG, format='%(asctime)s - %(levelname)s: %(message)s')


CIGAR_OP = {
    0: 'M',
    1: 'I',
    2: 'D',
    3: 'N',
    4: 'S',
    5: 'H',
    6: 'P',
    7: '=',
    8: 'X',
    9: 'B'
}

CONSUMES_QUERY = {
    0: True,
    1: True,
    2: False,
    3: False,
    4: True,
    5: False,
    6: False,
    7: True,
    8: True
}

CONSUMES_REFERENCE = {
    0: True,
    1: False,
    2: True,
    3: True,
    4: False,
    5: False,
    6: False,
    7: True,
    8: True
}

REV_CIGAR_OP = {v: k for k, v in CIGAR_OP.items()}

UNSUPPORTED_CIGAR_OPS = ['B'] # in the pysam library, but not defined in the BAM spec?

COMPLEMENTS = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A', 'N': 'N'}

translation_table = str.maketrans('ACGTNacgtn', ''.join([COMPLEMENTS[x.upper()] for x in 'ACGTNacgtn']))


def reverse_complement(s):
    return s.translate(translation_table)[::-1]


def remove_hardclipping(read):
    cigartuples = read.cigartuples
    is_hardclipped = False
    if CIGAR_OP[cigartuples[0][0]] == 'H':
        cigartuples = cigartuples[1:]
        is_hardclipped = True
    if CIGAR_OP[cigartuples[-1][0]] == 'H':
        cigartuples = cigartuples[:-1]
        is_hardclipped = True
    if is_hardclipped:
        read.cigartuples = cigartuples
        read.cigarstring = ''.join(['{}{}'.format(i[1], CIGAR_OP[i[0]]) for i in cigartuples])
    return read


def trim_softclipping(read, max_allowed, end='5prime'):
    assert(end in ['5prime', '3prime'])
    assert(isinstance(max_allowed, int))
    assert(max_allowed >= 0)
    cigartuples = read.cigartuples
    cigar_index = 0 if (read.is_forward and end == '5prime') or (not read.is_forward and end == '3prime') else -1
    cigar_op, cigar_length = cigartuples[cigar_index]
    cigar_op = CIGAR_OP[cigar_op]
    if cigar_op in UNSUPPORTED_CIGAR_OPS:
        raise ValueError('Do not support CIGAR operation: {} (read: {})'.format(cigar_op, read.query_name))
    if cigar_op == 'S':
        if cigar_length > max_allowed:
            trim_length = cigar_length - max_allowed
            new_cigartuples = cigartuples.copy()
            if max_allowed > 0:
                new_cigartuples[cigar_index] = (REV_CIGAR_OP['S'], max_allowed)
            else:
                new_cigartuples.pop(cigar_index)
            read.cigartuples = new_cigartuples
            read.cigarstring = ''.join(['{}{}'.format(i[1], CIGAR_OP[i[0]]) for i in new_cigartuples])
            if (read.is_forward and end == '5prime') or (not read.is_forward and end == '3prime'):
                q = read.query_qualities[trim_length:] # Assigning to query_sequence will invalidate any quality scores. Thus, to in-place edit the sequence and quality scores, copies of the quality scores need to be taken
                read.query_sequence = read.query_sequence[trim_length:]
                read.query_qualities = q
            else:
                q = read.query_qualities[:-trim_length] # Assigning to query_sequence will invalidate any quality scores. Thus, to in-place edit the sequence and quality scores, copies of the quality scores need to be taken
                read.query_sequence = read.query_sequence[:-trim_length]
                read.query_qualities = q
    return read


def amount_of_softclipping(read, end='5prime'):
    assert(end in ['5prime', '3prime'])
    cigartuples = read.cigartuples
    cigar_index = 0 if (read.is_forward and end == '5prime') or (not read.is_forward and end == '3prime') else -1
    cigar_op, cigar_length = cigartuples[cigar_index]
    cigar_op = CIGAR_OP[cigar_op]
    if cigar_op in UNSUPPORTED_CIGAR_OPS:
        raise ValueError('Do not support CIGAR operation: {} (read: {})'.format(cigar_op, read.query_name))
    if cigar_op == 'S':
        return cigar_length
    else:
        return 0


def trim(read, final_length, trim_from_end='3prime'):
    # TODO: could add a hardclipping CIGAR op to denote the trimming. But as this is really only meant for SCAFE, not bothering with that for now.
    assert(trim_from_end in ['5prime', '3prime'])
    assert(isinstance(final_length, int))
    assert(final_length >= 0)
    if len(read.query_sequence) <= final_length:
        return read
    new_cigartuples = []
    new_pos = read.reference_start
    if (trim_from_end == '3prime' and read.is_reverse) or (trim_from_end == '5prime' and not read.is_reverse): # trimming from the left when visualizing the BAM file
        new_qualities = read.query_qualities[-final_length:]
        new_sequence = read.query_sequence[-final_length:]
        cumulative_length = 0
        query_cigar_complete = False

        for i in read.cigartuples[::-1]:
            cigar_op, cigar_length = i
            if query_cigar_complete:
                if CONSUMES_REFERENCE[cigar_op]:
                    new_pos += cigar_length
            else:
                if CONSUMES_QUERY[cigar_op]:
                    if cumulative_length + cigar_length < final_length:
                        new_cigartuples.append(i)
                        cumulative_length += cigar_length
                    elif cumulative_length + cigar_length == final_length:
                        new_cigartuples.append(i)
                        query_cigar_complete = True
                    else:
                        new_cigartuples.append((cigar_op, final_length - cumulative_length))
                        query_cigar_complete = True
                        new_pos += (cigar_length - (final_length - cumulative_length))
                else:
                    new_cigartuples.append(i)
        new_cigartuples.reverse()

    else: # trimming from the right. No need to update mapping position
        new_qualities = read.query_qualities[:final_length]
        new_sequence = read.query_sequence[:final_length]
        cumulative_length = 0
        for i in read.cigartuples:
            cigar_op, cigar_length = i
            if CONSUMES_QUERY[cigar_op]:
                if cumulative_length + cigar_length < final_length:
                    new_cigartuples.append(i)
                    cumulative_length += cigar_length
                elif cumulative_length + cigar_length == final_length:
                    new_cigartuples.append(i)
                    break
                else:
                    new_cigartuples.append((cigar_op, final_length - cumulative_length))
                    break
            else:
                new_cigartuples.append(i)

    read.query_sequence = new_sequence
    read.query_qualities = new_qualities
    read.cigartuples = new_cigartuples
    read.cigarstring = ''.join(['{}{}'.format(i[1], CIGAR_OP[i[0]]) for i in new_cigartuples])
    read.reference_start = new_pos # todo: any need to update reference_end? I believe no because that's not actually written to the BAM file, just inferred by pysam I think.
    
    return read


def trim_n_bases(read, n, trim_from_end='3prime'):
    return trim(read, len(read.query_sequence) - n, trim_from_end)


parser = argparse.ArgumentParser(description="Preprocess 3' ONT-sequenced RNA-seq data for SCAFE.")
parser.add_argument('-i', '--input-bam', type=str, required=True, help='Input BAM file')
parser.add_argument('-o', '--output-bam', type=str, required=True, help='Output BAM file')
parser.add_argument('-m', '--max-softclipping', type=int, default=50, help="Maximum number of softclipped bases allowed at the 5' end of the read (reads with more than this will be skipped) (default: 50)")
parser.add_argument('-t', '--trim-to', type=int, default=100, help='Trim reads to this maximum length (default: 100 bp)')
parser.add_argument('-k', '--keep-tags', type=str, nargs='+', default=['CB', 'CR', 'CY', 'UB', 'UR', 'UY', 'GN', 'TR'], help="Tags to keep in the output BAM file (default: CB CR CY UB UR UY GN TR); should correspond to tags that won't be invalidated by trimming")
parser.add_argument('-g', '--ggg-mismatches-allowed', type=int, default=1, help="Number of mismatches allowed in the expected 5' GGG sequence (reads with more than this will be skipped) (default: 1)")
args = parser.parse_args()

BAM = args.input_bam
OUT_BAM = args.output_bam
MAX_SOFTCLIPPING = args.max_softclipping
TRIM_TO = args.trim_to
KEEP_TAGS = args.keep_tags
GGG_MISMATCHES_ALLOWED = args.ggg_mismatches_allowed

# BAM = '/scratch/scjp_root/scjp0/porchard/2023-HSM-ONT/work/rnaseq-b27ac1d/output/9266-VD-1/tagged.bam'
# MAX_SOFTCLIPPING = 20 # allow up to this many softclipped bases at the 5' end of the read
# TRIM_TO = 100 # trim reads to this maximum length.
# KEEP_TAGS = ['CB', 'CR', 'CY', 'UB', 'UR', 'UY', 'GN', 'TR']
# GGG_MISMATCHES_ALLOWED = 1

assert(TRIM_TO > MAX_SOFTCLIPPING) # at the end should also verify that "matches" remain in the CIGAR string

skipped_softclipping = 0
skipped_no_GGG = 0
skipped_unmapped = 0

total_reads = 0
with pysam.AlignmentFile(BAM, 'rb') as f:
    with pysam.AlignmentFile(OUT_BAM, 'wb', header=f.header) as out:
        for read in f.fetch(until_eof=True):
            total_reads += 1
            if total_reads % 1000000 == 0:
                logging.info('Processed {:,} reads so far...'.format(total_reads))
            if read.is_unmapped:
                skipped_unmapped += 1
                continue
            read = remove_hardclipping(read)
            if amount_of_softclipping(read) > MAX_SOFTCLIPPING:
                skipped_softclipping += 1
                continue
            # check for GGG at the beginning of the read
            sequence = read.query_sequence if read.is_forward else reverse_complement(read.query_sequence)
            first_bases = sequence[0:3]
            matches = sum([first_bases[i] == 'G' for i in range(3)])
            if matches < (3 - GGG_MISMATCHES_ALLOWED):
                skipped_no_GGG += 1
                continue
            # trim the GGG
            read = trim_n_bases(read, 3, trim_from_end='5prime')
            # trim the read to the desired length
            read = trim(read, TRIM_TO, trim_from_end='3prime')
            # trimming will invalidate tags, so need to remove them
            for (tag, tag_value) in read.get_tags():
                if tag not in KEEP_TAGS:
                    read.set_tag(tag, None)
            out.write(read)



logging.info('Processed {:,} reads total.'.format(total_reads))
logging.info('Skipped {:,} reads ({}%) total. Breakdown below.'.format(skipped_unmapped + skipped_softclipping + skipped_no_GGG, round((skipped_unmapped + skipped_softclipping + skipped_no_GGG)/total_reads*100, 3)))
logging.info('Skipped {:,} reads ({}%) due to being unmapped'.format(skipped_unmapped, round(skipped_unmapped/total_reads*100, 3)))
logging.info('Skipped {:,} reads ({}%) due to too much softclipping'.format(skipped_softclipping, round(skipped_softclipping/total_reads*100, 3)))
logging.info('Skipped {:,} reads ({}%) due to no GGG at the beginning (allowing for {} mismatch(es))'.format(skipped_no_GGG, round(skipped_no_GGG/total_reads*100, 3), GGG_MISMATCHES_ALLOWED))

