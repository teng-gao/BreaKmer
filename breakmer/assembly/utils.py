#! /usr/bin/local/python
# -*- coding: utf-8 -*-

'''
BreaKmer assembly utils module
'''

import re

__author__ = "Ryan Abo"
__copyright__ = "Copyright 2015, Ryan Abo"
__email__ = "ryanabo@gmail.com"
__license__ = "MIT"


class KmerUnit(object):

    '''Class to track value associated with a particular kmer sequence.

    Attributes:
    '''

    def __init__(self, seq, counts, kmer_seq_set, kmer_len):
        self.seq = seq
        self.counts = counts
        self.kmer_seq_set = kmer_seq_set
        self.kmer_len = kmer_len


def find_reads(kmer_seq, read_items, used_reads, order='for'):

    '''Return a list of tuples containing information from reads with the kmer sequence.

    First search all the read sequences for the given kmer sequence. Then,
    filter out used reads and order them according to position of the kmer
    sequence in the read sequence.

    Args:
        kmer_seq (str):    String of kmer sequence.
        read_items (list): List of fq_recs (key, value) tuples.
        used_reads (set):  Set of read IDs that have been previously used.
        order (str):       String indicating how the list of the identified reads
                           should be ordered.
    Returns:
        kmer_reads: List of tuples containing:
                    1. read object,
                    2. start position of kmer match in read seq
                    3. Boolean that a match was found.
                    4. Length of the read sequence.
                    5. Number of reads with this sequence.
    '''

    kmer_reads = []
    mapped_reads = filter(lambda x: x[2], map(read_search, [kmer_seq]*len(read_items), read_items))
    ids = map(lambda x: x[0].id, mapped_reads)
    filt_ids = set(ids) - set(used_reads)
    matched_reads = filter(lambda x: (x[0].id in filt_ids), mapped_reads)

    if order == 'rev':
        kmer_reads = sorted(matched_reads, key=lambda z: (-z[1], -z[3]))
    else: 
        kmer_reads = sorted(matched_reads, key=lambda z: (z[1], -z[3]))
    return kmer_reads

def read_search(kmer_seq, read_items):

    '''Return a tuple containing information regarding the alignment of the kmer sequence
    in a sequence read.

    This uses regex searching function re.search to determine if the kmerSeq
    is contained in the read sequence. If so, then it returns a 5 element
    tuple about information regarding this alignment. If no match, then return
    a 3 element tuple with None values.

    Args:
        kmer_seq (str):     Kmer sequence.
        read_items (list):  Fastq records in (key, value) tuples.
    Returns:
        searchResult: Tuple of result information.
                       1. read object,
                       2. start position of kmer match in read seq
                       3. Boolean that a match was found.
                       4. Length of the read sequence.
                       5. Number of reads with this sequence.
    '''

    search_result = (None, None, None)
    read_seq, fq_reads = read_items
    found_kmer_seq = re.search(kmer_seq, read_seq)
    if found_kmer_seq:
        search_result = (fq_reads[0], found_kmer_seq.start(), True, len(fq_reads[0].seq), len(fq_reads))
    return search_result