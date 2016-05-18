#! /usr/bin/local/python
# -*- coding: utf-8 -*-

import re
import logging
from collections import OrderedDict


__author__ = "Ryan Abo"
__copyright__ = "Copyright 2015, Ryan Abo"
__email__ = "ryanabo@gmail.com"
__license__ = "MIT"


def find_reads(mer, read_items, used_reads, order='for' ):

    '''Return a list of tuples containing information from reads with the kmer sequence.

    First search all the read sequences for the given kmer sequence. Then,
    filter out used reads and order them according to position of the kmer
    sequence in the read sequence.

    Args:
        kmerSeq: String of kmer sequence.
        readItems: List of fq_recs (key, value) tuples.
        usedReads: Set of read IDs that have been previously used.
        order: String indicating how the list of the identified reads
               should be ordered.
    Returns:
        kmerReads: List of tuples containing:
                    1. read object,
                    2. start position of kmer match in read seq
                    3. Boolean that a match was found.
                    4. Length of the read sequence.
                    5. Number of reads with this sequence.
    '''

    mer_reads = [] 
    mr = filter(lambda x: x[2], map(read_search, [mer]*len(read_items), read_items))
    ids = map(lambda x: x[0].id, mr)
    filt_ids = set(ids) - set(used_reads)
    matched_reads = filter(lambda x: (x[0].id in filt_ids), mr)

    if order == 'rev': 
        mer_reads = sorted(matched_reads, key=lambda z: (-z[1], -z[3])) 
    else: 
        mer_reads = sorted(matched_reads, key=lambda z: (z[1], -z[3])) 
    return mer_reads

def read_search(mer, read_values):

    '''Return a tuple containing information regarding the alignment of the kmerSeq
    in a sequence read.

    This uses regex searching function re.search to determine if the kmerSeq
    is contained in the read sequence. If so, then it returns a 5 element
    tuple about information regarding this alignment. If no match, then return
    a 3 element tuple with None values.

    Args:
        kmerSeq: String of kmer sequence.
        readItems: List of fq_recs (key, value) tuples.
    Returns:
        searchResult: Tuple of result information.
                       1. read object,
                       2. start position of kmer match in read seq
                       3. Boolean that a match was found.
                       4. Length of the read sequence.
                       5. Number of reads with this sequence.
    '''

    r = (None, None, None) 
    seq, reads = read_values
    x = re.search(mer, seq) 
    if x:
        r = (reads[0], x.start(), True, len(reads[0].seq), len(reads)) 
    return r