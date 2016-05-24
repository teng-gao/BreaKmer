#! /usr/bin/local/python
# -*- coding: utf-8 -*-

'''
BreaKmer assembler module

Assembles contigs
'''


from collections import OrderedDict
import breakmer.utils as utils
import breakmer.assembly.utils as assembly_utils
import breakmer.assembly.contig as contig


__author__ = "Ryan Abo"
__copyright__ = "Copyright 2015, Ryan Abo"
__email__ = "ryanabo@gmail.com"
__license__ = "MIT"


def init_assembly(kmers, fq_recs, kmer_len, rc_thresh, read_len):

    '''Entry function for assemblying a contiguous sequence from
    a pool of sample only kmers and the reads that contain them.

    A kmer tracker object is instantiated containing all the kmer seqs and
    their associated counts. These are sorted by

    Args:
      kmers (dict):      Sample-only kmers stored as key = kmer string, value = frequency observed in reads
      fq_recs (dict):    Sequence values stored as keys and a list of fq_read objects as values.
      kmer_len (int):    Number of base pairs in kmer string.
      rc_thresh (int):   The minimum number of reads supporting a contig needed to keep the contig.
      read_len (int):    Number of base pairs in the sequence reads.
    Returns:
      contigs (list):    Contig objects.
    '''

    logging_name = 'breakmer.assembly.assembler'
    contigs = []

    # Return if there are no kmers to build a contig.
    if len(kmers) == 0:
        utils.log(logging_name, 'info', 'No kmers to build contigs, returning.')
        return contigs

    # Store kmers in KmerTracker object.
    kmer_tracker = KmerTracker()
    for kmer_str in kmers:
        kmer_tracker.add_kmer(kmer_str, kmers[kmer_str])

    # While there are kmers to analyze continue to build contigs.
    contig_buffer = ContigBuffer()
    # Sort all of the kmers by frequency (descending) and store in order
    kmer_tracker.set_all_kmer_values()
    # Check if there are any kmers left to seed the build process
    while kmer_tracker.has_kmers():
        # Update the set of kmers to consider
        kmer_tracker.update_kmer_set()
        # Ge kmer seed for new contig
        kmer_seq, kmer_count = kmer_tracker.get_kmer()
        # Only analyze kmers that exist in 2 or more reads
        if kmer_count < 2: 
            continue

        utils.log(logging_name, 'info', 'Initiating kmer %s, found in %d reads' % (kmer_seq, kmer_count))
        setup_contigs(kmer_seq, fq_recs, kmer_len, kmer_tracker, contig_buffer)

        while len(contig_buffer.contigs) > 0:
            assembly_contig = contig_buffer.get_contig()
            assembly_contig.grow(fq_recs, kmer_tracker, kmer_len, contig_buffer)
            if assembly_contig.get_total_read_support() < int(rc_thresh) or assembly_contig.get_contig_len() <= read_len:
                utils.log(logging_name, 'info', 'Contig did not meet the read count threshold %d, with %d or contig length (%d) < read_len (%d)' % (rc_thresh, len(assembly_contig.reads), len(assembly_contig.seq.value), read_len))
            else: 
                utils.log(logging_name, 'info', 'Adding contig to buffer')
                contigs.append(assembly_contig)
        contig_buffer.remove_kmers(kmer_tracker)
        contig_buffer.remove_reads(fq_recs)
    return contigs


def setup_contigs(kmer_seq, fq_recs, kmer_len, kmer_tracker, contig_buffer):

    '''Create a contig instance starting with a seed kmer and associated reads.

    First find the reads containing the kmerSeq value, iterate through reads and
    either create a new contig or add to existing contig.

    Args:
        kmer_seq (str):               Kmer sequence.
        fq_recs (dict):               Read sequences stored as keys with associated list of fq_read objects.
        kmer_len (int):               Number of base pairs in kmer sequence.
        kmer_tracker (KmerTracker):   Object containing all the kmer values.
        contig_buffer (ContigBuffer): Object to track the buffered contig objects.
    Returns:
        None
    '''

    contig_unit = None

    # Find all reads with kmer sequence passed in.
    # kmer_reads contains a list of tuples.
    #   1. fq_read object defined in breakmer.utils.py
    #   2. Starting position of the kmer match in the read sequence
    #   3. Boolean that a match was found.
    #   4. Length of the read sequence.
    #   5. Number of reads with this sequence.
    kmer_reads = assembly_utils.find_reads(kmer_seq, fq_recs.items(), set())
    contig_buffer.add_used_kmer(kmer_seq)
    kmer_unit = assembly_utils.KmerUnit(kmer_seq, kmer_tracker.get_count(kmer_seq), kmer_tracker.kmer_seqs, kmer_len)
    for read_values in kmer_reads:
        read, kmer_pos, kmer_found, read_len, nreads = read_values
        read_align_values = {'read': read,
                             'align_pos': kmer_pos,
                             'nreads': nreads}
        contig_buffer.add_used_read(read.id)
        if not contig_unit: 
            contig_unit = contig.Contig(kmer_unit, read_align_values)
            contig_buffer.add_contig(read, contig_unit)
        else:
            contig_unit.check_read(kmer_seq, kmer_tracker.get_count(kmer_seq), read, kmer_pos, nreads, kmer_tracker.kmer_seqs, 'setup')
    if contig_unit: 
        contig_unit.finalize(fq_recs, kmer_tracker, contig_buffer, 'setup')

# def same_reads(seq1, seq2):
#   same = False
#   aln = olc.nw(seq1, seq2)
#   if aln[3] == 0 and aln[5] == 0 and aln[6] > 0.95*(len(seq1)):
#     same = True
#   return same

# Check if seq2 is a subseq of seq1 
# def subseq(seq1, seq2): 
#   aln = olc.nw(seq2, seq1)
#   seq2_sub = (False,None)
#   if aln[2] == len(seq2) and aln[3] == 0 and aln[6] >= (0.85*(len(seq2))): 
#     if len(seq2) < len(seq1):
#       seq2_sub = (True,None)
#     else:
#       seq2_sub = (True,aln[6]) 
#   else:
#     seq2_sub = (False,aln[6])
#   return seq2_sub


# def sim_seqs(seq1, b_read): 
#   sim = False
#   if not b_read.redundant: 
#     seq2 = b_read.read.seq
#     if same_reads(seq1, seq2) or subseq(seq2, seq1): 
#       sim = True 
#   return sim


class KmerTracker(object):

    '''Wrapper class for storing the kmer objects. Useful for adding
    and extracting kmers.

    Attributes:
        kmers (list):                  Tuples containing kmer count, kmer sequence, kmer object.
        ordered_kmers (OrderedDict):   Object containing kmer seq as key and kmer count as value.
                                       The top values are the most frequent kmer values.
        kmer_seqs (set):               Kmer sequences that exist in ordered_kmers.
    '''

    def __init__(self):
        self.kmers = []
        self.ordered_kmers = OrderedDict()
        self.kmer_seqs = set()

    def add_kmer(self, kmer_seq, count):

        '''Add a kmer object to the list. Stores a tuple with kmer count and kmer sequence string.
        This allows easy sorting.

        Args:
            kmer_seq (str): Kmer sequence.
            count (int):    Number of reads containing kmer sequence.
        Returns:
            None
        '''

        if len(set(kmer_seq)) > 1:
            self.kmers.append((int(count), kmer_seq))

    def set_all_kmer_values(self):

        '''Sort the kmer list by number of reads (descending) first and then
        by sequence value and store them in an ordered dictionary.

        Args:
            None
        Returns:
            None
        '''

        # Sort the kmer sequences by their frequency in read sequences.
        for kmer_tuple in sorted(self.kmers, key=lambda x: (int(x[0]), x[1]), reverse=True):
            self.ordered_kmers[kmer_tuple[1]] = kmer_tuple[0]

    def has_kmers(self):

        '''Check if there are any kmers left in the dictionary.

        Args:
            None
        Return:
            True if there are items in the dictionary and the counts of those items are > 1.
            False if there are no items in the dictionary or the counts of those items are <= 1.
        '''

        return (len(self.ordered_kmers) > 0) and (max(self.ordered_kmers.values()) > 1)

    def update_kmer_set(self):

        '''Update the set of kmer values. The ordered_kmers dictionary
        dynamically changes as kmers are taken out.

        Args:
            None
        Returns:
            None
        '''

        self.kmer_seqs = set(self.ordered_kmers.keys())

    def get_kmer(self):
        """Return the first kmer in the ordered dictionary"""
        return self.ordered_kmers.items()[0]

    def get_count(self, kmerSeq):
        """Return the number of reads the kmer_seq is within."""
        return self.ordered_kmers[kmerSeq]

    def remove_kmer(self, kmerSeq):
        """Delete the record associated with kmer sequence."""
        del self.ordered_kmers[kmerSeq]


class ContigBuffer(object):

    '''A class to track the used kmers and reads and their relation to contigs.

    Attributes:
        used_kmers (set):       Kmer sequences that have been used to build contigs.
        used_reads (set):       Read IDs that have been used to build contigs.
        contigs (OrderedDict):  Object to track the sequence reads and the contigs in which they are involved.
    '''

    def __init__(self):
        self.used_kmers = set()
        self.used_reads = set()
        self.contigs = OrderedDict()

    def add_contig(self, read, contig):

        '''Add read to contigs dict with contig object it is connected to.
        Set key to read ID and value to the contig object. Set the read used to True.

        Args:
            read:   fq_read object
            contig: Contig object.
        Return:
            None
        '''

        # Tie a contig to the seed read ID and store in dictionary.
        if read.id not in self.contigs and not read.used:
            self.contigs[read.id] = contig
            read.used = True

    def remove_contig(self, read_id):

        """Remove read ID from contigs dictionary.
        Args:
            read_id: String of read ID.
        Return: None
        """

        if read_id in self.contigs:
            del self.contigs[read_id]

    def get_contig(self):

        '''Return the contig associated with the first record the contigs dictionary.
        Delete the entry.

        Args:
            None
        Returns:
            Contig object to grow or complete.
        '''

        read_id = self.contigs.keys()[0]
        contig = self.contigs[read_id]
        del self.contigs[read_id]
        return contig

    def add_used_read(self, read_id):

        '''Add read ID to used set.

        Args:
            read_id (str): read sequence ID.
        Returns:
            None
        '''

        self.used_reads.add(read_id)

    def add_used_kmer(self, kmer_seq):

        '''Add kmer sequence to used set.

        Args:
            kmer_seq (str): Kmer sequence.
        Return: None
        '''

        self.used_kmers.add(kmer_seq)

    def remove_kmers(self, kmer_tracker):

        """Remove used kmer sequences from the kmer tracking object and reset used
        kmer set.
        Args:
            kmer_tracker (KmerTracker): Object to 
        Returns:
            None
        """

        map(kmer_tracker.remove_kmer, list(self.used_kmers))
        self.used_kmers = set()

    def remove_reads(self, fq_reads):
        """Remove the used reads from the fq_reads dictionary.
        Args:
            fqReads: Dictionary of fq_reads.
        Return: None
        """
        del_used = filter(lambda x: x in fq_reads, list(self.used_reads))
        map(fq_reads.__delitem__, del_used)
        self.used_reads = set()
