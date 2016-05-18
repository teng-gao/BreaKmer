#! /usr/bin/local/python

import re
import logging
from collections import OrderedDict
import breakmer.utils as utils
import breakmer.assembly.utils as assembly_utils
import breakmer.assembly.contig as contig

def setup_contigs(mer, fq_recs, kmer_len, akmers, buff):

    '''
    '''

    logger = logging.getLogger('breakmer.assembly')
    ct = None
    mer_reads = assembly_utils.find_reads(mer, fq_recs.items(), set())
    buff.add_used_mer(mer)
    for read_vals in mer_reads:
        read, mer_pos, bool, rlen, nreads = read_vals
        buff.add_used_read(read.id)
        if not ct: 
            ct = contig.AssemblyContig(mer, read, mer_pos, nreads, kmer_len)
            buff.add_contig(read, ct)
        else:
            ct.check_read(mer, akmers.get_count(mer), read, mer_pos, nreads, akmers.smers_set, 'setup')
    if ct: 
#    print 'Setup done'
        ct.finalize(fq_recs, akmers, buff, 'setup')


def init_assembly(mers, fq_recs, kmer_len, rc_thresh, read_len):

    '''Entry function for assemblying a contiguous sequence from
    a pool of sample only kmers and the reads that contain them.

    A kmer tracker object is instantiated containing all the kmer seqs and
    their associated counts. These are sorted by

    Args:
        kmers:      Dictionary of kmers only in the sample key = kmer, value = count in reads
        fqRecs:     Dictionary with sequence values as keys and a list of fq_read objects.
        kmerLen:    Integer of kmer size.
        rcThresh:   Integer representing the minimum readcount threshold for keeping a contig.
        readLen:    Integer of the read length.
    Returns:
        contigs:    List of contig objects.
    '''

    logger = logging.getLogger('breakmer.assembly')
    kmer_clusters = []
    if len(mers) == 0: 
        return kmer_clusters

    ks = kmers()
    for mer in mers: 
        ks.add_kmer(mer, mers[mer])

    buff = contig.buffer()
    akmers = ks.get_all_kmer_values()
    while akmers.has_mers():
        akmers.update_smer_set()
        mer, count = akmers.mers.items()[0]
        if count < 2: 
            continue
        logger.info('Initiating kmer %s, found in %d reads'%(mer,count))
        setup_contigs(mer, fq_recs, kmer_len, akmers, buff)
        while len(buff.contigs) > 0:
            ct = buff.get_contig()
            ct.grow(fq_recs, ks, akmers, kmer_len, buff)
            if ct.get_total_read_support() < int(rc_thresh) or ct.get_contig_len() <= read_len:
                logger.info('Contig did not meet the read count threshold %d, with %d or contig length (%d) < read_len (%d)'%(rc_thresh, len(ct.reads), len(ct.aseq.seq), read_len))
            else: 
                logger.info('Adding contig to buffer')
                kmer_clusters.append(ct)
        buff.remove_kmers(akmers)
        buff.remove_reads(fq_recs)
    return kmer_clusters


def same_reads(seq1, seq2):
  same = False
  aln = olc.nw(seq1, seq2)
  if aln[3] == 0 and aln[5] == 0 and aln[6] > 0.95*(len(seq1)):
    same = True
  return same


# Check if seq2 is a subseq of seq1 
def subseq(seq1, seq2): 
  aln = olc.nw(seq2, seq1)
  seq2_sub = (False,None)
  if aln[2] == len(seq2) and aln[3] == 0 and aln[6] >= (0.85*(len(seq2))): 
    if len(seq2) < len(seq1):
      seq2_sub = (True,None)
    else:
      seq2_sub = (True,aln[6]) 
  else:
    seq2_sub = (False,aln[6])
#  print seq2_sub, aln
  return seq2_sub

 

def sim_seqs(seq1, b_read): 
  sim = False
  if not b_read.redundant: 
    seq2 = b_read.read.seq
    if same_reads(seq1, seq2) or subseq(seq2, seq1): 
      sim = True 
  return sim


class kmer:
  def __init__(self,val,count):
    self.value = val
    self.count = count
    self.done = False

class kmers: 
  def __init__(self):
    self.s = []
  #*******************************************
  def add_kmer(self,mer,count):
    if len(set(mer)) > 1:
      self.s.append((int(count), mer, kmer(mer,count)))
  #*******************************************
  def get_all_kmer_values(self):
    mers_sorted = sorted(self.s, key=lambda x: (int(x[0]), x[1]), reverse=True)
    am = akmers()
    for mer in mers_sorted:
      am.add_mer(mer[1], mer[0])
    return am
  #*******************************************

class akmers:
  def __init__(self):
    self.mers = OrderedDict()
    self.smers_set = set()
  #*******************************************
  def get_mers(self):
    return set(self.mers.keys())
  #*******************************************
  def update_smer_set(self):
    self.smers_set = set(self.mers.keys())
  #*******************************************
  def add_mer(self, mer, count):
    self.mers[mer] = count
   #*******************************************
  def remove_mer(self, mer):
    del self.mers[mer]
  #*******************************************
  def has_mers(self):
    if len(self.mers) > 0 and max(self.mers.values()) > 1:
      return True
    else:
      return False
  #*******************************************
  def get_count(self, mer):
    return self.mers[mer]
  #*******************************************



