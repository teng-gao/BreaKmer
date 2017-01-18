#! /usr/bin/python
# -*- coding: utf-8 -*-


'''
BreaKmer target module
'''


import sys
import os
import subprocess
import shutil
import pysam
import breakmer.utils as utils
import breakmer.assembly.assembler as assembler
from collections import deque
import breakmer.caller.sv_caller2 as sv_caller

__author__ = "Ryan Abo"
__copyright__ = "Copyright 2015, Ryan Abo"
__email__ = "ryanabo@gmail.com"
__license__ = "MIT"


def pe_meta(aread):
    '''
    '''

    # First check if read is from a proper paired-end mapping --> <--
    proper_map = False
    overlap_reads = False
    proper_map1 = ((aread.flag == 83) or (aread.flag == 147)) and (aread.tlen < 0)
    proper_map2 = ((aread.flag == 99) or (aread.flag == 163)) and (aread.tlen > 0)
    if proper_map1 or proper_map2:
        proper_map = True
        if abs(aread.tlen) < (2 * len(aread.seq)):
            overlap_reads = True
    return proper_map, overlap_reads


def add_discordant_pe(aread, read_d, bamfile):

    '''
    '''

    qname = aread.qname
    # Keep discordant read pairs where the map quality is > 0, the paired reads are mapped to different chroms or > 1000 bp apart, and
    # the mate is mapped.
    if aread.mapq > 0 and ((aread.rnext != -1 and aread.tid != aread.rnext) or abs(aread.tlen) > 1000) and not aread.mate_is_unmapped:
        mate_refid = bamfile.getrname(aread.rnext)  # Grab the paired read
        # mate_read = None
        # try:    
        #     mate_read = bamfile.mate(aread)
        # except:
        #     print 'Skipping read'
        #     pass

        # if mate_read is not None:
        #     if mate_read.mapq > 0:
        if mate_refid not in read_d['disc']:
            read_d['disc'][mate_refid] = []
        read_d['disc'][mate_refid].append((aread.pos, aread.pnext))  # Store the read position and the mate position

    if aread.mapq > 0 and not aread.mate_is_unmapped and aread.tid == aread.mrnm:
        if aread.is_read1:
            read_positions = None
            if aread.is_reverse and aread.mate_is_reverse:
                # reverse -- reverse, samflag 115 (note: only considering read1, read2 samflag 179)
                read_positions = (aread.pos, aread.mpos, 0, 0, qname)
                if aread.mpos < aread.pos:
                    read_positions = (aread.mpos, aread.pos, 0, 0, qname)
                read_d['inv_reads'].append(read_positions)
            elif not aread.is_reverse and not aread.mate_is_reverse:
                # forward -- forward = samflag 67 (note: only considering read1, read2 samflag 131)
                read_positions = (aread.pos, aread.mpos, 1, 1, qname)
                if aread.mpos < aread.pos:
                    read_positions = (aread.mpos, aread.pos, 1, 1, qname)
                read_d['inv_reads'].append(read_positions)
            elif aread.is_reverse and not aread.mate_is_reverse and aread.pos < aread.mpos:
                # reverse -- forward = samflag 83 with positive insert (read2 samflag 163 with + insert size)
                read_positions = (aread.pos, aread.mpos, 0, 1, aread.qname)
                read_d['td_reads'].append(read_positions)
            elif not aread.is_reverse and aread.mate_is_reverse and aread.mpos < aread.pos:
                # reverse -- forward = samflag 99 with - insert (read2 samflag 147 with - insert)
                read_positions = (aread.mpos, aread.pos, 1, 0, qname)
                read_d['td_reads'].append(read_positions)
            if read_positions:
                read_d['other'].append(read_positions)

class KmerHash(object):
    '''
    '''

    def __init__(self, chrom, start, end):
        self.kmers = set()
        self.limit = 1000000
        self.chrom = chrom
        self.start = int(start)
        self.end = int(end)

    def add(self, kmer):
        self.kmers.add(kmer)
        if len(self.kmers) > self.limit:
            return True
        return False

    def set_window(self, pos):
        self.end = pos

    def reset_kmers(self):
        self.kmers = set()


class ClusterManager(object):
    def __init__(self, insert_size):
        # Store by partner chromosome and then discordance flag
        self.clusters = {}
        self.clip_positions_dict = None
        self.insert_size = insert_size

    def check_cluster_event(self):
        '''
        33 = + -
        0 = + NA
        1 = + +
        49 = - -
        16 = - NA
        17 = - +
        51 = - - same chrom inv
        3 = + + same chrom inv
        35 = + - same chrom (dist > ins thresh)
        19 = - + same chrom (dist > ins thresh)
        outtie = 
        '''

        for chrom in self.clusters:
            for dr_flag in self.clusters[chrom]:
                for rc in self.clusters[chrom][dr_flag]['closed']:
                    if rc.dr == 1:
                        continue
                    # Translocations
                    if dr_flag in [33, 1, 49, 17]:
                        print dr_flag, 'trl'
                    # Insertions
                    elif dr_flag in [0, 16]:
                        print dr_flag, 'ins'
                    # Inversions
                    elif dr_flag in [51, 3]:
                        print dr_flag, 'inv'
                    # Tandem dups
                    elif dr_flag == "outtie":
                        print dr_flag, 'td'
                    # Deletions
                    elif dr_flag in [35, 19]:
                        sc_dir = '+' if dr_flag == 35 else '-'
                        for cp in self.clip_positions_dict:
                            print 'Checking clip position', cp
                            if sc_dir in self.clip_positions_dict[cp]:
                                dc_cluster_position = rc.last_pos if dr_flag == 35 else rc.start
                                print 'Clip position dist from cluster', dc_cluster_position, rc.last_pos, rc.start, sc_dir
                                if abs(int(cp) - int(dc_cluster_position)) <= self.insert_size:
                                    print 'Cluster has softclip support', cp, dr_flag, sc_dir, [x for x in self.clip_positions_dict[cp][sc_dir]]


    def check_read(self, read, outties=False):
        '''
        '''

        # Set unmapped mates to the same chrom for easy checking later
        next_ref_name = read.reference_name if read.next_reference_id==-1 else read.next_reference_name
        if next_ref_name not in self.clusters:
            self.clusters[next_ref_name] = {}

        dr_flag = 'outtie' if outties else read.flag
        if not outties and read.is_paired:
            dr_flag -= 64 if read.is_read1 else 128

        if dr_flag not in self.clusters[next_ref_name]:
            self.clusters[next_ref_name][dr_flag] = {'active': deque([]), 'closed': []}
        # if read.next_reference_id not in self.clusters[dr_flag]:
        #     self.clusters[dr_flag][read.next_reference_id] = {'active': deque([]), 'closed': []}

        # Check active clusters
        clusters = self.clusters[next_ref_name][dr_flag]['active']
        add_to_cluster = False
        pop_indices = []
        for i,cluster in enumerate(clusters):
            # print '\tChecking window', win.start, win.end, read.reference_start
            if cluster.end < read.reference_start:
                pop_indices.append(i)
            else:
                # Add read to window
                add_to_cluster = True
                # print read
                cluster.add_read(read)
                # print '\tAdded read to cluster', dr_flag, read.next_reference_id, cluster.dr
                break
        if not add_to_cluster:
            self.clusters[next_ref_name][dr_flag]['active'].append(ReadCluster(read))
        for i in pop_indices:
            # Remove cluster from active and add to closed
            add_cluster = clusters.popleft()
            if add_cluster.dr > 1:
                self.clusters[next_ref_name][dr_flag]['closed'].append(add_cluster)
            self.clusters[next_ref_name][dr_flag]['active'] = clusters

    def check_clusters(self, brkpts):
        '''
        '''

        # Check the clusters in the target based on dr_flag for strand and type.
        # Strand + trl dr_flag = 
        # DR_FLAG 
        # 0 = not paired mapped forward
        # 1 = paired, forward, mapped
        # 
        # Look in disc_reads.clusters for chrom2 first, then look for specific DR_FLAG
        # 33 = + -
        # 0 = + NA
        # 1 = + +
        # 49 = - -
        # 16 = - NA
        # 17 = - +
        # 51 = - - same chrom inv
        # 3 = + + same chrom inv
        # 35 = + - same chrom
        # 19 = - + same chrom

        bp = sorted(brkpts, key=lambda brkpt: brkpt[0], reverse=True) 
        # print bp
        ndiscordant_reads = 0
        for i, v1 in enumerate(bp):
            in_target, chrom, pos, strand = v1
            # print in_target, chrom, pos, strand
            if in_target == 0:
                break
            for j, v2 in enumerate(bp[i+1:]):
                in_target2, chrom2, pos2, strand2 = v2
                # print in_target2, chrom2, pos2, strand2
                if chrom2 in self.clusters:
                    dr_flags = self.clusters[chrom2]
                    dr_flag = None
                    if (strand == '+') and (strand2 == '-') and (chrom != chrom2):
                        # Translocation event DR_FLAG 33
                        dr_flag = 33
                    elif (strand == '+') and (strand2 == '+') and (chrom != chrom2):
                        # Translocation event DR_FLAG 33
                        dr_flag = 1
                    elif (strand == '-') and (strand2 == '-') and (chrom != chrom2):
                        dr_flag = 49
                    elif (strand == '-') and (strand2 == '+') and (chrom != chrom2):
                        dr_flag = 17
                    if dr_flag in dr_flags:
                        for readcluster in dr_flags[dr_flag]['closed']:
                            # Strand is '+' so look upstream, last_pos can't be past breakpoint
                            if abs(pos - readcluster.last_pos) < self.insert_size:
                                ndiscordant_reads += readcluster.dr
                                # print readcluster.start, readcluster.last_pos, pos, readcluster.dr
        return ndiscordant_reads

    def check_clusters_td(self, chrom, brkpt1, brkpt2):
        '''

        There are two distinct possibilities where a contig contains a breakpoint from an inversion event. These are enumerated below:

        Sample sequence is observed as this:
                                             1A  1B +
        ----------||---------- <DUP> ----------||---------- <DUP> ----------||----------
        ----------||---------- <DUP> ----------||---------- <DUP> ----------||----------
                                             2B  2A - 

        Realignment to reference genome will be observed as this:
        1.         1A+        1B+  Look for DR_FLAG - + = "outtie" with 1B < DR1, DR2, < 1A 
            ----------(down)||----------

        2.  (rc)   2A-      2B+ Look for DR_FLAG - - = 51 with 2A < DR1 < 2B and 2B < DR2
            ----------(up)||----------
        '''

        ndiscordant_reads = 0
        if brkpt1[1] == '+' and brkpt2[1] == '+':
            # The DR_FLAG should be - - and between brkpt1 and brkpt2 51
            if chrom in self.clusters:
                if "outtie" in self.clusters[chrom]:
                    for readcluster in self.clusters[chrom]["outtie"]['closed']:
                        # Strand is '+' so look upstream, last_pos can't be past breakpoint
                        if brkpt1[0] > brkpt2[0]:
                            # Event #1
                            left_pos = brkpt2[0] - 10
                            right_pos = brkpt1[0] + 10
                            if (left_pos < readcluster.start < right_pos) and (left_pos < readcluster.last_pos < right_pos):
                                ndiscordant_reads += readcluster.dr
        elif brkpt1[1] == '-' and brkpt2[1] == '-':
            # The DR_FLAG should be - - and between brkpt1 and brkpt2 51
            if chrom in self.clusters:
                if "outtie" in self.clusters[chrom]:
                    for readcluster in self.clusters[chrom]["outtie"]['closed']:
                        # Strand is '-' so look upstream, last_pos can't be past breakpoint
                        if brkpt1[0] > brkpt2[0]:
                            # Event #2
                            left_pos = brkpt2[0] - 10
                            right_pos = brkpt1[0] + 10
                            if (left_pos < readcluster.start < right_pos) and (left_pos < readcluster.last_pos < right_pos):
                                ndiscordant_reads += readcluster.dr
        return ndiscordant_reads


    def check_clusters_inv(self, chrom, brkpt1, brkpt2):
        '''

        There are four distinct possibilities where a contig contains a breakpoint from an inversion event. These are enumerated below:

        Sample sequence is observed as this:
                1A  1B +                     2A  2B +
        ----------||---------- <INV> ----------||---------- 
        ----------||---------- <INV> ----------||---------- 
                3B  3A -                     4B  4A - 

        Realignment to reference genome will be observed as this:
        1.         1A+        1B+  Look for DR_FLAG + + = 3 with DR1 < 1A and 1A < DR2 < 1B 
            ----------(down)||----------

        2.  (rc)   2A-      2B+ Look for DR_FLAG - - = 51 with 2A < DR1 < 2B and 2B < DR2
            ----------(up)||----------

        3.         3A+        3B- (reverse complement) Look for DR_FLAG + + = 3 with DR1 < 3B and 3B < DR2 < 3A
            ----------(down)||----------

        4.  (rc)   4A-        4B+ Look for DR_FLAG - - = 51 with 4A > DR1 and 4B < DR2 < 4A
            ----------(down)||----------

        '''

        # for key in self.clusters.keys():
        #     print key, self.clusters[key].keys()
        #     for kk in self.clusters[key].keys():
        #         if kk == 49:
        #             print self.clusters[key][kk]['closed'][0].reads

        ndiscordant_reads = 0
        if brkpt1[1] == '-' and brkpt2[1] == '+':
            # The DR_FLAG should be - - and between brkpt1 and brkpt2 51
            if chrom in self.clusters:
                if 51 in self.clusters[chrom]:
                    for readcluster in self.clusters[chrom][51]['closed']:
                        # Strand is '-' so look upstream, last_pos can't be past breakpoint
                        if brkpt1[0] > brkpt2[0]:
                            # Event #4
                            left_pos = brkpt2[0] - 10
                            right_pos = brkpt1[0] + 10
                        else:
                            # Event #2
                            left_pos = brkpt1[0] - 10
                            right_pos = brkpt2[0] + 10
                        if (left_pos < readcluster.start < right_pos) and (left_pos < readcluster.last_pos < right_pos):
                            ndiscordant_reads += readcluster.dr
        elif brkpt1[1] == '+' and brkpt2[1] == '-':
            # The DR_FLAG should be - - and between brkpt1 and brkpt2 51
            if chrom in self.clusters:
                if 3 in self.clusters[chrom]:
                    for readcluster in self.clusters[chrom][3]['closed']:
                        # Strand is '-' so look upstream, last_pos can't be past breakpoint
                        if brkpt1[0] > brkpt2[0]:
                            # Event #3
                            left_pos = brkpt2[0] - 10
                            right_pos = brkpt1[0] + 10
                        else:
                            # Event #1
                            left_pos = brkpt1[0] - 10
                            right_pos = brkpt2[0] + 10
                        if (left_pos < readcluster.start < right_pos) and (left_pos < readcluster.last_pos < right_pos):
                            ndiscordant_reads += readcluster.dr
        return ndiscordant_reads

class ReadCluster(object):
    def __init__(self, read, insert_size=1000):
        self.reads = [read.query_name]
        self.start = read.reference_start
        self.end = self.start + insert_size + 50
        self.last_pos = read.reference_start + len(read.query_sequence)
        self.dr = 1
        self.sr = 0
        self.n = 0
        self.both = 0

    def print_values(self):
        '''
        '''

        print self.start, self.end, self.dr

    def add_read(self, read, read_type='dr'):
        '''
        '''

        # print read.query_name
        if read.query_name in self.reads:
            return
        self.reads.append(read.query_name)
        if read_type == 'dr':
            self.dr += 1
            self.last_pos = read.reference_start + len(read.query_sequence)
        elif read_type == 'sr':
            self.sr += 1
        elif read_type == 'both':
            self.both += 1

    def close_win(self):
        '''
        '''

        print 'Window closing', self.start, self.end, len(self.reads)


class AssemblyBatch(object):
    '''
    '''

    def __init__(self, name, index, bamfile, data_path):
        self.index = index
        self.fq = None
        self.bam = None
        self.bam_sorted = None
        self.assembled_fq = None
        self.nreads = 0
        self.reads = []
        self.start = None
        self.last_pos = None
        self.contig_fn = None
        # self.batch_mers = set()
        self.setup_files(name, bamfile, data_path)

    def setup_files(self, name, bamfile, data_path):
        '''
        '''

        self.fq = os.path.join(data_path, "%s_sv_reads_%d.fastq" % (name, self.index))
        self.bam = pysam.AlignmentFile(os.path.join(data_path, "%s_sv_reads_%d.bam" % (name, self.index)), 'wb', template=bamfile)
        self.bam_sorted = os.path.join(data_path, "%s_sv_reads_%d.sorted.bam" % (name, self.index))

    def check_read(self, read):
        '''
        '''

        if self.start is None:
            self.start = read.reference_start
            self.add_read(read)
            return True
        elif abs(self.last_pos - read.reference_start) > len(read.query_sequence): # and self.nreads > 10000:
            # Reject read and make it go into a new batch
            # Close the files.
            # print 'Closing batch'
            self.close_batch()
            return False
        else:
            self.add_read(read)
            return True

    def close_batch(self):
        '''
        '''
        # print 'Close batch'
        if self.nreads > 1:
            # Write files
            fq_f = open(self.fq, 'w')
            # print self.reads
            for read in self.reads:
                fq_f.write("@" + read.qname + "\n" + read.seq + "\n+\n" + read.qual + "\n")
                self.bam.write(read)
            # print 'Closing fq', self.fq
            fq_f.close()
            self.bam.close()
    #     utils.log(self.logging_name, 'info', 'Sorting bam file %s to %s' % (self.files['sv_bam'], self.files['sv_bam_sorted']))
            pysam.sort(self.bam.filename, self.bam_sorted.replace('.bam', ''))
    #     utils.log(self.logging_name, 'info', 'Indexing sorted bam file %s' % self.files['sv_bam_sorted'])
            pysam.index(self.bam_sorted)
        else:
            self.bam.close()
    # def print_values(self):
    #     '''
    #     '''
    #     print 'Assembly reads', self.index, self.fq, self.start, self.last_pos, self.nreads

    def add_read(self, read):
        '''
        '''

        # for kmer in utils.sliding_window(read.seq, 15):
        #     self.batch_mers.add(kmer)
        self.nreads += 1
        self.reads.append(read)
        self.last_pos = read.reference_start + len(read.query_sequence)


class Contig(object):
    def __init__(self, file_path, contig_id, seq, nreads, sv_khash, refmers, kmer_size):
        self.contig_id = contig_id
        self.seq = seq
        self.kmers = {}
        self.nreads = nreads
        self.file_path = file_path
        self.contig_fa_fn = None
        self.sv_khash = sv_khash
        self.refmers = refmers
        self.kmer_size = kmer_size
        self.all_kmers = {}
        self.breakpoint_coverages = {}
        self.exact_brkpt_coverages = []
        # self.kmer_coverage_list = [0] * len(self.seq)
        self.setup()

    def setup(self):
        '''
        '''

        # kiter = 0
        for kmer in utils.get_kmer_set(self.seq, self.kmer_size):
            # kiter = self.seq.find(kmer)
            # if kmer not in self.all_kmers:
            #     self.all_kmers[kmer] = []
            rc = utils.reverse_complement(kmer)
            # self.all_kmers[kmer].append((kmer, kmer in self.sv_khash, kmer in self.refmers, rc in self.sv_khash, rc in self.refmers))
            if kmer in self.sv_khash:
                # print 'kiter', kiter, kmer, self.sv_khash[kmer]
                self.kmers[kmer] = self.sv_khash[kmer]
                # self.kmer_coverage_list[kiter:(kiter+len(kmer))] = map(lambda x: x + self.sv_khash[kmer], self.kmer_coverage_list[kiter:(kiter+len(kmer))])
            elif rc in self.sv_khash:
                # print 'kiter', kiter, kmer, self.sv_khash[rc]
                self.kmers[kmer] = self.sv_khash[rc]
                # self.kmer_coverage_list[kiter:(kiter+len(kmer))] = map(lambda x: x + self.sv_khash[rc], self.kmer_coverage_list[kiter:(kiter+len(kmer))])
            # kiter += 1
        # print 'CONTIG SETUP', self.seq, self.kmer_coverage_list

        # Trim ends of contig
        # liter = 0
        # riter = 0
        # while self.kmer_coverage_list[liter] < 1:
        #     liter += 1
        # while self.kmer_coverage_list[::-1][riter] < 1:
        #     riter += 1
        # print 'Old contig seq', self.seq, len(self.seq)
        # right_trim_index = len(self.seq) - riter
        # self.seq = self.seq[liter:right_trim_index]
        # print 'Trimmed contig', self.seq, liter, riter, right_trim_index
        # print 'New contig seq', self.seq
        # self.kmers [k for k in  if k in sv_khash or utils.reverse_complement(k) in sv_khash]

        self.create_fa()

    def create_fa(self, ):
        '''
        '''
        # print 'Contig kmers', self.seq
        # for kmer in self.kmers:
        #     print kmer, self.sv_khash[kmer]
        self.contig_fa_fn = os.path.join(self.file_path, self.contig_id + ".fa")

        half2_idx = ((len(self.seq)/2) - 20)
        seq_half1 = self.seq[0:((len(self.seq)/2) + 20)]
        seq_half2 = self.seq[((len(self.seq)/2) - 20): len(self.seq)]
        fa_f = open(self.contig_fa_fn, 'w')
        fa_f.write(">" + self.contig_id + "\n" + self.seq + "\n")
        fa_f.write(">" + self.contig_id + ":1_" + str(half2_idx)  + "\n" + seq_half1 + "\n")
        fa_f.write(">" + self.contig_id + ":2_" + str(half2_idx) + "\n" + seq_half2)
        fa_f.close()

    def get_brkpt_coverage(self, contig_pos, store=False):
        '''
        '''

        # print 'Calling coverage', contig_pos, store, self.breakpoint_coverages
        if not store and contig_pos in self.breakpoint_coverages:
            # print 'Get brkpt coverage stored', contig_pos, self.breakpoint_coverages[contig_pos]
            return self.breakpoint_coverages[contig_pos]

        start = max(0, contig_pos - self.kmer_size - 1)
        end = min(len(self.seq), contig_pos + (self.kmer_size - 1))
        # print 'contig seq chunk', self.seq[start:end]
        bp_kmers = utils.get_kmer_set(self.seq[start:end], self.kmer_size)
        # print 'Kmers', [kmer for kmer in kmers]
        # print 'Finding breakpoint coverage in contig', self.seq, contig_pos, bp_kmers
        # for kmer in bp_kmers:
        #     print kmer, kmer in self.refmers, kmer in self.kmers, kmer in self.sv_khash, self.all_kmers[kmer]
        bp_cov = [self.kmers[kmer] for kmer in bp_kmers if kmer in self.kmers]
        cov = 0
        if len(bp_cov) != 0:
            cov = max([self.kmers[kmer] for kmer in bp_kmers if kmer in self.kmers])
        
        self.breakpoint_coverages[contig_pos] = cov
        # for read in reads:
        #     readmers = utils.get_kmer_set(read.seq, 15)
        #     jaccard_bp_kmers = float(len(set(bp_kmers).intersection(set(readmers)))) / float(len(set(bp_kmers).union(set(readmers))))
        #     jaccard_contig = float(len(set(self.kmers).intersection(set(readmers)))) / float(len(set(self.kmers).union(set(readmers))))
        #     if jaccard_bp_kmers > 0 and jaccard_contig > 0:
        #         print jaccard_bp_kmers, jaccard_contig
        #     if any(k in read.seq for k in bp_kmers):
        #         print 'Found read', read.seq
        #         cov += 1
        if store:
            self.exact_brkpt_coverages.append(cov)

        # print 'New brkpt coverage', bp_kmers, contig_pos, bp_cov
        # for bp_kmer in bp_kmers:
        #     print bp_kmer, bp_kmer in self.kmers
        # for kmer in self.kmers:
        #     print kmer, utils.reverse_complement(kmer), self.kmers[kmer]
        return cov

    def valid(self):
        '''
        '''

        # for k in self.sv_khash:
        #     print k, self.sv_khash[k]
        # print self.seq, self.kmers, [(kmer, self.kmers[kmer]) for kmer in self.kmers]
        if len(self.kmers) < 5 or max([self.kmers[kmer] for kmer in self.kmers]) < 2:
            return False
        else:
            return True


class TargetManager(object):

    '''TargetManager class handles all the high level information relating to a target.
    The analysis is peformed at the target level, so this class contains all the information
    necessary to perform an independent analysis.

    Attributes:
        params (ParamManager):      Parameters for breakmer analysis.
        logging_name (str):         Module name for logging file purposes.
        name (str):                 Target name specified in the input bed file.
        chrom (str):                Chromosome ID as specified in the input bed file.
        start (int):                Genomic position for the target region (minimum value among all intervals).
        end (int):                  Genomic position for the target region (maximum value among all intervals).
        paths (dict):               Contains the analysis paths for this target.
        files (dict):               Dicionary containing paths to file names needed for analysis.
        read_len (int):             Length of a single read.
        variation (Variation):      Stores data for variants identified within the target.
        regionBuffer (int):         Base pairs to add or subtract from the target region end and start locations.

    '''

    def __init__(self, intervals, params):
        '''
        '''

        self.params = params
        self.name = None
        self.chrom = None
        self.start = None
        self.end = None
        self.paths = {}
        self.files = {}
        self.disc_reads = None
        self.sv_reads = []
        self.cleaned_read_recs = None
        self.kmer_clusters = []
        self.kmers = {}
        self.contigs = []
        self.results = []
        self.formatted_results = []
        self.svs = {'trl':[0, '-'], 'indel':[0, ''], 'rearrangement':[0, '']}
        self.target_intervals = intervals
        self.repeat_mask = None
        self.logging_name = 'breakmer.processor.target'
        self.call_manager = sv_caller.SVCallManager(params)
        self.normal_hash = None
        self.dr_clusters = ClusterManager(params.get_param('insertsize_thresh'))
        self.assembly_batches = []
        self.sv_khash = {}
        self.ref_khash = {}
        self.setup()

    def setup(self):
        '''Setup the target object with the input params.

        Define the location (chrom, start, end), file paths, directory paths, and name.

        Args:
            None
        Returns:
            None
        '''

        for value in self.target_intervals:
            if not self.name:
                self.name = value[3]
            if not self.chrom:
                self.chrom = value[0]
            if not self.start:
                self.start = int(value[1])
            if not self.end:
                self.end = int(value[2])
            if int(value[1]) < self.start:
                self.start = int(value[1])
            if int(value[2]) > self.end:
                self.end = int(value[2])

        '''
        Create the proper paths for the target analysis.

        Each target analyzed has a set of directories associated with it.

        targets/
            <target name>/
                data/
                contigs/
                kmers/
        There is separate directory for each target in the output directory.

        output/
            <target name>/
        '''

        self.add_path('base', os.path.join(self.params.paths['targets'], self.name))
        self.add_path('ref_data', os.path.join(self.params.paths['ref_data'], self.name))
        self.add_path('data', os.path.join(self.paths['base'], 'data'))
        self.add_path('contigs', os.path.join(self.paths['base'], 'contigs'))
        self.add_path('kmers', os.path.join(self.paths['base'], 'kmers'))
        self.add_path('output', os.path.join(self.params.paths['output'], self.name))

        # Set reference paths
        if 'keep_repeat_regions' in self.params.opts:
            if not self.params.opts['keep_repeat_regions']:
                if 'repeat_mask_file' not in self.params.opts:
                    utils.log(self.logging_name, 'error', 'Keep repeat regions option is false, but no repeat mask bed file provided. All repeat region variants will be reported.')
                    self.params.opts['keep_repeat_regions'] = True
                else:
                    self.files['rep_mask_fn'] = os.path.join(self.paths['ref_data'], self.name+'_rep_mask.bed')

        '''
        Each target has reference files associated with it.

        <ref_data_dir>/
            <target_name>/
                <target_name>_forward_refseq.fa
                <target_name>_reverse_refseq.fa
                <target_name>_forward_refseq.fa_dump
                <target_name>_reverse_refseq.fa_dump
        '''
        self.files['target_ref_fn'] = [os.path.join(self.paths['ref_data'], self.name + '_forward_refseq.fa'), os.path.join(self.paths['ref_data'], self.name + '_reverse_refseq.fa')]

        ref_fa_marker_f = open(os.path.join(self.paths['ref_data'], '.reference_fasta'), 'w')
        ref_fa_marker_f.write(self.params.opts['reference_fasta'])
        ref_fa_marker_f.close()

        self.files['ref_kmer_dump_fn'] = [os.path.join(self.paths['ref_data'], self.name + '_forward_refseq.fa_dump'), os.path.join(self.paths['ref_data'], self.name + '_reverse_refseq.fa_dump')]

    def hash_normal_reads(self):
        '''
        '''

        bamfile = pysam.AlignmentFile(self.params.opts['normal_bam_file'], 'rb')
        buffer_size = int(self.params.get_param('buffer_size'))
        kmer_size = int(self.params.get_param('kmer_size'))

        aligned_reads = None
        if self.normal_hash is None:
            self.normal_hash = KmerHash(self.chrom, self.start - buffer_size, self.end + buffer_size)
            utils.log(self.logging_name, 'debug', 'Fetching bam file reads from %s, %s %d %d' % (self.params.opts['normal_bam_file'], self.chrom, self.start - buffer_size, self.end + buffer_size))
            aligned_reads = bamfile.fetch(self.chrom, self.start - buffer_size, self.end + buffer_size)
        else:
            utils.log(self.logging_name, 'debug', 'Fetching bam file reads from %s, %s %d %d' % (self.params.opts['normal_bam_file'], self.chrom, self.start - buffer_size, self.end + buffer_size))
            aligned_reads = bamfile.fetch(self.chrom, self.normal_hash.end, self.end + buffer_size)
            self.normal_hash.start = self.normal_hash.end
            self.normal_hash.end = self.end
            self.normal_hash.reset_kmers()

        hit_limit = False
        for read in aligned_reads:
            if not utils.check_read_mismapping(read):
                continue
            for kmer in utils.sliding_window(read.query_sequence, kmer_size):
                hit_limit = self.normal_hash.add(kmer)
            if hit_limit:
                print 'Hit limit'
                self.normal_hash.set_window(read.reference_start)
                break
        print 'Normal kmers', len(self.normal_hash.kmers)

    def get_sv_reads(self):
        '''
        '''

        if 'normal_bam_file' in self.params.opts:
            print 'Hashing normal read kmers'
            self.hash_normal_reads()

            # self.extract_bam_reads('norm')
            # self.clean_reads('norm')
        self.extract_bam_reads('sv')
        # print self.dr_clusters.clusters

        # check = True
        # if not self.clean_reads('sv'):
        #     shutil.rmtree(self.paths['output'])
        #     check = False
        return True

    def setup_read_extraction_files(self, sample_type):
        '''
        '''

        self.files['%s_fq' % sample_type] = os.path.join(self.paths['data'], self.name + "_sv_reads.fastq")
        self.files['%s_sc_unmapped_fa' % sample_type] = os.path.join(self.paths['data'], self.name + "_sv_sc_seqs.fa")
        if sample_type == 'sv':
          self.files['sv_bam'] = os.path.join(self.paths['data'], self.name + "_sv_reads.bam")
          self.files['sv_bam_sorted'] = os.path.join(self.paths['data'], self.name + "_sv_reads.sorted.bam")

    # def new_extract_reads(self, reads):
    #     '''
    #     '''

    #     self.counter += 1
    #     for read in reads:

    #         if read.is_duplicate or read.is_secondary or read.is_qcfail:
    #             return

    #         edits = False
    #         clips = False
    #         all_matches = True
    #         if read.has_tag("MD"):
    #             all_matches = read.get_tag('MD') == str(len(read.query_sequence))
    #         # elif read.has_tag("nM"):
    #         #     all_matches = read.get_tag("nM") == 0
    #         # elif read.has_tag('NM'):
    #         #     all_matches = read.get_tag("NM") < 2
    #         #     edits = read.get_tag("NM") >= 2

    #         perfect_align = True if (len(read.cigartuples) == 1) and (all_matches) and (read.flag in [83, 163, 99, 147]) else False
    #         if perfect_align:
    #             continue
    #         elif read.mapq == 0:
    #             continue

    #         clip_indices = []
    #         if read.cigartuples is not None:
    #             for index, cigar in enumerate(read.cigartuples):
    #                 if cigar[0] == 4:
    #                     clip_indices.append((index, cigar[1]))
    #                     sindex = index if index == 0 else len(read.query_sequence) - cigar[1]
    #                     eindex = cigar[1] + 1 if index == 0 else len(read.query_sequence)
    #                     seq = read.query_sequence[sindex:eindex]
    #                     qual = read.query_qualities[sindex:eindex]
    #                     if max(qual) < 3:
    #                         continue
    #                     clips = True

    #         if clips:
    #             add_to_win = False
    #             pop_indices = []
    #             for i,win in enumerate(self.clip_windows):
    #                 if win.end < read.reference_start:
    #                     win.close_window(self.inbam)
    #                     pop_indices.append(i)
    #                 else:
    #                     add_to_win = True
    #                     win.add_read(read)
    #                     break
    #             if not add_to_win:
    #                 self.clip_windows.append(clip_window(read))
    #             for i in pop_indices:
    #                 self.clip_windows.popleft()

    def close_dr_clusters(self, clip_positions_dict):
        '''
        '''

        rm_flags = []
        for dr_flag in self.dr_clusters.clusters.keys():
            keep_flag = False
            rm_chroms = []
            for chrom in self.dr_clusters.clusters[dr_flag].keys():
                keep_chrom = False
                for i,rc in enumerate(self.dr_clusters.clusters[dr_flag][chrom]['active']):
                    if rc.dr > 1: 
                        self.dr_clusters.clusters[dr_flag][chrom]['closed'].append(rc)
                del self.dr_clusters.clusters[dr_flag][chrom]['active']
                # print dr_flag, chrom, 'C:', len(self.dr_clusters.clusters[dr_flag][chrom]['closed'])
                if len(self.dr_clusters.clusters[dr_flag][chrom]['closed']) > 0:
                    # print 'Closed clusters', len(self.dr_clusters.clusters[dr_flag][chrom]['closed'])
                    # for i,rc in enumerate(self.dr_clusters.clusters[dr_flag][chrom]['closed']):
                    #     print '\t', i, 'Closed', dr_flag, chrom, rc.dr
                    keep_chrom = True
                if not keep_chrom:
                    # print 'Removing chrom', dr_flag, chrom
                    rm_chroms.append(chrom)
            for rm_chrom in rm_chroms:
                del self.dr_clusters.clusters[dr_flag][rm_chrom]
            # print 'Remove chroms', rm_chroms
            if len(rm_chroms) == self.dr_clusters.clusters.keys():
                rm_flags.append(dr_flag)

        # print 'Remove flags', rm_flags
        for rm_flag in rm_flags:
            del self.dr_clusters.clusters[rm_flag]

        pos_keys = clip_positions_dict.keys()
        for pos in pos_keys:
            clip_dirs = clip_positions_dict[pos].keys()
            for clip_dir in clip_dirs:
                if len(clip_positions_dict[pos][clip_dir]) < 2:
                    del clip_positions_dict[pos][clip_dir]
            if len(clip_positions_dict[pos]) == 0:
                del clip_positions_dict[pos]
        self.dr_clusters.clip_positions_dict = clip_positions_dict

    def extract_bam_reads(self, sample_type):
        '''
        '''
        self.setup_read_extraction_files(sample_type)
        bam_type = 'sample'
        utils.log(self.logging_name, 'info', 'Extracting bam reads from %s to %s' % (self.params.opts['%s_bam_file' % bam_type], self.files['sv_fq']))
        bamfile = pysam.AlignmentFile(self.params.opts['%s_bam_file' % bam_type], 'rb')

        buffer_size = int(self.params.get_param('buffer_size'))
        kmer_size = int(self.params.get_param('kmer_size'))

        utils.log(self.logging_name, 'debug', 'Fetching bam file reads from %s, %s %d %d' % (self.params.opts['%s_bam_file' % bam_type], self.chrom, self.start - buffer_size, self.end + buffer_size))
        aligned_reads = bamfile.fetch(self.chrom, self.start - buffer_size, self.end + buffer_size)

        sv_fq = open(self.files['sv_fq'], 'w')
        batch_index = 1
        
        clip_pos = {}
        self.assembly_batches.append(AssemblyBatch(self.name, batch_index, bamfile, self.paths['data']))
        for aligned_read in aligned_reads:
            sr_read = False
            if aligned_read.is_duplicate or aligned_read.is_secondary or aligned_read.is_qcfail:
                continue

            clips = False
            # all_matches = True
            # if aligned_read.has_tag("MD"):
            #     all_matches = aligned_read.get_tag('MD') == str(len(aligned_read.query_sequence))
            # elif aligned_read.has_tag("nM"):
            #     all_matches = aligned_read.get_tag("nM") == 0
            # elif aligned_read.has_tag('NM'):
            #     all_matches = aligned_read.get_tag("NM") < 2
            #     edits = aligned_read.get_tag("NM") >= 2

            perfect_align = True if utils.skip_read(aligned_read, self.params.get_param("insertsize_thresh")) else False
            if perfect_align:
                ref_kmer = aligned_read.query_sequence[0:kmer_size]
                if ref_kmer not in self.ref_khash:
                    self.ref_khash[ref_kmer] = 0
                self.ref_khash[ref_kmer] += 1
                # for kmer in utils.sliding_window(aligned_read.query_sequence, kmer_size):
                #     if kmer not in ref_khash:
                #         ref_khash[kmer] = 0
                #     ref_khash[kmer] += 1
                continue
            elif aligned_read.mapq == 0:
                continue

            '''
        align_pos_index = read.reference_start
        align_positions = [align_pos_index]
        if read.cigartuples is not None:
            for index, cigar in enumerate(read.cigartuples):
                if index == 0:
                    # Check if first cigar is a softclip, if so adjust read start and 
                    if cigar[0] == 4:
                        align_pos_index = align_pos_index - cigar[1]
                        align_positions[0] = align_pos_index
                        regions += self.get_regions(align_pos_index, kmer_size)
                    else:
                        regions += self.get_regions(align_pos_index, kmer_size)
                if cigar[0] == 0:  # Match
                    align_pos_index += cigar[1]
                elif cigar[0] == 3:  # Gap
                    align_pos_index += cigar[1]
                    regions += self.get_regions(align_pos_index, kmer_size)
                    align_positions.append(align_pos_index)
            '''

            # clip_indices = []
            read_kmers = set()
            alignment_pos = aligned_read.reference_start
            if aligned_read.cigartuples is not None:
                for index, cigar in enumerate(aligned_read.cigartuples):
                    if cigar[0] in [0, 3, 2]:
                        alignment_pos += cigar[1]
                        continue
                    elif cigar[0] != 4:
                        continue
                    # clip_indices.append((index, cigar[1]))
                    sindex = index if index == 0 else len(aligned_read.query_sequence) - cigar[1]
                    eindex = cigar[1] if index == 0 else len(aligned_read.query_sequence)
                    seq = aligned_read.query_sequence[sindex:eindex]
                    qual = aligned_read.query_qualities[sindex:eindex]

                    if max(qual) < 3:
                        continue

                    '''
                    Typical clipped sequence at the beginning of a read - trim the first N bases with low qual - keep new start index
                    seq  ACTGCGTGCGTGC 
                    qual 2222222225566


                    '''
                    idx = 0
                    if index == 0:
                        while qual[idx] < 3 and idx < (len(qual) - 1):
                            idx += 1
                            sindex += 1
                    else:
                        while qual[::-1][idx] < 3 and idx < (len(qual) - 1):
                            # print idx, qual[::-1][idx]
                            idx += 1
                            eindex -= 1

                    # print '1. Clip sequence indices', sindex, eindex, '\n'

                    if alignment_pos not in clip_pos:
                        clip_pos[alignment_pos] = {}
                    clip_dir = '-' if (not aligned_read.is_reverse and index == 0) or (aligned_read.is_reverse and index == 0) else '+'
                    if clip_dir not in clip_pos[alignment_pos]:
                        clip_pos[alignment_pos][clip_dir] = []
                    clip_type = 'e' if (not aligned_read.is_reverse and index == 0) or (aligned_read.is_reverse and index > 0) else 'i'
                    clip_pos[alignment_pos][clip_dir].append((clip_type, aligned_read.query_name))

                    # print '2. Clip qual, post', cigar, seq, len(seq), qual, sr_read, '\n'
                    sindex = sindex if index == 0 else max(0, sindex - kmer_size)
                    eindex = min(len(aligned_read.query_sequence), eindex + 1 + kmer_size) if index == 0 else eindex
                    seq = aligned_read.query_sequence[sindex:eindex]
                    qual = aligned_read.query_qualities[sindex:eindex]

                    # print '3. Clipped seq', sindex, eindex, seq, len(seq)
                    # xx = utils.reverse_complement('AGGAGAGTAATGGGAGTTCTGCAACACATAAG')
                    # if aligned_read.seq.find("AGGAGAGTAATGGGAGTTCTGCAACACATAAG") > -1 or aligned_read.seq.find(xx) > -1:
                    #     print 'FOUND READ BLAH', max(qual), seq, sindex, eindex, seq.find('AGGAGAGTAATGGGAGTTCTGCAACACATAAG'), seq.find(xx)
                    #     print aligned_read

                    for kmer in utils.get_kmer_set(str(seq), kmer_size):
                        # print 'SR Read kmer', kmer
                        # if kmer == 'TGCTAGTGGGAATGTAAAATGGTGCAGCCACT' or utils.reverse_complement('TGCTAGTGGGAATGTAAAATGGTGCAGCCACT') == kmer:
                        #      print '\tFOUND', kmer, kmer in self.ref_khash, self.sv_khash.get(kmer)
                        #      print aligned_read
                        #      sys.exit()
                        if kmer not in self.ref_khash:
                            if kmer not in self.sv_khash:
                                self.sv_khash[kmer] = 0
                            self.sv_khash[kmer] += 1
                            read_kmers.add(kmer)
                            # print kmer, sv_khash[kmer]
                        elif kmer in self.sv_khash:
                            # print 'Remove kmer', kmer
                            del self.sv_khash[kmer]
                    # clips = True

                    # print '4. Read added', aligned_read, seq
                    sr_read = True
                # elif cigar[0] in [1,2]:
                #     sr_read = True

            # print 'Checking if aligned read is soft clipped', aligned_read.query_name, sr_read
            if sr_read:
                # print 'Adding read', sr_read, aligned_read
                if self.normal_hash is not None:
                    if aligned_read.reference_start > self.normal_hash.end:
                        self.hash_normal_reads()
                    # read_kmers = set()
                    # for kmer in utils.sliding_window(aligned_read.query_sequence, kmer_size):
                    #     read_kmers.add(kmer)
                    if len(read_kmers.difference(self.normal_hash.kmers)) == 0:
                    # if len(read_kmers) == 0:
                        continue

                self.sv_reads.append(aligned_read)
                if not self.assembly_batches[-1].check_read(aligned_read):
                    # print 'Closing assembly batch', batch_index
                    batch_index += 1
                    self.assembly_batches.append(AssemblyBatch(self.name, batch_index, bamfile, self.paths['data']))
                    self.assembly_batches[-1].add_read(aligned_read)
                sv_fq.write("@" + aligned_read.qname + "\n" + aligned_read.seq + "\n+\n" + aligned_read.qual + "\n")

            dr_read = False
            if aligned_read.flag not in [83, 163, 99, 147]:
                dr_read = True
                self.dr_clusters.check_read(aligned_read)
            elif utils.check_outties(aligned_read):
                dr_read = True
                self.dr_clusters.check_read(aligned_read, True)
            elif abs(aligned_read.template_length) > self.params.get_param("insertsize_thresh"):
                dr_read = True
                # print 'INSERT SIZE READ', aligned_read.query_name, aligned_read.template_length, aligned_read.reference_start
                self.dr_clusters.check_read(aligned_read)

            # if dr_read:
            #     try:
            #         mate = bamfile.mate(aligned_read)
            #         print 'Getting mate read', aligned_read, mate
            #     except ValueError:
            #         continue
        # for c in clip_pos:
        #     for cc in clip_pos[c]:
        #         print c, cc, clip_pos[c][cc]
        # Close out dr_clusters
        self.assembly_batches[-1].close_batch()
        self.close_dr_clusters(clip_pos)

        # for dr_flag in self.dr_clusters.clusters.keys():
        #     for chrom in self.dr_clusters.clusters[dr_flag].keys():
        #         for i,rc in enumerate(self.dr_clusters.clusters[dr_flag][chrom]['closed']):
        #             print dr_flag, chrom, rc.dr, rc.start, rc.end

        sv_fq.close()

        # if 'TGCTAGTGGGAATGTAAAATGGTGCAGCCACT' in self.sv_khash:
        #     print 'BLAHB LADADFADF'
        #     sys.exit()
        # print len(sv_khash)
        # print sv_khash
        # sys.exit()

        # for aligned_read in aligned_reads:

        #     if aligned_read.is_duplicate or aligned_read.is_qcfail:  # Skip duplicates and failures
        #         continue
        #     if aligned_read.is_unmapped:  # Store unmapped reads
        #         read_d['unmapped'][aligned_read.qname] = aligned_read
        #         continue

        #     if aligned_read.mate_is_unmapped or aligned_read.rnext == -1:  # Indicate that mate is unmapped
        #         aligned_read.mate_is_unmapped = True

        #     proper_map = False
        #     overlap_reads = False

        #     # These two functions can operate on the first read of the pair.
        #     # Check if fragment hasn't been checked yet and that the mate is mapped.
        #     if aligned_read.qname not in pair_indices and not aligned_read.mate_is_unmapped:
        #         add_discordant_pe(aligned_read, read_d, bamfile)
        #         proper_map, overlap_reads = pe_meta(aligned_read)
        #     valid_reads.append((aligned_read, proper_map, overlap_reads))

        #     if aligned_read.qname not in pair_indices and not aligned_read.mate_is_unmapped:
        #         pair_indices[aligned_read.qname] = {}
        #     if aligned_read.qname in pair_indices:
        #         pair_indices[aligned_read.qname][int(aligned_read.is_read1)] = len(valid_reads) - 1

        #     # If read is mapped and mate is unmapped
        #     if (aligned_read.pos >= self.start and aligned_read.pos <= self.end) and aligned_read.mapq > 0 and aligned_read.mate_is_unmapped:
        #         read_d['unmapped_keep'].append(aligned_read.qname)
        # # pair_indices, valid_reads = process_reads(areads, read_d, bamfile)  # Deprecated

        # # for aread, proper_map, overlap_reads in valid_reads:  # Deprecated
        #     # Only take soft-clips from outer regions of properly mapped reads, take all others
        #     if (aligned_read.cigar is None) or (len(aligned_read.cigar) <= 1):  # cigar is a list of tuples 
        #         continue

        # # if aligned_read.cigar and len(aligned_read.cigar) > 1:  
        #     trim_coords = utils.trim_coords(aligned_read.qual, 3)  # Identify the read positions with qual > 2
        #     clip_coords = utils.get_clip_coords(aligned_read.qual, aligned_read.cigar)

        #     # Only keep reads that have a soft clip in sequence that has not been trimmed 
        #     # due to low quality sequence.           
        #     # if clip_coords[0] > trim_coords[0] or clip_coords[1] < trim_coords[1]:  # Deprecated
        #     if clip_coords[0] <= trim_coords[0] and clip_coords[1] >= trim_coords[1]:
        #         continue

        #     sc_seq = {'clipped':[], 'buffered':[]}
        #     new_clip_coords = [0, 0]
        #     start_coord, end_coord = clip_coords
        #     add_sc = [False, False]
        #     indel_only = False
        #     start_sc = start_coord > 0
        #     end_sc = end_coord < len(aligned_read.qual)
        #     seq = aligned_read.seq

        #     if start_sc and end_sc:
        #         add_sc = [True, True]
        #     else:
        #         if start_sc: 
        #             add_sc[0] = True
        #             new_clip_coords = [0, start_coord]
        #             if overlap_reads and aligned_read.is_reverse: 
        #                 mate_seq = valid_reads[pair_indices[aligned_read.qname][int(aligned_read.is_read1)]][0].seq
        #                 add_sc[0] = self.check_pair_overlap(mate_seq, aligned_read, [0, start_coord], 'back')
        #             if proper_map:
        #                 indel_only = aligned_read.is_reverse
        #         elif end_sc: 
        #             new_clip_coords = [end_coord, len(seq)]
        #             add_sc[1] = True
        #             if overlap_reads and not aligned_read.is_reverse: 
        #                 mate_seq = valid_reads[pair_indices[aligned_read.qname][int(aligned_read.is_read1)]][0].seq
        #                 add_sc[1] = self.check_pair_overlap(mate_seq, aligned_read, [end_coord, len(seq)], 'front')
        #             if proper_map:
        #                 indel_only = (indel_only and False) if aligned_read.is_reverse else (indel_only and True)
        #     final_add = add_sc[0] or add_sc[1]
        #     if add_sc[0]:
        #         sc_seq['buffered'].append(aligned_read.seq[0:(start_coord + kmer_size)])
        #         sc_seq['clipped'].append(aligned_read.seq[0:start_coord])
        #     if add_sc[1]:
        #         sc_seq['buffered'].append(seq[(end_coord - kmer_size):len(seq)])
        #         sc_seq['clipped'].append(seq[end_coord:len(seq)])
        #     if final_add:
        #         read_d['sv'][utils.get_seq_readname(aligned_read)] = (aligned_read, sc_seq, new_clip_coords, indel_only)
        # # end for loop

        # sv_fq = open(self.files['sv_fq'], 'w')
        # sv_sc_fa = open(self.files['sv_sc_unmapped_fa'], 'w')

        # for qname in read_d['unmapped_keep']:
        #     if qname in read_d['unmapped']:
        #         read = read_d['unmapped'][qname]
        #         read_d['sv'][utils.get_seq_readname(read)] = (read, None, None, False)
        #         sv_sc_fa.write(">" + read.qname + "\n" + str(read.seq) + "\n")

        # if not self.sv_reads:
        #     self.sv_reads = {}
        # self.sv_reads[sample_type] = {}
        # for qname in read_d['sv']:
        #     aligned_read, sc_seq, clip_coords, indel_only = read_d['sv'][qname]
        #     self.sv_reads[sample_type][qname] = read_d['sv'][qname]
        #     if sample_type == 'sv': 
        #         sv_bam.write(aligned_read)
        #     lout = utils.fq_line(aligned_read, indel_only, int(self.params.get_param('kmer_size')), True)
        #     if lout is not None:
        #         sv_fq.write(lout)
        #     if sc_seq:
        #         for clip_seq in sc_seq['buffered']: 
        #             sv_sc_fa.write(">" + qname + "\n" + clip_seq + "\n")
        # self.disc_reads = {'disc':read_d['disc'], 'inv':read_d['inv_reads'], 'td':read_d['td_reads'], 'other':read_d['other']}
        # sv_fq.close()
        # sv_sc_fa.close()
        bamfile.close()

        # if sample_type == 'sv':
        #     sv_bam.close()
        #     utils.log(self.logging_name, 'info', 'Sorting bam file %s to %s' % (self.files['sv_bam'], self.files['sv_bam_sorted']))
        #     pysam.sort(self.files['sv_bam'], self.files['sv_bam_sorted'].replace('.bam', ''))
        #     utils.log(self.logging_name, 'info', 'Indexing sorted bam file %s' % self.files['sv_bam_sorted'])
        #     pysam.index(self.files['sv_bam_sorted'])

    def clean_reads(self, sample_type):
        '''
        '''

        # Run cleaning program
        cutadapt = self.params.get_param('cutadapt')
        cutadapt_config = self.params.get_param('cutadapt_config_file')
        utils.log(self.logging_name, 'info', 'Cleaning reads using %s with configuration file %s' % (cutadapt, cutadapt_config))

        self.files['%s_cleaned_fq' % sample_type] = os.path.join(self.paths['data'], self.name + "_%s_reads_cleaned.fastq" % sample_type)

        utils.log(self.logging_name, 'info', 'Writing clean reads to %s' % self.files['%s_cleaned_fq' % sample_type])
        cutadapt_parameters = utils.stringify(cutadapt_config)
        cutadapt_cmd = '%s %s %s %s > %s' % (sys.executable, cutadapt, cutadapt_parameters, self.files['%s_fq' % sample_type], self.files['%s_cleaned_fq' % sample_type])
        utils.log(self.logging_name, 'debug', 'Cutadapt system command %s' % cutadapt_cmd)
        cutadapt_proc = subprocess.Popen(cutadapt_cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True)
        output, errors = cutadapt_proc.communicate()
        utils.log(self.logging_name, 'debug', 'Clean reads output %s' % output)
        utils.log(self.logging_name, 'debug', 'Clean reads errors %s' % errors)

        # Use these for pulling out reads after finding sample-only kmers.
        # Filter the cleaned reads to make sure soft clips were not adapters, re-write fastq
        # if not self.cleaned_read_recs:
        #     self.cleaned_read_recs = {}
        # self.cleaned_read_recs[sample_type] = None
        # self.files['%s_cleaned_fq' % sample_type], self.cleaned_read_recs[sample_type] = utils.get_fastq_reads(self.files['%s_cleaned_fq' % sample_type], self.sv_reads[sample_type])
        # self.sv_reads[sample_type] = None
        # check = True
        # if len(self.cleaned_read_recs[sample_type]) == 0:
        #     check = False

        check = True
        utils.log(self.logging_name, 'info', 'Check there are cleaned reads %r' % check)
        return check

    def run_assembly(self):
        '''
        '''

        # Call to fermi-lite
        # fq = self.files['sv_cleaned_fq']

        # print 'Run assembly'
        for ab in self.assembly_batches:
            if len(ab.reads) > 1:
                contig_fn = os.path.join(self.paths['contigs'], "%s_contigs.%d.fastq" % (self.name, ab.index))
                ab.contig_fn = contig_fn
                utils.run_fermi(self.params.get_param('fml-asm'), ab.fq, contig_fn)
                self.parse_fermi_out(ab.index, contig_fn)

    def parse_fermi_out(self, batch_idx, contig_fn):
        '''
        '''

        # print 'FERMI OUT'
        for header,seq,qual in utils.FastqFile(contig_fn):
            contig_values = header.split()
            nreads = int(contig_values[1])
            contig_id = str(batch_idx) + "_" + contig_values[0].lstrip('@').split(':')[0]
            # contig_seq = seq
            contig = Contig(self.paths['contigs'], self.name + '_contig-' + contig_id, seq, nreads, self.sv_khash, self.ref_khash, self.params.get_param('kmer_size'))
            if contig.valid():
                # print contig.seq, contig.valid()
                self.contigs.append(contig)


    def compare_kmers(self):

        '''
        '''

        self.kmers['ref'] = {}
        jellyfish = self.params.get_param('jellyfish')
        kmer_size = int(self.params.get_param('kmer_size'))

        # for i in range(len(self.files['target_ref_fn'])):
        #     utils.log(self.logging_name, 'info', 'Indexing kmers for reference sequence %s' % self.files['target_ref_fn'][i])
        #     self.kmers['ref'] = utils.load_kmers(utils.run_jellyfish(self.files['target_ref_fn'][i], jellyfish, kmer_size), self.kmers['ref'])

        # utils.log(self.logging_name, 'info', 'Indexing kmers for sample sequence %s' % self.files['sv_cleaned_fq'])
        self.kmers['case'] = {}
        self.kmers['case'] = utils.load_kmers(utils.run_jellyfish(self.files['sv_cleaned_fq'], jellyfish, kmer_size), self.kmers['case'])
        # self.kmers['case_sc'] = {}
        # self.kmers['case_sc'] = utils.load_kmers(utils.run_jellyfish(self.files['sv_sc_unmapped_fa'], jellyfish, kmer_size), self.kmers['case_sc'])
        sc_mers = set(self.kmers['case'].keys()) & set(self.kmers['case_sc'])
        sample_only_mers = list(sc_mers.difference(set(self.kmers['ref'].keys())))

        if 'normal_bam_file' in self.params.opts:
            norm_kmers = {}
            norm_kmers = utils.load_kmers(utils.run_jellyfish(self.files['norm_cleaned_fq'], jellyfish, kmer_size), norm_kmers)
            sample_only_mers = set(sample_only_mers).difference(set(norm_kmers.keys()))


        # sample_only_mers = list(sample_only_mers)

        # # Write case only kmers out to file.
        # self.files['sample_kmers'] = os.path.join(self.paths['kmers'], self.name + "_sample_kmers.out")
        # sample_kmer_fout = open(self.files['sample_kmers'], 'w')

        # self.kmers['case_only'] = {}
        # for mer in sample_only_mers:
        #     sample_kmer_fout.write("\t".join([str(x) for x in [mer, str(self.kmers['case'][mer])]]) + "\n")
        #     self.kmers['case_only'][mer] = self.kmers['case'][mer]
        # sample_kmer_fout.close()

        # self.kmers['ref'] = {}
        # self.kmers['case'] = {}
        # self.kmers['case_sc'] = {}

        # utils.log(self.logging_name, 'info', 'Writing %d sample-only kmers to file %s' % (len(self.kmers['case_only']), self.files['sample_kmers']))
        # self.files['kmer_clusters'] = os.path.join(self.paths['kmers'], self.name + "_sample_kmers_merged.out")
        # utils.log(self.logging_name, 'info', 'Writing kmer clusters to file %s' % self.files['kmer_clusters'])
        
        # self.contigs = assembler.init_assembly(self.kmers['case_only'], self.cleaned_read_recs['sv'], kmer_size, int(self.params.get_param('trl_sr_thresh')), self.params.get_param('read_len'))
        # self.cleaned_read_recs = None
        # self.kmers['case_only'] = {}
        # self.finalize_contigs()

    def finalize_contigs(self):
        '''
        '''

        utils.log(self.logging_name, 'info', 'Finalizing %d assembled contigs' % len(self.contigs))
        for contig_iter, assembled_contig in enumerate(self.contigs):
            utils.log(self.logging_name, 'info', 'Finalizing contig %s' % assembled_contig.seq.value)
            contig_id = self.name + '-contig' + str(contig_iter + 1)
            assembled_contig.write_contig_values(contig_id, self.files['kmer_clusters'], self.paths['contigs'])

    def resolve_sv(self):
        '''
        '''

        utils.log(self.logging_name, 'info', 'Resolving structural variants from %d kmer clusters' % len(self.contigs))
        self.results = self.call_manager.resolve_sv_calls(self.contigs, self.files['target_ref_fn'][0], self.get_values(), self.dr_clusters, self.sv_reads)

        # print self.results
        # sys.exit()
        # contig_iter = 1
        # utils.log(self.logging_name, 'info', 'Resolving structural variants from %d kmer clusters' % len(self.contigs))
        # for assembled_contig in self.contigs:
        #     utils.log(self.logging_name, 'info', 'Assessing contig %s' % assembled_contig.seq.value)
        #     contig_id = 'contig' + str(contig_iter)
        #     ctig = contig.TargetContig(self, contig_id, assembled_contig)
        #     ctig.query_ref(self.files['target_ref_fn'][0], self.get_values())
        #     ctig.make_calls(self.get_values(), self.disc_reads, self.repeat_mask)

        #     if ctig.has_result():
        #         ctig.write_result(self.paths['output'])
        #         ctig.write_bam(self.files['sv_bam_sorted'], self.paths['output'])
        #         self.results.append(ctig.result)
        #     else:
        #         utils.log(self.logging_name, 'info', '%s has no structural variant result.' % ctig.id)
        #     contig_iter += 1

    def get_values(self):
        '''
        '''

        return (self.chrom, self.start, self.end, self.name, self.target_intervals)

    def has_results(self):
        '''
        '''

        return len(self.results) > 0

    def add_path(self, key, path):
        '''Utility function to create all the output directories.

        Args:
            key (str):  String value to store the file path value.
            path (str): File path value.
        Returns:
            None
        Raises:
            None
        '''

        utils.log(self.logging_name, 'info', 'Creating %s %s path (%s)' % (self.name, key, path))
        self.paths[key] = path
        if not os.path.exists(self.paths[key]):
            os.makedirs(self.paths[key])

    def set_ref_data(self):

        '''
        '''

        # Write rmask bed file if needed.
        # if not self.params.opts['keep_repeat_regions'] and 'repeat_mask_file' in self.params.opts: 
        #   self.logger.info('Extracting repeat mask regions for target gene %s.' % self.name)
        #   self.repeat_mask = setup_rmask(self.get_values(), self.paths['ref_data'], self.params.opts['repeat_mask_file'])
         
        # Write reference fasta file if needed.
        for target_refseq_fn in self.files['target_ref_fn']:
            direction = "forward"
            if target_refseq_fn.find("forward") == -1:
                direction = "reverse"
            utils.log(self.logging_name, 'info', 'Extracting refseq sequence and writing %s' % target_refseq_fn)
            utils.extract_refseq_fa(self.get_values(), self.paths['ref_data'], self.params.get_param('reference_fasta'), direction, target_refseq_fn, self.params.get_param('buffer_size'))

    # def setup_rmask(self,marker_fn):

    #     '''
    #     '''

    #     # Iterate through genes in target list and find repeats in those genes.
    #     self.repeat_mask = [] 
    #     if not os.path.isfile(marker_fn):
    #       out_fn = self.files['rep_mask_fn']
    #       fout = open(out_fn,'w')
    #       f = open(self.params.opts['repeat_mask_file'],'rU')
    #       flines = f.readlines()
    #       for line in flines:
    #         line = line.strip()
    #         rchr,rbp1,rbp2,rname = line.split("\t")[0:4]
    #         rchr = rchr.replace('chr','')
    #         if rchr == self.chrom:
    #           if int(rbp1) >= self.start and int(rbp2) <= self.end: 
    #             fout.write("\t".join([str(x) for x in [rchr,int(rbp1),int(rbp2),rname]])+"\n")
    #             self.repeat_mask.append((rchr,int(rbp1),int(rbp2),rname))
    #       f.close()
    #       fout.close()
    #       cmd = 'touch %s'%marker_fn
    #       p = subprocess.Popen(cmd,stdout=subprocess.PIPE,stderr=subprocess.PIPE,shell=True)
    #       output, errors = p.communicate()  
    #       self.logger.info('Completed writing repeat mask file %s, touching marker file %s'%(out_fn,marker_fn))
    #     else:
    #       rep_f = open(self.files['rep_mask_fn'],'rU')
    #       rep_flines = rep_f.readlines()
    #       for line in rep_flines:
    #         line = line.strip()
    #         rchr,rbp1,rbp2,rname = line.split()
    #         self.repeat_mask.append((rchr,int(rbp1),int(rbp2),rname))
    #       rep_f.close()

    # def add_discordant_pe(self, aread, read_d, bamfile):
    #     qname = aread.qname
    #     # Keep discordant read pairs
    #     if aread.mapq > 0 and ((aread.rnext!=-1 and aread.tid != aread.rnext) or abs(aread.tlen) > 1000) and not aread.mate_is_unmapped:
    #       mate_refid = bamfile.getrname(aread.rnext)
    #       mate_read = bamfile.mate(aread)
    #       if mate_read.mapq > 0: 
    #         if mate_refid not in read_d['disc']: read_d['disc'][mate_refid] = []
    #         read_d['disc'][mate_refid].append((aread.pos, aread.pnext))
         
    #     if aread.mapq > 0 and not aread.mate_is_unmapped and aread.tid == aread.mrnm:
    #       if aread.is_read1:
    #         read_positions = None
    #         if aread.is_reverse and aread.mate_is_reverse:
    #           # reverse -- reverse, samflag 115 (note: only considering read1, read2 samflag 179)
    #           read_positions = (aread.pos, aread.mpos, 0, 0, qname)
    #           if aread.mpos < aread.pos: read_positions = (aread.mpos, aread.pos, 0, 0, qname)
    #           read_d['inv_reads'].append(read_positions)
    #         elif not aread.is_reverse and not aread.mate_is_reverse:
    #           # forward -- forward = samflag 67 (note: only considering read1, read2 samflag 131)
    #           read_positions = (aread.pos, aread.mpos, 1, 1, qname) 
    #           if aread.mpos < aread.pos: read_positions = (aread.mpos, aread.pos, 1, 1, qname)
    #           read_d['inv_reads'].append(read_positions)
    #         elif aread.is_reverse and not aread.mate_is_reverse and aread.pos < aread.mpos:
    #           # reverse -- forward = samflag 83 with positive insert (read2 samflag 163 with + insert size)
    #           read_positions = (aread.pos, aread.mpos, 0, 1, aread.qname)
    #           read_d['td_reads'].append(read_positions)
    #         elif not aread.is_reverse and aread.mate_is_reverse and aread.mpos < aread.pos:
    #           # reverse -- forward = samflag 99 with - insert (read2 samflag 147 with - insert)
    #           read_positions = (aread.mpos, aread.pos, 1, 0, qname)
    #           read_d['td_reads'].append(read_positions)
    #         if read_positions: read_d['other'].append(read_positions)

    # def pe_meta(self, aread):

    #     '''
    #     '''

    #     # First check if read is from a proper paired-end mapping --> <--    
    #     proper_map = False
    #     overlap_reads = False
    #     if ( ((aread.flag==83) or (aread.flag==147)) and (aread.isize<0) ) or (((aread.flag==99) or (aread.flag==163)) and (aread.isize>0)):
    #       proper_map = True
    #       if abs(aread.isize) < 2*len(aread.seq):
    #         overlap_reads = True   
    #     return proper_map, overlap_reads

    def check_overlap(self, dir, mseq, sc_seq):
        '''
        '''

        if dir == 'back':
            return mseq.find(sc_seq) != (len(mseq)-len(sc_seq))
        else: return mseq.find(sc_seq) != 0


    def check_pair_overlap(self, mate_seq, read, coords, trim_dir):
        '''
        '''

        nmisses = 0
        add_sc = True
        sc_seq = read.seq[coords[0]:coords[1]]
        sc_len = coords[1] - coords[0]
        
        if abs(read.isize) < len(read.seq):
            # Adapter seq
            if abs(len(read.seq) - (abs(read.isize)+1)) >= sc_len: 
                add_sc = False 
    #           print 'Adapter seq', sc_len, abs(read.isize), abs(len(read.seq) - abs(read.isize)), add_sc
        else:
            # abs((2*len(read.seq) - (abs(read.isize)+1)) - sc_len) < 5: add_sc_len_check = False
            while self.check_overlap(trim_dir, mate_seq, sc_seq) and nmisses < 5 and len(sc_seq) > 0:
                if trim_dir == 'back':
                    sc_seq = sc_seq[0:(len(sc_seq)-1)]
                else:
                    sc_seq = sc_seq[1:len(sc_seq)]
                nmisses += 1
        #      print 'Done checking', sc_seq, nmisses
            add_sc = (len(sc_seq) == 0) or (nmisses == 5)
    #      if trim_dir == 'back':
    #        q = read.qual
    #        read.seq = read.seq[coords[1]:len(q)]
    #        read.qual = q[coords[1]:len(q)]
    #      else:
    #        indx = read.seq.find(sc_seq)
    #        q = read.qual
    #        read.seq = read.seq[0:coords[0]]
    #        read.qual = q[0:coords[0]]
    #    print 'Checked read pair overlap', read.qname, read.seq
    #    print 'Using mate seq check', add_sc, sc_seq, mate_seq
        return add_sc #, read

    # def write_results(self):

    #     '''
    #     '''

    #     result_files = {}
    #     for res in self.results:
    #         tag = res[6]
    #         if tag.find('rearrangement') > -1: 
    #             tag = 'rearrangement'
    #         if tag not in result_files:  
    #             header = "\t".join(['genes', 'target_breakpoints', 'align_cigar', 'mismatches', 'strands', 'rep_overlap_segment_len', 'sv_type', 'split_read_count', 'nkmers', 'disc_read_count', 'breakpoint_coverages', 'contig_id', 'contig_seq']) + "\n"
    #             res_fn = os.path.join(self.paths['output'], self.name + "_" + tag + "_svs.out")
    #             utils.log(self.logging_name, 'info', 'Writing %s results to file %s' % (tag, res_fn))
    #             result_files[tag] = open(res_fn, 'w')
    #             if not self.params.opts['no_output_header']:
    #               result_files[tag].write(header)
    #         result_files[tag].write("\t".join([str(x) for x in res]) + "\n")
    #     for f in result_files:
    #         result_files[f].close()

    def write_results(self):
        '''
        '''

        res_fn = os.path.join(self.paths['output'], self.name + "_svs.out")
        result_file = open(res_fn, 'w')
        header = "\t".join(['genes', 'target_breakpoints', 'mismatches', 'strands', 'total_matching', 'sv_type', 'sv_subtype', 'split_read_count', 'disc_read_count', 'breakpoint_coverages', 'contig_id', 'contig_seq']) + "\n"
        result_file.write(header)

        for res in self.results:
            utils.log(self.logging_name, 'info', 'Writing results to file: %s' % res_fn)
            formatted_result_str = res.get_output_string()
            result_file.write(formatted_result_str)
            self.formatted_results.append(formatted_result_str)
        result_file.close()

    def get_sv_counts(self):
        '''
        '''

        total = 0
        rearr_genes = []
        for res in self.results:
          tag = res[6]
          if tag.find('rearrangement') > -1: 
            tag = 'rearrangement'
          if tag == 'rearrangment':
            genes = res[0].split(",")
            genes.sort()
            rearr_genes.append(";".join(genes))
          else: 
            self.svs[tag][0] += 1
            total += 1
        if len(set(rearr_genes)) > 0: 
          total += len(set(rearr_genes))
          self.svs[tag][0] = len(set(rearr_genes))
          self.svs[tag][1] = ",".join(list(set(rearr_genes)))
        return total

    def get_summary(self):
        '''
        '''

        header = ['Target','N_contigs', 'Total_variants']
        total = self.get_sv_counts()
        str_out = self.name + '\t' + str(len(self.contigs)) + '\t' + str(total) + '\t'
        keys = self.svs.keys()
        keys.sort()
        header += ['N_'+str(x) for x in keys]
        rearrs = '-'
        for t in keys:
          if t == 'rearrangment':
            rearrs = self.svs[t][1]
          str_out += str(self.svs[t][0]) +'\t'
        header.append('Rearrangements')
        str_out += rearrs
        return "\t".join(header), str_out

    def rm_output_dir(self):
        '''
        '''

        shutil.rmtree(self.paths['output'])
