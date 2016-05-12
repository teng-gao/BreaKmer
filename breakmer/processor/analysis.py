#! /usr/bin/python
# -*- coding: utf-8 -*-

import os
import sys
import logging
import time
import math
import subprocess
import multiprocessing
import breakmer.processor.target as target
import breakmer.utils as utils
import sv_caller
import pysam
from itertools import izip, islice, repeat, izip_longest

__author__ = "Ryan Abo"
__copyright__ = "Copyright 2015, Ryan Abo"
__email__ = "ryanabo@gmail.com"
__license__ = "MIT"


def process_reads(areads, read_d, bamfile):

    '''
    '''

    pair_indices = {}
    valid_reads = []

    for aread in areads:
        skip = False
        if aread.mate_is_unmapped or aread.rnext == -1: # Indicate that mate is unmapped
            aread.mate_is_unmapped = True
        if aread.is_duplicate or aread.is_qcfail: # Skip duplicates and failures
            skip = True
        if aread.is_unmapped: # Store unmapped reads
            read_d['unmapped'][aread.qname] = aread
            skip = True

        # If read is unmapped or duplicate or qcfail, then don't store
        if not skip:
            proper_map = False
            overlap_reads = False
            # These two functions can opeate on the first read of the pair.
            # Check if fragment hasn't been checked yet and that the mate is mapped.
            if aread.qname not in pair_indices and not aread.mate_is_unmapped:
                add_discordant_pe(aread, read_d, bamfile)
                proper_map, overlap_reads = pe_meta(aread)
            valid_reads.append((aread, proper_map, overlap_reads))

            if aread.qname not in pair_indices and not aread.mate_is_unmapped:
                pair_indices[aread.qname] = {}
            if aread.qname in pair_indices:
                pair_indices[aread.qname][int(aread.is_read1)] = len(valid_reads)-1
    return pair_indices, valid_reads


def pe_meta(aread):
    # First check if read is from a proper paired-end mapping --> <--    
    proper_map = False
    overlap_reads = False
    if (((aread.flag == 83) or (aread.flag == 147)) and (aread.tlen < 0)) or (((aread.flag == 99) or (aread.flag == 163)) and (aread.tlen > 0)):
        proper_map = True
        if abs(aread.tlen) < (2 * len(aread.seq)):
            overlap_reads = True
    return proper_map, overlap_reads


def add_discordant_pe(aread, read_d, bamfile):
  qname = aread.qname
  # Keep discordant read pairs where the map quality is > 0, the paired reads are mapped to different chroms or > 1000 bp apart, and
  # the mate is mapped.
  if aread.mapq > 0 and ((aread.rnext != -1 and aread.tid != aread.rnext) or abs(aread.tlen) > 1000) and not aread.mate_is_unmapped:
    mate_refid = bamfile.getrname(aread.rnext) # Grab the paired read
    mate_read = bamfile.mate(aread)
    if mate_read.mapq > 0:
      if mate_refid not in read_d['disc']:
        read_d['disc'][mate_refid] = []
      read_d['disc'][mate_refid].append((aread.pos, aread.pnext)) # Store the read position and the mate position

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


class RunTracker(object):

    """Class to manage the running of all the target region analyses.
    The params object is passed in with all the input information.
    The run() function creates the target region objects from the
    param inputs and then starts the analysis for each target.

    Args:
        params: ParamManager object.
    Returns:
        None
    """

    def __init__(self, params):

        '''
        '''

        self.params = params
        self.logging_name = 'breakmer.processor.analysis'
        self.results = []
        self.targets = {}
        self.summary = {}
        self.summary_header = ''

    def preset_ref_data(self):

        '''
        '''

        nprocs = int(self.params.get_param('nprocs'))
        ngroups = nprocs
        ntargets = len(self.params.targets)
        ntargets_per_group = ntargets/nprocs
        modval = math.fmod(ntargets, nprocs)
        if modval > 0:
            ngroups += 1
        p = multiprocessing.Pool(nprocs)
        trgt_groups = []
        trgt_group = []

        target_names = self.targets.keys()
        target_names.sort()
        for trgt_name in target_names:
            trgt = self.targets[trgt_name]
            trgt_vals = [trgt.chrom, trgt.start, trgt.end, trgt.name, trgt.target_intervals]
            if len(trgt_group) == ntargets_per_group:
                trgt_groups.append(trgt_group)
                trgt_group = []
            trgt_group.append(trgt_vals)

        if len(trgt_group) < ntargets_per_group:
            trgt_groups[-1].extend(trgt_group)
        else:
            trgt_groups.append(trgt_group)

        mask_fn = None
        if 'keep_repeat_regions' in self.params.opts:
            if not self.params.opts['keep_repeat_regions']:
                if 'repeat_mask_file' not in self.params.opts:
                    utils.log(self.logging_name, 'error', 'Keep repeat regions option is false, but no repeat mask bed file provided. All repeat region variants will be reported.')
                    self.params.opts['keep_repeat_regions'] = True
                else:
                    mask_fn = self.params.opts['repeat_mask_file']

        ref_fa_fn = self.params.opts['reference_fasta']
        altref_fa_fns = None
        if 'alternate_reference_fastas' in self.params.opts:
            altref_fa_fns = ",".join(self.params.opts['alternate_reference_fastas'])
        ref_data_dir = self.params.opts['reference_data_dir']
        jfish_path = self.params.opts['jellyfish']
        blat_path = self.params.opts['blat']
        ref_params = [mask_fn, ref_fa_fn, altref_fa_fns, ref_data_dir, jfish_path, blat_path, self.params.get_param('kmer_size')]
        setup_params = izip(trgt_groups, repeat(ref_params))
        p.map(utils.setup_ref_data, setup_params)

    def create_targets(self):

        '''
        '''

        targets = self.params.targets.keys()
        targets.sort()
        for trgt_name in targets:
            self.targets[trgt_name] = target.target(self.params.targets[trgt_name], self.params)
        # return targets

    def run(self):

        '''
        '''

        start_time = time.clock()  # Track the run time.
        self.create_targets()

        # if self.params.opts['preset_ref_data']:
        if self.params.fnc_cmd == 'prepare_reference_data':
            utils.log(self.logging_name, 'info', 'Creating all reference data.')
            self.preset_ref_data()
            print 'Reference data created!'
            return

        utils.start_blat_server(self.params)
        if self.params.fnc_cmd == 'start_blat_server':
            print 'Server started!'
            return

        trgt_lst = self.params.targets.keys()
        trgt_lst.sort()
        for trgt_name in trgt_lst:
            trgt = self.targets[trgt_name]
            utils.log(self.logging_name, 'info', 'Analyzing %s' % trgt.name)

            # Write reference sequence fasta for gene if it doesn't exist.
            if not self.params.opts['preset_ref_data']:
                trgt.set_ref_data()

            if not trgt.get_sv_reads():
                continue

            trgt.compare_kmers() # Get reference and case kmers
            trgt.resolve_sv() # Build contigs and blat them against the reference genome
            self.summary_header, trgt_summary = trgt.get_summary()
            self.summary[trgt.name] = trgt_summary
            utils.log(self.logging_name, 'info', '%s summary\n%s\n%s' % (trgt.name, self.summary_header, trgt_summary))
            if trgt.has_results():
                trgt.write_results()
                self.results.extend(trgt.results)
            else:
                trgt.rm_output_dir()
         
        self.write_output()
        time_to_complete = time.clock() - start_time
        utils.log(self.logging_name, 'info', 'Analysis complete, %s' % str(time_to_complete))

        # Stop gfServer
        if not self.params.get_param('keep_blat_server'):
            utils.log(self.logging_name, 'info', 'Stopping blat server on port %d' % self.params.get_param('blat_port'))
            gfserver_stop_cmd = '%s stop localhost %d' % (self.params.opts['gfserver'], self.params.opts['blat_port'])
        os.system(gfserver_stop_cmd)

        print 'Analysis complete.'

    def write_output(self):

        '''
        '''

        result_files = {}
        for res in self.results:
            tag = res[6]
            if tag not in result_files:  
                header = "\t".join(['genes', 'target_breakpoints', 'align_cigar', 'mismatches', 'strands', 'rep_overlap_segment_len', 'sv_type', 'split_read_count', 'nkmers', 'disc_read_count', 'breakpoint_coverages', 'contig_id', 'contig_seq']) + "\n"

                res_fn = os.path.join(self.params.paths['output'], self.params.opts['analysis_name'] + "_" + tag + "_svs.out")

                utils.log(self.logging_name, 'info', 'Writing %s output file %s' % (tag, res_fn))
                result_files[tag] = open(res_fn, 'w')
                if not self.params.opts['no_output_header']:
                    result_files[tag].write(header)
                result_files[tag].write("\t".join([str(x) for x in res]) + "\n")

        for result_filename in result_files:
            result_files[result_filename].close()

        summary_fn = os.path.join(self.params.paths['output'], self.params.opts['analysis_name'] + "_summary.out")
        summary_f = open(summary_fn, 'w')
        utils.log(self.logging_name, 'info', 'Writing summary file to %s' % summary_fn)
        summary_f.write(self.summary_header+"\n")

        keys = self.summary.keys()
        keys.sort()
        for gene in keys:
            summary_f.write(self.summary[gene]+"\n")
        summary_f.close()

##### OLD ##########

class contig(object):

    '''
    '''

    def __init__(self, parent_target, contig_id, assembly):
        self.params = parent_target.params
        self.id = contig_id
        self.query_region = parent_target.get_values()
        self.path = os.path.join(parent_target.paths['contigs'], contig_id)
        self.assembly_fq_fn = os.path.join(parent_target.paths['contigs'], contig_id, contig_id + ".fq")
        self.reads = assembly.reads
        self.kmers = assembly.kmers
        self.contig_seq = assembly.get_contig_seq()
        self.contig_rcounts = assembly.get_contig_counts()
        self.contig_kmer_locs = assembly.get_kmer_locs()
        self.result = None
        self.contig_fa_fn = None
        self.query_res_fn = None
        self.logger = logging.getLogger('root')
        self.setup(parent_target.files['kmer_clusters'])

    def setup(self, cluster_fn):
        self.logger.info('Setting up contig path %s'%self.path)
        if not os.path.exists(self.path): os.makedirs(self.path)
        self.write_cluster_file(cluster_fn)
        self.write_read_fq()
        self.write_contig_fa()

    def write_cluster_file(self, cluster_fn):
        cluster_f = open(cluster_fn, 'w')
        cluster_f.write(self.id + " " + str(len(self.kmers)) + "\n")
        cluster_f.write(",".join([x[0] for x in self.kmers]) + "\n")
        cluster_f.write(",".join([x.id for x in self.reads]) + "\n\n")
        cluster_f.close()

    def write_read_fq(self):
        assembly_fq = open(self.assembly_fq_fn, 'w')
        self.logger.info('Writing reads containing kmers to fastq %s'%self.assembly_fq_fn)
        for read in self.reads:
            assembly_fq.write(read.id+"\n"+read.seq+"\n+\n"+read.qual+"\n")
        assembly_fq.close()

    def write_contig_fa(self):
        self.contig_fa_fn = os.path.join(self.path, self.id+".fa")
        self.logger.info('Writing contig fasta file for blatting %s'%self.contig_fa_fn)
        blat_f = open(self.contig_fa_fn, 'w')
        blat_f.write(">contig1"+"\n"+self.contig_seq)
        blat_f.close()

    def has_result(self):
        if self.result:
            return True
        else:
            return False

    def write_result(self, output_path):
        res_fn = os.path.join(self.path,self.id + "_svs.out")
        self.logger.info('Writing %s result file %s'%(self.id,res_fn))
        if self.result: 
            res_f = open(res_fn,'w')
            res_f.write("\t".join([str(x) for x in self.result]))
            res_f.close()
            shutil.copyfile(res_fn, os.path.join(output_path, self.id+"_svs.out"))

    def write_bam(self, bam_in, path):
        bam_out_fn = os.path.join(path,self.id+"_reads.bam")
        self.logger.info('Writing contig reads bam file %s'%bam_out_fn)
        bam_out_sorted_fn = os.path.join(path,self.id+"_reads.sorted.bam")
        bamf = Samfile(bam_in,'rb')
        bam_out_f = Samfile(bam_out_fn,"wb",template=bamf)
        for bam_read in bamf.fetch():
          for read in self.reads:
            rid, idx = read.id.lstrip("@").split("/")
            ridx, indel_only_read = idx.split("_")
            if (bam_read.qname == rid) and ((ridx=='2' and bam_read.is_read2) or (ridx=='1' and bam_read.is_read1)):
              bam_out_f.write(bam_read)
        bamf.close()
        bam_out_f.close()
        self.logger.info('Sorting bam file %s to %s'%(bam_out_fn,bam_out_sorted_fn))
        pysam.sort(bam_out_fn,bam_out_sorted_fn.replace('.bam',''))
        self.logger.info('Indexing bam file %s'%bam_out_sorted_fn)
        pysam.index(bam_out_sorted_fn)

    def query_ref(self, target_ref_fn, query_region):
        if self.contig_fa_fn:
          self.run_blat(target_ref_fn, 'target') # Run blat against target reference sequence first for speed.
          if not self.query_res_fn:
            self.logger.info('No blat results file %s, no calls for %s.'%(self.query_res_fn, self.id))
            return
          if not self.check_target_blat(query_region):
            # Blat against whole genome reference fasta
            self.run_blat(self.params.opts['reference_fasta'], 'all')

    def run_blat(self, db, name):
        self.query_res_fn = os.path.join(self.path,'blat_res.'+name+'.psl')
        if not os.path.isfile(self.query_res_fn):
          if name == 'all': 
            self.logger.info('Running blat %s, storing results in %s'%(self.params.opts['gfclient'],self.query_res_fn))
            cmd = '%s -t=dna -q=dna -out=psl -minScore=20 -nohead localhost %d %s %s %s'%(self.params.opts['gfclient'],self.params.opts['blat_port'], self.params.opts['reference_fasta_dir'], self.contig_fa_fn, self.query_res_fn)
          else:
            self.logger.info('Running blat %s, storing results in %s'%(self.params.opts['blat'],self.query_res_fn))
            cmd = '%s -t=dna -q=dna -out=psl -minScore=20 -stepSize=10 -minMatch=2 -repeats=lower -noHead %s %s %s'%(self.params.opts['blat'], db, self.contig_fa_fn, self.query_res_fn)
    #        cmd = '%s -t=dna -q=dna -out=psl -minScore=20 -stepSize=1 -minMatch=1 -repeats=lower -noHead %s %s %s'%(self.params.opts['blat'], db, self.contig_fa_fn, self.query_res_fn)

          self.logger.info('Blat system command %s'%cmd)
          p = subprocess.Popen(cmd,stdout=subprocess.PIPE,stderr=subprocess.PIPE,shell=True)
          output, errors = p.communicate()
          self.logger.info('Blat output %s'%output)
          if errors != '': self.logger.info('Blat errors %s'%errors)
        else: self.logger.info('Blat already run, results file %s exists, continuing'%self.query_res_fn)

    def check_target_blat(self, query_region ):
        meta_dict = {'offset': query_region[1]-200, 'tname': query_region[0].replace('chr',''), 'query_res_fn': self.query_res_fn, 'sbam': self.params.opts['sample_bam_file']}
        am = sv_caller.align_manager(meta_dict)
        hit, self.query_res_fn = am.check_target_results()
        return hit

    def make_calls(self, query_region, disc_reads, rep_mask):
        meta_dict = {'params': self.params, 'repeat_mask': rep_mask, 'query_region': query_region, 'query_res_fn': self.query_res_fn, 'disc_reads': disc_reads, 'contig_vals': (self.contig_seq,self.contig_rcounts,self.id,self.reads,len(self.kmers), self.contig_kmer_locs), 'sbam': self.params.opts['sample_bam_file']}
        am = sv_caller.align_manager(meta_dict)
        self.result = am.get_result()

'''
    if not self.query_res_fn: 
      self.logger.info('No blat results file %s, no calls for %s.'%(self.query_res_fn, self.id))
      return

    self.logger.info('Making variant calls from blat results %s'%self.query_res_fn)
    blat_f = open(self.query_res_fn,'r') # no header blat result psl
    clipped_queries = [] 
    se = None
    qsize = None
    blat_results = []
    hit_freq = []
    for line in blat_f.readlines():
      line = line.strip()
      res_d['blat_values'] = line.split("\t")
      br = blat_res(res_d)
      if not qsize: 
        qsize = br.get_size('query')
        hit_freq = [0]*qsize
      for i in range(br.get_coords('query')[0],br.get_coords('query')[1]): hit_freq[i] += 1 
      blat_results.append((br.get_nmatch_total(),br))

    blat_results_sorted = sorted(blat_results, key=lambda blat_results: blat_results[0]) 
    blat_results_sorted.reverse()
    for i in range(len(blat_results_sorted)):
      nmatch, br = blat_results_sorted[i]
      mean_freq = float(sum(hit_freq[br.get_coords('query')[0]:br.get_coords('query')[1]])) / float(len(hit_freq[br.get_coords('query')[0]:br.get_coords('query')[1]]))
      if br.get_size('query') == (br.get_coords('query')[1]-br.get_coords('query')[0]):
        if br.valid and mean_freq < 2 and br.in_target and not br.in_repeat:
          se = sv_event(br,query_region,self.contig_seq,self.contig_rcounts,self.id)
          self.logger.debug('Top hit contains whole query sequence, indel variant')
          break
        else: return None
      elif len(blat_results_sorted) == 1 and br.in_target:
        if not br.in_repeat:
          se = sv_event(br,query_region,self.contig_seq,self.contig_rcounts,self.id)
          self.logger.debug('One blat result within target, indel variant')  
        else: return None
      elif (mean_freq < 4 and ((br.get_nmatch_total()<30 and not br.in_repeat) or br.get_nmatch_total()>=30)) or (br.in_target and mean_freq < 10): 
        clipped_queries.append((br.get_coords('query')[0],br.get_coords('query')[1],br,i))
        
    gaps = [(0,qsize)]
    if len(clipped_queries) > 1:
      self.logger.debug('Iterating through %d clipped blat results.'%len(clipped_queries))
      merged_clip = [0,None]
      for i in range(len(clipped_queries)):
        qs, qe, br, idx = clipped_queries[i]
        segment_size = qe-qs
        self.logger.debug('Blat result with start %d, end %d, chrom %s'%(qs,qe,br.get_name('hit')))
        new_gaps = []
        for gap in gaps:
          gs, ge = gap
          self.logger.debug('Gap coords %d, %d'%(gs, ge))
          if (qs >= gs and qs <= ge) or (qe <= ge and qe >= gs):
            ngap = []
            if qs > gs: 
              if (qs-1-gs) > 10: 
                ngap.append((gs,qs-1))
            if qe < ge: 
              if (ge-qe+1) > 10:
                ngap.append((qe+1,ge))
            if i==0: 
              se = sv_event(br,query_region,self.contig_seq,self.contig_rcounts,self.id)
              new_gaps.extend(ngap)
            else:
              # Calc % of segment overlaps with gap
              over_perc = round((float(min(qe,ge)-max(qs,gs)) / float(segment_size)) * 100)
              # Check overlap with other aligned segments
              ov_right = 0
              if qe > ge: ov_right = abs(qe-ge)
              ov_left = 0
              if qs < gs: ov_left = abs(qs-gs)
              max_seg_overlap = max(ov_right,ov_left)
              self.logger.debug('Blat query segment overlaps gap by %f'%over_perc)
              self.logger.debug('Max segment overlap %f'%max_seg_overlap)
              self.logger.debug('Event in target %r and blat result in target %r'%(se.in_target, br.in_target))
              if over_perc >= 50 and max_seg_overlap < 15 and (se.in_target or br.in_target): 
                self.logger.debug('Adding blat result to event')
                new_gaps.extend(ngap)
                se.add(br)
          else:
              new_gaps.append(gap)
          self.logger.debug('New gap coords %s'%(",".join([str(x) for x in new_gaps])))
        gaps = new_gaps
        if se.qlen > merged_clip[0] and se.in_target: merged_clip = [se.qlen,se]
      se = merged_clip[1]
            
    blat_f.close()
    if se:
      self.result = se.get_result(query_region, [len(self.reads),len(self.kmers)], disc_reads, self.params, rep_mask)
'''
