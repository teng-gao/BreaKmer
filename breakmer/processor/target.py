#! /usr/bin/python
# -*- coding: utf-8 -*-

import os
import pysam
import subprocess
import logging
import breakmer.utils as utils

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


class target(object):

    '''

    '''

    def __init__(self, intervals, params):

        '''
        '''

        self.name = None
        self.params = params
        self.chrom = None
        self.start = None
        self.end = None
        self.paths = {}
        self.files = {}
        self.disc_reads = None
        self.sv_reads = None
        self.cleaned_read_recs = None
        self.read_len = 0
        self.kmer_clusters = []
        self.kmers = {}
        self.results = []
        self.svs = {'trl':[0, '-'], 'indel':[0, ''], 'rearrangement':[0, '']}
        self.logger = logging.getLogger('breakmer.processor.target')
        self.target_intervals = intervals
        self.repeat_mask = None
        self.logging_name = 'breakmer.processor.target'
        self.setup()

    def setup(self):

        '''
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
                    self.logger.error('Keep repeat regions option is false, but no repeat mask bed file provided. All repeat region variants will be reported.')
                    self.params.opts['keep_repeat_regions'] = True
                else:
                    self.files['rep_mask_fn'] = os.path.join(self.paths['ref_data'], self.name+'_rep_mask.bed')

        self.files['target_ref_fn'] = [os.path.join(self.paths['ref_data'], self.name+'_forward_refseq.fa'), os.path.join(self.paths['ref_data'], self.name+'_reverse_refseq.fa')]

        ref_fa_marker_f = open(os.path.join(self.paths['ref_data'], '.reference_fasta'), 'w')
        ref_fa_marker_f.write(self.params.opts['reference_fasta'])
        ref_fa_marker_f.close()

        self.files['ref_kmer_dump_fn'] = [os.path.join(self.paths['ref_data'], self.name+'_forward_refseq.fa_dump'), os.path.join(self.paths['ref_data'], self.name+'_reverse_refseq.fa_dump')]

    def extract_bam_reads(self, type):
        self.setup_read_extraction_files(type)

        bam_type = 'sample'
        if type == 'norm':
            bam_type = 'normal'

        utils.log(self.logging_name, 'info', 'Extracting bam reads from %s to %s'%(self.params.opts['%s_bam_file' % bam_type], self.files['sv_fq']))
        
        bamfile = pysam.Samfile(self.params.opts['%s_bam_file' % bam_type], 'rb')
        if type == 'sv':
            sv_bam = pysam.Samfile(self.files['sv_bam'], 'wb', template=bamfile)

        read_d = {'unmapped':{},
                  'disc':{},
                  'sv':{},
                  'unmapped_keep':[],
                  'inv_reads':[],
                  'td_reads':[],
                  'other':[]
                }
        utils.log(self.logging_name, 'debug', 'Fetching bam file reads from %s, %s %d %d' % (self.params.opts['%s_bam_file' % bam_type], self.chrom, self.start - 200, self.end + 200))
        areads = bamfile.fetch(self.chrom, self.start - 200, self.end + 200)
        kmer_size = int(self.params.get_param('kmer_size'))
        pair_indices, valid_reads = process_reads(areads, read_d, bamfile)

        for aread, proper_map, overlap_reads in valid_reads:
    #    for aread in areads:
    #      if aread.mate_is_unmapped or aread.rnext == -1:
    #        aread.mate_is_unmapped = True
    #      qname = aread.qname
    #      if aread.is_duplicate or aread.is_qcfail: 
    #        continue
    #      elif aread.is_unmapped:
    #        read_d['unmapped'][aread.qname] = aread
    #        continue

    #      self.add_discordant_pe(aread, read_d, bamfile)
    #      proper_map, overlap_reads = self.pe_meta(aread)

          # Only take soft-clips from outer regions of properly mapped reads, take all others
          # Cigar is a list of tuples 
          if aread.cigar and len(aread.cigar) > 1: 
            tc_coords = trim_coords(aread.qual, 3)
            sc_coords = [0,len(aread.qual)]
            coords = [0,0]
            for i in range(len(aread.cigar)):
              code,clen = aread.cigar[i]
              if not code == 2 and not code == 4: coords[1] += clen
              if code == 4:
                if i == 0: 
                  coords[0] = clen
                  coords[1] += clen
              sc_coords = coords
            # Only keep reads that have a soft clip in sequence that has not been trimmed 
            # due to low quality sequence.           
            sc_seq = {'clipped':[], 'buffered':[]}
            if sc_coords[0] > tc_coords[0] or sc_coords[1] < tc_coords[1]: 
              clip_coords = [0,0]
              s,e = sc_coords
              add_sc = [False, False]
              indel_only = False
              start_sc = s > 0
              end_sc = e < len(aread.qual)
              seq = aread.seq
              ll = len(seq)
              if start_sc and end_sc:
                add_sc = [True, True]
              else:
                if start_sc: 
                  add_sc[0] = True
                  clip_coords = [0,s]
                  if overlap_reads and aread.is_reverse: 
                    mate_seq = valid_reads[pair_indices[aread.qname][int(aread.is_read1)]][0].seq
                    add_sc[0] = self.check_pair_overlap(mate_seq, aread, [0,s], 'back')
                  if proper_map:
                    if aread.is_reverse: indel_only = True
                    else: indel_only = False
                elif end_sc: 
                  clip_coords = [e,ll]
                  add_sc[1] = True
                  if overlap_reads and not aread.is_reverse: 
                    mate_seq = valid_reads[pair_indices[aread.qname][int(aread.is_read1)]][0].seq
                    add_sc[1] = self.check_pair_overlap(mate_seq, aread, [e,ll], 'front')
                  if proper_map:
                    if aread.is_reverse: indel_only = indel_only and False
                    else: indel_only = indel_only and True
              final_add = add_sc[0] or add_sc[1]
              if add_sc[0]:
                sc_seq['buffered'].append(aread.seq[0:(s+kmer_size)])
                sc_seq['clipped'].append(aread.seq[0:s])
              if add_sc[1]:
                sc_seq['buffered'].append(seq[(e-kmer_size):ll])
                sc_seq['clipped'].append(seq[e:ll])
              if final_add:
                read_d['sv'][get_seq_readname(aread)] = (aread,sc_seq,clip_coords,indel_only)

          # If read is mapped and mate is unmapped
          if (aread.pos >= self.start and aread.pos <= self.end) and aread.mapq > 0 and aread.mate_is_unmapped: 
            read_d['unmapped_keep'].append(aread.qname)

        sv_fq = open(self.files['sv_fq'],'w')
        sv_sc_fa = open(self.files['sv_sc_unmapped_fa'],'w')

        for qname in read_d['unmapped_keep']:
          if qname in read_d['unmapped']:
            read = read_d['unmapped'][qname]
            read_d['sv'][get_seq_readname(read)] = (read,None,None,False)
            lout = ">" + read.qname + "\n" + str(read.seq)
            sv_sc_fa.write(lout+"\n") 

        if not self.sv_reads:
          self.sv_reads = {}
        self.sv_reads[type] = {}
        for qname in read_d['sv']:
          aread, sc_seq, cc, indel_only = read_d['sv'][qname]
          self.sv_reads[type][qname] = read_d['sv'][qname]
          if type == 'sv': 
            sv_bam.write(aread)
          lout = fq_line(aread, indel_only, self.params.get_kmer_size(), True)
          if lout: sv_fq.write(lout)
          if sc_seq:
            for sc in sc_seq['buffered']: 
              sv_sc_fa.write(">"+qname+"\n"+sc+"\n")
        self.disc_reads = {'disc':read_d['disc'], 'inv':read_d['inv_reads'], 'td':read_d['td_reads'], 'other':read_d['other']}
        sv_fq.close()
        sv_sc_fa.close()
        bamfile.close()

        if type == 'sv':
          sv_bam.close()
          self.logger.info('Sorting bam file %s to %s'%(self.files['sv_bam'],self.files['sv_bam_sorted']))
          sort(self.files['sv_bam'],self.files['sv_bam_sorted'].replace('.bam',''))
          self.logger.info('Indexing sorted bam file %s'%self.files['sv_bam_sorted'])
          index(self.files['sv_bam_sorted'])

    def get_values(self):
        return (self.chrom, self.start, self.end, self.name, self.target_intervals)

    def rm_output_dir(self):
        shutil.rmtree(self.paths['output'])

    def has_results(self):

        '''
        '''

        if len(self.results) > 0:
            return True
        else:
            return False

    def add_path(self,key,path):

        '''
        '''

        self.logger.info('Creating %s %s path (%s)' % (self.name, key, path))
        self.paths[key] = path
        if not os.path.exists(self.paths[key]):
            os.makedirs(self.paths[key])

    def setup_rmask(self,marker_fn):

        '''
        '''

        # Iterate through genes in target list and find repeats in those genes.
        self.repeat_mask = [] 
        if not os.path.isfile(marker_fn):
          out_fn = self.files['rep_mask_fn']
          fout = open(out_fn,'w')
          f = open(self.params.opts['repeat_mask_file'],'rU')
          flines = f.readlines()
          for line in flines:
            line = line.strip()
            rchr,rbp1,rbp2,rname = line.split("\t")[0:4]
            rchr = rchr.replace('chr','')
            if rchr == self.chrom:
              if int(rbp1) >= self.start and int(rbp2) <= self.end: 
                fout.write("\t".join([str(x) for x in [rchr,int(rbp1),int(rbp2),rname]])+"\n")
                self.repeat_mask.append((rchr,int(rbp1),int(rbp2),rname))
          f.close()
          fout.close()
          cmd = 'touch %s'%marker_fn
          p = subprocess.Popen(cmd,stdout=subprocess.PIPE,stderr=subprocess.PIPE,shell=True)
          output, errors = p.communicate()  
          self.logger.info('Completed writing repeat mask file %s, touching marker file %s'%(out_fn,marker_fn))
        else:
          rep_f = open(self.files['rep_mask_fn'],'rU')
          rep_flines = rep_f.readlines()
          for line in rep_flines:
            line = line.strip()
            rchr,rbp1,rbp2,rname = line.split()
            self.repeat_mask.append((rchr,int(rbp1),int(rbp2),rname))
          rep_f.close()

    def set_ref_data(self):

        '''
        '''

        # Write rmask bed file if needed.
        if not self.params.opts['keep_repeat_regions'] and 'repeat_mask_file' in self.params.opts: 
          self.logger.info('Extracting repeat mask regions for target gene %s.' % self.name)
          self.repeat_mask = setup_rmask(self.get_values(), self.paths['ref_data'], self.params.opts['repeat_mask_file'])
         
        # Write reference fasta file if needed.
        for i in range(len(self.files['target_ref_fn'])):
          fn = self.files['target_ref_fn'][i]
          direction = "forward"
          if fn.find("forward") == -1: direction = "reverse"
          self.logger.info('Extracting refseq sequence and writing %s'%fn)
          extract_refseq_fa(self.get_values(), self.paths['ref_data'], self.params.opts['reference_fasta'], direction, fn)

        # Write alternate reference files.
        if 'target_altref_fn' in self.files:
          if not create_ref_test_fa(os.path.join(self.paths['ref_data'], self.name + "_forward_refseq.fa"), os.path.join(self.paths['ref_data'], self.name + "_start_end_refseq.fa")):
            return

          altref_fns = []
          alt_iter = 1
          for i in range(len(self.files['target_altref_fn'])):
            for j in range(len(self.files['target_altref_fn'][i])):
               fn = self.files['target_altref_fn'][i][j]
               marker_fn = get_marker_fn(fn) 
               if not os.path.isfile(marker_fn):
                 altref_fns.append((self.params.opts['alternate_reference_fastas'][i], fn, alt_iter))
            alt_iter += 1
          if len(altref_fns) > 0:
            create_ref_test_fa(os.path.join(self.paths['ref_data'], self.name + "_forward_refseq.fa"), os.path.join(self.paths['ref_data'], self.name + "_start_end_refseq.fa"))
            for i in range(len(altref_fns)): #range(len(self.files['target_altref_fn'])):
              alt_gene_coords = get_altref_genecoords(self.params.opts['blat'], altref_fns[i][0], os.path.join(self.paths['ref_data'], self.name + "_start_end_refseq.fa"), self.chrom, os.path.join(self.paths['ref_data'], self.name + '_altref_blat_' + str(altref_fns[i][2]) + '.psl'))
              if not alt_gene_coords[2]:
                self.logger.info("No sequence for target gene in %s, no reference kmers extracted."%altref_fns[i][0])
                continue
              gene_vals = (self.chrom, alt_gene_coords[0][1], alt_gene_coords[1][1], self.name, self.target_intervals)
              fn = altref_fns[i][1] #self.files['target_altref_fn'][i][j]
              direction = "forward"
              if fn.find("forward") == -1: direction = "reverse"
              self.logger.info('Extracting alternate refseq sequence and writing %s'%fn)
              extract_refseq_fa(gene_vals, self.paths['ref_data'], altref_fns[i][0], direction, fn)
            # Clean up BLAT files!
    #        os.remove(os.path.join(self.paths['ref_data'], self.name + "_start_end_refseq.fa")) 

    def add_discordant_pe(self, aread, read_d, bamfile):
        qname = aread.qname
        # Keep discordant read pairs
        if aread.mapq > 0 and ((aread.rnext!=-1 and aread.tid != aread.rnext) or abs(aread.tlen) > 1000) and not aread.mate_is_unmapped:
          mate_refid = bamfile.getrname(aread.rnext)
          mate_read = bamfile.mate(aread)
          if mate_read.mapq > 0: 
            if mate_refid not in read_d['disc']: read_d['disc'][mate_refid] = []
            read_d['disc'][mate_refid].append((aread.pos, aread.pnext))
         
        if aread.mapq > 0 and not aread.mate_is_unmapped and aread.tid == aread.mrnm:
          if aread.is_read1:
            read_positions = None
            if aread.is_reverse and aread.mate_is_reverse:
              # reverse -- reverse, samflag 115 (note: only considering read1, read2 samflag 179)
              read_positions = (aread.pos, aread.mpos, 0, 0, qname)
              if aread.mpos < aread.pos: read_positions = (aread.mpos, aread.pos, 0, 0, qname)
              read_d['inv_reads'].append(read_positions)
            elif not aread.is_reverse and not aread.mate_is_reverse:
              # forward -- forward = samflag 67 (note: only considering read1, read2 samflag 131)
              read_positions = (aread.pos, aread.mpos, 1, 1, qname) 
              if aread.mpos < aread.pos: read_positions = (aread.mpos, aread.pos, 1, 1, qname)
              read_d['inv_reads'].append(read_positions)
            elif aread.is_reverse and not aread.mate_is_reverse and aread.pos < aread.mpos:
              # reverse -- forward = samflag 83 with positive insert (read2 samflag 163 with + insert size)
              read_positions = (aread.pos, aread.mpos, 0, 1, aread.qname)
              read_d['td_reads'].append(read_positions)
            elif not aread.is_reverse and aread.mate_is_reverse and aread.mpos < aread.pos:
              # reverse -- forward = samflag 99 with - insert (read2 samflag 147 with - insert)
              read_positions = (aread.mpos, aread.pos, 1, 0, qname)
              read_d['td_reads'].append(read_positions)
            if read_positions: read_d['other'].append(read_positions)

    def pe_meta(self, aread):
        # First check if read is from a proper paired-end mapping --> <--    
        proper_map = False
        overlap_reads = False
        if ( ((aread.flag==83) or (aread.flag==147)) and (aread.isize<0) ) or (((aread.flag==99) or (aread.flag==163)) and (aread.isize>0)):
          proper_map = True
          if abs(aread.isize) < 2*len(aread.seq):
            overlap_reads = True   
        return proper_map, overlap_reads

    def get_sv_reads(self):
        self.extract_bam_reads('sv')
        if 'normal_bam_file' in self.params.opts:
          self.extract_bam_reads('norm')
          self.clean_reads('norm')
        
        check = True
        if not self.clean_reads('sv'):
          self.rm_output_dir()
          check = False
        return check

    def setup_read_extraction_files(self, type):
        self.files['%s_fq'%type] = os.path.join(self.paths['data'],self.name + "_sv_reads.fastq")
        self.files['%s_sc_unmapped_fa'%type] = os.path.join(self.paths['data'],self.name + "_sv_sc_seqs.fa")
        if type == 'sv':
          self.files['sv_bam'] = os.path.join(self.paths['data'],self.name + "_sv_reads.bam")
          self.files['sv_bam_sorted'] = os.path.join(self.paths['data'],self.name + "_sv_reads.sorted.bam")

    def check_overlap(self, dir, mseq, sc_seq):
    #    print dir, sc_seq, mseq, mseq.find(sc_seq)
        if dir == 'back': return mseq.find(sc_seq) != (len(mseq)-len(sc_seq))
        else: return mseq.find(sc_seq) != 0

  # Move to utils
    def check_pair_overlap(self, mate_seq, read, coords, trim_dir):
        nmisses = 0
        add_sc = True
        sc_seq = read.seq[coords[0]:coords[1]]
        sc_len = coords[1] - coords[0]
        
        if abs(read.isize) < len(read.seq):
          # Adapter seq
          if abs(len(read.seq) - (abs(read.isize)+1)) >= sc_len: 
            add_sc = False 
    #      print 'Adapter seq', sc_len, abs(read.isize), abs(len(read.seq) - abs(read.isize)), add_sc
        else:
          # abs((2*len(read.seq) - (abs(read.isize)+1)) - sc_len) < 5: add_sc_len_check = False
          while self.check_overlap(trim_dir, mate_seq, sc_seq) and nmisses < 5 and len(sc_seq) > 0:
            if trim_dir == 'back': sc_seq = sc_seq[0:(len(sc_seq)-1)]
            else: sc_seq = sc_seq[1:len(sc_seq)]
            nmisses += 1
    #      print 'Done checking', sc_seq, nmisses
          if len(sc_seq) == 0 or nmisses == 5:
            add_sc = True
          else:
            add_sc = False
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

    def clean_reads(self, type):
        # Run cleaning program
        cutadapt = self.params.opts['cutadapt']
        cutadapt_config = self.params.opts['cutadapt_config_file']
        self.logger.info('Cleaning reads using %s with configuration file %s'%(cutadapt,cutadapt_config))
        self.files['%s_cleaned_fq'%type] = os.path.join(self.paths['data'],self.name + "_%s_reads_cleaned.fastq"%type)
        self.logger.info('Writing clean reads to %s'%self.files['%s_cleaned_fq'%type])
        cutadapt_parameters = stringify(cutadapt_config)
        cmd = '%s %s %s %s > %s'%(sys.executable, cutadapt, cutadapt_parameters, self.files['%s_fq'%type], self.files['%s_cleaned_fq'%type])
        self.logger.debug('Cutadapt system command %s'%cmd)
        p = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True)
        output, errors = p.communicate()
        self.logger.debug('Clean reads output %s'%output)
        self.logger.debug('Clean reads errors %s'%errors)

        # Use these for pulling out reads after finding sample-only kmers.
        # Filter the cleaned reads to make sure soft clips were not adapters, re-write fastq
        if not self.cleaned_read_recs:
          self.cleaned_read_recs = {}
        self.cleaned_read_recs[type] = None
        self.files['%s_cleaned_fq'%type], self.cleaned_read_recs[type], self.read_len = get_fastq_reads(self.files['%s_cleaned_fq'%type], self.sv_reads[type])
        self.sv_reads[type] = None
        check = True
        if len(self.cleaned_read_recs[type]) == 0: check = False
        self.logger.info('Check there are cleaned reads %r'%check)
        return check 

    def compare_kmers(self):
        self.kmers['ref'] = {}
        jellyfish = self.params.opts['jellyfish']
        kmer_size = self.params.get_kmer_size()
        for i in range(len(self.files['target_ref_fn'])): 
          self.logger.info('Indexing kmers for reference sequence %s'%self.files['target_ref_fn'][i])
          self.kmers['ref'] = load_kmers(run_jellyfish(self.files['target_ref_fn'][i], jellyfish, kmer_size), self.kmers['ref'])

        if 'target_altref_fn' in self.files:
          for i in range(len(self.files['target_altref_fn'])):
            for j in range(len(self.files['target_altref_fn'][i])):
              self.logger.info('Indexing kmers for reference sequence %s'%self.files['target_altref_fn'][i])
              self.kmers['ref'] = load_kmers(run_jellyfish(self.files['target_altref_fn'][i][j], jellyfish, kmer_size), self.kmers['ref'])

        self.logger.info('Indexing kmers for sample sequence %s'%self.files['sv_cleaned_fq'])
        self.kmers['case'] = {} 
        self.kmers['case'] = load_kmers(run_jellyfish(self.files['sv_cleaned_fq'],jellyfish,kmer_size),self.kmers['case'])
        self.kmers['case_sc'] = {}
        self.kmers['case_sc'] = load_kmers(run_jellyfish(self.files['sv_sc_unmapped_fa'],jellyfish,kmer_size),self.kmers['case_sc'])
        sc_mers = set(self.kmers['case'].keys()) & set(self.kmers['case_sc']) 
        sample_only_mers = list(sc_mers.difference(set(self.kmers['ref'].keys())))

        if 'normal_bam_file' in self.params.opts:
          norm_kmers = {}
          norm_kmers = load_kmers(run_jellyfish(self.files['norm_cleaned_fq'],jellyfish,kmer_size),norm_kmers)
          sample_only_mers = set(sample_only_mers).difference(set(norm_kmers.keys()))

        sample_only_mers = list(sample_only_mers)

        # Write case only kmers out to file.
        self.files['sample_kmers'] = os.path.join(self.paths['kmers'],self.name + "_sample_kmers.out")
        sample_kmer_fout = open(self.files['sample_kmers'],'w')
        kmer_counter = 1
        self.kmers['case_only'] = {}
        for mer in sample_only_mers:
          sample_kmer_fout.write("\t".join([str(x) for x in [mer,str(self.kmers['case'][mer])]])+"\n")
          self.kmers['case_only'][mer] = self.kmers['case'][mer]
        sample_kmer_fout.close()

        self.kmers['ref'] = {}
        self.kmers['case'] = {}
        self.kmers['case_sc'] = {}

        self.logger.info('Writing %d sample-only kmers to file %s'%(len(self.kmers['case_only']),self.files['sample_kmers']))
        self.files['kmer_clusters'] = os.path.join(self.paths['kmers'],self.name + "_sample_kmers_merged.out")
        self.logger.info('Writing kmer clusters to file %s'%self.files['kmer_clusters'])
        
        self.kmers['clusters'] = init_assembly(self.kmers['case_only'], self.cleaned_read_recs['sv'], kmer_size, self.params.get_sr_thresh('min'), self.read_len)
        self.cleaned_read_recs = None
        self.kmers['case_only'] = {}

    def resolve_sv(self):
        iter = 1
        self.logger.info('Resolving structural variants from %d kmer clusters'%len(self.kmers['clusters']))
        for kc in self.kmers['clusters']:
          self.logger.info('Assessing contig %s'%kc.aseq.seq)
          contig_id = 'contig' + str(iter)
          ctig = contig(self, contig_id, kc)
          ctig.query_ref(self.files['target_ref_fn'][0], self.get_values()) 
          ctig.make_calls(self.get_values(), self.disc_reads, self.repeat_mask)

          if ctig.has_result():
            ctig.write_result(self.paths['output'])
            ctig.write_bam(self.files['sv_bam_sorted'],self.paths['output'])
            self.results.append(ctig.result)
          else: 
            self.logger.info('%s has no structural variant result.'%ctig.id)
          iter += 1

    def write_results(self):
        result_files = {}
        for res in self.results:
          tag = res[6]
          if tag.find('rearrangement') > -1: 
            tag = 'rearrangement'
          if tag not in result_files:  
            header = "\t".join(['genes', 'target_breakpoints', 'align_cigar', 'mismatches', 'strands', 'rep_overlap_segment_len', 'sv_type', 'split_read_count', 'nkmers', 'disc_read_count', 'breakpoint_coverages', 'contig_id', 'contig_seq']) + "\n"
            res_fn = os.path.join(self.paths['output'],self.name+"_"+tag+"_svs.out")
            self.logger.info('Writing %s results to file %s'%(tag,res_fn))
            result_files[tag] = open(res_fn,'w')
            if not self.params.opts['no_output_header']:
              result_files[tag].write(header)
          result_files[tag].write("\t".join([str(x) for x in res]) + "\n")
        for f in result_files: result_files[f].close()

    def get_sv_counts(self):
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
        header = ['Target','N_contigs', 'Total_variants']
        total = self.get_sv_counts()
        str_out = self.name + '\t' + str(len(self.kmers['clusters'])) + '\t' + str(total) + '\t'
        keys = self.svs.keys()
        keys.sort()
        header += ['N_'+str(x) for x in keys]
        rearrs = '-'
        for t in keys:
          if t == 'rearrangment': rearrs = self.svs[t][1]
          str_out += str(self.svs[t][0]) +'\t'
        header.append('Rearrangements')
        str_out += rearrs
        return "\t".join(header), str_out
