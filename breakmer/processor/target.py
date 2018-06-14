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
import breakmer.caller.sv_caller2 as sv_caller
import pdb

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
        self.sv_reads = None
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

    def get_sv_reads(self):

        '''
        '''

        self.extract_bam_reads('sv')
        if 'normal_bam_file' in self.params.opts:
            self.extract_bam_reads('norm')
            self.clean_reads('norm')

        check = True
        if not self.clean_reads('sv'):
            shutil.rmtree(self.paths['output'])
            check = False
        return check

    def setup_read_extraction_files(self, sample_type):

        '''
        '''

        self.files['%s_fq' % sample_type] = os.path.join(self.paths['data'], self.name + "_sv_reads.fastq")
        self.files['%s_sc_unmapped_fa' % sample_type] = os.path.join(self.paths['data'], self.name + "_sv_sc_seqs.fa")
        if sample_type == 'sv':
          self.files['sv_bam'] = os.path.join(self.paths['data'], self.name + "_sv_reads.bam")
          self.files['sv_bam_sorted'] = os.path.join(self.paths['data'], self.name + "_sv_reads.sorted.bam")

    def extract_bam_reads(self, sample_type):

        '''
        '''
        self.setup_read_extraction_files(sample_type)

        bam_type = 'sample'
        if sample_type == 'norm':
            bam_type = 'normal'

        utils.log(self.logging_name, 'info', 'Extracting bam reads from %s to %s' % (self.params.opts['%s_bam_file' % bam_type], self.files['sv_fq']))

        bamfile = pysam.Samfile(self.params.opts['%s_bam_file' % bam_type], 'rb')
        if sample_type == 'sv':
            sv_bam = pysam.Samfile(self.files['sv_bam'], 'wb', template=bamfile)

        read_d = {'unmapped':{},
                  'disc':{},
                  'sv':{},
                  'unmapped_keep':[],
                  'inv_reads':[],
                  'td_reads':[],
                  'other':[]
                 }

        buffer_size = int(self.params.get_param('buffer_size'))
        kmer_size = int(self.params.get_param('kmer_size'))

        utils.log(self.logging_name, 'debug', 'Fetching bam file reads from %s, %s %d %d' % (self.params.opts['%s_bam_file' % bam_type], self.chrom, self.start - buffer_size, self.end + buffer_size))
        aligned_reads = bamfile.fetch(self.chrom, self.start - buffer_size, self.end + buffer_size)

        pair_indices = {}
        valid_reads = []

        for aligned_read in aligned_reads:

            if aligned_read.is_duplicate or aligned_read.is_qcfail:  # Skip duplicates and failures
                continue
            if aligned_read.is_unmapped:  # Store unmapped reads
                read_d['unmapped'][aligned_read.qname] = aligned_read
                continue

            if aligned_read.mate_is_unmapped or aligned_read.rnext == -1:  # Indicate that mate is unmapped
                aligned_read.mate_is_unmapped = True

            proper_map = False
            overlap_reads = False

            # These two functions can operate on the first read of the pair.
            # Check if fragment hasn't been checked yet and that the mate is mapped.
            if aligned_read.qname not in pair_indices and not aligned_read.mate_is_unmapped:
                add_discordant_pe(aligned_read, read_d, bamfile)
                proper_map, overlap_reads = pe_meta(aligned_read)
            valid_reads.append((aligned_read, proper_map, overlap_reads))

            if aligned_read.qname not in pair_indices and not aligned_read.mate_is_unmapped:
                pair_indices[aligned_read.qname] = {}
            if aligned_read.qname in pair_indices:
                pair_indices[aligned_read.qname][int(aligned_read.is_read1)] = len(valid_reads) - 1

            # If read is mapped and mate is unmapped
            if (aligned_read.pos >= self.start and aligned_read.pos <= self.end) and aligned_read.mapq > 0 and aligned_read.mate_is_unmapped:
                read_d['unmapped_keep'].append(aligned_read.qname)
        # pair_indices, valid_reads = process_reads(areads, read_d, bamfile)  # Deprecated

        # for aread, proper_map, overlap_reads in valid_reads:  # Deprecated
            # Only take soft-clips from outer regions of properly mapped reads, take all others
            if (aligned_read.cigar is None) or (len(aligned_read.cigar) <= 1):  # cigar is a list of tuples
                continue

        # if aligned_read.cigar and len(aligned_read.cigar) > 1:
            trim_coords = utils.trim_coords(aligned_read.qual, 3)  # Identify the read positions with qual > 2
            clip_coords = utils.get_clip_coords(aligned_read.qual, aligned_read.cigar)

            # Only keep reads that have a soft clip in sequence that has not been trimmed
            # due to low quality sequence.
            # if clip_coords[0] > trim_coords[0] or clip_coords[1] < trim_coords[1]:  # Deprecated
            if clip_coords[0] <= trim_coords[0] and clip_coords[1] >= trim_coords[1]:
                continue

            sc_seq = {'clipped':[], 'buffered':[]}
            new_clip_coords = [0, 0]
            start_coord, end_coord = clip_coords
            add_sc = [False, False]
            indel_only = False
            start_sc = start_coord > 0
            end_sc = end_coord < len(aligned_read.qual)
            seq = aligned_read.seq

            if start_sc and end_sc:
                add_sc = [True, True]
            else:
                if start_sc:
                    add_sc[0] = True
                    new_clip_coords = [0, start_coord]
                    if overlap_reads and aligned_read.is_reverse:
                        mate_seq = valid_reads[pair_indices[aligned_read.qname][int(aligned_read.is_read1)]][0].seq
                        add_sc[0] = self.check_pair_overlap(mate_seq, aligned_read, [0, start_coord], 'back')
                    if proper_map:
                        indel_only = aligned_read.is_reverse
                elif end_sc:
                    new_clip_coords = [end_coord, len(seq)]
                    add_sc[1] = True
                    if overlap_reads and not aligned_read.is_reverse:
                        mate_seq = valid_reads[pair_indices[aligned_read.qname][int(aligned_read.is_read1)]][0].seq
                        add_sc[1] = self.check_pair_overlap(mate_seq, aligned_read, [end_coord, len(seq)], 'front')
                    if proper_map:
                        indel_only = (indel_only and False) if aligned_read.is_reverse else (indel_only and True)
            final_add = add_sc[0] or add_sc[1]
            if add_sc[0]:
                sc_seq['buffered'].append(aligned_read.seq[0:(start_coord + kmer_size)])
                sc_seq['clipped'].append(aligned_read.seq[0:start_coord])
            if add_sc[1]:
                sc_seq['buffered'].append(seq[(end_coord - kmer_size):len(seq)])
                sc_seq['clipped'].append(seq[end_coord:len(seq)])
            if final_add:
                read_d['sv'][utils.get_seq_readname(aligned_read)] = (aligned_read, sc_seq, new_clip_coords, indel_only)
        # end for loop

        sv_fq = open(self.files['sv_fq'], 'w')
        sv_sc_fa = open(self.files['sv_sc_unmapped_fa'], 'w')

        for qname in read_d['unmapped_keep']:
            if qname in read_d['unmapped']:
                read = read_d['unmapped'][qname]
                read_d['sv'][utils.get_seq_readname(read)] = (read, None, None, False)
                sv_sc_fa.write(">" + read.qname + "\n" + str(read.seq) + "\n")

        if not self.sv_reads:
            self.sv_reads = {}
        self.sv_reads[sample_type] = {}
        for qname in read_d['sv']:
            aligned_read, sc_seq, clip_coords, indel_only = read_d['sv'][qname]
            self.sv_reads[sample_type][qname] = read_d['sv'][qname]
            if sample_type == 'sv':
                sv_bam.write(aligned_read)
            lout = utils.fq_line(aligned_read, indel_only, int(self.params.get_param('kmer_size')), True)
            if lout is not None:
                sv_fq.write(lout)
            if sc_seq:
                for clip_seq in sc_seq['buffered']:
                    sv_sc_fa.write(">" + qname + "\n" + clip_seq + "\n")
        self.disc_reads = {'disc':read_d['disc'], 'inv':read_d['inv_reads'], 'td':read_d['td_reads'], 'other':read_d['other']}
        sv_fq.close()
        sv_sc_fa.close()
        bamfile.close()

        if sample_type == 'sv':
            sv_bam.close()
            utils.log(self.logging_name, 'info', 'Sorting bam file %s to %s' % (self.files['sv_bam'], self.files['sv_bam_sorted']))
            pysam.sort("-o", self.files['sv_bam_sorted'], self.files['sv_bam'])
            utils.log(self.logging_name, 'info', 'Indexing sorted bam file %s' % self.files['sv_bam_sorted'])
            pysam.index(self.files['sv_bam_sorted'])

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
        if not self.cleaned_read_recs:
            self.cleaned_read_recs = {}
        self.cleaned_read_recs[sample_type] = None
        self.files['%s_cleaned_fq' % sample_type], self.cleaned_read_recs[sample_type] = utils.get_fastq_reads(self.files['%s_cleaned_fq' % sample_type], self.sv_reads[sample_type])
        self.sv_reads[sample_type] = None
        check = True
        if len(self.cleaned_read_recs[sample_type]) == 0:
            check = False

        utils.log(self.logging_name, 'info', 'Check there are cleaned reads %r' % check)
        return check

    def compare_kmers(self):

        '''
        '''

        self.kmers['ref'] = {}
        jellyfish = self.params.get_param('jellyfish')
        kmer_size = int(self.params.get_param('kmer_size'))

        for i in range(len(self.files['target_ref_fn'])):
            utils.log(self.logging_name, 'info', 'Indexing kmers for reference sequence %s' % self.files['target_ref_fn'][i])
            self.kmers['ref'] = utils.load_kmers(utils.run_jellyfish(self.files['target_ref_fn'][i], jellyfish, kmer_size), self.kmers['ref'])

        # if 'target_altref_fn' in self.files:
        #     for i in range(len(self.files['target_altref_fn'])):
        #         for j in range(len(self.files['target_altref_fn'][i])):
        #             utils.log(self.logging_name, 'info', 'Indexing kmers for reference sequence %s' % self.files['target_altref_fn'][i])
        #             self.kmers['ref'] = utils.load_kmers(utils.run_jellyfish(self.files['target_altref_fn'][i][j], jellyfish, kmer_size), self.kmers['ref'])

        utils.log(self.logging_name, 'info', 'Indexing kmers for sample sequence %s' % self.files['sv_cleaned_fq'])
        self.kmers['case'] = {}
        self.kmers['case'] = utils.load_kmers(utils.run_jellyfish(self.files['sv_cleaned_fq'], jellyfish, kmer_size), self.kmers['case'])
        self.kmers['case_sc'] = {}
        self.kmers['case_sc'] = utils.load_kmers(utils.run_jellyfish(self.files['sv_sc_unmapped_fa'], jellyfish, kmer_size), self.kmers['case_sc'])
        sc_mers = set(self.kmers['case'].keys()) & set(self.kmers['case_sc'])
        sample_only_mers = list(sc_mers.difference(set(self.kmers['ref'].keys())))

        if 'normal_bam_file' in self.params.opts:
            norm_kmers = {}
            norm_kmers = utils.load_kmers(utils.run_jellyfish(self.files['norm_cleaned_fq'], jellyfish, kmer_size), norm_kmers)
            sample_only_mers = set(sample_only_mers).difference(set(norm_kmers.keys()))

        sample_only_mers = list(sample_only_mers)

        # Write case only kmers out to file.
        self.files['sample_kmers'] = os.path.join(self.paths['kmers'], self.name + "_sample_kmers.out")
        sample_kmer_fout = open(self.files['sample_kmers'], 'w')

        self.kmers['case_only'] = {}
        for mer in sample_only_mers:
            sample_kmer_fout.write("\t".join([str(x) for x in [mer, str(self.kmers['case'][mer])]]) + "\n")
            self.kmers['case_only'][mer] = self.kmers['case'][mer]
        sample_kmer_fout.close()

        self.kmers['ref'] = {}
        self.kmers['case'] = {}
        self.kmers['case_sc'] = {}

        utils.log(self.logging_name, 'info', 'Writing %d sample-only kmers to file %s' % (len(self.kmers['case_only']), self.files['sample_kmers']))
        self.files['kmer_clusters'] = os.path.join(self.paths['kmers'], self.name + "_sample_kmers_merged.out")
        utils.log(self.logging_name, 'info', 'Writing kmer clusters to file %s' % self.files['kmer_clusters'])

        self.contigs = assembler.init_assembly(self.kmers['case_only'], self.cleaned_read_recs['sv'], kmer_size, int(self.params.get_param('trl_sr_thresh')), self.params.get_param('read_len'))
        self.cleaned_read_recs = None
        self.kmers['case_only'] = {}
        self.finalize_contigs()

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
        self.results = self.call_manager.resolve_sv_calls(self.contigs, self.files['target_ref_fn'][0], self.get_values(), self.disc_reads)
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
