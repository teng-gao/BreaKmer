#! /usr/bin/python
# -*- coding: utf-8 -*-

'''
BreaKmer processor.analysis module
'''

import os
import sys
import time
import math
import multiprocessing
from itertools import izip, repeat
import breakmer.processor.target as target
import breakmer.utils as utils
import pdb


__author__ = "Ryan Abo"
__copyright__ = "Copyright 2015, Ryan Abo"
__email__ = "ryanabo@gmail.com"
__license__ = "MIT"


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
        proc_pool = multiprocessing.Pool(nprocs)
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
        ref_data_dir = self.params.opts['reference_data_dir']
        jfish_path = self.params.opts['jellyfish']
        blat_path = self.params.opts['blat']
        ref_params = [mask_fn, ref_fa_fn, ref_data_dir, jfish_path, blat_path, self.params.get_param('kmer_size'), int(self.params.get_param('buffer_size'))]
        setup_params = izip(trgt_groups, repeat(ref_params))
        proc_pool.map(utils.setup_ref_data, setup_params)

    def create_targets(self):

        '''
        '''

        targets = self.params.targets.keys()
        targets.sort()
        for trgt_name in targets:
            self.targets[trgt_name] = target.TargetManager(self.params.targets[trgt_name], self.params)
        # return targets

    def run(self):

        '''
        '''

        start_time = time.clock()  # Track the run time.

        if self.params.fnc_cmd == 'profile_data':
            utils.log(self.logging_name, 'info', 'Profiling the coverage of the sample bam file.')
            utils.profile_data(self.params.get_param('sample_bam_file'))
            return

        self.create_targets()

        if self.params.fnc_cmd == 'prepare_reference_data':
            utils.log(self.logging_name, 'info', 'Creating all reference data.')
            self.preset_ref_data()
            print 'Reference data created!'
            return

        utils.start_blat_server(self.params)
        if self.params.fnc_cmd == 'start_blat_server':
            print 'Server started on host %s and port %s' % (self.params.get_param('blat_hostname'), self.params.get_param('blat_port'))
            return

        trgt_lst = self.params.targets.keys()
        trgt_lst.sort()
        for trgt_name in trgt_lst:
            trgt = self.targets[trgt_name]
            utils.log(self.logging_name, 'info', 'Analyzing %s' % trgt.name)

            # Write reference sequence fasta for gene if it doesn't exist.
            # if not self.params.opts['preset_ref_data']:
            trgt.set_ref_data()

            if not trgt.get_sv_reads():
                continue

            trgt.compare_kmers()  # Get reference and case kmers
            trgt.resolve_sv()  # Build contigs and blat them against the reference genome
            # self.summary_header, trgt_summary = trgt.get_summary()
            # self.summary[trgt.name] = trgt_summary

            # utils.log(self.logging_name, 'info', '%s summary\n%s\n%s' % (trgt.name, self.summary_header, trgt_summary))
            if trgt.has_results():
                trgt.write_results()
                self.results.extend(trgt.formatted_results)
            else:
                trgt.rm_output_dir()

        self.write_output()
        time_to_complete = time.clock() - start_time
        utils.log(self.logging_name, 'info', 'Analysis complete, %s' % str(time_to_complete))

        if not self.params.get_param('keep_blat_server'):
            # Stop gfServer
            utils.log(self.logging_name, 'info', 'Stopping blat server on port %d' % self.params.get_param('blat_port'))
            gfserver_stop_cmd = '%s stop localhost %d' % (self.params.get_param('gfserver'), self.params.get_param('blat_port'))
            os.system(gfserver_stop_cmd)
        print 'Analysis complete.'

    def write_output(self):

        '''
        '''

        res_fn = os.path.join(self.params.paths['output'], self.params.opts['analysis_name'] + "_svs.out")
        result_file = open(res_fn, 'w')
        header = "\t".join(['genes', 'target_breakpoints', 'mismatches', 'strands', 'total_matching', 'sv_type', 'sv_subtype', 'split_read_count', 'disc_read_count', 'breakpoint_coverages', 'contig_id', 'contig_seq']) + "\n"
        result_file.write(header)

        for result_str in self.results:
            utils.log(self.logging_name, 'info', 'Writing results to file: %s' % res_fn)
            result_file.write(result_str)
        result_file.close()

        # result_files = {}
        # for res in self.results:
        #     tag = res[6]
        #     if tag not in result_files:
        #         header = "\t".join(['genes', 'target_breakpoints', 'align_cigar', 'mismatches', 'strands', 'rep_overlap_segment_len', 'sv_type', 'split_read_count', 'nkmers', 'disc_read_count', 'breakpoint_coverages', 'contig_id', 'contig_seq']) + "\n"

        #         res_fn = os.path.join(self.params.paths['output'], self.params.opts['analysis_name'] + "_" + tag + "_svs.out")

        #         utils.log(self.logging_name, 'info', 'Writing %s output file %s' % (tag, res_fn))
        #         result_files[tag] = open(res_fn, 'w')
        #         if not self.params.opts['no_output_header']:
        #             result_files[tag].write(header)
        #         result_files[tag].write("\t".join([str(x) for x in res]) + "\n")

        # for result_filename in result_files:
        #     result_files[result_filename].close()

        # summary_fn = os.path.join(self.params.paths['output'], self.params.opts['analysis_name'] + "_summary.out")
        # summary_f = open(summary_fn, 'w')
        # utils.log(self.logging_name, 'info', 'Writing summary file to %s' % summary_fn)
        # summary_f.write(self.summary_header+"\n")

        # keys = self.summary.keys()
        # keys.sort()
        # for gene in keys:
        #     summary_f.write(self.summary[gene] + "\n")
        # summary_f.close()
