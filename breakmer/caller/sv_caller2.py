#! /usr/bin/local/python
# -*- coding: utf-8 -*-


'''
BreaKmer sv_caller module
'''


import breakmer.utils as utils
import breakmer.realignment.realign as realigner
import breakmer.results.results as results
import pdb


__author__ = "Ryan Abo"
__copyright__ = "Copyright 2015, Ryan Abo"
__email__ = "ryanabo@gmail.com"
__license__ = "MIT"



class SVCallManager(object):

    '''
    '''

    def __init__(self, params):
        self.params = params
        self.logging_name = 'breakmer.caller.sv_caller'
        self.align_manager = realigner.RealignManager(params)
        self.filter_manager = FilterManager(params)

    def resolve_sv_calls(self, contigs, target_ref_fn, target_region_values, disc_reads):

        '''
        '''
        sv_results = []
        for assembled_contig in contigs:
            utils.log(self.logging_name, 'info', 'Assessing contig %s' % assembled_contig.seq.value)
            # pdb.set_trace()
            realignment_set = self.align_manager.realignment(assembled_contig, target_ref_fn, target_region_values)
            sv_result = results.SVResult(self.make_call(assembled_contig, target_region_values, realignment_set), self.params.get_param('sample_bam_file'), assembled_contig, target_region_values, disc_reads)
            if not sv_result.filter:
                self.filter_manager.filter_result(sv_result)
            if not sv_result.filter:
                sv_results.append(sv_result)
        return sv_results

    def make_call(self, contig, region_values, realignment_result_set):

        '''
        '''

        sv_event = None
        if not realignment_result_set.has_results:
            utils.log(self.logging_name, 'info', 'No blat results file %s, no calls for %s.' % (self.align_manager.query_res_fn, contig.contig_id))
        else:
            utils.log(self.logging_name, 'info', 'Making variant calls from blat results %s' % self.align_manager.query_res_fn)
            if realignment_result_set.has_indel():
                sv_event = realignment_result_set.sv_event  # get_indel_result()
            elif realignment_result_set.check_svs():
                sv_event = realignment_result_set.sv_event  #bm.get_svs_result()
        return sv_event

# class SVCaller(object):

#     '''

#     '''

#     def __init__(self):
#         self.logging_name = 'breakmer.caller.sv_caller'

#     def get_result(self, realignment_result_set):

#         '''
#         '''

#         contig_id = None
#         if not realignment_result_set.has_results:
#             utils.log(self.logging_name, 'info', 'No blat results file %s, no calls for %s.' % (self.query_res_fn, contig_id))
#         else:
#             print 'checking results'
#             utils.log(self.logging_name, 'info', 'Making variant calls from blat results %s' % self.query_res_fn)
#             if realignment_result_set.has_indel():
#                 sv_event = realignment_result_set.get_indel_result()
#             elif realignment_result_set.check_svs():
#                 sv_event = self.bm.get_svs_result()
#         return self.result

class FilterManager(object):

    '''
    '''

    def __init__(self, params):
        self.params = params
        self.logging_name = 'breakmer.caller.sv_caller'

    def filter_result(self, sv_result):

        '''
        '''

        filter_result = False
        if sv_result.sv_type == 'indel':
            filter_result = self.filter_indel(sv_result)
        else:
            filter_result = self.filter_svs(sv_result)
        sv_result.filter = filter_result

    def filter_indel(self, sv_result):

        '''
        '''

        indel_size_thresh = int(self.params.get_param('indel_size'))
        indel_segment = sv_result.sv_event.realignments[0]
        utils.log(self.logging_name, 'info', 'Blat result spans query (%r) or only one blat result (%r) and blat result in target (%r)' % (indel_segment.spans_query(), (len(sv_result.sv_event.realignments) == 1), indel_segment.in_target))

        keep_result = indel_segment.valid and (indel_segment.mean_cov < 2) and indel_segment.in_target and (indel_segment.indel_maxevent_size[0] >= indel_size_thresh) # and (not indel_segment.rep_man.breakpoint_in_rep[0] and not indel_segment.rep_man.breakpoint_in_rep[1])
        utils.log(self.logging_name, 'debug', 'Keep blat result %r' % keep_result)

        filter_result = False
        if not keep_result:
            filter_result = True
            utils.log(self.logging_name, 'debug', 'Indel failed checking criteria: in annotated gene: %r, mean query coverage < 2: %r, in target: %r, in repeat: %r, indel size < %d: %r' % (indel_segment.valid, indel_segment.mean_cov, indel_segment.in_target, False, indel_size_thresh, indel_segment.indel_maxevent_size[0] < indel_size_thresh))
            return filter_result

        flank_match_thresh = True
        for flank_match in indel_segment.indel_flank_match:
            fm_perc = round((float(flank_match)/float(indel_segment.get_size('query')))*100, 2)
            if fm_perc < 10.0:
                flank_match_thresh = False
            utils.log(self.logging_name, 'info', 'Indel result has matching flanking sequence of largest indel event of %d (%d%% of query)'%(flank_match, fm_perc))

        utils.log(self.logging_name, 'info', 'Indel result has matching flanking sequence of largest indel event (10 %% of query) on both sides (%r)' % flank_match_thresh)
        in_ff, span_ff = utils.filter_by_feature(indel_segment.get_brkpt_locs(), sv_result.query_region, self.params.get_param('keep_intron_vars'))

        brkpt_cov = [sv_result.contig.get_contig_counts().get_counts(x, x, 'indel') for x in indel_segment.query_brkpts]
        low_cov = min(brkpt_cov) < self.params.get_param('indel_sr_thresh')
        if not in_ff and not low_cov and flank_match_thresh:
            utils.log(self.logging_name, 'debug', 'Top hit contains whole query sequence, indel variant')
        else:
            filter_result = True
            utils.log(self.logging_name, 'debug', 'Indel in intron (%r) or low coverage at breakpoints (%r) or minimum segment size < 20 (%r), filtering out.' % (in_ff, low_cov, min(indel_segment.query_blocksizes)) )
        return filter_result


    def filter_svs(self, sv_result):


        '''
        '''

        filter_result = False
        nmissing_query_cov = len(filter(lambda y: y, map(lambda x: x==0, sv_result.query_cov)))
        if nmissing_query_cov >= self.params.get_param('trl_minseg_len'):
            filter_result = True

        if sv_result.values['sv_subtype'] == 'trl':
            if max(sv_result.contig.get_contig_counts().others) >= self.params.get_param('trl_sr_thresh'):
                filter_result = self.filter_trl(sv_result, 0.0)
        # if not self.multiple_genes(brkpts['chrs'], brkpts['r'], res_values['anno_genes']):
        #     # brkpt_counts, brkpt_kmers, brkpt_rep_filt = self.get_brkpt_counts_filt(brkpts, 'rearr')

        #     rearr_type, disc_read_support = self.define_rearr(brkpts['r'], res_values['strands'], brkpts['tcoords'], disc_reads)
        #     # if not self.filter_rearr(query_region, params, brkpts['r'], brkpt_counts, brkpt_kmers, rearr_type, disc_read_support):
        else:
            filter_result = self.filter_rearr(sv_result)
        #         # result = self.format_result(res_values)
        # elif max(self.contig_rcounts.others) >= params.get_param('trl_sr_thresh'):
        # else:
            # if max(sv_result.contig.get_contig_counts().others) >= self.params.get_param('trl_sr_thresh'):
        #     brkpt_counts, brkpt_kmers, brkpt_rep_filt = self.get_brkpt_counts_filt(brkpts, 'trl')
        #     disc_read_count = self.check_disc_reads(brkpts['t'], query_region, disc_reads['disc'])
                # filter_result = self.filter_trl(valid_rearrangement, query_region, params, brkpt_counts, brkpt_kmers, disc_read_count, res_values['anno_genes'], max_repeat, brkpt_rep_filt)

    def filter_rearr(self, sv_result): # query_region, params, brkpts, brkpt_counts, brkpt_kmers, rearr_type, disc_read_count):

        '''
        '''

        in_ff, span_ff = utils.filter_by_feature(sv_result.breakpoint_values['ref_pos'], sv_result.query_region, self.params.get_param('keep_intron_vars'))
        match_sorted_realignments = sorted(sv_result.sv_event.realignments, key=lambda x: x.get_nmatch_total())
        top_realigned_segment = match_sorted_realignments[0]

        check1 = (min(sv_result.breakpoint_values['counts']['n']) < self.params.get_param('rearr_sr_thresh'))
        check2 = top_realigned_segment.get_nmatch_total() < self.params.get_param('rearr_minseg_len')
        check3 = (in_ff and span_ff)
        check4 = (sv_result.values['disc_read_count'] < 1)
        check5 = (sv_result.values['sv_subtype'] == 'NA')
        check6 = (min(sv_result.breakpoint_values['kmers']) == 0)
        filter_result = check1 or check2 or check3 or check4 or check5 or check6
        utils.log(self.logging_name, 'info', 'Check filter for rearrangement')
        utils.log(self.logging_name, 'info', 'Filter by feature for being in exon (%r) or spanning exon (%r)' % (in_ff, span_ff))
        utils.log(self.logging_name, 'info', 'Split read threshold %d, breakpoint read counts %d' % (min(sv_result.breakpoint_values['counts']['n']), self.params.get_param('rearr_minseg_len')))
        utils.log(self.logging_name, 'info', 'Minimum segment length observed (%d) meets threshold (%d)' % (top_realigned_segment.get_nmatch_total(), self.params.get_param('rearr_minseg_len')))
        utils.log(self.logging_name, 'info', 'Minimum discordant read pairs for rearrangement (%d)' % (sv_result.values['disc_read_count']))
        return filter_result

    def filter_trl(self, sv_result, max_repeat): # br_valid, query_region, params, brkpt_counts, brkpt_kmers, disc_read_count, anno_genes, max_repeat, rep_filt):

        '''
        '''

        filter_result = sv_result.valid_rearrangement[1] or (max(sv_result.breakpoint_values['counts']['d']) < self.params.get_param('trl_sr_thresh')) #or not br_valid[0]
        utils.log(self.logging_name, 'debug', 'Check translocation filter')
        utils.log(self.logging_name, 'debug', 'All blat result segments are within annotated or pre-specified regions %r' % sv_result.valid_rearrangement[0])
        utils.log(self.logging_name, 'debug', 'All blat result segments are within simple repeat regions that cover > 75.0 percent of the segment %r' % sv_result.valid_rearrangement[1])
        utils.log(self.logging_name, 'debug', 'The maximum read count support around breakpoints %d meets split read threshold %d' % (max(sv_result.breakpoint_values['counts']['d']), self.params.get_param('trl_sr_thresh')))
        utils.log(self.logging_name, 'debug', 'The minimum number of kmers at breakpoints %d' % min(sv_result.breakpoint_values['kmers']))
        utils.log(self.logging_name, 'debug', 'The maximum repeat overlap by a blat result: %f' % max_repeat)

        if not filter_result:
            utils.log(self.logging_name, 'debug', 'Filter %r, checking discordant read counts %d' % (filter_result, sv_result.values['disc_read_count']))
            if sv_result.values['disc_read_count'] < 2:
                match_sorted_realignments = sorted(sv_result.sv_event.realignments, key=lambda x: x.get_nmatch_total())
                top_realigned_segment = match_sorted_realignments[0]
                if (top_realigned_segment.get_nmatch_total() < self.params.get_param('trl_min_seg_len')) or (min(sv_result.breakpoint_values['counts']['n']) < self.params.get_param('trl_sr_thresh')) or (min(sv_result.breakpoint_values['kmers']) == 0) or sv_result.breakpoint_values['rep_filter']:
                    utils.log(self.logging_name, 'debug', 'Shortest segment is < %d bp with %d discordant reads. Filtering.' % (self.params.get_param('trl_minseg_len'), sv_result.values['disc_read_count']))
                    utils.log(self.logging_name, 'debug', 'The minimum read count support for breakpoints %d meets split read threshold %d'%(min(sv_result.breakpoint_values['counts']['n']), self.params.get_param('trl_sr_thresh')))
                    utils.log(self.logging_name, 'debug', 'The minimum number of kmers at breakpoints %d' % min(sv_result.breakpoint_values['kmers']))
                    filter_result = True
                elif sv_result.values['disc_read_count'] == 0:
                    # Check a number of metrics for shortest blat segment
                    br_qs = top_realigned_segment.qstart()
                    br_qe = top_realigned_segment.qend()
                    low_complexity = self.minseq_complexity(sv_result.contig.seq.value[br_qs:br_qe],3) < 25.0 # Complexity of blat segment
                    missing_qcov = self.missing_query_coverage() > 5.0
                    short = top_realigned_segment.get_nmatch_total() <= round(float(len(sv_result.contig.seq.value))/float(4.0))
                    utils.log(self.logging_name, 'debug', 'Checking length of shortest sequence, considered too short %r, %d, %f' % (short, top_realigned_segment.get_nmatch_total(), round(float(len(sv_result.contig.seq.value))/float(4.0))) )
                    overlap = max(top_realigned_segment.seg_overlap) > 5
                    gaps_exist = max(top_realigned_segment.gaps['query'][0], top_realigned_segment.gaps['hit'][0]) > 0
                    low_uniqueness = self.check_uniqueness(match_sorted_realignments)
                    intergenic_regions = 'intergenic' in sv_result.values['anno_genes']
                    read_strand_bias = self.check_read_strands(sv_result)
                    check_values = [low_complexity, missing_qcov, short, overlap, gaps_exist, low_uniqueness, read_strand_bias, intergenic_regions]
                    utils.log(self.logging_name, 'debug', 'Discordant read count of 0 checks %s' % (",".join([str(x) for x in check_values])))
                    num_checks = 0
                    for check in check_values:
                        if check:
                            num_checks += 1
                    if num_checks > 1:
                        utils.log(self.logging_name, 'info', 'Two or more filter checks, setting filtering to true for contig')
                        filter_result = True
        return filter_result

    def missing_query_coverage(self, sv_result):

        '''
        '''

        missing_cov = 0
        for i in sv_result.query_cov:
          if i == 0:
            missing_cov += 1
          else:
            break

        for i in reversed(sv_result.query_cov):
          if i == 0:
            missing_cov += 1
          else:
            break

        perc_missing = round((float(missing_cov)/float(len(sv_result.contig.seq.value)))*100, 4)
        utils.log(self.logging_name, 'debug', 'Calculated %f missing coverage of blat query sequence at beginning and end' % perc_missing)
        return perc_missing

    def check_read_strands(self, sv_result):

        '''
        '''

        same_strand = False
        strands = []
        for read in sv_result.contig.reads:
            strand = read.id.split("/")[1]
            strands.append(strand)
        if len(set(strands)) == 1:
            same_strand = True
        utils.log(self.logging_name, 'debug', 'Checking read strands for contig reads %s' % (",".join([read.id for read in sv_result.contig.reads])))
        utils.log(self.logging_name, 'debug', 'Reads are on same strand: %r' % same_strand)
        return same_strand


    def check_uniqueness(self, realigned_segments):

        '''
        '''

        low_unique = False
        for br_vals in realigned_segments:
            if not br_vals[0].in_target:
                if br_vals[0].mean_cov > 4:
                    low_unique = True
            else:
                if br_vals[0].mean_cov > 10:
                    low_unique = True
        return low_unique

    def minseq_complexity(self, seq, N):

        '''
        '''

        utils.log(self.logging_name, 'debug', 'Checking sequence complexity of blat result segment %s using %d-mers' % (seq, N))
        nmers = {}
        total_possible = len(seq) - 2
        for i in range(len(seq) - (N - 1)):
            nmers[str(seq[i:i+N]).upper()] = True
        complexity = round((float(len(nmers))/float(total_possible))*100,4)
        utils.log(self.logging_name, 'debug', 'Complexity measure %f, based on %d unique %d-mers observed out of a total of %d %d-mers possible' % (complexity, len(nmers), N, total_possible, N))
        return complexity
