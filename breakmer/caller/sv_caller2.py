#! /usr/bin/local/python
# -*- coding: utf-8 -*-


'''
BreaKmer sv_caller module
'''


import breakmer.utils as utils
import breakmer.realignment.realign as realigner
import breakmer.results.results as results


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
            realignment_set = self.align_manager.realignment(assembled_contig, target_ref_fn, target_region_values)
            sv_result = results.SVResult(self.make_call(assembled_contig, target_region_values, realignment_set), assembled_contig, target_region_values, disc_reads)
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

        if sv_result.sv_type == 'indel':
            self.filter_indel(sv_result)
        else:
            self.filter_svs(sv_result)

    def filter_indel(self, sv_result):

        '''
        '''

        x = 1
        # utils.log(self.logging_name, 'info', 'Blat result spans query (%r) or only one blat result (%r) and blat result in target (%r)'%(br.spans_query(), (len(self.realignments) == 1), br.in_target))
        # keep_br = br.valid and br.mean_cov < 2 and br.in_target and (br.indel_maxevent_size[0] >= indel_size_thresh) and (not br.rep_man.breakpoint_in_rep[0] and not br.rep_man.breakpoint_in_rep[1])
        # utils.log(self.logging_name, 'debug', 'Keep blat result %r' % keep_br)

        # if keep_br:
        #     brkpt_cov = [self.meta_dict['contig_vals'][1].get_counts(x, x, 'indel') for x in br.query_brkpts]
        #     low_cov = min(brkpt_cov) < self.params.get_param('indel_sr_thresh')
        #     flank_match_thresh = True

        #     for flank_match in br.indel_flank_match:
        #         fm_perc = round((float(flank_match)/float(br.get_size('query')))*100, 2)
        #         if fm_perc < 10.0: 
        #             flank_match_thresh = False
        #         utils.log(self.logging_name, 'info', 'Indel result has matching flanking sequence of largest indel event of %d (%d of query)'%(flank_match, fm_perc))

        #     utils.log(self.logging_name, 'info', 'Indel result has matching flanking sequence of largest indel event (10 perc of query) on both sides (%r)' % flank_match_thresh)
        #     in_ff, span_ff = utils.filter_by_feature(br.get_brkpt_locs(), self.meta_dict['query_region'], self.meta_dict['params'].opts['keep_intron_vars'])

        #     if not in_ff and not low_cov and flank_match_thresh:
        #         self.sv_event = results.SVEvent(br) #, self.meta_dict['query_region'], self.meta_dict['contig_vals'], self.meta_dict['sbam'])
        #         utils.log(self.logging_name, 'debug', 'Top hit contains whole query sequence, indel variant')
        #     else: 
        #         utils.log(self.logging_name, 'debug', 'Indel in intron (%r) or low coverage at breakpoints (%r) or minimum segment size < 20 (%r), filtering out.' % (in_ff, low_cov, min(br.query_blocksizes)) )
        # else: 
        #     utils.log(self.logging_name, 'debug', 'Indel failed checking criteria: in annotated gene: %r, mean query coverage < 2: %r, in target: %r, in repeat: %r, indel size < %d: %r' % (br.valid, br.mean_cov, br.in_target, ",".join([str(x) for x in br.rep_man.breakpoint_in_rep]), indel_size_thresh, br.indel_maxevent_size[0] < indel_size_thresh))

    def filter_svs(self, sv_result):


        '''
        '''

        x = 1
        # nmissing_query_cov = len(filter(lambda y: y, map(lambda x: x==0, self.sv_event.query_cov)))
            # if nmissing_query_cov < self.params.get_param('trl_minseg_len'): 
            #     valid = True

        # if not self.multiple_genes(brkpts['chrs'], brkpts['r'], res_values['anno_genes']):
        #     # brkpt_counts, brkpt_kmers, brkpt_rep_filt = self.get_brkpt_counts_filt(brkpts, 'rearr')

        #     rearr_type, disc_read_support = self.define_rearr(brkpts['r'], res_values['strands'], brkpts['tcoords'], disc_reads)
        #     # if not self.filter_rearr(query_region, params, brkpts['r'], brkpt_counts, brkpt_kmers, rearr_type, disc_read_support):

        #         # result = self.format_result(res_values)
        # elif max(self.contig_rcounts.others) >= params.get_param('trl_sr_thresh'):
        #     brkpt_counts, brkpt_kmers, brkpt_rep_filt = self.get_brkpt_counts_filt(brkpts, 'trl')
        #     disc_read_count = self.check_disc_reads(brkpts['t'], query_region, disc_reads['disc'])
        #     if not self.filter_trl(valid_rearrangement, query_region, params, brkpt_counts, brkpt_kmers, disc_read_count, res_values['anno_genes'], max_repeat, brkpt_rep_filt):
        #         print 'hello'

    #   def get_brkpt_counts_filt(self, brkpts, sv_type):
    # #    print 'Contig seq', self.contig_seq
    # #    print 'Breakpoint simple repeat filter', brkpts['f'] 
    #     avg_comp, comp_vec = utils.calc_contig_complexity(self.contig_seq) 
    # #    print 'Contig avg complexity', avg_comp
    # #    print 'Contig complexity vec', comp_vec
    #     brkpt_rep_filt = False
    #     brkpt_counts = {'n':[],'d':[],'b':[]}
    #     brkpt_kmers = []
    #     for qb in brkpts['q'][1]:
    #       left_idx = qb[0] - min(qb[1],5)
    #       right_idx = qb[0] + min(qb[2],5)
    # #      print self.contig_rcounts.others
    # #      print self.contig_rcounts.indel_only
    # #      print qb[0], left_idx, right_idx
    #       bc = self.contig_rcounts.get_counts(left_idx, right_idx, sv_type)
    #       brkpt_counts['n'].append(min(bc))
    #       brkpt_counts['d'].append(min(self.contig_rcounts.get_counts((qb[0]-1), (qb[0]+1), sv_type)))
    # #      print 'Others counts', self.contig_rcounts.others, qb[0]
    # #      print 'Indel only counts', self.contig_rcounts.indel_only, qb[0]
    #       brkpt_counts['b'].append(self.contig_rcounts.get_counts(qb[0], qb[0], sv_type))
    #       brkpt_kmers.append(self.contig_kmer_locs[qb[0]])
    # #      print 'Breakpoint in contig', qb[0]
    #       brkpt_rep_filt = brkpt_rep_filt or (comp_vec[qb[0]] < (avg_comp/2))
    # #      print 'Breakpoint rep filter', brkpt_rep_filt, comp_vec[qb[0]]
    #       utils.log(self.logging_name, 'debug', 'Read count around breakpoint %d: %s'%(qb[0],",".join([str(x) for x in bc])))
    #     utils.log(self.logging_name, 'debug', 'Kmer count around breakpoints %s'%(",".join([str(x) for x in brkpt_kmers])))
    #     brkpt_rep_filt = brkpt_rep_filt or (len(filter(lambda x: x, brkpts['f'])) > 0)
    #     return brkpt_counts, brkpt_kmers, brkpt_rep_filt

    #   def filter_rearr(self, query_region, params, brkpts, brkpt_counts, brkpt_kmers, rearr_type, disc_read_count):
    #     print 'Breakpoint counts', brkpt_counts
    #     print params.opts
    #     in_ff, span_ff = utils.filter_by_feature(brkpts, query_region, params.get_param('keep_intron_vars'))
    #     filter = (min(brkpt_counts['n']) < params.get_param('rearr_sr_thresh')) or self.br_sorted[0][1] < params.get_param('rearr_minseg_len') or (in_ff and span_ff) or (disc_read_count < 1) or (rearr_type == 'rearrangement') or (min(brkpt_kmers) == 0)
    #     utils.log(self.logging_name, 'info' ,'Check filter for rearrangement')
    #     utils.log(self.logging_name, 'info', 'Filter by feature for being in exon (%r) or spanning exon (%r)'%(in_ff, span_ff))
    #     utils.log(self.logging_name, 'info', 'Split read threshold %d, breakpoint read counts %d'%(min(brkpt_counts['n']), params.get_param('rearr_minseg_len')))
    #     utils.log(self.logging_name, 'info', 'Minimum segment length observed (%d) meets threshold (%d)'%(self.br_sorted[0][1], params.get_param('rearr_minseg_len')))
    #     utils.log(self.logging_name, 'info', 'Minimum discordant read pairs for rearrangement (%d)'%(disc_read_count))
    #     return filter

    #   def filter_trl(self, br_valid, query_region, params, brkpt_counts, brkpt_kmers, disc_read_count, anno_genes, max_repeat, rep_filt):
    #     filter = br_valid[1] or (max(brkpt_counts['d']) < params.get_param('trl_sr_thresh')) #or not br_valid[0]
    #     utils.log(self.logging_name, 'debug', 'Check translocation filter')
    #     utils.log(self.logging_name, 'debug', 'All blat result segments are within annotated or pre-specified regions %r' % br_valid[0])
    #     utils.log(self.logging_name, 'debug', 'All blat result segments are within simple repeat regions that cover > 75.0 percent of the segment %r' % br_valid[1])
    #     utils.log(self.logging_name, 'debug', 'The maximum read count support around breakpoints %d meets split read threshold %d' % (max(brkpt_counts['d']), params.get_param('trl_sr_thresh')))
    #     utils.log(self.logging_name, 'debug', 'The minimum number of kmers at breakpoints %d' % min(brkpt_kmers))
    #     utils.log(self.logging_name, 'debug', 'The maximum repeat overlap by a blat result: %f' % max_repeat)
    #     if not filter:
    #       utils.log(self.logging_name, 'debug', 'Filter %r, checking discordant read counts %d'%(filter, disc_read_count)) 
    #       if disc_read_count < 2:
    # #        print 'Filter due to repeat', rep_filt
    #         if (self.br_sorted[0][1] < params.get_param('trl_min_seg_len')) or (min(brkpt_counts['n']) < params.get_param('trl_sr_thresh')) or (min(brkpt_kmers)==0) or rep_filt:
    #           utils.log(self.logging_name, 'debug', 'Shortest segment is < %d bp with %d discordant reads. Filtering.'%(params.get_param('trl_minseg_len'), disc_read_count))
    #           utils.log(self.logging_name, 'debug', 'The minimum read count support for breakpoints %d meets split read threshold %d'%(min(brkpt_counts['n']),params.get_param('trl_sr_thresh')))
    #           utils.log(self.logging_name, 'debug', 'The minimum number of kmers at breakpoints %d' % min(brkpt_kmers))
    #           filter = True
    #         elif disc_read_count == 0: 
    #           # Check a number of metrics for shortest blat segment
    #           br_qs = self.br_sorted[0][0].qstart()
    #           br_qe = self.br_sorted[0][0].qend()
    #           low_complexity = self.minseq_complexity(self.contig_seq[br_qs:br_qe],3) < 25.0 # Complexity of blat segment
    #           missing_qcov = self.missing_query_coverage() > 5.0
    #           short = self.br_sorted[0][1] <= round(float(len(self.contig_seq))/float(4.0))
    #           utils.log(self.logging_name, 'debug', 'Checking length of shortest sequence, considered too short %r, %d, %f'%(short, self.br_sorted[0][1], round(float(len(self.contig_seq))/float(4.0))) )
    #           overlap = max(self.br_sorted[0][0].seg_overlap) > 5
    #           gaps_exist = max(self.br_sorted[0][0].gaps['query'][0], self.br_sorted[0][0].gaps['hit'][0]) > 0
    #           low_uniqueness = self.check_uniqueness()
    #           intergenic_regions = 'intergenic' in anno_genes
    #           read_strand_bias = self.check_read_strands()
    #           check_values = [low_complexity, missing_qcov, short, overlap, gaps_exist, low_uniqueness, read_strand_bias, intergenic_regions]
    #           utils.log(self.logging_name, 'debug', 'Discordant read count of 0 checks %s' % (",".join([str(x) for x in check_values])))
    #           num_checks = 0
    #           for check in check_values:
    #             if check: 
    #               num_checks += 1
    #           if num_checks > 1: 
    #             utils.log(self.logging_name, 'info', 'Two or more filter checks, setting filtering to true for contig')
    #             filter = True
    #     return filter