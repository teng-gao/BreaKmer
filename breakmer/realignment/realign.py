#! /usr/bin/local/python
# -*- coding: utf-8 -*-

'''
BreaKmer realign module
'''

import os
import math
import breakmer.utils as utils
import breakmer.results.results as results
import pdb


__author__ = "Ryan Abo"
__copyright__ = "Copyright 2015, Ryan Abo"
__email__ = "ryanabo@gmail.com"
__license__ = "MIT"


class RealignedSegment(object):

    '''
    '''

    def __init__(self, scope, params, query_region, realignment_values):
        self.logging_name = 'breakmer.realignment.realign'
        self.values = realignment_values
        self.matches = {}
        self.gaps = {}
        self.vals = {}
        self.fragments = {}
        self.strand = ''
        self.genes = ''
        self.in_target = False
        # self.rep_man = blat_repeat_manager()
        self.in_repeat = False
        self.repeat_overlap = 0.0
        self.repeat_coords = None
        self.filter_reps_edges = [False, False]
        self.valid = True
        self.interval = None
        self.breakpts = []
        self.query_brkpts = []
        self.query_blocksizes = []
        self.indel_sizes = []
        self.indel_maxevent_size = [0, '']
        self.mean_cov = 0.0
        self.percent_identity = 0.0
        self.seg_overlap = [0, 0]
        self.cigar = ''
        self.indel_flank_match = [0, 0]
        self.ngaps = 0
        self.score = 0
        self.parse_result(scope, params, query_region)
        # self.report_blat_hit()

    def parse_result(self, scope, params, query_region):

        '''
        '''
        # pdb.set_trace()
        self.matches['match'] = int(self.values[0])
        self.matches['mis'] = int(self.values[1])
        self.matches['rep'] = int(self.values[2])
        self.gaps['hit'] = [int(self.values[6]), int(self.values[7])]
        self.gaps['query'] = [int(self.values[4]), int(self.values[5])]

        tname = self.values[13].replace('chr', '')
        toffset = 0
        if scope == 'target':
            self.values[13] = query_region[0].replace('chr', '')
            toffset = query_region[1] - int(params.get_param('buffer_size')) # 45003745 - 100
        tname = self.values[13]

        tcoords = [toffset + int(self.values[15]), toffset + int(self.values[16])] # 45003745 - 100 + q_start/q_end of alignment to B2M ref seq
        self.values[15] = tcoords[0]
        self.values[16] = tcoords[1]

        self.vals['hit'] = {'name': tname,
                            'size': int(self.values[14]),
                            'coords': tcoords
                           }
        self.vals['query'] = {'name': self.values[9],
                              'size':int(self.values[10]),
                              'coords': [int(self.values[11]), int(self.values[12])]
                              }

        self.strand = self.values[8]

        self.query_blocksizes = [int(x) for x in self.values[18].rstrip(",").split(",")]
        tstarts = [toffset + int(x) for x in self.values[20].rstrip(",").split(",")]
        self.values[20] = ",".join([str(x) for x in tstarts]) + ","
        qstarts = [int(x) for x in self.values[19].rstrip(",").split(",")]
        self.fragments['hit'] = []
        self.fragments['query'] = []
        for qstart, tstart, blocksize in zip(qstarts, tstarts, self.query_blocksizes):
            self.fragments['hit'].append((tstart, tstart + blocksize))
            self.fragments['query'].append((qstart, qstart + blocksize))
        self.fragments['count'] = len(self.query_blocksizes)

        self.set_gene_anno(params.gene_annotations, query_region)
        #   if 'repeat_mask' in d:
        #     self.set_repeat(d['repeat_mask'],d['params'].repeat_mask)

        self.set_indel_locs()
        self.percent_identity = 100.0 - self.calcMilliBad()

        self.ngaps = self.get_ngap_total()
        self.score = self.get_nmatch_total() + (float(self.get_nmatch_total())/float(self.get_size('query')))


    def set_segment_overlap(self, right, left):

        '''
        '''

        self.seg_overlap = [left, right]

    # def set_repeat(self, target_rep_mask, all_rep_mask):

    #     '''
    #     '''

    #     self.rep_man = blat_repeat_manager()
    #     if self.matches['rep'] > 0:
    #         self.in_repeat = True
    #     if target_rep_mask and all_rep_mask:
    #         # Check rep_mask if it exists.
    #         rmask = target_rep_mask
    #         if not self.in_target:
    #             rmask = None
    #             if self.vals['hit']['name'] in all_rep_mask:
    #                 rmask = all_rep_mask[self.vals['hit']['name']]
    #         if rmask:
    #             self.rep_man.setup(self.get_coords('hit'), rmask)
    #             self.in_repeat, self.repeat_overlap, self.repeat_coords, self.filter_reps_edges = self.rep_man.other_values

    def has_gaps(self):
        return min(self.gaps['hit']) > 0 or min(self.gaps['query']) > 0

    def get_coords(self, type):
        return self.vals[type]['coords']

    def qstart(self):
        return self.get_coords('query')[0]

    def qend(self):
        return self.get_coords('query')[1]

    def tstart(self):
        return self.get_coords('hit')[0]

    def tend(self):
        return self.get_coords('hit')[1]

    def get_name(self, type):
        return self.vals[type]['name']

    def get_size(self, type):
        return self.vals[type]['size']

    def get_query_span(self):
        return self.get_coords('query')[1] - self.get_coords('query')[0]

    def get_query_coverage(self):
        return round((float(self.get_query_span())/float(self.get_size('query')))*100,2)

    def spans_query(self):
        return self.get_size('query') == (self.get_coords('query')[1]- self.get_coords('query')[0])

    def get_ngap_total(self):
        return self.gaps['hit'][1] + self.gaps['query'][1]

    def get_num_gaps(self):
        return self.gaps['hit'][0] + self.gaps['query'][0]

    def get_nmatch_total(self):
        return self.matches['match'] + self.matches['rep']

    def get_nmatches(self, type):
        return self.matches[type]

    def sum_indel_flank_matches(self, flank_str):

        '''
        '''

        m_indxs = []
        match_sum = 0
        for i in range(len(flank_str)):
            if flank_str[i] == "M":
                m_indxs.append(i)
        for indx in m_indxs:
            nmatch = ''
            windx = indx-1
            while windx > -1 and utils.is_number(flank_str[windx]):
                nmatch = flank_str[windx] + nmatch
                windx = windx - 1
            match_sum += int(nmatch)
        return match_sum

    def set_indel_flank_matches(self):

        '''
        '''

        if self.indel_maxevent_size[0] > 0:
            csplit = self.cigar.split(str(self.indel_maxevent_size[0]) + self.indel_maxevent_size[1])
            lflank = csplit[0]
            self.indel_flank_match[0] += self.sum_indel_flank_matches(lflank)
            rflank = csplit[-1]
            self.indel_flank_match[1] += self.sum_indel_flank_matches(rflank)

    def set_indel_locs(self):

        '''
        '''

        # chrm = 'chr'+str(self.get_name('hit'))
        for i in range(self.fragments['count']-1):
            if i==0 and self.fragments['query'][i][0] > 0:
                self.cigar = str(self.fragments['query'][i][0]) + "S"
            qend1 = int(self.fragments['query'][i][1])
            qstart2 = int(self.fragments['query'][i+1][0])
            tend1 = int(self.fragments['hit'][i][1])
            tstart2 = int(self.fragments['hit'][i+1][0])
            ins_bp = qstart2 - qend1
            del_bp = tstart2 - tend1
            bp1 = tend1
            bp2 = tstart2
            self.cigar += str(self.query_blocksizes[i]) + "M"
            if ins_bp > 0:
                self.breakpts.append([bp1])
                self.indel_sizes.append("I"+str(ins_bp))
                self.add_query_brkpt(qend1)
                self.add_query_brkpt(qstart2)
                self.cigar += str(ins_bp) + "I"
                if ins_bp > self.indel_maxevent_size[0]:
                    self.indel_maxevent_size = [ins_bp, "I"]
            if del_bp > 0:
                self.breakpts.append([bp1, bp2])
                self.indel_sizes.append("D"+str(del_bp))
                self.add_query_brkpt(qend1)
                self.cigar += str(del_bp) + "D"
                if del_bp > self.indel_maxevent_size[0]:
                    self.indel_maxevent_size = [del_bp, "D"]

        self.cigar += str(self.query_blocksizes[-1]) + "M"
        end_clipped = self.get_size('query') - self.get_coords('query')[1]
        if end_clipped > 0:
            self.cigar += str(end_clipped) + "S"

        self.set_indel_flank_matches()

        if self.strand == "-":
            for i in range(len(self.query_brkpts)):
                self.query_brkpts[i] = self.get_size('query') - self.query_brkpts[i]

    def add_query_brkpt(self, brkpt):

        '''
        '''

        if brkpt not in self.query_brkpts:
            self.query_brkpts.append(brkpt)

    def get_brkpt_str(self, with_sizes=False):

        '''
        '''

        brkpt_out = []
        bp_str = []
        chrm = str(self.get_name('hit'))
        if len(self.breakpts) > 0:
            for b, s in zip(self.breakpts, self.indel_sizes):
                if len(b) > 1:
                    bb = "-".join([str(x) for x in b])
                else:
                    bb = str(b[0])
                bstr = chrm + ":" + bb
                if with_sizes:
                    bstr += " " + "(" + s + ")"
                bp_str.append(bstr)
            brkpt_out.append(",".join(bp_str))
        return ",".join(brkpt_out)

    def get_brkpt_locs(self):

        '''
        '''

        brkpt_locs = []
        for b in self.breakpts:
            brkpt_locs.extend(b)
        return brkpt_locs

    def get_gene_anno(self):

        '''
        '''

        return self.genes

    def get_blat_output(self):

        '''
        '''

        return "\t".join([str(x) for x in self.values])

    def get_len(self):

        '''
        '''

        return self.qend - self.qstart

    def set_gene_anno(self, annotations, query_region):

        '''
        '''

        br_start = self.get_coords('hit')[0]
        br_end = self.get_coords('hit')[1]
        start_in = br_start >= (query_region[1]-200) and br_start <= (query_region[2]+200)
        end_in = br_end <= (query_region[2]+200) and br_end >= (query_region[1]-200)
        if query_region[0] == self.get_name('hit') and (start_in or end_in):
            self.in_target = True
            self.genes = query_region[3]
        else:
            ann_genes = []
            chrom = self.get_name('hit')
            pos = self.get_coords('hit')
            if chrom.find('chr') == -1:
                chrom = 'chr'+str(chrom)
            for g in annotations.genes:
                gs = annotations.genes[g][1]
                ge = annotations.genes[g][2]
                if chrom == annotations.genes[g][0]:
                    if int(pos[0]) >= gs and int(pos[0]) <= ge:
                        ann_genes.append(g)
                        break
            if len(ann_genes) == 0:
                ann_genes = ['intergenic']
                self.valid = False
            self.genes = ",".join(ann_genes)

    def calcMilliBad(self):

        '''
        '''

        bad = 0.0
        qAliSize = self.qend() - self.qstart()
        tAliSize = self.tend() - self.tstart()
        aliSize = min(qAliSize, tAliSize)
        if aliSize <= 0:
          return 0.0
        sizeDif = qAliSize - tAliSize
        if sizeDif < 0:
          sizeDif = 0
        insertFactor = self.gaps['query'][0]
        total = self.matches['match'] + self.matches['rep'] + self.matches['mis']
        if total != 0:
          bad = (1000 * (self.matches['mis'] + insertFactor + round(3*math.log(1+sizeDif)))) / total
        return bad*0.1

    def report_blat_hit(self):

        '''
        '''

        print 'Blat hit'
        print 'Blat values', self.values
        print 'Matches', self.matches
        print 'Gaps', self.gaps
        print 'Fragments', self.fragments
        print 'Percent identity', self.percent_identity


class RealignResultSet(object):

    '''
    '''

    def __init__(self, params, realign_result_fn, query_region, scope):
        self.params = params
        self.realign_result_fn = realign_result_fn
        self.query_region = query_region
        self.scope = scope
        self.qsize = 0
        self.hit_freq = []
        self.nmismatches = 0
        self.ngaps = 0
        self.has_results = True
        # self.blat_results = []
        self.realignments = []
        # self.clipped_qs = []
        # self.se = None
        self.sv_event = None
        self.logging_name = 'breakmer.realign.realign'
        self.set_values()

    def set_values(self):

        '''
        '''
        # pdb.set_trace()
        self.check_results_exist()
        if not self.has_results:
            return

        with open(self.realign_result_fn, 'rU') as realign_result_file:
            for line in realign_result_file:
                realigned_segment = RealignedSegment(self.scope, self.params, self.query_region, line.strip().split('\t'))
                # score_raw = br.get_nmatch_total()
                # ngaps = br.get_ngap_total()
                # in_target = 1 if br.in_target else 0
                # score_frac = float(score_raw)/float(br.get_size('query'))
                # score = score_raw + score_frac
                # perc_ident = br.perc_ident
                self.realignments.append(realigned_segment)
                # self.blat_results.append((score, ngaps, in_target, br, perc_ident))
                self.add_summary_metrics(realigned_segment)
        self.realignments = sorted(self.realignments, key=lambda x: (-x.score, -x.percent_identity, x.ngaps))
        # self.blat_results = sorted(self.blat_results, key=lambda blat_results: (-blat_results[0], -blat_results[4], blat_results[1]) )

    def add_summary_metrics(self, realigned_segment):

        '''
        '''

        self.nmismatches += realigned_segment.get_nmatches('mis')
        self.ngaps += realigned_segment.get_num_gaps()
        if not self.qsize:
            self.qsize = realigned_segment.get_size('query')
            self.hit_freq = [0]*self.qsize
        for i in range(realigned_segment.qstart(), realigned_segment.qend()):
            self.hit_freq[i] += 1

    def check_results_exist(self):

        '''
        '''

        if self.realign_result_fn is None:
            self.has_results = False
        elif len(open(self.realign_result_fn, 'rU').readlines()) == 0:
            self.has_results = False

    def write_mod_result_file(self, filename):

        '''
        '''

        blat_ff = open(filename, 'w')
        for realignment in self.realignments:
            blat_ff.write(realignment.get_blat_output() + '\n')
        blat_ff.close()

    def target_hit(self):
    #single_res = (len(self.blat_results)==1) and (self.get_query_coverage() >= 90.0)
    #mult_res1 = (len(self.blat_results)>1 and (self.get_query_coverage()>=95.0) and (self.nmismatches<5) and (self.ngaps<3))
    #mult_res2 = self.blat_results[0][3].spans_query() and (self.blat_results[0][3].get_nmatches('mis')<5) and (self.blat_results[0][3].get_num_gaps()<3)
        # Single hit with 90% identity.
        indel_hit = self.realignments[0].spans_query() or (len(self.realignments) == 1 and self.get_query_coverage() >= 90.0)
    #    single_res = (len(self.blat_results)==1) and (self.get_query_coverage() >= 90.0)
    #    mult_res1 = False
    #    mult_res2 = False
    # Or multiple hits
    #    if len(self.blat_results) > 1:
    #      mult_res1 = self.blat_results[0][3].get_query_coverage() >= 90.0
    #      print 'blat percent ident', self.blat_results[0][3].perc_ident
    #      qcov = [0]*self.qsize
    #      for br in self.blat_results:
        # Calculate the total query coverage for blat results with 95% identity
    #        print 'Multiple hit', br[3].perc_ident
    #        if br[3].perc_ident >= 95.0:
    #          qcov[br[3].qstart():br[3].qend()] = map(lambda x: x+1, qcov[br[3].qstart():br[3].qend()])
    #      nhits = 0
    #      for i in qcov:
    #        if i>0: nhits += 1
    #      total_query_cov = round((float(nhits)/float(self.qsize))*100,2)
    #      mult_res2 = total_query_cov >= 95.0
    #    sys.exit()
        utils.log(self.logging_name, 'debug', 'Checking if query is a target hit or not %r' % indel_hit)
        return indel_hit #single_res or mult_res1 or mult_res2

    def has_indel(self):

        '''
        '''

        has_indel = False
        for i, realigned_segment in enumerate(self.realignments):
            # br.mean_cov = self.get_mean_cov(br.qstart(), br.qend())
            #keep_clipped = (mean_cov<4 and ((br.get_nmatch_total()<30 and not br.in_repeat) or br.get_nmatch_total()>=30))
            #keep_clipped = keep_clipped or (br.in_target and mean_cov<10)
            #print nmatch, ngaps, br.mean_cov
            # if i == 0 and self.check_segment_indel(realigned_segment):
            if i == 0:
                if realigned_segment.has_gaps() and (realigned_segment.spans_query() or (len(self.realignments) == 1 and realigned_segment.in_target)):
                    has_indel = True
                    utils.log(self.logging_name, 'info', 'Contig has indel, returning %r' % has_indel)
                    self.sv_event = results.SVEvent(realigned_segment)
                    self.sv_event.set_event_type('indel')
                    return has_indel
            # else:
            #     utils.log(self.logging_name, 'debug', 'Storing clipped blat result start %d, end %d' % (realigned_segment.qstart(), realigned_segment.qend()))
                # self.clipped_qs.append((br.qstart(), br.qend(), br, i))
        utils.log(self.logging_name, 'info', 'Contig does not have indel, return %r' % has_indel)
        return has_indel

    def get_query_coverage(self):

        '''
        '''

        nhits = 0
        for i in self.hit_freq:
            if i > 0:
                nhits += 1
        return round((float(nhits)/float(self.qsize))*100, 2)

    def get_mean_cov(self, start_coord, end_coord):

        '''
        '''

        return float(sum(self.hit_freq[start_coord:end_coord])) / float(len(self.hit_freq[start_coord:end_coord]))

    # def check_segment_indel(self, br):

    #     '''
    #     '''

    #     indel = False
    #     # indel_size_thresh = int(self.params.get_param('indel_size'))
    #     utils.log(self.logging_name, 'info', 'Checking if blat result contains an indel variant')

    #     if br.spans_query() or (len(self.realignments) == 1 and br.in_target):
            # utils.log(self.logging_name, 'info', 'Blat result spans query (%r) or only one blat result (%r) and blat result in target (%r)'%(br.spans_query(), (len(self.realignments) == 1), br.in_target))
            # indel = True
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
        # return indel

    # def get_indel_result(self):

    #     '''
    #     '''

    #     return self.sv_event
        # if self.se:
        #     return self.se.get_indel_result()
        # else:
        #     return None

    # def get_svs_result(self):

    #     '''
    #     '''

    #     if self.se:
    #       return self.se.get_svs_result(self.meta_dict['query_region'], self.meta_dict['params'], self.meta_dict['disc_reads'])
    #     else:
    #       return None

    def check_svs(self):

        '''
        '''

        utils.log(self.logging_name, 'info', 'Checking for SVs')
        gaps = [(0, self.qsize)]

        # if len(self.clipped_qs) > 1:
        utils.log(self.logging_name, 'debug', 'Iterating through %d clipped blat results.' % len(self.realignments))
        merged_clip = [0, None]

        for i, realigned_segment in enumerate(self.realignments):
            # qs, qe, blat_result, idx = clipped_qs
            utils.log(self.logging_name, 'debug', 'Blat result with start %d, end %d, chrom %s' % (realigned_segment.qstart(), realigned_segment.qend(), realigned_segment.get_name('hit')))
            gaps = self.iter_gaps(gaps, realigned_segment, i)
            if self.sv_event.qlen > merged_clip[0]:
                merged_clip = [self.sv_event.qlen, self.sv_event]
        self.sv_event = merged_clip[1]
        # else:
        #   utils.log(self.logging_name, 'info', 'There are no more than 1 clipped blat results, not continuing with SVs calling.')

        # print 'Check valid', self.se_valid()
        valid_result = self.sv_event is not None and (len(self.sv_event.realignments) > 1) and self.sv_event.in_target
        if valid_result:
            self.sv_event.set_event_type('rearrangement')
        return valid_result

    def check_add_br(self, qs, qe, gs, ge, br):

        '''
        '''

        utils.log(self.logging_name, 'info', 'Checking to add blat result with start %d, end %d' % (qs, qe))
        add = False
        over_perc = round((float(min(qe, ge) - max(qs, gs)) / float(qe - qs)) * 100)  # Calc % of segment overlaps with gap
        ov_right = 0 # Check overlap with other aligned segments
        if qe > ge:
            ov_right = abs(qe - ge)
        ov_left = 0
        if qs < gs:
            ov_left = abs(qs - gs)
        br.set_segment_overlap(ov_left, ov_right)
        max_seg_overlap = max(ov_right, ov_left)
        utils.log(self.logging_name, 'debug', 'Blat query segment overlaps gap by %f' % over_perc)
        utils.log(self.logging_name, 'debug', 'Max segment overlap %f' % max_seg_overlap)
        utils.log(self.logging_name, 'debug', 'Event in target %r and blat result in target %r' % (self.sv_event.in_target, br.in_target))
        if over_perc >= 50 and (max_seg_overlap < 15 or (br.in_target and self.sv_event.in_target)): # and (self.se.in_target or br.in_target):
            add = True
        utils.log(self.logging_name, 'debug', 'Add blat result to SV event %r' % add)
        return add

    def iter_gaps(self, gaps, realigned_segment, segment_idx):

        '''
        '''

        new_gaps = []
        # qs, qe, br, idx = cq
        segment_start = realigned_segment.qstart()
        segment_end = realigned_segment.qend()
        hit = False
        for gap in gaps:
            gap_start, gap_end = gap
            utils.log(self.logging_name, 'debug', 'Gap coords %d, %d' % (gap_start, gap_end))

            # Options to consider when to add segment to reconstruction
            start_within_gap = (segment_start >= gap_start and segment_start <= gap_end)
            end_within_gap = (segment_end <= gap_end and segment_end >= gap_start)
            gap_edge_dist_start = (segment_start <= gap_start) and ((gap_start - segment_start) < 15)
            gap_edge_dist_end = (segment_end >= gap_end) and ((segment_end - gap_end) < 15)
            opt1 = (gap_edge_dist_start and (end_within_gap or gap_edge_dist_end))
            opt2 = (gap_edge_dist_end and (start_within_gap or gap_edge_dist_start))

            if start_within_gap or end_within_gap or opt1 or opt2:
                ngap = []
                if segment_start > gap_start:
                    if (segment_start-1 - gap_start) > 10:
                        ngap.append((gap_start, segment_start-1))
                if segment_end < gap_end:
                    if (gap_end - segment_end+1) > 10:
                        ngap.append((segment_end + 1, gap_end))
                if segment_idx == 0:
                    utils.log(self.logging_name, 'debug', 'Creating SV event from blat result with start %d, end %d'%(segment_start, segment_end))
                    self.sv_event = results.SVEvent(realigned_segment) #, self.meta_dict['query_region'], self.meta_dict['contig_vals'], self.meta_dict['sbam'])
                    new_gaps.extend(ngap)
                    hit = True
                elif self.check_add_br(segment_start, segment_end, gap_start, gap_end, realigned_segment):
                    utils.log(self.logging_name, 'debug', 'Adding blat result to event')
                    new_gaps.extend(ngap)
                    self.sv_event.add(realigned_segment)
                    hit = True
                else:
                    new_gaps.append(gap)
            else:
                new_gaps.append(gap)
            utils.log(self.logging_name, 'debug', 'New gap coords %s'%(",".join([str(x) for x in new_gaps])))

        if not hit:
            self.sv_event.check_previous_add(realigned_segment)
        return new_gaps

    # def se_valid(self):

    #     '''
    #     '''

    #     return self.sv_event is not None and (len(self.sv_event.realignments) > 1) and self.sv_event.in_target
            # nmissing_query_cov = len(filter(lambda y: y, map(lambda x: x==0, self.sv_event.query_cov)))
            # if nmissing_query_cov < self.params.get_param('trl_minseg_len'):
            #     valid = True
        # return valid


class RealignManager(object):

    '''
    '''

    def __init__(self, params):
        self.logging_name = 'breakmer.realign.realign'
        self.params = params
        self.query_res_fn = None
        self.realign_results = None

    def realignment(self, contig, target_ref_fa_fn, target_region_values):

        '''
        '''
        # pdb.set_trace()

        if contig.contig_fa_fn is None:
            return

        self.query_res_fn = os.path.join(contig.file_path, 'blat_res.target.psl')
        realign_dict = {'binary': self.params.get_param('blat'),
                        'database': target_ref_fa_fn
                       }
        utils.run_blat(realign_dict, self.query_res_fn, contig.contig_fa_fn, 'target') # Run blat against target reference sequence first for speed.

        if self.query_res_fn is None:
            utils.log(self.logging_name, 'info', 'No blat results file %s, no calls for %s.' % (self.query_res_fn, contig.contig_id))
            return

        self.realign_results = RealignResultSet(self.params, self.query_res_fn, target_region_values, 'target')

        if not self.check_target_blat():
            # Blat against whole genome reference fasta
            # for some reason gfClient resets the reference path
            # realign_dict = {'binary': self.params.get_param('gfclient'), 'blat_port': self.params.get_param('blat_port'), 'database': self.params.get_param('reference_fasta_dir')}
            realign_dict = {'binary': self.params.get_param('gfclient'), 'blat_port': self.params.get_param('blat_port'), 'database': '/'}
            self.query_res_fn = os.path.join(contig.file_path, 'blat_res.genome.psl')
            utils.run_blat(realign_dict, self.query_res_fn, contig.contig_fa_fn, 'genome')
            self.realign_results = RealignResultSet(self.params, self.query_res_fn, target_region_values, 'genome')
        return self.realign_results

    def check_target_blat(self):

        '''
        '''

        # am = sv_caller.align_manager(meta_dict)
        # hit, self.query_res_fn = am.check_target_results()

        hit = False
        utils.log(self.logging_name, 'info', 'Checking if target blat contains most of query or if whole genome needs to queried.')

        if not self.realign_results.has_results:
            self.query_res_fn = None
            hit = True
        else:
            self.realign_results.write_mod_result_file(self.query_res_fn + '.mod')
            #print 'Checking target results'
            #print len(self.bm.blat_results)
            #print self.bm.get_query_coverage()
            #print self.bm.nmismatches
            #print self.bm.ngaps
            if self.realign_results.target_hit():
                hit = True
                utils.log(self.logging_name, 'debug', 'Top hit contains whole query sequence, indel variant')
                self.query_res_fn += '.mod'
    #      else:
    #        if max(contig_counts.others) < params.get_sr_thresh('trl'):
    #          hit = True
    #          self.query_res_fn = None
        return hit # , self.query_res_fn
