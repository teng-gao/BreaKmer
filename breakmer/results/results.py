#! /usr/bin/local/python
# -*- coding: utf-8 -*-

'''
BreaKmer results module
'''


import pysam
import breakmer.utils as utils


__author__ = "Ryan Abo"
__copyright__ = "Copyright 2015, Ryan Abo"
__email__ = "ryanabo@gmail.com"
__license__ = "MIT"


class SVResult(object):

    '''
    '''

    def __init__(self, sv_event, sample_bam_file, contig, target_region_values, disc_reads):
        self.logging_name = 'breakmer.results.results'
        self.sample_bam = sample_bam_file
        self.sv_event = sv_event
        self.contig = contig
        self.disc_reads = disc_reads
        self.filter = False
        self.query_region = target_region_values
        self.sv_type = None
        self.rearr_desc = None
        self.query_cov = []
        self.breakpoint_values = None
        self.values = {'anno_genes': None,
                       'target_breakpoints': None,
                       'align_cigar': '',
                       'mismatches': 0,
                       'strands': None,
                       'total_matching': None,
                       'sv_type': None,
                       'sv_subtype': 'NA',
                       'split_read_count': None,
                       'disc_read_count': 0,
                       'contig_id': contig.contig_id,
                       'contig_seq': contig.seq.value,
                       'sv_subtype': None,
                       'breakpoint_coverages': 0
                      }
        self.set_values()

    def set_values(self):

        '''
        '''
        if self.sv_event is None:
            self.filter = True
            return

        self.sv_type = self.sv_event.get_event_type()
        self.query_cov = [0] * len(self.contig.seq.value)
        for realigned_segment in self.sv_event.realignments:
            for i in range(realigned_segment.get_coords('query')[0], realigned_segment.get_coords('query')[1]):
                i = min(i, len(self.query_cov)-1)
                self.query_cov[i] += 1

        if self.sv_type == 'indel':
            self.set_indel_values()
        else:
            self.set_rearrangement_values()

    def set_indel_values(self):

        '''
        '''

        realigned_segment = self.sv_event.realignments[0]
        self.values['anno_genes'] = realigned_segment.get_gene_anno()
        self.values['total_matching'] = str(realigned_segment.get_nmatch_total())
        self.values['mismatches'] = realigned_segment.get_nmatches('mis')
        self.values['strands'] = realigned_segment.strand
        self.values['target_breakpoints'] = realigned_segment.get_brkpt_str(True)
        self.values['align_cigar'] = realigned_segment.cigar
        self.values['sv_type'] = 'indel'
        self.values['split_read_count'] = ",".join([str(self.contig.get_contig_counts().get_counts(x, x, 'indel')) for x in realigned_segment.query_brkpts])
        self.values['breakpoint_coverages'] = self.get_brkpt_coverages()

    def set_rearrangement_values(self):

        '''
        '''

  # def get_svs_result(self, query_region, params, disc_reads):
        utils.log(self.logging_name, 'info', 'Resolving a rearrangement call from the realignment results')
        coord_sorted_realignments = sorted(self.sv_event.realignments, key=lambda x: x.get_coords('query')[0])

        brkpts = {'t': {'in_target': None, 'other': None},
                  'formatted': [],
                  'r': [],
                  'q': [[0, 0], []],
                  'chrs': [],
                  'brkpt_str': [],
                  'tcoords': [],
                  'qcoords': [],
                  'f': []
                 }

        res_values = {'target_breakpoints': [],
                      'align_cigar': [],
                      'sv_type': '',
                      'sv_subtype': 'NA',
                      'strands': [],
                      'mismatches': [],
                      'repeat_matching': [],
                      'anno_genes': [],
                      'disc_read_count': 0,
                      'total_matching': []
                     }

        valid_rearrangement = [True, True]
        max_repeat = 0.0

        for i, realigned_segment in enumerate(coord_sorted_realignments):
            valid_rearrangement[0] = valid_rearrangement[0] and realigned_segment.valid
            # valid_rearrangement[1] = valid_rearrangement[1] and (realigned_segment.rep_man.simple_rep_overlap > 75.0)
            max_repeat = max(max_repeat, realigned_segment.repeat_overlap)
            # res_values['repeat_matching'].append(":".join([str(realigned_segment.repeat_overlap), str(realigned_segment.get_nmatch_total()), str(round(realigned_segment.mean_cov, 3))]))
            res_values['total_matching'].append(str(realigned_segment.get_nmatch_total()))
            res_values['anno_genes'].append(realigned_segment.get_gene_anno())
            res_values['align_cigar'].append(realigned_segment.cigar)
            res_values['strands'].append(realigned_segment.strand)
            res_values['mismatches'].append(realigned_segment.get_nmatches('mis'))
            brkpts = self.get_brkpt_info(realigned_segment, brkpts, i, i==(len(self.sv_event.realignments)-1))

        # result = None
        match_sorted_realignments = sorted(self.sv_event.realignments, key=lambda x: x.get_nmatch_total())
        # self.br_sorted = sorted(self.br_sorted, key=lambda br: br[1])

        if not self.multiple_genes(brkpts['chrs'], brkpts['r'], res_values['anno_genes']):
            brkpt_counts, brkpt_kmers, brkpt_rep_filter = self.get_brkpt_counts_filt(brkpts, 'rearr')
            self.breakpoint_values = {'kmers': brkpt_kmers, 'counts': brkpt_counts, 'rep_filter': brkpt_rep_filter, 'ref_pos': brkpts['r']}
            sv_type, sv_subtype, disc_read_support = self.define_rearr(self.disc_reads, brkpts['r'], res_values['strands'], brkpts['tcoords'], brkpts['qcoords'])
            # if not self.filter_rearr(query_region, params, brkpts['r'], brkpt_counts, brkpt_kmers, rearr_type, disc_read_support):
            res_values['sv_type'] = sv_type
            res_values['sv_subtype'] = sv_subtype
            res_values['disc_read_count'] = disc_read_support
            res_values['anno_genes'] = list(set(res_values['anno_genes']))
            res_values['target_breakpoints'] = brkpts['brkpt_str']
            res_values['split_read_count'] = brkpt_counts['b']
                # result = self.format_result(res_values)
        else:
            brkpt_counts, brkpt_kmers, brkpt_rep_filter = self.get_brkpt_counts_filt(brkpts, 'trl')
            self.breakpoint_values = {'kmers': brkpt_kmers, 'counts': brkpt_counts, 'rep_filter': brkpt_rep_filter, 'ref_pos': brkpts['r']}
            # if not self.filter_trl(valid_rearrangement, query_region, params, brkpt_counts, brkpt_kmers, disc_read_count, res_values['anno_genes'], max_repeat, brkpt_rep_filt):
            res_values['disc_read_count'] = self.check_disc_reads(brkpts['t'], self.disc_reads['disc'])
            res_values['sv_type'] = ['rearrangement']
            res_values['sv_subtype'] = ['trl']
            res_values['target_breakpoints'] = brkpts['brkpt_str']
            res_values['split_read_count'] = brkpt_counts['b']

        for key in res_values.keys():
            if key in self.values:
                self.values[key] = res_values[key]
        self.values['breakpoint_coverages'] = self.get_brkpt_coverages()
        self.valid_rearrangement = valid_rearrangement

    def get_brkpt_info(self, br, brkpt_d, i, last_iter):

        '''
        '''

        ts, te = br.get_coords('hit')
        qs, qe = br.get_coords('query')
        target_key = 'in_target' if br.in_target else 'other'
        brkpt_d['chrs'].append(br.get_name('hit'))
        brkpt_d['tcoords'].append((ts, te))
        brkpt_d['qcoords'].append((qs, qe))
        tbrkpt = []
        filt_rep_start = None
        if i == 0:
            brkpt_d['q'][0] = [max(0, qs - 1), qe]
            brkpt_d['q'][1].append([qe, qe-brkpt_d['q'][0][0], None])
            tbrkpt = [te]
            filt_rep_start = br.filter_reps_edges[0]
            if br.strand == '-':
                tbrkpt = [ts]
                filt_rep_start = br.filter_reps_edges[0]
        elif last_iter: 
            brkpt_d['q'][1][-1][2] = qe - brkpt_d['q'][1][-1][0]
            brkpt_d['q'][1].append([qs, qs-brkpt_d['q'][0][0], qe-qs])
            tbrkpt = [ts]
            filt_rep_start = br.filter_reps_edges[0]
            if br.strand == '-':
                tbrkpt = [te] 
                filt_rep_start = br.filter_reps_edges[1]
        else:
            brkpt_d['q'][1][-1][2] = qe - brkpt_d['q'][1][-1][1]
            brkpt_d['q'][1].append([qs, qs-brkpt_d['q'][0][0], qe-qs])
            brkpt_d['q'][1].append([qe, qe-qs, None])
            brkpt_d['q'][0] = [qs, qe]
            tbrkpt = [ts, te]
            if br.strand == "-":
                filt_rep_start = br.filter_reps_edges[1]
                tbrkpt = [te, ts]

        brkpt_d['brkpt_str'].append(str(br.get_name('hit')) + ":" + "-".join([str(x) for x in tbrkpt]))
        brkpt_d['r'].extend(tbrkpt)
        brkpt_d['f'].append(filt_rep_start)
        brkpt_d['t'][target_key] = (br.get_name('hit'), tbrkpt[0])
        brkpt_d['formatted'].append(str(br.get_name('hit')) + ":" + "-".join([str(x) for x in tbrkpt]))
        return brkpt_d

    def get_brkpt_coverages(self):

        '''
        '''

        brkpts = []
        tbp = self.values['target_breakpoints']
        if not isinstance(tbp, list):
            if tbp.find("(") > -1:
                tbp = [self.values['target_breakpoints'].split()[0]]
        # tbp = tbp.split(',')
        for bp in tbp:
            chrom, locs = bp.split(':')
            chrom = chrom.replace('chr', '')
            ll = locs.split('-')
            if len(ll) > 1: 
                brkpts.append((chrom, int(ll[0]), int(ll[0])+1))
                brkpts.append((chrom, int(ll[1]), int(ll[1])+1))
            else:
                brkpts.append((chrom, int(ll[0]), int(ll[0])+1))
          
        bamfile = pysam.Samfile(self.sample_bam,'rb')

        covs = [0]*len(brkpts)
        bp_index = 0
        for bp in brkpts:
            cov = 0
            c,s,e = bp
            try:
                areads = bamfile.fetch(str(c), s, e)
            except:
                # print 'Passing on brkpt', c, s, e
                continue


            for aread in areads:
                if aread.is_duplicate or aread.is_qcfail or aread.is_unmapped or aread.mapq < 10: 
                    continue
                cov += 1
            covs[bp_index] = cov
            bp_index += 1
        return ",".join([str(x) for x in covs])

    def multiple_genes(self, chrs, brkpts, anno_genes):

        '''
        '''

        mult_genes = True
        if len(set(anno_genes)) == 1:
            utils.log(self.logging_name, 'debug', 'One annotated gene among SVs breakpoints: %s' % ",".join(anno_genes))
            mult_genes = False
        elif self.dup_gene_names(anno_genes) and len(set(chrs)) == 1 and ((max(brkpts) - min(brkpts)) < 10000):
            utils.log(self.logging_name, 'debug', 'Anno genes are not the same, but similar and breakpoints are less than 10Kb from each other %s' % ",".join(anno_genes))
            mult_genes = False
        utils.log(self.logging_name, 'debug', 'Test whether SVs breakpoints are in multiple genes %r' % mult_genes)
        return mult_genes

    def dup_gene_names(self, anno_genes):

        '''
        '''

        found_dup = False
        for i in range(len(anno_genes)-1):
            gene1 = anno_genes[i]
            for gene2 in anno_genes[(i+1):]:
                if (gene1.find(gene2) > -1) or (gene2.find(gene1) > -1):
                    found_dup = True
        return found_dup

    def which_rearr(self, disc_reads, brkpts, strands, tcoords, qcoords):

        '''
        '''

        rearr_values = {'disc_read_count': None, 'sv_type': 'rearrangement', 'sv_subtype': None, 'hit': False}
        if not self.check_overlap(tcoords[0], tcoords[1]):
            utils.log(self.logging_name, 'debug', 'Checking rearrangement svType, strand1 %s, strand2 %s, breakpt1 %d, breakpt %d' % (strands[0], strands[1], brkpts[0], brkpts[1]))
            if strands[0] != strands[1]:
                # Inversion
                # Get discordantly mapped read-pairs
                utils.log(self.logging_name, 'debug', 'Inversion event identified.')
                rearr_values['hit'] = True
                rearr_values['sv_subtype'] = 'inversion'
                rearr_values['disc_read_count'] = 0
                brkpt1 = min(brkpts)
                brkpt2 = max(brkpts)
                bp_buffer = 20
                for read_pair in disc_reads['inv']:
                    r1p, r2p, r1s, r2s, qname = read_pair
                    if r1s == 1 and r2s == 1:
                        if (r1p <= (brkpt1 + bp_buffer)) and (r2p <= (brkpt2 + bp_buffer) and r2p >= (brkpt1 - bp_buffer)):
                            rearr_values['disc_read_count'] += 1
                    else:
                        if (r1p <= (brkpt2 + bp_buffer) and r1p >= (brkpt1 - bp_buffer)) and r2p >= (brkpt2 - bp_buffer):
                            rearr_values['disc_read_count'] += 1 
            elif strands[0] == strands[1]:
                tgap = brkpts[1] - brkpts[0]
                qgap = qcoords[1][0] - qcoords[0][1]
                brkpt1 = min(brkpts)
                brkpt2 = max(brkpts)
                bp_buffer = 20
                if tgap < 0:
                    utils.log(self.logging_name, 'debug', 'Tandem duplication event identified.')
                    rearr_values['hit'] = True
                    rearr_values['sv_subtype'] = 'tandem_dup'
                    rearr_values['disc_read_count'] = 0
                    for read_pair in disc_reads['td']:
                        r1p, r2p, r1s, r2s, qname = read_pair
                        disc_read_check = (r1p >= (brkpt1 - bp_buffer) and r1p <= (brkpt2 + bp_buffer)) and (r2p <= (brkpt2 + bp_buffer) and r2p >= (brkpt1 - bp_buffer))
                        if disc_read_check:
                            rearr_values['disc_read_count'] += 1
                elif tgap > qgap:
                    # Gapped deletion from Blast result
                    utils.log(self.logging_name, 'debug', 'Deletion event identified.')
                    rearr_values['hit'] = True
                    rearr_values['sv_type'] = 'indel'
                    rearr_values['indelSize'] = 'D' + str(tgap)
                else:
                    # Gapped insertion from Blast result
                    utils.log(self.logging_name, 'debug', 'Insertion event identified.')
                    rearr_values['hit'] = True
                    rearr_values['sv_type'] = 'indel'
                    rearr_values['indelSize'] = 'I' + str(qgap)
        return rearr_values

    def define_rearr(self, disc_reads, brkpts, strands, tcoords, qcoords):

        '''
        '''

        sv_type = 'rearrangement'
        sv_subtype = None
        rs = 0
        hit = False

        # if len(strands) < 3:
        #   if not self.check_overlap(tcoords[0], tcoords[1]):
        #     utils.log(self.logging_name, 'debug', 'Checking rearrangement type, strand1 %s, strand2 %s, breakpt1 %d, breakpt %d' % (strands[0], strands[1], brkpts[0], brkpts[1]))
        #     if (strands[0] != strands[1]):  # and (brkpts[0] < brkpts[1]):
        #       # Inversion
        #       # Get discordantly mapped read-pairs
        #       hit = True
        #       sv_type = 'inversion'
        #       for read_pair in disc_reads['inv']:
        #         r1p, r2p, r1s, r2s, qname = read_pair
        #         if r1s == 1 and r2s == 1:
        #           if (r1p <= brkpts[0]) and (r2p <= brkpts[1] and r2p >= brkpts[0]):
        #             rs += 1
        #         else:
        #           if (r1p <= brkpts[1] and r1p >= brkpts[0]) and r2p >= brkpts[1]:
        #             rs += 1 
        #     elif (strands[0] == "+" and strands[1] == "+") and (brkpts[0] > brkpts[1]):
        #       utils.log(self.logging_name, 'debug', 'HIT TANDEM DUP')
        #       hit = True
        #       sv_type = 'tandem_dup'
        #       # Tandem dup
        #       for read_pair in disc_reads['td']:
        #         r1p, r2p, r1s, r2s, qname = read_pair 
        #         if (r1p <= brkpts[0] and r1p >= brkpts[1]) and ():
        #             rs += 1

        rearr_hits = {}
        coord_sorted_realignments = sorted(self.sv_event.realignments, key=lambda x: x.get_coords('query')[0])
        for i in range(1, len(coord_sorted_realignments)):
            vals = self.which_rearr(disc_reads, brkpts[(i - 1):(i + 1)], strands[(i - 1):(i + 1)], tcoords[(i - 1):(i + 1)], qcoords[(i - 1):(i + 1)])
            if vals['hit']:
                if vals['sv_type'] not in rearr_hits:
                    rearr_hits[vals['sv_type']] = []
                rearr_hits[vals['sv_type']].append(vals)

        if 'rearrangement' not in rearr_hits:
            utils.log(self.logging_name, 'debug', 'Error in realignment parsing. Indel found without rearrangement event.')

        rearr_hit = False
        for rearr in rearr_hits:
            for i, rr in enumerate(rearr_hits[rearr]):
                if rearr == 'rearrangement':
                    if not rearr_hit:
                        sv_subtype = rearr_hits[rearr][i]['sv_subtype']
                        rs = int(rearr_hits[rearr][i]['disc_read_count'])
                        rearr_hit = True
                    else:
                        sv_subtype = None
                        if self.rearr_desc is None:
                            self.rearr_desc = [sv_subtype]
                        self.rearr_desc.append(rearr_hits[rearr][i]['sv_subtype'])
                else:
                    if self.rearr_desc is None:
                        self.rearr_desc = []
                    self.rearr_desc.append(rearr_hits[rearr][i]['indelSize'])

        if sv_subtype is None:
            utils.log(self.logging_name, 'debug', 'Not inversion or tandem dup, checking for odd read pairs around breakpoints')
            utils.log(self.logging_name, 'debug', 'Not inversion or tandem dup, checking for odd read pairs around breakpoints')
            rs = [0] * len(brkpts)
            for i in range(len(brkpts)): 
                b = brkpts[i]
                for read_pair in disc_reads['other']:
                    r1p, r2p, r1s, r2s, qname = read_pair
                    if abs(r1p-b) <= 300 or abs(r2p-b) <= 300:
                        utils.log(self.logging_name, 'debug', 'Adding read support from read %s, with strands %d, %d and positions %d, %d for breakpoint at %d' % (qname, r1s, r2s, r1p, r2p, b))
                        rs[i] += 1
            rs = max(rs)

        return sv_type, sv_subtype, rs

    def check_overlap(self, coord1, coord2):

        '''
        '''

        contained = False
        if coord1[0] >= coord2[0] and coord1[1] <= coord2[1]:
          contained = True
        elif coord2[0] >= coord1[0] and coord2[1] <= coord1[1]:
          contained = True
        return contained

    def check_disc_reads(self, brkpts, disc_reads):

        '''
        '''

        disc_read_count = 0
        if brkpts['other'][0] in disc_reads:
            for p1, p2 in disc_reads[brkpts['other'][0]]:
                d1 = abs(p1 - brkpts['in_target'][1])
                d2 = abs(p2 - brkpts['other'][1])
                if d1 <= 1000 and d2 <= 1000:
                    disc_read_count += 1
        return disc_read_count

    def get_brkpt_counts_filt(self, brkpts, sv_type):

        '''

        '''

        avg_comp, comp_vec = utils.calc_contig_complexity(self.contig.seq.value)
        brkpt_rep_filt = False
        brkpt_counts = {'n':[],'d':[],'b':[]}
        brkpt_kmers = []
        for qb in brkpts['q'][1]:
            left_idx = qb[0] - min(qb[1],5)
            right_idx = qb[0] + min(qb[2],5)
            bc = self.contig.get_contig_counts().get_counts(left_idx, right_idx, sv_type)
            brkpt_counts['n'].append(min(bc))
            brkpt_counts['d'].append(min(self.contig.get_contig_counts().get_counts((qb[0]-1), (qb[0]+1), sv_type)))
            brkpt_counts['b'].append(self.contig.get_contig_counts().get_counts(qb[0], qb[0], sv_type))
            brkpt_kmers.append(self.contig.get_kmer_locs()[qb[0]])
            brkpt_rep_filt = brkpt_rep_filt or (comp_vec[qb[0]] < (avg_comp/2))
            utils.log(self.logging_name, 'debug', 'Read count around breakpoint %d: %s'%(qb[0],",".join([str(x) for x in bc])))
        utils.log(self.logging_name, 'debug', 'Kmer count around breakpoints %s'%(",".join([str(x) for x in brkpt_kmers])))
        brkpt_rep_filt = brkpt_rep_filt or (len(filter(lambda x: x, brkpts['f'])) > 0)
        return brkpt_counts, brkpt_kmers, brkpt_rep_filt

    def get_output_string(self):

        '''
        '''
        output_list_keys = ['anno_genes', 'target_breakpoints', 'mismatches', 'strands', 'total_matching', 'sv_type', 'sv_subtype', 'split_read_count', 'disc_read_count', 'breakpoint_coverages', 'contig_id', 'contig_seq']
        output_list = []
        for key in output_list_keys:
            if not isinstance(self.values[key], list):
                output_list.append(str(self.values[key]))
            else:
                value_list = self.values[key]
                formatted_str = ",".join([str(x) for x in value_list])
                output_list.append(formatted_str)
        return '\t'.join(output_list) + '\n'

    # def format_result(self, values):

    #     '''
    #     '''

    #     res_lst = []
    #     if values: 
    #       for v in values:
    #         if not isinstance(values[v], list): values[v] = [values[v]]
    #         self.result_values[v] = ",".join([str(x) for x in values[v]])
    #     if self.result_values['sv_subtype']: self.result_values['sv_type'] += '_' + self.result_values['sv_subtype']
    #     return self.get_values()

        # return self.format_result(None)

class SVEvent(object):

    '''
    '''

    def __init__(self, realigned_segment):
        self.logging_name = 'breakmer.results.results'
        # self.events = []
        # self.blat_res = []
        # self.br_sorted = []
        # self.sample_bam = sample_bam
        self.qlen = 0
        self.nmatch = 0
        self.in_target = False
        # self.query_region = query_region
        # self.contig_seq = contig_vals[0]
        # self.contig_rcounts = contig_vals[1]
        # self.contig_id = contig_vals[2]
        # self.contig_reads = contig_vals[3]
        # self.nkmers = contig_vals[4]
        # self.contig_kmer_locs = contig_vals[5]
        # self.logging_name = 'breakmer.sv_caller'
        self.valid = True
        self.in_rep = False
        self.query_size = None
        self.realignments = []
        # self.query_cov = [0]*len(self.contig_seq)
        self.homology = {'non':[False, None], 'in':[False, None]}
        self.event_type = None
        # self.result_values = SVResult()
        self.add(realigned_segment)

    def add(self, realigned_segment):

        '''
        '''

        # self.realignments.append((realigned_segment.get_coords('query')[0], realigned_segment))
        self.realignments.append(realigned_segment)
        # for i in range(realigned_segment.get_coords('query')[0], realigned_segment.get_coords('query')[1]):
        #     self.query_cov[i] += 1

        if not self.query_size:
            self.query_size = realigned_segment.get_size('query')

        self.qlen += realigned_segment.get_query_span()
        self.nmatch += realigned_segment.get_nmatch_total()
        self.in_target = self.in_target or realigned_segment.in_target
        self.in_rep = self.in_rep and (realigned_segment.repeat_overlap > 75.0)
        self.valid = self.valid and realigned_segment.valid
        # self.br_sorted.append((realigned_segment, realigned_segment.get_nmatch_total()))

    def check_previous_add(self, new_segment):

        '''
        '''

        ncoords = new_segment.get_coords('query')
        # prev_br, prev_nmatch = self.br_sorted[-1]
        previous_segment = self.realignments[-1]
        previous_coords = previous_segment.get_coords('query')
        if ncoords[0] == previous_coords[0] and ncoords[1] == previous_coords[1]:
            # new_segment_nmatch = new_segment.get_nmatch_total()
            if abs(previous_segment.get_nmatch_total() - new_segment.get_nmatch_total()) < 10:
                if not previous_segment.in_target and new_segment.in_target:
                    self.realignments[-1] = new_segment
                    # self.blat_res[-1] = (ncoords[0], br)
                    self.in_target = True

    # def set_result_value(self, key, value):
    #     self.result_values[key] = value

    # def format_result(self, values):

    #     '''
    #     '''

    #     res_lst = []
    #     if values: 
    #       for v in values:
    #         if not isinstance(values[v], list): values[v] = [values[v]]
    #         self.result_values[v] = ",".join([str(x) for x in values[v]])
    #     if self.result_values['sv_subtype']: self.result_values['sv_type'] += '_' + self.result_values['sv_subtype']
    #     return self.get_values()

    def get_event_type(self):

        '''
        '''

        return self.event_type

    def set_event_type(self, event_type):

        '''
        '''

        self.event_type = event_type

#   def get_brkpt_coverages(self):
#     brkpts = []
#     tbp = self.result_values['target_breakpoints']
#     if self.result_values['target_breakpoints'].find("(") > -1:
#       tbp = self.result_values['target_breakpoints'].split()[0]
#     tbp = tbp.split(',')
#     for bp in tbp:
#       chrom,locs = bp.split(':')
#       chrom = chrom.replace('chr','')
#       ll = locs.split('-')
#       if len(ll) > 1: 
#         brkpts.append((chrom,int(ll[0]),int(ll[0])+1))
#         brkpts.append((chrom,int(ll[1]),int(ll[1])+1))
#       else:
#         brkpts.append((chrom, int(ll[0]), int(ll[0])+1))
      
#     bamfile = pysam.Samfile(self.sample_bam,'rb')

#     covs = [0]*len(brkpts)
#     bp_index = 0
#     for bp in brkpts:
#       cov = 0
#       c,s,e = bp
#       areads = bamfile.fetch(str(c), s, e)

#       for aread in areads:
#         if aread.is_duplicate or aread.is_qcfail or aread.is_unmapped or aread.mapq < 10: 
#           continue
#         cov += 1
#       covs[bp_index] = cov
#       bp_index += 1
#     return ",".join([str(x) for x in covs])

#   def get_values(self):
#     lst = ['anno_genes', 'target_breakpoints', 'align_cigar', 'mismatches', 'strands', 'repeat_matching', 'sv_type', 'split_read_count', 'nkmers', 'disc_read_count', 'breakpoint_coverages', 'contig_id', 'contig_seq']
#     out_lst = []
#     for l in lst:
#       if self.result_values[l] == 'trl':
#         self.result_values[l] = 'rearrangement'=
#       out_lst.append(str(self.result_values[l]))
#       if l == 'target_breakpoints':
#         self.result_values['breakpoint_coverages'] = self.get_brkpt_coverages()
#     return out_lst

#   def get_indel_result(self):
#     br = self.blat_res[0][1]
#     self.set_result_value('anno_genes', br.get_gene_anno()) 
#     self.set_result_value('repeat_matching', '0.0:'+str(br.get_nmatch_total()))
#     self.set_result_value('mismatches', br.get_nmatches('mis'))
#     self.set_result_value('strands', br.strand)
#     self.set_result_value('target_breakpoints', br.get_brkpt_str(True))
#     self.set_result_value('align_cigar', br.cigar)
#     self.set_result_value('sv_type', 'indel')
#     self.set_result_value('split_read_count',",".join([str(self.contig_rcounts.get_counts(x,x,'indel')) for x in br.query_brkpts]))
#     return self.format_result(None)

#   def get_brkpt_info(self, br, brkpt_d, i, last_iter):

#     '''
#     '''

#     ts, te = br.get_coords('hit')
#     qs, qe = br.get_coords('query')
#     target_key = 'in_target' if br.in_target else 'other'
#     brkpt_d['chrs'].append(br.get_name('hit'))
#     brkpt_d['tcoords'].append((ts,te))
#     tbrkpt = []
#     filt_rep_start = None
#     if i == 0:
#       brkpt_d['q'][0] = [max(0,qs-1),qe]
#       brkpt_d['q'][1].append([qe, qe-brkpt_d['q'][0][0], None])
#       tbrkpt = [te]
#       filt_rep_start = br.filter_reps_edges[0]
#       if br.strand == '-':
#         tbrkpt = [ts]
#         filt_rep_start = br.filter_reps_edges[0]
#     elif last_iter: 
#       brkpt_d['q'][1][-1][2] = qe - brkpt_d['q'][1][-1][0]
#       brkpt_d['q'][1].append([qs,qs-brkpt_d['q'][0][0],qe-qs])
#       tbrkpt = [ts]
#       filt_rep_start = br.filter_reps_edges[0]
#       if br.strand == '-':
#         tbrkpt = [te] 
#         filt_rep_start = br.filter_reps_edges[1]
#     else:
#       brkpt_d['q'][1][-1][2] = qe - brkpt_d['q'][1][-1][1]
#       brkpt_d['q'][1].append([qs,qs-brkpt_d['q'][0][0],qe-qs])
#       brkpt_d['q'][1].append([qe, qe-qs, None])
#       brkpt_d['q'][0] = [qs, qe]
#       tbrkpt = [ts, te]
#       if br.strand == "-":
#         filt_rep_start = br.filter_reps_edges[1]
#         tbrkpt = [te, ts]

#     brkpt_d['brkpt_str'].append('chr'+str(br.get_name('hit')) + ":" + "-".join([str(x) for x in tbrkpt]))
#     brkpt_d['r'].extend(tbrkpt)
#     brkpt_d['f'].append(filt_rep_start)
#     brkpt_d['t'][target_key] = (br.get_name('hit'),tbrkpt[0])
#     brkpt_d['formatted'].append( 'chr'+str(br.get_name('hit')) + ":" + "-".join([str(x) for x in tbrkpt]))
#     return brkpt_d

#   def get_svs_result(self, query_region, params, disc_reads):
#     utils.log(self.logging_name, 'info', 'Resolving SVs call from blat results')
#     blat_res = self.blat_res
#     blat_res_sorted = sorted(blat_res, key=lambda blat_res: blat_res[0]) 
#     brkpts = {'t':{'in_target':None, 'other':None }, 'formatted':[], 'r':[], 'q': [[0,0],[]], 'chrs':[], 'brkpt_str':[], 'tcoords':[], 'f':[]}
#     res_values = {'target_breakpoints':[], 'align_cigar':[], 'sv_type':'', 'strands':[], 'mismatches':[], 'repeat_matching':[], 'anno_genes':[], 'disc_read_count':0 }
#     br_valid = [True, True]
#     max_repeat = 0.0

#     for i in range(len(blat_res_sorted)):
#       br = blat_res_sorted[i][1]
#       br_valid[0] = br_valid[0] and br.valid
#       br_valid[1] = br_valid[1] and (br.rep_man.simple_rep_overlap > 75.0)
#       max_repeat = max(max_repeat, br.repeat_overlap)
#       res_values['repeat_matching'].append(":".join([str(br.repeat_overlap), str(br.get_nmatch_total()), str(round(br.mean_cov,3))]))
#       res_values['anno_genes'].append(br.get_gene_anno())
#       res_values['align_cigar'].append(br.cigar)
#       res_values['strands'].append(br.strand)
#       res_values['mismatches'].append(br.get_nmatches('mis'))
#       brkpts = self.get_brkpt_info(br, brkpts, i, i==(len(blat_res_sorted)-1))

#     result = None
#     self.br_sorted = sorted(self.br_sorted, key=lambda br: br[1])
#     if not self.multiple_genes(brkpts['chrs'], brkpts['r'], res_values['anno_genes']): 
#       brkpt_counts, brkpt_kmers, brkpt_rep_filt = self.get_brkpt_counts_filt(brkpts, 'rearr') 
#       rearr_type, disc_read_support = self.define_rearr(brkpts['r'], res_values['strands'], brkpts['tcoords'], disc_reads) 
#       if not self.filter_rearr(query_region, params, brkpts['r'], brkpt_counts, brkpt_kmers, rearr_type, disc_read_support):
#         res_values['sv_type'] = 'rearrangement'
#         if rearr_type != 'rearrangement':
#           res_values['sv_subtype'] = rearr_type
#         res_values['disc_read_count'] = disc_read_support
#         res_values['anno_genes'] = list(set(res_values['anno_genes']))
#         res_values['target_breakpoints'] = brkpts['brkpt_str']
#         res_values['split_read_count'] = brkpt_counts['b']
#         result = self.format_result(res_values)
#     elif max(self.contig_rcounts.others) >= params.get_param('trl_sr_thresh'):
#       brkpt_counts, brkpt_kmers, brkpt_rep_filt = self.get_brkpt_counts_filt(brkpts, 'trl') 
#       disc_read_count = self.check_disc_reads(brkpts['t'], query_region, disc_reads['disc'])
#       if not self.filter_trl(br_valid, query_region, params, brkpt_counts, brkpt_kmers, disc_read_count, res_values['anno_genes'], max_repeat, brkpt_rep_filt): 
#         res_values['disc_read_count'] = disc_read_count
#         res_values['sv_type'] = ['trl']
#         res_values['target_breakpoints'] = brkpts['brkpt_str']
#         res_values['split_read_count'] = brkpt_counts['b']
#         result = self.format_result(res_values)
#     return result

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

#   def call_rearr(self):
#       rearr_type, disc_read_support = self.define_rearr(brkpts['r'], res_values['strands'], brkpts['tcoords'], disc_reads) 
#       if not self.filter_rearr(query_region, params, brkpts['r'], brkpt_counts, brkpt_kmers, rearr_type, disc_read_support):
#         res_values['sv_type'] = 'rearrangement'
#         if rearr_type != 'rearrangement': res_values['sv_subtype'] = rearr_type
#         res_values['disc_read_count'] = disc_read_support
#         res_values['anno_genes'] = list(set(res_values['anno_genes']))
#         res_values['target_breakpoints'] = brkpts['brkpt_str']
#         res_values['split_read_count'] = brkpt_counts['b']
#         result = self.format_result(res_values)

#   def call_trl(self):
#       disc_read_count = self.check_disc_reads(brkpts['t'], query_region, disc_reads['disc'])
#       if not self.filter_trl(br_valid, query_region, params, brkpt_counts, brkpt_kmers, disc_read_count, res_values['anno_genes'], max_repeat, brkpt_rep_filt): 
#         res_values['disc_read_count'] = disc_read_count
#         res_values['sv_type'] = ['trl']
#         res_values['target_breakpoints'] = brkpts['brkpt_str']
#         res_values['split_read_count'] = brkpt_counts['b']
#         result = self.format_result(res_values)

#   def check_overlap(self, coord1, coord2 ):
#     contained = False
#     if coord1[0] >= coord2[0] and coord1[1] <= coord2[1]:
#       contained = True
#     elif coord2[0] >= coord1[0] and coord2[1] <= coord1[1]:
#       contained = True
#     return contained 

#   def define_rearr(self, brkpts, strands, tcoords, disc_reads):
#     type = 'rearrangement'
#     rs = 0
#     hit = False
#     if len(strands) < 3:
#       if not self.check_overlap(tcoords[0], tcoords[1]):
#         utils.log(self.logging_name, 'debug', 'Checking rearrangement type, strand1 %s, strand2 %s, breakpt1 %d, breakpt %d'%(strands[0], strands[1], brkpts[0], brkpts[1]))
#         if (strands[0] != strands[1]) and (brkpts[0] < brkpts[1]):
#           # Inversion
#           # Get discordantly mapped read-pairs
#           utils.log(self.logging_name, 'debug', 'HIT INVERSION')
#           hit = True
#           type = 'inversion'
#           for read_pair in disc_reads['inv']:
#             r1p, r2p, r1s, r2s, qname = read_pair
#             if r1s == 1 and r2s == 1:
#               if (r1p <= brkpts[0]) and (r2p <= brkpts[1] and r2p >= brkpts[0]):
#                 rs += 1
#             else:
#               if (r1p <= brkpts[1] and r1p >= brkpts[0]) and r2p >= brkpts[1]:
#                 rs += 1 
#         elif (strands[0] == "+" and strands[1] == "+") and (brkpts[0] > brkpts[1]):
#           utils.log(self.logging_name, 'debug', 'HIT TANDEM DUP')
#           hit = True
#           type = 'tandem_dup'
#           # Tandem dup
#           for read_pair in disc_reads['td']:
#             r1p, r2p, r1s, r2s, qname = read_pair 
#             if (r1p <= brkpts[0] and r1p >= brkpts[1]) and (): rs += 1
#     if not hit:
#       utils.log(self.logging_name, 'debug', 'Not inversion or tandem dup, checking for odd read pairs around breakpoints')
#       rs = [0] * len(brkpts)
#       for i in range(len(brkpts)): 
#         b = brkpts[i]
#         for read_pair in disc_reads['other']:
#           r1p, r2p, r1s, r2s, qname = read_pair
#           if abs(r1p-b) <= 300 or abs(r2p-b) <= 300:
#             utils.log(self.logging_name, 'debug', 'Adding read support from read %s, with strands %d, %d and positions %d, %d for breakpoint at %d' % (qname, r1s, r2s, r1p, r2p, b))
#             rs[i] += 1
#       rs = max(rs)
#     return type, rs   

