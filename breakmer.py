#! /usr/bin/local/python
# -*- coding: utf-8 -*-

import sys
import argparse
import breakmer.params as params
import sv_processor

__author__ = "Ryan Abo"
__copyright__ = "Copyright 2015, Ryan Abo"
__email__ = "ryanabo@gmail.com"
__license__ = "MIT"

'''
Main script that initiates the BreaKmer analysis or auxiliary functions to setup BreaKmer for analysis.
'''


def parse_config_file(arguments) :
    '''
    '''

    param_opts = {}
    # config_f = open(pArgs.config_fn, 'rU')
    # flines = config_f.readlines()

    for line in open(arguments.config_fn, 'rU'):
        line = line.strip()
        if line == '' or line.find('#') > -1:  # Allow for blank lines and comments
            continue
        linesplit = line.split("=")
        if len(linesplit) == 1:  # Make sure the lines in the configuration file are set properly.
            err_msg = 'Config line', line, ' not set correctly. Exiting.'
            print err_msg
            # utils.log(self.loggingName, 'error', err_msg)
            sys.exit(2)
        else:
            key, value = linesplit
            param_opts[key] = value

    # Store all the arguments into the self.opts dictionary.
    for opt in vars(arguments):
        param_opts[opt] = vars(arguments)[opt]
        # if (self.get_param(opt) is not None) and (vars(arguments)[opt] is None):
        #     utils.log(self.loggingName, 'info', 'Parameter %s is set in config file and not on the command line. Using config file value %s.' % (opt, self.get_param(opt)))
        # else:
        #     self.set_param(opt, vars(arguments)[opt])

    # for line in flines :
    #     line = line.strip()
    #     linesplit = line.split("=")
    #     if len(linesplit) == 1 : 
    #         print 'Config line', line, ' not set correctly. Exiting.'
    #         sys.exit(2)
    #     else :
    #         k, v = linesplit
    #         param_opts[k] = v

    # for opt in vars(pArgs) : 
    #     param_opts[opt] = vars(pArgs)[opt]
    return param_opts

PARSER = argparse.ArgumentParser(description='Program to identify structural variants within targeted locations.', usage='%(prog)s [options]', add_help=True)
PARSER.add_argument('-l', '--log_level', dest='log_level', default='DEBUG', help='Log level [default: %default]')
PARSER.add_argument('-a', '--keep_repeat_regions', dest='keep_repeat_regions', default=False, action='store_true', help='Keep indels in repeat regions. No repeat mask bed file required if set. [default: %default]')
PARSER.add_argument('-p', '--preset_ref_data', dest='preset_ref_data', default=False, action="store_true", help='Preset all the reference data for all the targets before running analysis. [default: %default]')
PARSER.add_argument('--indel_size', dest='indel_size', default=15, type=int, help='Indel size filter [default: %default]')
PARSER.add_argument('--trl_sr_thresh', dest='trl_sr_thresh', default=2, type=int, help='Split read support threshold for translocations [default: %default]')
PARSER.add_argument('--indel_sr_thresh', dest='indel_sr_thresh', default=5, type=int, help='Split read support threshold for indels [default: %default]')
PARSER.add_argument('--rearr_sr_thresh', dest='rearr_sr_thresh', default=3, type=int, help='Split read support threshold for rearrangements [default: %default]')
PARSER.add_argument('-g', '--gene_list', dest='gene_list', default=None, help='Gene list to consider for analysis [default: %default]')
PARSER.add_argument('-k', '--keep_intron_vars', dest='keep_intron_vars', default=False, action='store_true', help='Keep intronic indels or rearrangements [default: %default]')
PARSER.add_argument('-v', '--var_filter', dest='var_filter', default='all', help='Variant types to report (all, indel, trl, rearrangment) [default: %default]')
PARSER.add_argument('-m', '--rearr_min_seg_len', dest='rearr_minseg_len', default=30, type=int, help='Threshold for minimum segment to be rearranged [default: %default]')
PARSER.add_argument('-n', '--trl_min_seg_len', dest='trl_minseg_len', default=25, type=int, help='Threshold for minimum length of translocation segment [default: %default]')
PARSER.add_argument('-t', '--align_thresh', dest='align_thresh', default=.90, type=int, help='Threshold for minimum read alignment for assembly [default: %default]')
PARSER.add_argument('-z', '--no_output_header', dest='no_output_header', default=False, action='store_true', help='Suppress output headers [default: %default]')
PARSER.add_argument('-c', '--config', dest='config_fn', default=None, required=True, help='The configuration filename that contains additional parameters. [default: %(default)s]')

# usage = '%prog [options] <config file name>'
# desc = """Script to identify structural variants within targeted locations."""
# parser = OptionParser(usage=usage,description=desc)
# parser.add_option('-l', '--log_level', dest='log_level', default='DEBUG', type='string', help='Log level [default: %default]')
# parser.add_option('-a', '--keep_repeat_regions', dest='keep_repeat_regions', default=False, action='store_true', help='Keep indels in repeat regions. No repeat mask bed file required if set. [default: %default]')
# parser.add_option('-p', '--preset_ref_data', dest='preset_ref_data', default=False, action="store_true", help='Preset all the reference data for all the targets before running analysis. [default: %default]')
# parser.add_option('-s', '--indel_size', dest='indel_size', default=15, type='int', help='Indel size filter [default: %default]')
# parser.add_option('-c', '--trl_sr_thresh', dest='trl_sr_thresh', default=2, type='int', help='Split read support threshold for translocations [default: %default]')
# parser.add_option('-d', '--indel_sr_thresh', dest='indel_sr_thresh', default=5, type='int', help='Split read support threshold for indels [default: %default]')
# parser.add_option('-r', '--rearr_sr_thresh', dest='rearr_sr_thresh', default=3, type='int', help='Split read support threshold for rearrangements [default: %default]')
# parser.add_option('-g', '--gene_list', dest='gene_list', default=None, type='string', help='Gene list to consider for analysis [default: %default]')
# parser.add_option('-k', '--keep_intron_vars', dest='keep_intron_vars', default=False, action='store_true', help='Keep intronic indels or rearrangements [default: %default]')
# parser.add_option('-v', '--var_filter', dest='var_filter', default='all', type='string', help='Variant types to report (all, indel, trl, rearrangment) [default: %default]')
# parser.add_option('-m', '--rearr_min_seg_len', dest='rearr_minseg_len', default=30, type='int', help='Threshold for minimum segment to be rearranged [default: %default]')
# parser.add_option('-n', '--trl_min_seg_len', dest='trl_minseg_len', default=25, type='int', help='Threshold for minimum length of translocation segment [default: %default]')
# parser.add_option('-t', '--align_thresh', dest='align_thresh', default=.90, type='int', help='Threshold for minimum read alignment for assembly [default: %default]')
# parser.add_option('-z', '--no_output_header', dest='no_output_header', default=False, action='store_true', help='Suppress output headers [default: %default]')

PARGS = PARSER.parse_args()

# if __name__ == '__main__' :
#   try :
#     opts, args = parser.parse_args(sys.argv[1:])
#     if len(args) != 1 :
#       parser.error('Requires input of configuration file name')
#       usage()
#     else :
#       config_fn = args[0]
#   except parser.error:
#     usage()
#     sys.exit(2)

# tic = time.clock()
PMANAGER = params.ParamManager(PARGS)
# CONFIGDICT = parse_config_file(PARGS)
# utils.setup_logger(CONFIGDICT, 'root')
# sys.exit(2)
RUNNER = sv_processor.runner(CONFIGDICT)
RUNNER.run()
