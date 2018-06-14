#! /usr/bin/local/python
# -*- coding: utf-8 -*-

'''
BreaKmer utils module

Function and classes used throughout the program - most of these are
general purpose functions.
'''

import os
import sys
import glob
import logging
import subprocess
import random
import time
import math
import pysam
from Bio import SeqIO
import pdb


__author__ = "Ryan Abo"
__copyright__ = "Copyright 2015, Ryan Abo"
__email__ = "ryanabo@gmail.com"
__license__ = "MIT"


def profile_data(sample_bam_file):

    '''
    '''

    means = [random.random() for i in range(5)]
    param = 0.01
    thresh = 10
    print means

    target_regions = []
    region = []
    in_region = False
    bam_file = pysam.AlignmentFile(sample_bam_file, 'rb')

    bam_header = bam_file.header
    # for seq in bam_header['SQ']:
    #     print seq['SN'], seq['LN']
    # sys.exit()

    print 'Pass 1'
    seq_indices = []
    sampling = 0
    sample_size = 50000000
    chrom_sample_size = 5000000
    while sampling < sample_size:
        seq_index = random.randint(0,23)
        while seq_index in seq_indices:
            seq_index = random.randint(0,23)
        print 'Evaluating %s sites on chrom' % chrom_sample_size, bam_header['SQ'][seq_index]['SN']
        seq_indices.append(seq_index)
        for i, pileup_col in enumerate(bam_file.pileup(bam_header['SQ'][seq_index]['SN'])):
            cov = pileup_col.n
            # If coverage is above a threshold, then determine the
            # coverage using the high quality mapped reads
            if cov > thresh:
                cov2 = 0
                for pileupread in pileup_col.pileups:
                    if pileupread.alignment.mapq > 20:
                        cov2 += 1
                cov = cov2
            closest_k = 0
            smallest_error = 99999
            for k in enumerate(means):
                error = abs(cov-k[1])
                if error < smallest_error:
                    smallest_error = error
                    closest_k = k[0]

            means[closest_k] = means[closest_k]*(1-param) + cov*(param)
            # print pileup_col.reference_name, pileup_col.pos, 'cov:' + str(cov), 'Index:' + str(closest_k), 'means:' + ','.join([str(x) for x in means])
            if i > chrom_sample_size:
                print means
            #     print means
            #     print pileup_col.reference_name, pileup_col.pos
                break
            sampling += 1

    print 'Pass 2'
    cluster_window = [0] * 50
    for i, pileup_col in enumerate(bam_file.pileup()):
        cov = pileup_col.n
        # If coverage is above a threshold, then determine the
        # coverage using the high quality mapped reads
        if cov > thresh:
            cov2 = 0
            for pileupread in pileup_col.pileups:
                if pileupread.alignment.mapq > 20:
                    cov2 += 1
            cov = cov2
        cluster_window.pop(0)
        cluster_window.append(cov)
        avg_cov = float(sum(cluster_window)) / float(len(cluster_window))
        max_cov = max(cluster_window)
        closest_k = 0
        smallest_error = 99999
        for k in enumerate(means):
            error = abs(max_cov-k[1])
            if error < smallest_error:
                smallest_error = error
                closest_k = k[0]

        # print pileup_col.reference_name, pileup_col.pos, 'cov:' + str(cov), 'Index:' + str(closest_k), 'means:' + ','.join([str(x) for x in means])

        if min(means) != means[closest_k]:
            if not in_region:
                in_region = True
                region = [closest_k, pileup_col.reference_name, pileup_col.pos]
                print 'Starting a new region', pileup_col.reference_name, pileup_col.pos, avg_cov
        elif in_region:
            in_region = False
            region.extend([pileup_col.pos, means[region[0]]])
            print 'Ending region', 'chr' + region[1] + ':' + str(region[2]) + '-' + str(region[3]), avg_cov
            target_regions.append(region)

        # means[closest_k] = means[closest_k]*(1-param) + cov*(param)
        # print pileup_col.reference_name, pileup_col.pos, 'cov:' + str(cov), 'Index:' + str(closest_k), 'means:' + ','.join([str(x) for x in means])
        # if i > 5000000:
        #     print means
        #     print target_regions
        # #     print means
        # #     print pileup_col.reference_name, pileup_col.pos
        #     break

    # means.sort()
    # print means
    fout = open('regions.bed', 'w')
    for target_region in target_regions:
        bed_out = '\t'.join([str(x) for x in [target_region[1], target_region[2], target_region[3]]])
        fout.write(bed_out + '\n')
    fout.close()

def profile_data2(sample_bam_file):

    '''Iterate over the bam file coverage information and provide basic
    sequence coverage metrics.

    Args:
        sample_bam_file (str):  The absolute path to the sample bam file.
                                It is assumed to be coordinate sorted with a
                                corresponding index file in the same location.
    Returns:
        None
    '''

    moving_avg = 0
    coverage_window = [0]*10
    bam_file = pysam.AlignmentFile(sample_bam_file, 'rb')
    total = 0
    target_positions = []
    in_target = False
    region = []
    thresh = 10
    region_max = 0
    clusters = []
    for i, pileup_col in enumerate(bam_file.pileup()):
        cov = pileup_col.n
        # If coverage is above a threshold, then determine the
        # coverage using the high quality mapped reads
        if cov > thresh:
            cov2 = 0
            for pileupread in pileup_col.pileups:
                if pileupread.alignment.mapq > 20:
                    cov2 += 1
            cov = cov2

        coverage_window.pop(0)
        coverage_window.append(cov)
        window_avg = float(sum(coverage_window)) / len(coverage_window)
        if in_target:
            if window_avg <= 10*moving_avg or cov < 2:
                # End the target region
                region.append(pileup_col.pos)
                region.append(region_max)
                target_positions.append(region)
                region = []
                region_max = 0
                in_target = False
                if len(clusters) < 1:
                    clusters.append(target_positions[-1][-1])
                else:
                    for cluster in clusters:
                        if abs(cluster - target_positions[-1][-1]) > (cluster*1.5):
                            clusters.append(target_positions[-1][-1])
                            print cluster, target_positions[-1][-1], abs(cluster - target_positions[-1][-1]), (cluster*1.5)
                    print 'clusters', clusters
                if target_positions[-1][-1] >= 100:
                    print len(target_positions), 'target positions', target_positions[-1]
            else:
                region_max = max(region_max, cov)
        elif window_avg > 10*moving_avg and cov > thresh:
            # Initiate a target region
            in_target = True
            region = [pileup_col.reference_name, pileup_col.pos]
            region_max = cov
        else:
            total += cov
            moving_avg = float(total) / float(i+1)
        #print pileup_col.reference_name, pileup_col.pos, pileup_col.n, cov, moving_avg, total, i, window_avg #, target_positions

def setup_logger(log_filename_path, name):
    """Creates the logger object and associated text file to use throughout
    the analysis.

    It first creates a log.txt file in the specified analyis directory as the
    FileHandler. The console handler is then formatted to report the time,
    name of the source, level of the message and the message.

    Args:
        log_filename_path:  Absolute path for the directory that will contain the
                            log file.
        name:               The name of the package to initially setup the logger object.
    Returns:
        None
    """

    output_path = os.path.abspath(log_filename_path)
    if not os.path.exists(output_path):
        os.makedirs(output_path)

    logger = logging.getLogger(name)
    logger.setLevel(logging.DEBUG)

    # FileHandler
    file_handle = logging.FileHandler(os.path.join(output_path, 'log.txt'), mode='w')
    file_handle.setLevel(logging.DEBUG)

    # ConsoleHandler
    console_handle = logging.StreamHandler()
    console_handle.setLevel(logging.ERROR)

    formatter = logging.Formatter(fmt='%(asctime)s - %(name)s - %(levelname)s - %(message)s', datefmt='%m/%d/%Y %I:%M:%S %p')
    file_handle.setFormatter(formatter)
    console_handle.setFormatter(formatter)

    logger.addHandler(file_handle)
    logger.addHandler(console_handle)


def log(name, level, msg):

    '''Write log message to the appropriate level.

    Args:
        name (str):  The logger name, typically the module name.
        level (str): The level of debugging classification.
        msg (str):   The message to log.
    Returns:
        None
    '''

    logger = logging.getLogger(name)
    if level == 'info':
        logger.info(msg)
    elif level == 'debug':
        logger.debug(msg)
    elif level == 'error':
        logger.error(msg)


def start_blat_server(params):

    """Fire up a blat server instance using a random port number and localhost.
    The required files to start a blat server are first checked and created, if
    necessary. These include a genome-wide reference fasta file and a 2bit
    file generated from that fasta file. The faToTwoBit program is used if the
    2bit file needs to be generated on the fly. The gfServer is started and
    we wait while the server is successfully started.

    Starting the blat server can take some time depending on the type of machine being
    used. The different functions that can be run do not all need to the blat server.

    Functions:
      - prepare_reference_data - do not start
      - start_blat_server - start and exit
      - run - start if one is not already running.

    Args:
        params (ParamManager): ParamManager object
    Return:
        None
    """

    logging_name = 'breakmer.utils'
    if params.fnc_cmd == 'prepare_reference_data':  # Do not start blat server for this function.
        return
    elif params.fnc_cmd == 'start_blat_server':
        port = params.get_param('blat_port')
        hostname = params.get_param('blat_hostname')
        params.set_param('blat_hostname', hostname)
        params.set_param('blat_port', port)
        # If no port is specified for this function, then randomly select a port between 8000-9500.
        if port is None:
            params.set_param('blat_port', random.randint(8000, 9500))
            log(logging_name, 'info', 'Starting blat server on port %d on host %s.' % (params.get_param('blat_port'), params.get_param('blat_hostname')))
    elif params.fnc_cmd == 'run':  # Start the blat server if it is not already running.
        if not params.get_param('start_blat_server'):  # Start blat server option is not set. Check that one is running, if not, start it.
            port = params.get_param('blat_port')
            hostname = params.get_param('blat_hostname')
            params.set_param('blat_hostname', hostname)
            if port is None:  # No port is specified for a server that should be running. It will start a new one on a random numbered port.
                log(logging_name, 'debug', 'BreaKmer set to run and start_blat_server is set to False, but no blat server port is specified. Setting blat port to random value and starting blat server.')
                params.set_param('blat_port', random.randint(8000, 9500))
            else:  # Blat server is already running in this instance. Check it to make sure with a test blat.
                params.set_param('blat_port', int(params.get_param('blat_port')))
                if check_blat_server(params):  # Both port and hostname are specified. Check that the server is running.
                    return
                else:
                    log(logging_name, 'debug', 'Blat server with port %d and hostname %s did not pass test query. Please check specifications.' % (params.get_param('blat_port'), params.get_param('blat_hostname')))

    reference_fasta_dir, reference_fasta_filename = os.path.split(params.get_param('reference_fasta'))
    params.set_param('reference_fasta_dir', reference_fasta_dir)
    ref_fasta_name, file_ext = os.path.splitext(reference_fasta_filename)

    # Check if 2bit file exists - if not then create it.
    params.set_param('blat_2bit', os.path.join(params.get_param('reference_fasta_dir'), ref_fasta_name + ".2bit"))
    if not os.path.exists(params.get_param('blat_2bit')):  # Create 2bit file to use for running the blat server.
        log(logging_name, 'info', 'Creating 2bit from %s reference fasta' % params.get_param('reference_fasta'))
        curdir = os.getcwd()
        os.chdir(params.get_param('reference_fasta_dir'))
        twobit_cmd = '%s %s %s' % (params.get_param('fatotwobit'), ref_fasta_name + file_ext, ref_fasta_name + ".2bit")
        twobit_process = subprocess.Popen(twobit_cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True)
        output, errors = twobit_process.communicate()
        if errors != '':
            sys.stderr.write('Creation of 2bit file for fasta file failed. Exiting')
            sys.exit(1)
        os.chdir(curdir)

    curdir = os.getcwd()
    os.chdir(params.get_param('reference_fasta_dir'))

    # Start gfServer, change dir to 2bit file, gfServer start localhost 8000 .2bit
    params.set_param('gfserver_log', os.path.join(params.paths['output'], 'gfserver_%d.log' % params.get_param('blat_port')))
    gfserver_cmd = '%s -canStop -log=%s start %s %d %s &' % (params.get_param('gfserver'), params.get_param('gfserver_log'), params.get_param('blat_hostname'), params.get_param('blat_port'), os.path.join(reference_fasta_dir, ref_fasta_name + ".2bit"))
    log(logging_name, 'info', "Starting gfServer %s" % gfserver_cmd)
    gfserver_process = subprocess.Popen(gfserver_cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True)
    start_time = time.time()
    gfserver_state = server_ready(params.get_param('gfserver_log'))
    while gfserver_state != 'ready':  # Wait for the blat server to initiate. Timeout if it has not started in 15 minutes.
        print("Trying to start blat server..")
        new_time = time.time()
        wait_time = new_time - start_time
        if (wait_time > 1000) or gfserver_state == 'error':
            log(logging_name, 'error', 'gfServer failed to start, please check the gfServer log file %s, exiting' % params.get_param('gfserver_log'))
            sys.exit(1)
        log(logging_name, 'info', 'Waiting for blat gfServer to load reference seq')
        time.sleep(60)
        gfserver_state = server_ready(params.get_param('gfserver_log'))
    log(logging_name, 'info', 'Server ready!')
    print("Blat server ready")
    os.chdir(curdir)


def check_blat_server(params):

    """Run a test query on the specified blat server to make sure it is running.

    Args:
        None
    Returns:
        server_success (boolean): Indicates whether the test ran without errors.
    Raises:
        None
    """

    logging_name = 'breakmer.utils'
    test_dir = os.path.join(params.paths['analysis'], 'blatserver_test')
    test_fasta_filename = os.path.join(test_dir, 'test.fa')
    if not os.path.exists(test_dir):
        os.makedirs(test_dir)

    test_fasta_file = open(test_fasta_filename, 'w')
    test_fasta_file.write('>test\nCCAAGGGAGACTTCAAGCAGAAAATCTTTAAGGGACCCTTGCATAGCCAGAAGTCCTTTTCAGGCTGATGTACATAAAATATTTAGTAGCCAGGACAGTAGAAGGACTGAAGAGTGAGAGGAGCTCCCAGGGCCTGGAAAGGCCACTTTGTAAGCTCATTCTTG')
    test_fasta_file.close()

    result_filename = os.path.join(test_dir, 'blatserver_test.psl')
    blat_cmd = '%s -t=dna -q=dna -out=psl -minScore=20 -nohead %s %d %s %s %s' % (params.get_param('gfclient'), params.get_param('blat_hostname'), params.get_param('blat_port'), params.get_param('reference_fasta_dir'), test_fasta_filename, result_filename)
    log(logging_name, 'info', 'Blat server test system command %s' % blat_cmd)
    blat_process = subprocess.Popen(blat_cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True)
    output, errors = blat_process.communicate()
    log(logging_name, 'info', 'Realignment output file %s' % result_filename)

    server_success = True
    if errors != '':
        server_success = False
        log(logging_name, 'info', 'Realignment errors %s' % errors)
    return server_success


def server_ready(out_file):

    '''Function called by start_blat_server() to poll the gfServer to check if
    the server has been started.

    Args:
        out_file (str):   The gfserver log file path
    Returns:
        ready (boolean): Indicator that the gfServer is ready for queries.
    '''

    logger = logging.getLogger('breakmer.utils')
    while not os.path.exists(out_file):
        logger.info('Waiting for log file %s' % out_file)
        time.sleep(10)
    ready = 'waiting'
    with open(out_file, 'r') as server_file:
        for line in server_file:
            if line.find('Server ready for queries') > -1:
                ready = 'ready'
            elif line.find('error:') > -1 or line.find('abort') > -1:
                ready = 'error'
    return ready


def run_jellyfish(fa_fn, jellyfish, kmer_size):

    '''
    '''

    logging_name = 'breakmer.utils'
    file_path = os.path.split(fa_fn)[0]
    file_base = os.path.basename(fa_fn)
    dump_fn = os.path.join(file_path, file_base + "_" + str(kmer_size) + "mers_dump")
    dump_marker_fn = get_marker_fn(dump_fn)
    if not os.path.isfile(dump_marker_fn):
        if not os.path.exists(fa_fn):
            log(logging_name, 'info', '%s does not exist.' % fa_fn)
            dump_fn = None
            return dump_fn

        version_cmd = '%s --version' % jellyfish
        version_process = subprocess.Popen(version_cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True)
        output, errors = version_process.communicate()

        jfish_version = int(output.split()[1].split('.')[0])
        log(logging_name, 'info', 'Using jellyfish version %d' % jfish_version)

        count_fn = os.path.join(file_path, file_base + "_" + str(kmer_size) + "mers_counts")

        log(logging_name, 'info', 'Running %s on file %s to determine kmers' % (jellyfish, fa_fn))
        jellyfish_cmd = '%s count -m %d -s %d -t %d -o %s %s' % (jellyfish, int(kmer_size), 100000000, 8, count_fn, fa_fn)
        log(logging_name, 'info', 'Jellyfish counts system command %s' % jellyfish_cmd)
        jellyfish_process = subprocess.Popen(jellyfish_cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True)
        output, errors = jellyfish_process.communicate()
        log(logging_name, 'info', 'Jellyfish count output %s and errors %s' % (output, errors))

        if jfish_version < 2:
            count_fn += '_0'
        dump_cmd = '%s dump -c -o %s %s'%(jellyfish, dump_fn, count_fn)
        log(logging_name, 'info', 'Jellyfish dump system command %s' % dump_cmd)
        dump_process = subprocess.Popen(dump_cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True)
        output, errors = dump_process.communicate()
        log(logging_name, 'info', 'Jellyfish dump output %s and errors %s' % (output, errors))
        marker_cmd = 'touch %s' % dump_marker_fn
        marker_process = subprocess.Popen(marker_cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True)
        output, errors = marker_process.communicate()
        log(logging_name, 'info', 'Completed jellyfish dump %s, touching marker file %s' % (dump_fn, dump_marker_fn))
        count_fns = glob.glob(os.path.join(file_path, "*mers_counts*"))

        for count_fn in count_fns:
            os.remove(count_fn)
    else:
        log(logging_name, 'info', 'Jellfish already run and kmers already generated for target.')
    return dump_fn


def extract_refseq_fa(gene_coords, ref_path, ref_fa, direction, target_fa_fn, buffer_size):

    '''
    '''

    # pdb.set_trace()

    logging_name = 'breakmer.utils'
    chrom, start_coord, end_coord, name, intervals = gene_coords
    marker_fn = get_marker_fn(target_fa_fn)

    if not os.path.isfile(marker_fn):
        # ref_d = SeqIO.to_dict(SeqIO.parse(ref_fa, 'fasta'))
        # parse the chromesome fa file instead cuz run out of memory. Quick and dirty fix, need to be unified later
        ref_d = SeqIO.to_dict(SeqIO.parse('%s/chr%s.fa' % (os.path.split(ref_fa)[0], chrom), 'fasta'))
        seq_str = ''
        # seq = ref_d[chrom].seq[(start_coord - buffer_size):(end_coord + buffer_size)]
        seq = ref_d['chr' + chrom].seq[(start_coord - buffer_size):(end_coord + buffer_size)]
        if direction == "reverse":
            seq_str = str(seq.reverse_complement())
        else:
            seq_str = str(seq)
        target_refseq_fasta_filename = open(target_fa_fn, 'w')
        target_refseq_fasta_filename.write(">" + name + "\n" + seq_str + "\n")
        target_refseq_fasta_filename.close()

        touch_cmd = 'touch %s' % marker_fn
        touch_process = subprocess.Popen(touch_cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True)
        output, errors = touch_process.communicate()
        log(logging_name, 'info', 'Completed writing refseq fasta file %s, touching marker file %s' % (target_fa_fn, marker_fn))
    else:
        log(logging_name, 'info', 'Refseq sequence fasta (%s) exists already' % target_fa_fn)

    return target_fa_fn


def setup_ref_data(setup_params):

    '''
    '''

    logger = logging.getLogger('breakmer.utils')
    genes = setup_params[0]
    rep_mask, ref_fa, ref_path, jfish_path, blat_path, kmer_size, buffer_size = setup_params[1]

    for gene in genes:
        chrom, bp1, bp2, name, intvs = gene
        gene_ref_path = os.path.join(ref_path, name)
        if rep_mask:
            logger.info('Extracting repeat mask regions for target gene %s.' % name)
            setup_rmask(gene, gene_ref_path, rep_mask)

        logger.info('Extracting refseq sequence for %s, %s:%d-%d' % (name, chrom, bp1, bp2))
        directions = ['forward', 'reverse']
        for direction in directions:
            target_fa_fn = os.path.join(gene_ref_path, name + '_' + direction + '_refseq.fa')
            ref_fn = extract_refseq_fa(gene, gene_ref_path, ref_fa, direction, target_fa_fn, buffer_size)
            run_jellyfish(ref_fn, jfish_path, kmer_size)

def trim_coords(qual_str, min_qual):

    '''
    '''

    start = seq_trim(qual_str, min_qual)
    if start == len(qual_str):
        return (0, 0, 0)
    else:
        end = len(qual_str) - seq_trim(qual_str[::-1], min_qual)
        lngth = end - start
        return (start, end, lngth)

def get_seq_readname(read):

    '''
    '''

    end = '1'
    if read.is_read2:
        end = '2'
    return read.qname + "/" + end

def get_fastq_reads(fn, sv_reads):

    '''
    '''

    read_len = 0
    filtered_fq_fn = fn.split(".fastq")[0] + "_filtered.fastq"
    filt_fq = open(filtered_fq_fn, 'w')
    fq_recs = {}

    for header,seq,qual in FastqFile(fn):
        qname_split = header.lstrip("@").split("_")
        indel_only = qname_split[-1]
        qname = "_".join(qname_split[0:len(qname_split)-1])
        add = False
        if qname in sv_reads:
            oseq, sc_seqs, clip_coords, indel_meta = sv_reads[qname]
            cleaned_seq = seq
            old_seq = oseq.seq
            add = True
            if str(cleaned_seq) != str(old_seq) and sc_seqs:
                sc_clips = sc_seqs['clipped']
                idx = old_seq.find(cleaned_seq)
                trimmed_seq = ''
                if idx == 0:
                    trimmed_seq = old_seq[len(cleaned_seq):len(old_seq)]
                else: trimmed_seq = old_seq[0:idx]
                sc_lens = 0
                for sc_seq in sc_clips:
                    sc_lens += len(sc_seq)
                    if trimmed_seq.find(sc_seq) > -1:
                        add = False
                if len(cleaned_seq) == (len(old_seq) - sc_lens):
                    for sc_seq in sc_clips:
                        if cleaned_seq.find(sc_seq) == -1:
                      # Don't add, just trimmed clipped portion.
                            add = False
        if add:
            filt_fq.write(header + "\n" + seq + "\n+\n" + qual + "\n")
            fr = fq_read(header, seq, qual, indel_meta)
            read_len = max(read_len, len(fr.seq))
            seq = fr.seq
            if seq not in fq_recs:
                fq_recs[seq] = []
            fq_recs[seq].append(fr)
    filt_fq.close()
    return filtered_fq_fn, fq_recs


def fq_line(read, indel_only, min_len, trim=True):
    '''
    '''

    add_val = '0'
    if indel_only:
        add_val = '1'
    lineout = None
    if trim:
        read = trim_qual(read, 5, min_len)
    if read is not None:
        lineout = "@" + get_seq_readname(read) + "_" + add_val + "\n" + read.seq + "\n+\n" + read.qual + "\n"
    return lineout


def trim_qual(read, min_qual, min_len):

    '''
    '''

    qual_str = read.qual
    start = seq_trim(qual_str, min_qual)
    if start == len(qual_str):
        return None
    else:
        end = len(qual_str) - seq_trim(qual_str[::-1], min_qual)
        lngth = end - start
        if lngth < min_len:
            return None
        nseq = read.seq[start:end]
        nqual = qual_str[start:end]
        read.seq = nseq
        read.qual = nqual
        return read


def seq_trim(qual_str, min_qual):

    '''
    '''

    counter = 0
    while (ord(qual_str[counter]) - 33) < min_qual:
        counter += 1
        if counter == len(qual_str):
            break
    return counter

def get_clip_coords(read_quals, read_cigar):

    '''
    '''

    clip_coords = [0, len(read_quals)]
    coords = [0, 0]
    for i, cigar_tuple in enumerate(read_cigar):
        code, clen = cigar_tuple
        if (code != 2) and (code != 4):
            coords[1] += clen
        if code == 4:
            if i == 0:
                coords[0] = clen
                coords[1] += clen
        clip_coords = coords
    return clip_coords

def is_number(value):

    '''
    '''

    try:
        float(value)
        return True
    except ValueError:
        pass

    try:
        import unicodedata
        unicodedata.numeric(value)
        return True
    except (TypeError, ValueError):
        pass
    return False


def test_cutadapt(fq_fn, cutadapt_binary, cutadapt_config):

    '''
    '''

    fq_clean = os.path.basename(fq_fn).split('.')[0] + "_cleaned.fq"
    fq_clean_fn = os.path.join(os.path.dirname(fq_fn), fq_clean)
    cutadapt_params = stringify(cutadapt_config)
    cutadapt_cmd = '%s %s %s %s > %s'%(sys.executable, cutadapt_binary, cutadapt_params, fq_fn, fq_clean_fn)
    cutadapt_process = subprocess.Popen(cutadapt_cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True)
    output, errors = cutadapt_process.communicate()
    return_code = cutadapt_process.returncode
    if return_code != 0:
        return (None, return_code)
    else:
        return (fq_clean_fn, return_code)

def test_jellyfish(jellyfish_binary, fa_fn, analysis_dir):

    '''
    '''

    version_cmd = '%s --version' % jellyfish_binary
    version_process = subprocess.Popen(version_cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True)
    output, errors = version_process.communicate()
    jellyfish_version = int(output.split()[1].split('.')[0])

    kmer_size = 15
    count_fn = os.path.join(analysis_dir, "test_jellyfish_counts")
    count_cmd = '%s count -m %d -s %d -t %d -o %s %s'%(jellyfish_binary, kmer_size, 100000000, 8, count_fn, fa_fn)
    count_process = subprocess.Popen(count_cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True)
    output, errors = count_process.communicate()
    if count_process.returncode != 0:
        return ("Jellyfish counts", count_process.returncode)

    if jellyfish_version < 2:
        count_fn += '_0'
    dump_fn = os.path.join(analysis_dir, "test_jellyfish_dump")

    dump_cmd = '%s dump -c -o %s %s'%(jellyfish_binary, dump_fn, count_fn)
    dump_process = subprocess.Popen(dump_cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True)
    output, errors = dump_process.communicate()
    if dump_process.returncode != 0:
        return ("Jellyfish dump", dump_process.returncode)
    return ("Jellyfish", 0)


def which(program):

    '''
    '''

    def is_exe(fpath):
        return os.path.isfile(fpath) and os.access(fpath, os.X_OK)

    fpath, fname = os.path.split(program)
    if fpath:
        if is_exe(program):
            return program
    else:
        for path in os.environ["PATH"].split(os.pathsep):
            path = path.strip('"')
            exe_file = os.path.join(path, program)
            if is_exe(exe_file):
                return exe_file
    return None


def mean(lst):

    '''Calculates the mean/average of a list of numbers
    '''

    summation = 0
    for value in lst:
        summation += value
    return summation / len(lst)



def median(lst):

    '''Returns the median value of a list of values.

    Args:
        lst (list):       List of numeric values
    Returns:
        median (numeric): Median value from list.
    '''

    lst = sorted(lst)
    if len(lst) < 1:
        return None
    if len(lst) % 2 == 1:
        return lst[((len(lst) + 1) / 2) - 1]
    else:
        return float(sum(lst[(len(lst) / 2) - 1:(len(lst) / 2) + 1])) / 2.0


def stddev(lst):

    '''Calculates the standard deviation from the list of input  values.

    Args:
        lst (list): List of numeric values to calculate standard deviation
    Returns:
        Standard deviation (float)
    Raises:
        None
    '''

    summation = 0
    avg = mean(lst)
    for value in lst:
        summation += pow((value - avg), 2)
    return math.sqrt(summation / len(lst) - 1)


def remove_outliers(lst):

    '''
    '''

    qnt1 = percentile(lst, 0.25)
    qnt2 = percentile(lst, 0.75)
    intv_cut = 1.5 * (qnt2 - qnt1)
    lst.sort()
    i = 0
    while lst[i] < (qnt1 - intv_cut):
        i += 1
    cut1 = i

    lst.sort(reverse=True)
    while lst[i] > (qnt2 + intv_cut):
        i += 1

    cut2 = len(lst) - i

    lst.sort()
    return lst[cut1:cut2]


def percentile(lst, percent, key=lambda x: x):

    '''Find the percentile of a list of values.
    Args:
        lst (list):       A list of numeric values.
        percent (float):  Percentage value ranging from 0.0 to 1.0.
        key (lambda):     A function to compute a value from each element of N.
    Returns:
        percentile (numeric): The percentile of the values
    '''

    if not lst:
        return None
    lst.sort()
    cut_val = (len(lst) - 1) * percent
    low_val = math.floor(cut_val)
    hi_val = math.ceil(cut_val)
    if low_val == hi_val:
        return key(lst[int(cut_val)])
    d0 = key(lst[int(low_val)]) * (hi_val - cut_val)
    d1 = key(lst[int(hi_val)]) * (cut_val - low_val)
    return d0 + d1

################################
# OLD
################################

class fq_read:
  def __init__(self, header, seq, qual, indel_only):
    self.id = header
    self.seq = str(seq)
    self.qual = str(qual)
    self.used = False
    self.dup = False
    self.indel_only = indel_only


class FastqFile(object):
  def __init__(self,f):
    if isinstance(f,str):
      f = open(f)
      self._f = f
  def __iter__(self):
    return self

  def next(self):
    header, seq, qual_header, qual = [self._f.next() for _ in range(4)]
    header = header.strip()
    # inst,lane,tile,x,y_end = header.split(':')
    seq = seq.strip()
    qual = qual.strip()
    '''
    bc = None
    y = y_end
    if y.find('/') > -1:
      y, end = y.split('/')
    if y.find('#') > -1:
      y, bc = y.split('#')
     header_dict = {'inst':inst,
                   'lane':int(lane),
                   'tile':int(tile),
                   'x':int(x),
                   'y':int(y),
                   'end':end,
                   'bc': bc}
    '''
    return (header, seq, qual)

def stringify(fn):
  # Turn file contents into a space delimited string
  str = []
  for line in open(fn,'rU').readlines():
    line = line.strip()
    str.append(line)
  return ' '.join(str)

def create_ref_test_fa(target_fa_in, test_fa_out):
  if not os.path.isfile(get_marker_fn(test_fa_out)):
    fa_in = open(target_fa_in, "rU")
    fa_out = open(test_fa_out, "w")

    record = SeqIO.read(fa_in, "fasta")
    ref_target_seq = str(record.seq)
    end = min(len(ref_target_seq), 1500)
    start = max(0,len(ref_target_seq)-1500)
    fa_out.write(">"+record.id + "_start\n" + ref_target_seq[0:end] + "\n>" + record.id + "_end\n" + ref_target_seq[start:len(ref_target_seq)] + "\n")
    fa_out.close()

    cmd = 'touch %s'%get_marker_fn(test_fa_out)
    p = subprocess.Popen(cmd,stdout=subprocess.PIPE,stderr=subprocess.PIPE,shell=True)
    output, errors = p.communicate()
    return True
  else:
    return False

def run_blat(realign_value_dict, result_fn, query_fn, scope):

    '''
    '''
    # pdb.set_trace()

    logging_name = 'breakmer.utils'
    if not os.path.isfile(result_fn):
        realign_cmd = None
        if scope == 'genome':
            log(logging_name, 'info', 'Running blat %s, storing results in %s' % (realign_value_dict['binary'], result_fn))
            realign_cmd = '%s -t=dna -q=dna -out=psl -minScore=20 -nohead localhost %d %s %s %s' % (realign_value_dict['binary'], realign_value_dict['blat_port'], realign_value_dict['database'], query_fn, result_fn)
        else:
            log(logging_name, 'info', 'Running blat %s, storing results in %s'%(realign_value_dict['binary'], result_fn))
            realign_cmd = '%s -t=dna -q=dna -out=psl -minScore=20 -stepSize=10 -minMatch=2 -repeats=lower -noHead %s %s %s' % (realign_value_dict['binary'], realign_value_dict['database'], query_fn, result_fn)
            # cmd = '%s -t=dna -q=dna -out=psl -minScore=20 -stepSize=1 -minMatch=1 -repeats=lower -noHead %s %s %s'%(self.params.opts['blat'], db, self.contig_fa_fn, self.query_res_fn)

        log(logging_name, 'info', 'Blat system command %s' % realign_cmd)
        realign_process = subprocess.Popen(realign_cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True)
        output, errors = realign_process.communicate()
        log(logging_name, 'info', 'Blat output %s' % output)
        if errors != '':
            log(logging_name, 'info', 'Blat errors %s' % errors)
    else:
        log(logging_name, 'info', 'Blat already run, results file %s exists, continuing' % result_fn)

# def get_altref_genecoords(blat_path, altref_fa, query_fa_fn, chr, out_fn):
#   altref_twobit = os.path.splitext(altref_fa)[0] + ".2bit"
#   blat_db = altref_twobit + ":" + str(chr)
#   cmd = "%s -noHead %s %s %s"%(blat_path, blat_db, query_fa_fn, out_fn)
#   p = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True)
#   output, errors = p.communicate()

#   coords = [[0,0],[0,0], False]
#   blat_res = open(out_fn, 'rU')
#   hits = [False, False]
#   for line in blat_res.readlines():
#     line = line.strip()
#     linesplit = line.split()
#     res_id = linesplit[9]
#     if res_id.find("start") > -1:
#       if coords[0][0] < int(linesplit[0]):
#         coords[0][0] = int(linesplit[0])
#         coords[0][1] = int(linesplit[15])
#         hits[0] = True
#     elif res_id.find("end") > -1:
#       if coords[1][0] < int(linesplit[0]):
#         coords[1][0] = int(linesplit[0])
#         coords[1][1] = int(linesplit[16])
#         hits[1] = True
#   coords[2] = hits[0] and hits[1]
#   blat_res.close()
# #  os.remove(out_fn)
#   return coords

def calc_contig_complexity(seq, N=3, w=6):
  cmers = []
  for i in range(len(seq)):
    s = max(0,i-w)
    e = min(len(seq),i+w)
    cmer = count_nmers(seq[s:e], N)
    n = len(cmer)
    nmod = float(n) / float(e-s)
    cmers.append(round(nmod,2))
  cmers_mean = sum(cmers) / len(cmers)
  return cmers_mean, cmers

def count_nmers(seq, N):
  nmers = {}
  total_possible = len(seq) - 2
  for i in range(len(seq) - (N - 1)):
    mer = str(seq[i:i+N]).upper()
    if mer not in nmers: nmers[mer] = 0
    nmers[mer] += 1
  return nmers


def filter_by_feature(brkpts, query_region, keep_intron_vars):
    in_filter = False
    span_filter = False
    if not keep_intron_vars:
        in_vals, span_vals = check_intervals(brkpts, query_region)
    if in_vals[0]:
      if 'exon' not in in_vals[1]:
        in_filter = True
    else:
      in_filter = True
    if span_vals[0]:
      if 'exon' not in span_vals[1]:
        span_filter = True
    else:
      span_filter = True
    return in_filter, span_filter


def check_intervals(breakpts, query_region):
  in_values = [False, [], []]
  span_values = [False, [], []]
  bp_in_interval = False
  contains_interval = False
  in_features = []
  in_interval = None
  for bp in breakpts:
    for interval in query_region[4]:
      if (int(bp) >= (interval[1]-20)) and (int(bp) <= (interval[2]+20)):
        in_values[1].append(interval[4])
        in_values[0] = True
        in_values[2].append(interval)
      if (interval[2] <= max(breakpts)) and (interval[1] >= min(breakpts)):
        span_values[0] = True
        span_values[1].append(interval[4])
        span_values[2].append(interval)
  return in_values, span_values


def check_repeat_regions(coords, repeat_locs):
  start, end = coords
  seg_len = float(end-start)
  in_repeat = False
  rep_overlap = 0.0
  rep_coords = []
  filter_reps_edges = [False, False]
  for rloc in repeat_locs:
    rchr, rbp1, rbp2, rname = rloc
    if (rbp1 >= start and rbp1 <= end) or (rbp2 >= start and rbp2 <= end) or (rbp1 <= start and rbp2 >= end):
      in_repeat = True
      rep_overlap += float(min(rbp2,end)-max(rbp1,start))
      rep_coords.append((rbp1,rbp2))
      # Simple or low complexity seq repeat for filtering
      if rname.find(")n") > -1 or rname.find("_rich") > -1:
        if (rbp1<=start and rbp2>=start): filter_reps_edges[0] = True
        elif (rbp1<=end and rbp2>=end): filter_reps_edges[1] = True
#      if rep_overlap >= seg_len:
#        break
  roverlap = round( (float(min(rep_overlap, seg_len)) / float(seg_len))*100, 2)
  #break
  return in_repeat, roverlap, rep_coords, filter_reps_edges


def get_marker_fn(fn):
  return os.path.join(os.path.split(fn)[0],"."+os.path.basename(fn))


# def get_fastq_reads_old(fn, sv_reads):
#   read_len = 0
#   fq_recs = {}
#   f = open(fn,'r')
# #  fq_recs = list(SeqIO.parse(f,'fastq'))
#   for header,seq,qual in FastqFile(fn):
#     qname = header.lstrip("@")
#     if qname in sv_reads:
#       oseq, sc_seqs, clip_coords = sv_reads[qname]
#       cleaned_seq = seq
#       old_seq = oseq.seq
#       add = True
#       if str(cleaned_seq) != str(old_seq) and sc_seqs:
#         idx = old_seq.find(cleaned_seq)
#         trimmed_seq = ''
#         if idx == 0:
#           trimmed_seq = old_seq[len(cleaned_seq):len(old_seq)]
#         else: trimmed_seq = old_seq[0:idx]
#         sc_lens = 0
#         for sc_seq in sc_seqs:
#           sc_lens += len(sc_seq)
#           if trimmed_seq.find(sc_seq) > -1:
#             add = False
#         if len(cleaned_seq) == (len(old_seq) - sc_lens):
#           for sc_seq in sc_seqs:
#             if cleaned_seq.find(sc_seq) == -1:
#               # Don't add, just trimmed clipped portion.
#               add = False
# #    else: print qname, 'not in sv reads'
#     if add:
#       fr = fq_read(header, seq, qual)
#       read_len = max(read_len, len(fr.seq))
#       fq_recs[fr.id] = fr
#   return fq_recs, read_len


def load_kmers(fns, kmers):
  if not fns:
    return kmers

  fns = fns.split(",")
  for fn in fns:
    f = open(fn,'rU')
    for line in f.readlines():
      line = line.strip()
      mer, count = line.split()
      if not mer in kmers: kmers[mer] = 0
      kmers[mer] += int(count)
  return kmers


# Store all repeats in a dict by chrom
def setup_rmask_all(rmask_fn):
  logger = logging.getLogger('root')

  rmask = {}
  f = open(rmask_fn,'rU')
  flines = f.readlines()
  for line in flines:
    line = line.strip()
    rchr,rbp1,rbp2,rname = line.split("\t")[0:4]
    rchr = rchr.replace('chr','')
    if rchr not in rmask: rmask[rchr] = []
    rmask[rchr].append((rchr,int(rbp1),int(rbp2),rname))
  return rmask


# Stores repeats for a specific gene
def setup_rmask(gene_coords, ref_path, rmask_fn):
  logger = logging.getLogger('root')
  chrom, s, e, name, intervals = gene_coords
  mask_out_fn = os.path.join(ref_path,name+'_rep_mask.bed')
  marker_fn = get_marker_fn(mask_out_fn)

  rmask = []
  if not os.path.isfile(marker_fn):
    fout = open(mask_out_fn,'w')
    f = open(rmask_fn,'rU')
    flines = f.readlines()
    for line in flines:
      line = line.strip()
      rchr,rbp1,rbp2,rname = line.split("\t")[0:4]
      rchr = rchr.replace('chr','')
      if rchr == chrom:
        if int(rbp1) >= int(s) and int(rbp2) <= int(e):
          fout.write("\t".join([str(x) for x in [rchr,int(rbp1),int(rbp2),rname]])+"\n")
          rmask.append((rchr,int(rbp1),int(rbp2),rname))
    f.close()
    fout.close()

    cmd = 'touch %s'%marker_fn
    p = subprocess.Popen(cmd,stdout=subprocess.PIPE,stderr=subprocess.PIPE,shell=True)
    output, errors = p.communicate()
    logger.info('Completed writing repeat mask file %s, touching marker file %s'%(mask_out_fn,marker_fn))
  else:
    rep_f = open(mask_out_fn,'rU')
    rep_flines = rep_f.readlines()
    for line in rep_flines:
      line = line.strip()
      rchr,rbp1,rbp2,rname = line.split()
      rmask.append((rchr,int(rbp1),int(rbp2),rname))
  return rmask


def get_overlap_index_nomm(a,b):
  i = 0
  while a[i:] != b[:len(a[i:])]:
    i += 1
  return i


def seq_complexity(seq, N):
  nmers = {}
  total_possible = len(seq) - 2
  for i in range(len(seq) - (N - 1)):
    nmers[str(seq[i:i+N]).upper()] = True
  complexity = round((float(len(nmers))/float(total_possible))*100,4)
  return complexity


def get_overlap_index_mm(a, b):
  i = 0
  nmismatch = [0,0]
  match = False
  while i < len(a) and not match:
    nmismatch = [0,0]
    c = 0
    match_len = min(len(a[i:]), len(b[:len(a[i:])]))
    for aa,bb in zip(a[i:],b[:len(a[i:])]):
      if aa != bb:
        nmismatch[0] += 1
        nmismatch[1] += 1
      else:
        nmismatch[0] = 0
      if nmismatch[0] > 1 or nmismatch[1] > 3:
        break
      c += 1
#    print c, match_len, i, c== match_len, a[i:], b[:len(a[i:])]
    if c == match_len:
      match = True
    i += 1
  return i-1


def get_read_kmers(seq,l,skmers):
  kmers = []
  i = 0
  while (i+l) <= len(seq):
    k = seq[i:i+l]
    #if k in skmers:
    kmers.append(k)
    i += 1
  return list(set(kmers)&set(skmers))


def get_overlap_index(a,b):
  i = 0
  nmismatch = 10
  while nmismatch > 1:
    nmismatch = 0
    for aa,bb in zip(a[i:],b[:len(a[i:])]):
      if aa != bb: nmismatch += 1
    i += 1
#  while a[i:] != b[:len(a[i:])]:
#    i += 1
  return i-1
