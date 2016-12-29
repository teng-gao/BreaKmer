#! /usr/bin/python
# -*- coding: utf-8 -*-

import os
import sys
import shutil
import pysam
import breakmer.utils as utils

__author__ = "Ryan Abo"
__copyright__ = "Copyright 2015, Ryan Abo"
__email__ = "ryanabo@gmail.com"
__license__ = "MIT"


class ParamManager(object):

    """ParamManager class stores all the input specifications provided to the program to run. These include
    file paths, thresholds, directories, etc...

    Attributes:
        opts (dict):                   Containing parameter options and input values as key-values.
        gene_annotations (Annotation): Tracks the annotation information.
        targets (dict):                Target region coordinates, key-values.
        paths (dict):                  Dictionary containing the top level directories for the analysis output.
        logging_name (str):            Logging string name object for logging messages.
    """

    def __init__(self, arguments):

        """Initialize the ParamManager class.

        Args:
            arguments (dict):   The argparse dictionary object from the command line parameters.
        Returns:
            None
        Raises:
            None
        """

        self.logging_name = 'breakmer.params'
        self.opts = {}
        self.gene_annotations = None
        self.filter = None
        self.targets = {}
        self.paths = {}
        self.fnc_cmd = arguments.fnc_cmd
        self.input_args = arguments
        self.set_params()

    def set_params(self):

        """Organize and format all input parameters into class variables to access
        later. Specific instances of parameters are checked and set. All the parameters that are
        set are logged. The target objects are set along with the paths.

        Args:
            None
        Returns:
            None
        Raises:
            None
        """

        log_msgs = self.parse_opts()  # Parse the config file and command line parameters into the self.opts dictionary.
        utils.setup_logger(self.get_param('analysis_dir', True), 'breakmer')  # Create logging object.
        # Log the parameter setup after the logging object is created.
        utils.log(self.logging_name, 'info', 'Setting up parameters')

        # Setup the logger first before checking required params - log if there is a missing parameter.
        self.check_required_params()

        # Log all parameters passed in, warn for poor paths
        for param_key, param_value in self.opts.items():
            utils.log(self.logging_name, 'info', '%s = %s' % (param_key, param_value))

        # Log parameter overwritten by configuration file input values.
        for log_msg in log_msgs:
            utils.log(self.logging_name, log_msg[0], log_msg[1])

        if self.get_param('targets_bed_file') is not None:
            self.set_targets()
            self.gene_annotations = Anno()
            self.gene_annotations.add_genes(self.get_param('gene_annotation_file'))

        if 'reference_data_dir' in self.opts:
            self.paths['ref_data'] = os.path.abspath(os.path.normpath(self.opts['reference_data_dir']))  # Path to target reference sequence fast files.
        if 'reference_fasta' in self.opts:
            self.set_param('reference_fasta_dir', os.path.split(self.opts['reference_fasta'])[0])  # Path to genome fasta file.

        # Setup directories
        self.paths['analysis'] = os.path.abspath(os.path.normpath(self.opts['analysis_dir']))
        self.paths['output'] = os.path.join(self.paths['analysis'], 'output')
        if 'targets_dir' in self.opts:
            self.paths['targets'] = os.path.abspath(os.path.normpath(self.opts['targets_dir']))
        else:
            self.paths['targets'] = os.path.join(self.paths['analysis'], 'targets')

        # Create all the paths.
        for path in self.paths:
            utils.log(self.logging_name, 'info', 'Creating %s directory (%s)' % (path, self.paths[path]))
            if not os.path.exists(self.paths[path]):
                os.makedirs(self.paths[path])

        # If starting the blat server then return.
        if self.fnc_cmd == 'start_blat_server':
            utils.log(self.logging_name, 'info', 'Starting the blat server.')
            return
        if self.fnc_cmd == 'profile_data':
            return

        self.check_binaries()  # Check if Jellyfish and Cutadapt work.
        # self.filter = resultfilter.ResultFilter(self.get_param('filterList'), self)  # Instantiate the filter class.
        self.set_insertsize_thresh()  # Set the expected insert size threshold from the properly mapped read

    def parse_opts(self):

        """Formats input parameters into self.opts dictionary. It first parses the configuration file and stores the key, values in the self.opts dictionary.

        It will exit with an error if the configuration file does not have lines in the proper format (i.e., key=value).
        It will also iterate through the command line paramaters and store the keys and values in the opts dictionary.
        A final check is performed for the required parameters depending on the parameters that have been passed.
        
        Args:
            None
        Returns:
            None
        Raises:
            None
        """

        if self.input_args.config_fn is not None:
            for line in open(self.input_args.config_fn, 'rU'):
                line = line.strip()
                if line == '' or line.find('#') > -1:  # Allow for blank lines and comments
                    continue
                linesplit = line.split("=")
                if len(linesplit) == 1:  # Make sure the lines in the configuration file are set properly.
                    err_msg = 'Configuration file line', line, ' is not set correctly. Exiting.'
                    print err_msg
                    # utils.log(self.logging_name, 'error', err_msg)
                    sys.exit(1)
                else:
                    key, value = linesplit
                    self.set_param(key, value)  # Store key-value in opts dictionary.

        log_msgs = []
        # Store all the arguments into the self.opts dictionary.
        for opt in vars(self.input_args):
            if (self.get_param(opt) is not None) and (vars(self.input_args)[opt] is None):
                # Save log message
                # utils.log(self.logging_name, 'info', 'Parameter %s is set in config file and not on the command line. Using config file value %s.' % (opt, self.get_param(opt)))
                log_msgs.append(('info', 'Parameter %s is None in input parameters and set in configuration file. Using configuration file value %s for analysis.' % (opt, self.get_param(opt))))
            else:
                self.set_param(opt, vars(self.input_args)[opt])
        return log_msgs

    def check_required_params(self):

        '''Iterate over the required parameter values and check if they exist.

        Args:
            None
        Returns:
            None
        '''

        # Check that the required parameters are set.
        required = ['analysis_name',
                    'targets_bed_file',
                    'sample_bam_file',
                    'analysis_dir',
                    'reference_data_dir',
                    'cutadapt_config_file',
                    'reference_fasta']
#                    'gene_annotation_file']
        if self.fnc_cmd == 'prepare_reference_data':
            required = ['reference_data_dir', 'reference_fasta', 'targets_bed_file']

        if self.fnc_cmd == 'profile_data':
            required = ['analysis_dir', 'sample_bam_file']

        for req in required:
            value = self.get_param(req, True)
            if req == 'sample_bam_file':
                # Check that the index is present
                bam_file = pysam.AlignmentFile(value)
                if not bam_file.check_index():
                    utils.log(self.logging_name, 'debug', 'Missing an index file for the input bam file. Exiting.')
                    sys.exit(1)

    def check_binaries(self):
        # Binaries required for operation: blat, gfserver, gfclient, fatotwobit, cutadapt, jellyfish
        binaries = ('blat',
                    'gfserver',
                    'gfclient',
                    'fatotwobit',
                    'cutadapt',
                    'jellyfish')

        for binary_name in binaries:
            binary_path = self.get_param(binary_name)
            binary_check = None
            if binary_path is not None:
                binary_check = utils.which(binary_path)
            else:
                binary_check = utils.which(binary_name)
                self.set_param(binary_name, binary_check)
            if not binary_check:
                print 'Missing path/executable for', binary_name
                sys.exit(2)
            utils.log(self.logging_name, 'debug', '%s path = %s' % (binary_name, binary_check))

        utils.log(self.logging_name, 'debug', 'All the required binaries have been check successfully!')

        # Test cutadapt and jellyfish binaries
        test_dir = os.path.join(self.paths['analysis'], 'binary_test')
        test_fq = os.path.join(test_dir, 'test.fq')
        if not os.path.exists(test_dir):
            os.makedirs(test_dir)

        fq_file = open(test_fq, 'w')
        fq_file.write("@H91H9ADXX140327:1:2102:19465:23489/2\nCACCCCCACTGAAAAAGATGAGTATGCCTGCCGTGTGAACCATGTGACTTTACAATCTGCATATTGGGATTGTCAGGGAATGTTCTTAAAGATC\n+\n69EEEFBAFBFABCCFFBEFFFDDEEHHDGH@FEFEFCAGGCDEEEBGEEBCGBCCGDFGCBBECFFEBDCDCEDEEEAABCCAEC@>>BB?@C\n@H91H9ADXX140327:2:2212:12198:89759/2\nTCTTGTACTACACTGAATTCACCCCCACTGAAAAAGATGAGTATGCCTGCCGTGTGAACCATGTGACTTTACAATCTGCATATTGGGATTGTCAGGGA\n+\nA@C>C;?AB@BBACDBCAABBDDCDDCDEFCDDDDEBBFCEABCGDBDEEF>@GBGCEDGEDGCGFECAACFEGDFFGFECB@DFGCBABFAECEB?=")
        fq_file.close()

        clean_fq, return_code = utils.test_cutadapt(test_fq, self.get_param('cutadapt'), self.get_param('cutadapt_config_file'))
        if clean_fq:
            utils.log(self.logging_name, 'info', 'Test cutadapt ran successfully')
            jfish_prgm, return_code = utils.test_jellyfish(self.opts['jellyfish'], clean_fq, test_dir)
            if return_code != 0:
                utils.log(self.logging_name, 'error', '%s unable to run successfully, exit code %s. Check installation and correct version.' % (jfish_prgm, str(return_code)))
                sys.exit(2)
            else:
                utils.log(self.logging_name, 'info', 'Test jellyfish ran successfully')
        else:
            utils.log(self.logging_name, 'error', 'Cutadapt failed to run, exit code %s. Check installation and version.' % str(return_code))
            sys.exit(2)
        shutil.rmtree(test_dir)

    def set_targets(self):

        """Parse the targets bed file and store them in a dictionary. Limit to a gene
        list if input.

        A list of genes can be passed in by the user to limit the analysis. This will
        limit which targets are stored in the dictionary as the target bed file is parsed.
        The target bed file is a tab-delimited text file that should have at minimum,
        four columns (chromosome, start, end, name) with an optional fourth column
        containing a coding feature (i.e., exon or intron). Each row is either a tiled
        region with sequencing coverage or it is just a region to analyze by BreaKmer.
        
        The name can be applied to multiple rows, and if multiple tiled regions are input
        with the same name they are aggregated together under the same key.
        
        Store the target information in the self.target dictionary with the name as the key
        and a list of tuples of interval genomic locations as the values.

        self.target[gene_name] = [(chrom, start_bp, end_bp, name, feature),...]

        Args:
            None
        Returns:
            None
        Raises:
            None
        """

        # Get the gene list file path if it exists.
        gene_list = self.get_param('gene_list')
        region_list = None
        if gene_list:
            region_list = []
            # Each line contains a gene name.
            for line in open(gene_list, 'r'):
                region_list.append(line.strip().upper())

        utils.log(self.logging_name, 'info', 'Parsing target list')

        # TODO: Check to make sure there aren't duplicate genes.
        with open(self.get_param('targets_bed_file'), 'rU') as target_bed_file:
            for target in target_bed_file:
                # Each target is formatted like a bed, chr bp1 bp2 name
                target = target.strip()
                targetsplit = target.split()
                chrm, bp1, bp2, name = targetsplit[0:4]
                if region_list:
                    if name.upper() not in region_list:
                        continue
                # Allow a fifth column containing indication of what type of region it is.
                # Typically exon/intron designation. This will be deprecated.
                feature = None if len(targetsplit) <= 4 else targetsplit[4]
                self.targets.setdefault(name.upper(), [])
                self.targets[name.upper()].append((chrm, int(bp1), int(bp2), name, feature))
        utils.log(self.logging_name, 'info', '%d targets' % len(self.targets))

    def set_insertsize_thresh(self):

        '''Store the insert sizes for a small number of "properly mapped" reads
        and determine an upperbound cutoff to use to determine discordantly mapped read
        pairs.

        Args:
            None
        Returns:
            None
        Raises:
            None
        '''

        nsample_reads = 100000
        bam_file = pysam.Samfile(self.get_param('sample_bam_file'), 'rb')
        test_reads = bam_file.fetch()
        insert_sizes = []
        read_iter = 0
        for read in test_reads:
            if read.is_duplicate or read.mapq == 0:
                continue
            proper_map = read.flag == 83 or read.flag == 99
            if read.is_read1 and proper_map:  # Sample the read and store the insert size to its partner.
                read_iter += 1
                insert_sizes.append(abs(read.tlen))
                if 'read_len' not in self.opts:  # Store the read length if it is not already stored.
                    self.set_param('read_len', read.rlen)
            if read_iter == nsample_reads:
                break
        insertsize_median = utils.median(insert_sizes)
        insertsize_sd = utils.stddev(utils.remove_outliers(insert_sizes))  # Calculate the standard deviation of the sample read pairs insert sizes.
        self.set_param('insertsize_thresh', insertsize_median + (5 * insertsize_sd))  # Set the threshold to be median + 5 standard deviations.

    def get_param(self, key, required=False):
        """Get the parameter value in the self.opts dictionary.
        If the parameer is required to be availale, then exit the program
        and throw an error.

        Args:
            key (str): The key in the opts dictionary to access the parameter value.
                       required: Boolean value to indicate if the key should be required to
                       be in the dictionary or not.
        Returns:
            value (int, str, boolean): The value of the parameter if it is found. If the parameter is
                                       required and not found the program will exit with error. If the parameter is
                                       not required and not found, it will return None.
        Raises:
            None
        """

        value = None
        if key in self.opts:
            value = self.opts[key]
        elif required:
            utils.log(self.logging_name, 'error', 'Missing required parameter %s, exiting.' % key)
            sys.exit(1)
        return value

    def set_param(self, key, value):
        """Set the parameter value in the self.opts dict.
        Args:
            key (str):              Dictionary key
            value (int/str/boolean):  Value to store
        Returns:
            None
        Raises:
            None
        """

        self.opts[key] = value


class Anno(object):

    '''
    '''

    def __init__(self):

        '''
        '''

        self.genes = {}
        self.logging_name = 'breakmer.params.anno'

    def add_genes(self, gene_fn):

        '''
        '''

        utils.log(self.logging_name, 'info', 'Adding gene annotations from %s' % gene_fn)
        gene_f = open(gene_fn, 'r')
        gene_flines = gene_f.readlines()
        for line in gene_flines[1:]:
            line = line.strip()
            linesplit = line.split()
            chrom = linesplit[2]
            start = int(linesplit[4])
            end = int(linesplit[5])
            geneid = linesplit[12]
            if geneid in self.genes:
                if start <= self.genes[geneid][1] and end >= self.genes[geneid][2]:
                    self.genes[geneid] = [chrom, start, end]
            else:
                self.genes[geneid] = [chrom, start, end]
        gene_f.close()

    def add_regions(self, regions_bed_fn):

        '''
        '''

        region_f = open(regions_bed_fn, 'rU')
        region_lines = region_f.readlines()
        for line in region_lines:
            line = line.strip()
            chrom, start, end, name = line.split()
            if name not in self.genes:
                self.genes[name] = [chrom, int(start), int(end)]
        utils.log(self.logging_name, 'info', 'Adding in %d other target regions' % len(region_lines))

    def set_gene(self, chrom, pos):

        '''
        '''

        ann_genes = []
        if chrom.find('chr') == -1:
            chrom = 'chr' + str(chrom)

        for g in self.genes:
            gs = self.genes[g][1]
            ge = self.genes[g][2]
            if chrom == self.genes[g][0]:
                if len(pos) == 1:
                    if int(pos[0]) >= gs and int(pos[0]) <= ge:
                        ann_genes.append(g)
                        break
            else:
                # Find genes between pos1 and pos2
                if (int(pos[0]) >= gs and int(pos[0]) <= ge) or (int(pos[1]) >= gs and int(pos[1]) <= ge):
                    ann_genes.append(g)
        if len(ann_genes) == 0:
            ann_genes = ['intergenic']
        return ",".join(ann_genes)
