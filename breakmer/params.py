#! /usr/bin/python
# -*- coding: utf-8 -*-

import os
import sys
import subprocess
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

        # Log all parameters passed in, warn for poor paths
        for param_key, param_value in self.opts.items():
            utils.log(self.logging_name, 'info', '%s = %s' % (param_key, param_value))
        # Log parameter overwritten by configuration file input values.
        for log_msg in log_msgs:
            utils.log(self.logging_name, log_msg[0], log_msg[1])

        self.set_targets()
        self.paths['ref_data'] = os.path.abspath(os.path.normpath(self.opts['reference_data_dir']))  # Path to target reference sequence fast files.
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

        # self.check_binaries()  # Check if Jellyfish and Cutadapt work.
        # self.filter = resultfilter.ResultFilter(self.get_param('filterList'), self)  # Instantiate the filter class.
        # self.set_insertsize_thresh()  # Set the expected insert size threshold from the properly mapped read

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

        # Check that the required parameters are set.
        required = ['analysis_name',
                    'targets_bed_file',
                    'sample_bam_file',
                    'analysis_dir',
                    'reference_data_dir',
                    'cutadapt_config_file',
                    'reference_fasta',
                    'gene_annotation_file']
        if self.fnc_cmd == 'prepare_reference_data':
            required = ['reference_data_dir', 'reference_fasta', 'targets_bed_file']

        for req in required:
            self.get_param(req, True)
        # Return log messages that need to be back logged once the logger is setup
        return log_msgs

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
            if value == '':
                value = None
        if required and value is None:
            utils.log(self.logging_name, 'error', 'Missing required parameter %s, exiting.' % key)
            sys.exit(1)
        sys.stderr.write('Warning: %s is set to None\n' % key)
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
