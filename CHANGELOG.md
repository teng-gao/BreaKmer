# Change Log
All notable changes to this project will be documented in this file.

## [0.0.4-beta.1] - 2016-06-07
    - Added an initial sample from the published, publicly available bam files from Abo et al. (2014) manuscript to the example_data directory.
    - Added a check for the index file for the input sample bam file.
    - Created new function for profiling the sample bam input file sequence coverage. This function will aid in the generation of a bed file when none exists or just for general purpose metrics of a bam file.

## [0.0.4-beta.1] - 2016-05-31
    - Bug fix, reimplemented function to remove target output directory if there is no target output results.

## [0.0.4-beta] - 2016-05-26
    - Breakpoint read coverages were reimplemented and now in the output file.
    - Basic filters were reimplemented.
    - Update README documentation.
    - v0.0.4-beta released!

## [0.0.4-beta] - 2016-05-24
    - Many small modifications and code rearrangements. Caller, realignment, and results modules were created to contain code for their respective functions. Realignment functionality and calling functionality has been modularized. Running blat function has been moved to utils module. Results have also been placed in a module and simplified.
    - Removed "nkmers" and "repeat_overlap" output columns.
    - Currently suspended the summary output file.
    - Decoupled the filtering process from the calling process. The filters need to be reimplemented.
    - Merged the two Contig classes into a single class.
    - Removed 'chr' concatenation to output chromosome breakpoint strings.

## [0.0.4-beta] - 2016-05-19
    - Removed parameter 'var_filter'. All variants will be reported and the user can filter beyond that.
    - Major bug fix in function to determine the blat hits that make up a contig sequence. This should allow for more sensitive calling where a contig sequence contains multiple breakpoint events.
    - Gene annotations added back as a required parameter, along with keep_intron_vars parameter. These will be transitioned in the next version.
    - Check binaries function added as an initial check before analysis begins.
    - Added estimation of paired-end insert size distribution based on a sample of well-mapped reads.

## [0.0.4-beta] - 2016-05-17
    - Added print statement for the blat server hostname and port used when the server is started.
    - Minor bug fix when starting blat server and checking existence of a blat server. The blat server is started from the directory that the reference file is located and ONLY the reference file name is used, not the full path. When the gfClient is run against the server the absolute path to the reference 2bit file is used.
    - Added a buffer size parameter option to indicate the number of base pairs to add to both sides of a target region for extraction aligned reads. The default is set to 100bp.
    - Moved sv_assembly module to assembly/assembler and assembly/contig modules. Also moved olc module to assembly/olc.
    - Moved sv_caller module into breakmer module directory.

## [0.0.4-beta] - 2016-05-12
    - Bug fixes for starting the blat server. It is now more robust in handling both .fasta and .fa file extensions, and reports if the gfServer aborts during the startup process.
    - Added faToTwoBit binary to the bin directory for direct usage. Note that the permissions should be set to executable.
    - Both start_blat_server and prepare_reference_data functions were tested successfully.

## [0.0.4-beta] - 2016-05-11
    - Begin transition from camel casing to underscores for local variable names.
    - Moved sv_processor.py code to processor/analysis.py module.
    - Moved functions to start blat server into utils module.
    - Adding parameter requirement checks while initializing the parameter object.
    - Adding functional run mode options - prepare_reference_data, run, start_blat_server.
    - Removing support for alternative reference sequences.
    - Added a bin directory containing cutadapt v1.9.1, jellyfish v1.1.11, blat, gfClient, gfServer

## [0.0.4-beta] - 2016-05-10
    - Changed from OptParser to argparser
    - Added initial module structure with breakmer directory and params.py and utils.py modules.
    - Added parameter manager class to parse and track input parameters and configuration file values.
    - Simplified and cleaned main script code in breakmer.py.
