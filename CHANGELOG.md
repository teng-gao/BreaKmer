# Change Log
All notable changes to this project will be documented in this file.

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
