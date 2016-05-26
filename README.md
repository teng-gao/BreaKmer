BreaKmer
========

A method to identify genomic structural variation in target regions/genes from reference-aligned high-throughput sequence data. It uses a “kmer” strategy to assemble misaligned sequence reads for predicting insertions, deletions, inversions, tandem duplications, and translocations at base-pair resolution.

For detailed documentation related to the usage of BreaKmer see the Documentation below or visit this website: http://ryanabo.github.io/BreaKmer

A google groups forum page is also available for posting issues: https://groups.google.com/forum/#!forum/breakmer

News
====

#### May 26, 2016
Preparing to deploy a new version release (0.0.4-beta). I am looking for beta testers to run their data through it or sharing some data for me to test.

More work will follow for benchmarking and testing, but an initial test on a publically available sample [SRR1304138]() from the published dataset with a BCR-ABL1 translocation had a run time of < 1 minute on a Mac when only analyzing BCR and ABL1 genes with ~300x depths.

#### May 24, 2016
Major updates to structure of the code, it is much more modularized and succinct. Improvements should be much easier now. I have merged the different output types into a single output file ("_svs.out") and taken out the "nkmers" and "repeat_overlap" columns. The contig IDs are now labelled with the target region name and there is a "sv_subtype" included in the subtype for the different type of rearrangement categories. There is also no concatentation of 'chr' onto the chromosome names for indels or rearrangements as before. There were a few other minor algorithmic changes to the processing of realignment calls and how the discordant reads are being extracted and matched to the variant calls. These will continue to be modified, as the ultimate goal is to call variants with splitread and discordant read evidence jointly rather than in serial. Currently, the code is stable and works on a few test cases; however, the basic filters that were in place need to be reimplemented.

#### May 19, 2016
All BreaKmer functionality is currently working on a single test sample (a translocation) with updates and changes. Work is progressing to clean the code and organize into modules and remove PEP8 incompatibilities.

#### May 17, 2016
Developments continue to scrub and improve old, "stable" code with new features. Modifications and testing for the "run" function have been the bulk of work lately.

#### May 12, 2016
Continuing to aim for a "runnable" and updated code base in the master branch. Both start_blat_server and prepare_reference_data functions were tested successfully. Once the "run" function works I will begin to develop in the development branch and merge back to master branch. I have also create a google groups forum page for posting issues (https://groups.google.com/forum/#!forum/breakmer).

#### May 11, 2016
The focus thus far has been on cleaning up the main script to have cleaner and more succinct, modularized code. The parameters are packaged into the params module and utility functions pushed to utils module. Algorithmically nothing will change in this next release. The main features that will be implemented are allowing BreaKmer to run in separate "function" modes - prepare reference data, start blat server, run, etc...

#### May 10, 2016
Due to issues in the latest release of BreaKmer, I am taking a step back and reverting to the original "stable" version (v0.0.3) to redevelop new features in a more robust manner. I will aim to document the regular progress and changes as I go and attempt to address the numerous bugs that have come about in the past year. Please submit any feature requests, questions, or comments on this page and I will do my best to address it. If I don't know about an issue, I can't fix it. For those using the development version (old master), you can find that code in the beta_version branch. Onwards to version 0.0.4!

Documentation
=============

Installation
----------

The following are required for installation:
- [Python 2.7](https://www.python.org/download/releases/2.7)

Download the python scripts and run the command:
```
python setup.py install
```
Use appropriate commands for installing locally:
```
python setup.py install --user
```

Dependencies that need to be downloaded and installed if not already:
- [Biopython 1.62](http://biopython.org/wiki/Main_Page)
- [Pysam 0.6](https://code.google.com/p/pysam/)

BreaKmer has currently been installed and tested on:
- 64bit OSX
- 64bit linux using CentOS release 5.5 and python2.7.2 
- 64bit linux using Ubuntu release 14.04 and python2.7.6

Usage
---------

####### List the available command line parameters.
```
python <PATH TO BREAKMER DIR>/breakmer.py -h
```

###### Setup the reference data prior to analysis.
  - This will extract the reference sequence data into the specified reference directory for each of the target regions defined in the input bed file.
```
python <PATH TO BREAKMER DIR>/breakmer.py setup_reference_data -c <PATH TO BREAKMER CONFIGURATION FILE>
```

###### Start blat server
  - Note that this is a convenience function. If the server is not started before using the 'run' function, a new server will be initiated prior to analysis.

Start the server on the 'localhost' with a random port number. The port that is choosen will be printed in the standard output, which should be specified as an input parameter when using the 'run' function.
```
python <PATH TO BREAKMER DIR>/breakmer.py start_blat_server -c <PATH TO BREAKMER CONFIGURATION FILE>
```

Start the server on a specified host with a random port number.
```
python <PATH TO BREAKMER DIR>/breakmer.py start_blat_server -c <PATH TO BREAKMER CONFIGURATION FILE> --hostname <HOSTNAME>
```

Start the server on a specified host and specified port number.
```
python <PATH TO BREAKMER DIR>/breakmer.py start_blat_server -c <PATH TO BREAKMER CONFIGURATION FILE> --hostname <HOSTNAME> -p <PORT NUMBER>
```

###### Run analysis
  - The reference files and blat server do not need to be setup before executing this command, but they may save time if the analysis needs to be rerun multiple times. It will first check to see if the reference sequence files for the targets have been extracted.
  - If the blat port and hostname are not specified in the command, it will start a new blat server.

Analyze all the target genes specified in the targets bed file.
```
python <PATH TO BREAKMER DIR>/breakmer.py -c <PATH TO BREAKMER CONFIGURATION FILE>
```

Analyze a subset of genes specified in a file.
```
python <PATH TO BREAKMER DIR>/breakmer.py -c <PATH TO BREAKMER CONFIGURATION FILE> -g <PATH TO TEXT FILE CONTAINING A LIST OF TARGET REGIONS TO ANALYZE>
```

Requirements
---------

### Programs
- BLAT standalone and server binaries ([blat, gfServer, gfClient, faToTwoBit](http://hgdownload.cse.ucsc.edu/admin/exe/)).
  - Re-alignment to reference sequence.
  - Versions tested :
    - standalone BLAT v35x1
    - gfServer v35x1
    - gfClient v35x1
- [Cutadapt](https://code.google.com/p/cutadapt/)
  - Trims adapter sequence from aligned reads.
  - v1.5 tested
- [Jellyfish](http://www.cbcb.umd.edu/software/jellyfish)
  - Generating kmers.
  - [v1.1.11](http://www.cbcb.umd.edu/software/jellyfish/jellyfish-1.1.11.tar.gz) tested
  - [v2.1.3](http://www.genome.umd.edu/jellyfish.html) tested

When these programs are installed, the paths to the binaries can be either be specified in the BreaKmer configuration file or put in the path (e.g. export PATH=$PATH:/path/to/binary).

### Sequence data
- Sample bam file
  - BreaKmer requires sequence reads that have been aligned to a reference sequence in a binary alignment format (BAM). The alignment program must soft-clip or trim reads that can be partially aligned to the reference sequence (e.g., bwa, bowtie, novoalign, mosaik). The partially aligned reads and unmapped reads with mapped mates (paired-end data) are used to build contigs with potential SV. This bam file is required as input in the configuration file as the "sample_bam_file".
    - The BAM files need to be sorted and indexed, with the indexed files in the same directory as the BAM file.
    - There are a number of aligners that soft-clip partially aligned sequences, bwa and Bowtie are two well-known tools:
      - [bwa](http://bio-bwa.sourceforge.net/)
      - [Bowtie](http://bowtie-bio.sourceforge.net/index.shtml)
    - It is also useful to mark duplicate reads prior to using BreaKmer, as these will inflate read counts for identified variants.
      - [Picard MarkDuplicates](http://picard.sourceforge.net/command-line-overview.shtml#MarkDuplicates)
- Normal bam file (Currently untested)
  - A matched normal with similarly targeted sequencing data as the sample can be input to help filter germline events. The sequences in the normal bam file are processed for each target in a similar manner to the sample sequences and used to further filter our kmers. This is an optional input in the configuration file as "normal_bam_file".

### Reference data
- Reference fasta file 
  - BreaKmer makes use of reference sequence data throughout the program. The genomic reference sequence used to align the short sequence reads is required as an input in the configuration file as the "reference_fasta".
  - Format: single fasta file with chromosome number/id as names (i.e. '>1', '>2', '>3')
  - Hg19 fasta files can be downloaded from UCSC Genome Browser(http://hgdownload.soe.ucsc.edu/goldenPath/hg19/bigZips/)
  - This file should be placed in a writeable directory. A 2bit file will be generated from this file to start the blat server.
- Reference genome annotation
  - This file is used by BreaKmer to annotate the locations of the breakpoints identified and the genes involved. This is required as input in the configuration file as the "gene_annotation_file".
  - Format: tab delimited file containing a row for each RefSeq transcript with multiple columns describing the coding coordinates of the transcript.
  - Hg19 RefSeq annotation file can be downloaded from UCSC Genome Browser [Table Browser](http://genome.ucsc.edu/cgi-bin/hgTables?command=start) using assembly: Feb. 2009 (GRCh37/hg19), group: Genes and Gene Predictions, track: RefSeq Genes, table: refGene

Configuration file
------------

- A template configuration file with all the options is available, breakmer.cfg.
- Below lists and describes the parameters that can be used in the configuration file. Do not keep commented text (i.e., #Required parameters) in the configuration file.
- Note that the paths to the six required program binaries (Cutadapt, Jellyfish, blat, gfServer, and gfClient) can be set in the configuration file
  or these binaries can be included in the users path (e.g., for linux users: export PATH=$PATH:/path/to/binary).
- Use full paths (e.g., /home/bob and not ~/)
```
# Required parameters
analysis_name=<sample_id, string value that all the output files will contain>
targets_bed_file=<path to bed file containing locations of target regions>
sample_bam_file=<path to sorted and indexed sample bam file>
analysis_dir=<path to analysis directory, location where analysis and final output files will be located>
reference_data_dir=<path to where reference files will/are stored for each of the targeted genes analyzed> 
cutadapt=<path to cutadapt binary v1.5, i.e. /usr/bin/cutadapt-1.5/bin/cutadapt> 
cutadapt_config_file=<path to cutadapt configuration file> 
jellyfish=<path to Jellyfish binary, i.e. /usr/bin/jellyfish>
blat=<path to blat binary, i.e. /usr/bin/blat>
gfserver=<path to gfServer binary, i.e. /usr/bin/gfServer>
gfclient=<path to gfClient binary, i.e. /usr/bin/gfClient>
fatotwobit=<path to faToTwoBit binary, i.e. /usr/bin/faToTwoBit>
reference_fasta=<path to whole genome reference fasta file, one file with all records>
gene_annotation_file=<path to gene annotation file, e.g., ucsc_hg19_refgene.txt>
kmer_size=15

# Optional parameters
normal_bam_file=<path to normal bam file, can be used to filter germline events with matched-normal sample>
```

Input file formats
-----------

- targets_bed_file = tab-delimited file with columns: chr, start, end, region_name, feature_name
  - Each "target" region can contain multiple subregions that are annotated by feature type (i.e., intron, exon). These feature types are used in the filtering steps with certain parameters. Note that the minimum coordinate and maximum coordinate +/- 200 bp are used as the boundaries for each region (e.g., 45003745-200,45008540+200 for B2M).
```
15      45003745        45003811        B2M     exon
15      45007621        45007922        B2M     exon
15      45008527        45008540        B2M     exon
```

- cutadapt_config_file = each row corresponds to a parameter for cutadapt (see cutadapt.cfg example file or cutadapt documentation)
  - The file provided is intended for data generated using the paired-end Illumina TruSeq library.
  - Many of the Illumina library sequences have been annotated [elsewhere](https://wikis.utexas.edu/display/GSAF/Illumina+-+all+flavors).
- reference_fasta = genome reference fasta formatted file containing all the chromosome reference sequences that were used to initially align the data.
  - This should be a single file containing all the sequences.
- gene_annotation = Annotation file containing the location of reference genes. These can be downloaded from UCSC Genome Browser, (i.e., [hg19 refGene table](http://hgdownload.cse.ucsc.edu/goldenPath/hg19/database/))
```
#bin    name    chrom   strand  txStart txEnd   cdsStart        cdsEnd  exonCount       exonStarts      exonEnds        score   name2   cdsStartStat    cdsEndStat      exonFrames
0       NM_032291       chr1    +       66999824        67210768        67000041        67208778        25      66999824,67091529,67098752,67101626,67105459,67108492,67109226,67126195,67
133212,67136677,67137626,67138963,67142686,67145360,67147551,67154830,67155872,67161116,67184976,67194946,67199430,67205017,67206340,67206954,67208755, 67000051,67091593,67098777,6710169
8,67105516,67108547,67109402,67126207,67133224,67136702,67137678,67139049,67142779,67145435,67148052,67154958,67155999,67161176,67185088,67195102,67199563,67205220,67206405,67207119,6721
0768,   0       SGIP1   cmpl    cmpl    0,1,2,0,0,0,1,0,0,0,1,2,1,1,1,1,0,1,1,2,2,0,2,1,1,
1       NM_032785       chr1    -       48998526        50489626        48999844        50489468        14      48998526,49000561,49005313,49052675,49056504,49100164,49119008,49128823,49
332862,49511255,49711441,50162984,50317067,50489434,    48999965,49000588,49005410,49052838,49056657,49100276,49119123,49128913,49332902,49511472,49711536,50163109,50317190,50489626,  0AGBL4    cmpl    cmpl    2,2,1,0,0,2,1,1,0,2,0,1,1,0,
1       NM_018090       chr1    +       16767166        16786584        16767256        16785385        8       16767166,16770126,16774364,16774554,16775587,16778332,16782312,16785336, 16767348,16770227,16774469,16774636,16775696,16778510,16782388,16786584, 0       NECAP2  cmpl    cmpl    0,2,1,1,2,0,1,2,
1       NM_052998       chr1    +       33546713        33585995        33547850        33585783        12      33546713,33546988,33547201,33547778,33549554,33557650,33558882,33560148,33
562307,33563667,33583502,33585644,      33546895,33547109,33547413,33547955,33549728,33557823,33559017,33560314,33562470,33563780,33583717,33585995,    0       ADC     cmpl    cmpl    -1
,-1,-1,0,0,0,2,2,0,1,0,2,
...
```

Output files and formats
-----------

- While the program is running a log file analysis_dir/log.txt will continually be updated with information regarding the status of the analysis.

### Logging
  - While the program is running a log file will (\<analysis\_directory\>/log.txt) continually be updated with information regarding the status of the analysis.

### Target data
  - For each region analyzed, a directory with the name of the region, as specified in the targets bed file, is created in a 'targets' directory (\<analysis\_dir\>/targets).
  - Each target directory contains 'data', 'contigs', and 'kmers' directories. 
    - data - Contains the extracted reference-aligned reads from this target region as well as the kmers created from these extracted reads.
    - contigs - Contains a directory for each contig that was created, which contains the contig sequence in fasta format, and the reads used to assemble the contig in fastq format. The BLAT results are also stored with formatted output if a SV was called.
    - kmers - Contains the file with all the sample-only kmer sequences and how many reads in which it was contained as well as a file with the kmers and read ids that were used in assembling each contig.

### Final output 
  - When the program completes, the final output files are directed into a directory labeled 'output' within the specified analysis directory. 
  - Output for all structural variants are output in a single tab-delimited file, labeled \<analysis_name\>\_svs.out
    - The columns are:
       - genes - Gene/Target names
       - target_breakpoints - Target genomic breakpoints (chr:pos) (I/D<size> for indels)
       - mismatches - Number of re-alignment mismatches
       - strands - Strand(s) of contig when realigned
       - total_matches - The number of alignment matches between the contig sequence segment and the reference sequence, based on the realignment.
       - sv_type - Type of variation detected
       - sv_subtype - Specific type of rearrangement (inversion, tandem duplication, translocation)
       - split_read_count - Number of assembled reads that cover the assembled contig at the inferred breakpoint
       - disc_read_count - Number of discordantly-mapped paired-end reads that support the event (Not applicable to indels).
       - breakpoint_coverages - Number of non-duplicated reads aligned at the inferred breakpoint locations.
       - contig_id - Contig ID
       - contig_seq - Contig sequence.
  - Each target gene in which a SV was detected has a separate output directory (\<analysis\_dir\>/output/\<target\_name\>) containing formatted output specific to the target and the related reference-aligned sequence reads for the contigs that contain the structural variants detected in BAM format.
=======

Documentation:
http://ryanabo.github.io/BreaKmer

