import sys
import logging
import os
import subprocess
def link_reads(sample_freq, read_ones):
    """Link files into temp directory"""
    for sample_name, value in sample_freq.items():
        if value > 1:
            sample_first_reads = [
                x for x in read_ones if sample_name in x]
            for sample_read_one in sample_first_reads:
                combined_file_one = \
                    sample_name + "_L999_R1_" + \
                    (sample_read_one.split("_")[-1]).split(".")[0] + \
                    ".fastq.gz"
                combined_file_two = \
                    sample_name + \
                    "_L999_R2_" + \
                    (sample_read_one.split("_")[-1]).split(".")[0] + \
                    ".fastq.gz"

                try:
                    subprocess.call(f"zcat -f {__directory__}/{sample_read_one} | \
                                    gzip >> {TMPDIR}/{combined_file_one}",
                                    shell=True, stderr=subprocess.STDOUT)
                except subprocess.CalledProcessError as subprocess_error:
                    logging.error(f"Preprocessing: [Merging lanes] Could not merge {sample_name}")
                    logging.error(f"ERROR: {subprocess_error}")
                    sys.exit(f"Could not merge read files for {sample_name}!")
                else:
                    logging.debug(f"Preprocessing: [Merging lanes] Merged reads for {sample_name}")
                sample_read_two = sample_read_one.replace("_R1_", "_R2_")

                try:
                    subprocess.call(f"zcat -f {__directory__}/{sample_read_two} | \
                                    gzip >> {TMPDIR}/{combined_file_two}",
                                    shell=True, stderr=subprocess.STDOUT)
                except subprocess.CalledProcessError as subprocess_error:
                    logging.error(f"Preprocessing: [Merging lanes] Could not merge {sample_name}")
                    logging.error(f"ERROR: {subprocess_error}")
                    sys.exit(f"Could not merge read files for sample {sample_name}!")
                else:
                    logging.debug(f"Preprocessing: [Merging lanes] Merged read for {sample_name}")
        else:
            full_sample_name = [x for x in read_ones if sample_name in x]
            sys_call_sample_one = f"ln -sL {__directory__}/{full_sample_name[0]} \
                                            {TMPDIR}/{full_sample_name[0]}"
            try:
                subprocess.call(sys_call_sample_one,
                                shell=True, stderr=subprocess.STDOUT)
            except subprocess.CalledProcessError as subprocess_error:
                logging.error(f"Preprocessing: [Linking reads] Could not link {sample_name} files")
                logging.error(f"ERROR: {subprocess_error}")
                sys.exit(f"Could not link read files for sample {sample_name}!")
            else:
                logging.debug(f"Preprocessing: [Linking reads] Linked reads for {sample_name}")

            sample_two = full_sample_name[0].replace("_R1_", "_R2_")
            sys_call_sample_two = f"ln -sL {__directory__}/{sample_two} {TMPDIR}/{sample_two}"
            try:
                subprocess.call(sys_call_sample_two,
                                shell=True, stderr=subprocess.STDOUT)
            except subprocess.CalledProcessError as subprocess_error:
                logging.error(f"Preprocessing: [Linking reads] Could not link {sample_name} files")
                logging.error(f"ERROR: {subprocess_error}")
                sys.exit(f"Could not link read files for sample {sample_name}!")
            else:
                logging.debug(f"Preprocessing: [Linking reads] Linked reads for {sample_name}")

def weight_profile(allele_count, weight_dict):
    """
    Function   : weight_profile
    Input      : allele count global var, weight factors
    Output/Desc: Normalizes __allele_count__ by weight factor
    """
    logging.debug("Post-processing: Adjusting counts based on marker size")
    weighted_dict = {}
    for sample in allele_count:
        weighted_dict[sample] = {}
        for loc in allele_count[sample]:
            weighted_dict[sample][loc] = {}
            for a_num in allele_count[sample][loc]:
                if a_num in weight_dict[loc]:
                    weight = allele_count[sample][loc][a_num] // weight_dict[loc][a_num]
                    weighted_dict[sample][loc][a_num] = weight
                else:
                    weighted_dict[sample][loc][a_num] = allele_count[sample][loc][a_num]

    return weighted_dict
HELP_TEXT_SMALL = """
To build a database:
abil_URDOcaller.py --buildDB -c <config file> [-k <int>] [-P|--prefix <database prefix>] [-a <log file path>]

To predict and call markers:
abil_URDOcaller.py --predict -c <config file> -1 <fwd read FASTQ> -2 <rev read FASTQ> [-d <input directory>] [-o <output file>] [-P | --prefix <database prefix>] [-r] -[x]

abil_URDOcaller.py --help for more detailed instructions
"""

HELP_TEXT = """
Readme for abil_URDOcaller
=============================================================================================
Usage
./abil_URDOcaller.py
[--buildDB]
[--predict]
[-1 filename_fastq1][--fastq1 filename_fastq1]
[-2 filename_fastq2][--fastq2 filename_fastq2]
[-d directory][--dir directory][--directory directory]
[-c][--config]
[-P][--prefix]
[-a]
[-k]
[-o output_filename][--output output_filename]
[-x][--overwrite]
[-r]
[-t][--threads]
[-v]
[-h][--help]
==============================================================================================

There are two steps to predicting ST using abil_URDOcaller.
1. Create DB : abil_URDOcaller.py --buildDB
2. Predict : abil_URDOcaller --predict

1. abil_URDOcaller.py --buildDB

Synopsis:
abil_URDOcaller.py --buildDB -c <config file> -k <kmer length(optional)> -P <DB prefix(optional)>
  config file : is a tab delimited file which has the information for typing scheme ie loci, its multifasta file and profile definition file.
    Format :
      [loci]
      locus1    locusFile1
      locus2    locusFile2
      [profile]
      profile   profile_file
  kmer length : is the kmer length for the db. Note, while processing this should be smaller than the read length.
    We suggest kmer lengths of 35, 66 depending on the read length.
  DB prefix(optional) : holds the information for DB files to be created and their location. This module creates 3 files with this prefix.
    You can use a folder structure with prefix to store your db at particular location.

Required arguments
--buildDB
  Identifier for build db module
-c,--config = <configuration file>
  Config file in the format described above.

Optional arguments
-k = <kmer length>
  Kmer size for which the db has to be formed(Default k = 35). Note the tool works best with kmer length in between 35 and 66
  for read lengths of 55 to 150 bp. Kmer size can be increased accordingly. It is advised to keep lower kmer sizes
  if the quality of reads is not very good.
-P,--prefix = <prefix>
  Prefix for db and log files to be created(Default = kmer). Also you can specify folder where you want the dbb to be created.
-a
        File location to write build log
-h,--help
  Prints the help manual for this application

 --------------------------------------------------------------------------------------------

2. abil_URDOcaller.py --predict

abil_URDOcaller --predict : can run in two modes
  1) single sample (default mode)
  2) multi-sample : run abil_URDOcaller for all the samples in a folder (for a particular specie)

Synopsis
abil_URDOcaller.py --predict -1 <fastq file> -2 <fastq file> -d <directory location> - -P <DB prefix(optional)> -k <kmer length(optional)> -o <output file> -x

Required arguments
--predict
  Identifier for predict module
-c,--config = <configuration file>
  Config file in the format described above.

Optional arguments
-1,--fastq1 = <fastq1_filename>
  Path to first fastq file for paired end sample and path to the fastq file for single end file.
  Should have extension fastq or fq.
-2,--fastq2 = <fastq2_filename>
  Path to second fastq file for paired end sample.
  Should have extension fastq or fq.
-d,--dir,--directory = <directory>
  Directory containing paired end read files for multi-sample prediction
-k = <kmer_length>
  Kmer length for which the db was created (Default k = 35).
takes one line. For paired end samples the 2 files should be tab separated on single line.
-o,--output = <output_filename>
  Prints the output to a file instead of stdio.
-P,--prefix = <prefix>
  Prefix using which the db was created (Defaults = kmer). The location of the db could also be provided.
-r
  A FASTQ file is generated in the current directory for each sample containing reads with kmer matches.
-R, --readsdir = < output directory >
  A FASTQ file is generated in the specified directory for each sample containing reads with kmer matches.
-t
  Integer number of threads to use to process samples
-v
  Prints the version of the software.
-x,--overwrite
  By default abil_URDOcaller appends the results to the output_filename if same name is used.
  This argument overwrites the previously specified output file.
-h,--help
  Prints the help manual for this application
"""
