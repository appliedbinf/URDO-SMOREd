"""
Heplper functions for SMORE'd
"""
import sys
import logging
import os
import subprocess
def link_reads(sample_freq, read_ones):
    """Link files into temp directory"""
    for sample_name, value in sample_freq.items():
        if value > 1:
            sample_first_reads = [
                x for x in read_ones if x.startswith(sample_name)]
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
            full_sample_name = [x for x in read_ones if x.startswith(sample_name)]
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

def read_mapping_file(amp2tax):
    """
    Function   : make_report
    Input      : results, output file, overwrite?
    Output     : Takes in Excel template, outputs one report per sample
    Description: Create Excel reports
    """
    tax_dict = {}
    amplicons = open(amp2tax, 'r').read()
    entries = [x for x in amplicons.split('\n') if len(x) != 0]
    for entry in entries:
        org, taxonomy = entry.rstrip().split('\t')
        if org not in tax_dict:
            tax_dict[org] = {'kingdom' : "", "phylum" : "", "class" : "",
                             "order" : "", "family" : "", "genus" : "",
                             "species" : "", "drug" : "", "type" : ""}
        taxonomy = taxonomy.split(';')
        for tax in taxonomy:
            if tax.startswith('sk:'):
                tax_dict[org]['kingdom'] = tax.split(":", 1)[1]
            elif tax.startswith('p:'):
                tax_dict[org]['phylum'] = tax.split(":", 1)[1]
            elif tax.startswith('c:'):
                tax_dict[org]['class'] = tax.split(":", 1)[1]
            elif tax.startswith('o:'):
                tax_dict[org]['order'] = tax.split(":", 1)[1]
            elif tax.startswith('f:'):
                tax_dict[org]['family'] = tax.split(":", 1)[1]
            elif tax.startswith('g:'):
                tax_dict[org]['genus'] = tax.split(":", 1)[1]
            elif tax.startswith('s:'):
                tax_dict[org]['species'] = tax.split(":", 1)[1]
            elif tax.startswith('d:'):
                tax_dict[org]['drug'] = tax.split(":", 1)[1]
            elif tax.startswith('t:'):
                tax_dict[org]['type'] = tax.split(":", 1)[1]
    return tax_dict

def make_report(report_data, read_count, output_filename, template_fp, tax_dict):
    """
    Function   : make_report
    Input      : results, output file, overwrite?
    Output     : Takes in Excel template, outputs one report per sample
    Description: Create Excel reports
    """
    from collections import Counter
    import openpyxl as pyxl
    sample_names = report_data.pop("sample")
    for sample in sample_names:
        template = pyxl.load_workbook(filename=template_fp)
        report = template.active
        report['B3'] = sample
        report['B4'] = subprocess.check_output('date "+%Y-%m-%d"',
                                               shell=True).decode('utf-8').rstrip()
        report['B5'] = read_count[sample_names.index(sample)]
        total_count = 0
        data = {}
        top_three_idx = 9
        start_row = 15 # Magic numbers are bad...fix for later
        bacteria_count = 0
        virus_count = 0
        other_count = 0
        for org in report_data:
            try:
                data[org] = int(report_data[org][sample])
                total_count += int(report_data[org][sample])
            except KeyError:
                data[org] = 0
        unclassified_counts = data.pop("Unclassified")
        report['B6'] = total_count
        for org in Counter(data).most_common(3):
            common_name, count = org
            try:
                genus = tax_dict[common_name]['genus']
                species = tax_dict[common_name]['species']
            except KeyError:
                genus = ""
                species = ""
            report[f"A{top_three_idx}"] = common_name
            report[f"B{top_three_idx}"] = genus
            report[f"C{top_three_idx}"] = species
            report[f"E{top_three_idx}"] = count
            report[f"F{top_three_idx}"] = "{0:.2f}".format(count/total_count * 100)
            top_three_idx += 1
        for org in sorted(data):
            print_row = 0
            if data[org] == 0:
                continue
            common_name = org
            if tax_dict[org]['kingdom'] == 'Bacteria':
                try:
                    genus = tax_dict[org]['genus']
                    species = tax_dict[org]['species']
                    resistance = tax_dict[org]['drug']
                    typ = tax_dict[org]['type']
                except KeyError:
                    genus = ""
                    species = org
                    resistance = ""
                    typ = ""
                print_row = start_row + bacteria_count
                report.insert_rows(print_row)
                report[f"A{print_row}"] = common_name
                report[f"B{print_row}"] = genus
                try:
                    species = species.split(" ", 1)[1]
                except IndexError:
                    pass
                report[f"C{print_row}"] = species
                for idx in ["A", "B", "C"]:
                    report[f"{idx}{print_row}"].font = pyxl.styles.Font(italic=True)
                report[f"D{print_row}"] = resistance
                report[f"E{print_row}"] = data[org]
                report[f"F{print_row}"] = "{0:.2f}".format(data[org]/total_count * 100)
                bacteria_count += 1

            elif tax_dict[org]['kingdom'] == 'Viruses':
                try:
                    family = tax_dict[org]['family']
                    species = tax_dict[org]['species']
                    resistance = tax_dict[org]['drug']
                    typ = tax_dict[org]['type']
                except KeyError:
                    family = ""
                    species = org
                    resistance = ""
                    typ = ""
                print_row = start_row + bacteria_count + 3 + virus_count
                report.insert_rows(print_row)
                report[f"A{print_row}"] = common_name
                report[f"B{print_row}"] = family
                report[f"C{print_row}"] = species
                for idx in ["A", "C"]:
                    report[f"{idx}{print_row}"].font = pyxl.styles.Font(italic=True)
                report[f"D{print_row}"] = resistance
                report[f"E{print_row}"] = data[org]
                report[f"F{print_row}"] = "{0:.2f}".format(data[org]/total_count * 100)
                virus_count += 1

            else:
                print_row = start_row + bacteria_count + 6 + virus_count + other_count
                report.insert_rows(print_row)
                report[f"A{print_row}"] = common_name
                report[f"E{print_row}"] = data[org]
                report[f"F{print_row}"] = "{0:.2f}".format(data[org]/total_count * 100)
                other_count += 1
        print_row = start_row + bacteria_count + 6 + virus_count + other_count
        report[f"A{print_row}"] = "Unclassified"
        report[f"A{print_row}"].font = pyxl.styles.Font(italic=False)
        report[f"E{print_row}"] = unclassified_counts
        report[f"F{print_row}"] = "{0:.2f}".format(unclassified_counts/total_count * 100)
        virus_header = start_row + bacteria_count+1
        report.merge_cells(f'A{virus_header}:F{virus_header}')
        report[f"A{virus_header}"].alignment = pyxl.styles.Alignment(horizontal="center")
        if output_filename is not None:
            output_path = os.path.dirname(os.path.realpath(output_filename))
            filename = os.path.join(output_path, f"{sample}.xlsx")
        else:
            filename = f'{sample}.xlsx'
        template.save(filename)

HELP_TEXT_SMALL = """
To build a database:
smored --buildDB -c <config file> [-k <int>] [-P|--prefix <database prefix>] [-a <log file path>]

To predict and call markers:
smored -c <config file> -1 <fwd read FASTQ> -2 <rev read FASTQ> [-d <input directory>] [-o <output file>] [-P <database prefix>] [-r] [--report] [-u] [-x]

smored --help for more detailed instructions
"""

HELP_TEXT = """
Readme for smored
=============================================================================================
Usage
./smored
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
[-R read/output/path]
[-r]
[--report]
[-U unclassified/read/path]
[-u]
[-t][--threads]
[-v]
[-h][--help]
==============================================================================================

There are two steps to predicting ST using smored.
1. Create DB : smored --buildDB
2. Predict : smored --predict

1. smored --buildDB

Synopsis:
smored --buildDB -c <config file> -k <kmer length(optional)> -P <DB prefix(optional)>
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

2. smored --predict

smored --predict : can run in two modes
  1) single sample (default mode)
  2) multi-sample : run smored for all the samples in a folder (for a particular specie)

Synopsis
smored --predict -1 <fastq file> -2 <fastq file> -d <directory location> - -P <DB prefix(optional)> -k <kmer length(optional)> -o <output file> -x

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
--report
  Generate per-sample reports in Excel format.  If an output file path is provided, per sample reports will
  be deposited in the same folder
-u
  A FASTQ file is generated in the current directory for each sample containing reads with no kmer matches (Unclassified reads).
-U, --unclassified = < output directory >
  A FASTQ file is generated in the specified directory for each sample containing reads wwith no kmer matches (Unclassified reads).
-t
  Integer number of threads to use to process samples
-v
  Prints the version of the software.
-x,--overwrite
  By default smored appends the results to the output_filename if same name is used.
  This argument overwrites the previously specified output file.
-h,--help
  Prints the help manual for this application
"""
