#!/usr/bin/env python
"""
The program has 3 basic modes:
    buildDB: for building the required databases
    predict, single sample: Predict markers and phenotypes from raw reads
    predict, multisample: Like above, but for multiple samples stored
        at a common location (paired end samples)
"""
import getopt
import sys
import logging
import os
import subprocess
import tempfile
import shutil
from itertools import islice
import operator
VERSION = """ abil_URDOcaller ALPHA.1 (updated : July 13, 2018) """
"""
abil_URDOcaller free for academic users and requires permission before any
commercial or government usage of any version of this code/algorithm.
If you are a commercial or governmental user, please contact abil@ihrc.com
for permissions.

abil_URDOcaller is licensed under a modified Creative Commons By-NC-SA v4
license, please see the LICENSE file for specific terms.

For additional terms and conditions for government employees, see
"For Government Employees" section
"""

#predict part starts here
############################################################
HELP_TEXT_SMALL = """
To build a database:
abil_URDOcaller.py --buildDB -c <config file> [-k <int>] [-P|--prefix <database prefix>] [-a <log file path>]

To predict and call markers:
abil_URDOcaller.py --predict -c <config file> -1 <fwd read FASTQ> -2 <rev read FASTQ> [-d < input directory] [-o <output file>] [-P | --prefix <database prefix>] [-r] -[x] 

abil_URDOcaller.py --help for more detailed instructions
"""

HELP_TEXT = ""
TMPDIR = tempfile.mkdtemp()
def batch_tool(fdir, kmer):
    """
    Function   : batch_tool
    Input      : Directory name, paired only, k value
    Output     : STs and allelic profiles for each FASTQ file
    Description: Processes all FASTQ files present in the input directory
    """
    if not fdir.endswith('/'):
        fdir += '/'
    freq_dict_samples = {}
    all_first_reads = [x for x in os.listdir(fdir) if "_R1_" in x]
    for first_read_name in all_first_reads:
        sample_name = "_".join(first_read_name.split("_")[:-3])
        if sample_name in freq_dict_samples:
            freq_dict_samples[sample_name] += 1
        else:
            freq_dict_samples[sample_name] = 1

    for sample_name, value in freq_dict_samples.items():
        if value > 1:
            sample_first_reads = [x for x in all_first_reads if sample_name in x]
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
                    subprocess.call(f"zcat -f {fdir}/{sample_read_one} | \
                                    gzip >> {TMPDIR}/{combined_file_one}",
                                    shell=True, stderr=subprocess.STDOUT)
                except subprocess.CalledProcessError as subprocess_error:
                    logging.error(f"Preprocessing: [Merging lanes] Could not merge read files for sample {sample_name}")
                    sys.exit(f"Could not merge read files for sample {sample_name}!")
                else:
                    logging.debug(f"Preprocessing: [Merging lanes] Merged read files for sample {sample_name}")
                sample_read_two = sample_read_one.replace("_R1_", "_R2_")

                try:
                    subprocess.call(f"zcat -f {fdir}/{sample_read_two} | \
                                    gzip >> {TMPDIR}/{combined_file_two}",
                                    shell=True, stderr=subprocess.STDOUT)
                except subprocess.CalledProcessError as subprocess_error:
                    logging.error(f"Preprocessing: [Merging lanes] Could not merge read files for sample {sample_name}\n{subprocess_error}")
                    sys.exit(f"Could not merge read files for sample {sample_name}!")
                else:
                    logging.debug(f"Preprocessing: [Merging lanes] Merged read files for sample {sample_name}")
        else:
            full_sample_name = [x for x in all_first_reads if sample_name in x]
            sys_call_string_sample_one = f"ln -sL {fdir}/{full_sample_name[0]} \
                                            {TMPDIR}/{full_sample_name[0]}"
            try:
                subprocess.call(sys_call_string_sample_one,
                                shell=True, stderr=subprocess.STDOUT)
            except subprocess.CalledProcessError as subprocess_error:
                logging.error(f"Preprocessing: [Linking reads] Could not link read files for sample {sample_name}\n{subprocess_error}")
                sys.exit(f"Could not link read files for sample {sample_name}!")
            else:
                logging.debug(f"Preprocessing: [Linking reads] Linked read files for sample {sample_name}")

            sample_two = full_sample_name[0].replace("_R1_", "_R2_")
            sys_call_string_sample_two = f"ln -sL {fdir}/{sample_two} {TMPDIR}/{sample_two}"
            try:
                subprocess.call(sys_call_string_sample_two,
                                shell=True, stderr=subprocess.STDOUT)
            except subprocess.CalledProcessError as subprocess_error:
                logging.error(f"Preprocessing: [Linking reads] Could not link read files for sample {sample_name}\n{subprocess_error}")
                sys.exit(f"Could not link read files for sample {sample_name}!")
            else:
                logging.debug(f"Preprocessing: [Linking reads] Linked read files for sample {sample_name}")


    file_list = [x for x in os.listdir(TMPDIR) if "_R1_" in x]
    for read_one in file_list:
        fastq1_processed = f"{TMPDIR}/{read_one}"
        fastq2_processed = fastq1_processed.replace("_R1_", "_R2_")
        kCount = single_sample_tool(fastq1_processed, fastq2_processed,  kmer, rawCounts)
    shutil.rmtree(TMPDIR)
    return kCount
def single_sample_tool(fastq1, fastq2, k, results):
    """
    Function   : single_sample_tool
    Input      : fastq file 1 and 2, paired or single, k value, output dictionary
    Output     : STs and allelic profiles for each FASTQ file
    Description: Processes both FASTQ files passed to the function
    """
    # pp.pprint(kmerDict)
    if __reads__:
        read_file_name = fastq1.split('/')[-1].split('.')[0][:-1] + '_reads.fq'
        global read_file
        read_file = open(read_file_name, 'w+')
    sName = fastq1.split('/')[-1].split('_')[0]
    msg = f"Preprocessing: working with {fastq1} and {fastq2}"
    logging.debug(msg)
    global alleleCount
    alleleCount = {}
    logging.debug(f"Preprocessing: [Merging reads] Merging {fastq1} and {fastq2}")
    vsearch_cmd = f"vsearch --fastq_mergepairs {fastq1} --reverse {fastq2} --fastqout {TMPDIR}/reads.fq 2>/dev/null 1>/dev/null"
    logging.debug(f"Preprocessing: [Merging reads] VSEARCH command\n\t{vsearch_cmd}")
    try:
        subprocess.call(vsearch_cmd, shell=True, stderr=subprocess.STDOUT)
    except subprocess.CalledProcessError as subprocess_error:
        logging.error(f"Preprocessing: [Merging reads] Could not merge read files for sample {sample_name}")
        sys.exit(f"Could not merge read files for sample {sample_name}!")
    read_processor(TMPDIR, k, sName)
    if kCount == {}:
        string = f"No k-mer matches were found for the sample {fastq1} and {fastq2}"
        string += f"\n\tProbable cause of the error:  low quality data/too many N's in the data"
        logging.error(f"ERROR: {string}")
        print(string)
        exit()
    if __reads__:
        read_file.close()
    return kCount
def read_processor(fastq, k, sName):
    """
    Function   : read_processor
    Input      : fastq file, k value
    Output     : Edits a global dictionary - results
    Description: Processes the single fastq file
    """
    fastq = "/".join([fastq, "reads.fq"])
    msg = f"Analysis: Begin processing merged reads ({fastq})"
    logging.debug(msg)
    if os.path.isfile(fastq):
        logging.debug(f"Analysis: Kmer counting using k={k}")
        finalProfile = {}
        if sName not in kCount:
            kCount[sName] = {}
        f = open(fastq)
        for lines in iter(lambda: tuple(islice(f, 4)), ()):
            if len(lines) < 4:
                dbg = "ERROR: Please verify the input FASTQ files are correct"
                logging.debug(dbg)
            try:
                if len(lines[1]) < k:
                    m1 = f"ERROR: Read ID: {[0][1:]} is length {len(lines[1])} and < {k}"
                    print(m1)
                    logging.debug(m1)
                    return 0
            except Exception:
                m2 = f"ERROR: Check fastq file {file}"
                logging.debug(m2)
                return 0
            start = int((len(lines[1])-k)//2)
            firstKmer = str(lines[1][:k])
            midKmer = str(lines[1][start:k+start])
            lastKmer = str(lines[1][-35:])
            if firstKmer in kmerDict[k] or midKmer in kmerDict[k] or lastKmer in kmerDict[k]:
                count_kmers(lines[1], k, sName)
                if __reads__: read_file.write('\n'.join('{}'.format(l) for l in lines))
    else:
        logging.error(f"ERROR: File does not exist: {fastq}")
def count_kmers(read, k, sName):
    """
    Function   : goodReads
    Input      : sequence read, k, step size
    Output     : Edits the count of global variable alleleCount
    Description: Increment the count for each k-mer match
    """
    alleleCount = {}
    n = 0
    line = read.rstrip()
    for n in range(len(line)-k-1):
        s = str(line[n:n+k])
        if s in kmerDict[k]:
            for probLoc in kmerDict[k][s]:
                if probLoc not in alleleCount:
                    alleleCount[probLoc] = {}
                a = kmerDict[k][s][probLoc]
                for allele in a:
                    allele = allele.rstrip()
                    if allele in alleleCount[probLoc]:
                        alleleCount[probLoc][allele] += 1
                    else:
                        alleleCount[probLoc][allele] = 1
        n += 1
    maxSupports = {}
    for allele in alleleCount:
        alleleNumber = max(alleleCount[allele].items(), key=operator.itemgetter(1))[0]
        alleleKcount = alleleCount[allele][max(alleleCount[allele].items(), key=operator.itemgetter(1))[0]]
        if allele not in maxSupports:
            maxSupports[allele] = {}
        maxSupports[allele][alleleNumber] = alleleKcount
        if allele not in kCount[sName]:
            kCount[sName][allele] = {}
        if alleleNumber not in kCount[sName][allele]:
            kCount[sName][allele][alleleNumber] = 1
        else:
            kCount[sName][allele][alleleNumber] += 1
def weight_profile(alleleCount, weightDict):
    """
    Function   : weight_profile
    Input      : allele count global var, weight factors
    Output/Desc: Normalizes alleleCount by weight factor
    """
    logging.debug("Post-processing: Adjusting counts based on marker size")
    weightedDict = {}
    for sample in alleleCount:
        weightedDict[sample] = {}
        for loc in alleleCount[sample]:
            weightedDict[sample][loc] = {}
            for aNum in alleleCount[sample][loc]:
                if aNum in weightDict[loc]:
                    weightedDict[sample][loc][aNum] = alleleCount[sample][loc][aNum] // weightDict[loc][aNum]
                else:
                    weightedDict[sample][loc][aNum] = alleleCount[sample][loc][aNum]

    return weightedDict

def load_module(k, dbPrefix):
    """
    Function   : load_module
    Input      : k value and prefix of the DB file
    Output     : Updates the DB dictionary variables
    Description: Used in loading the DB as set of variables by calling other functions
    """
    global dbFile
    dbFile = dbPrefix+'_'+str(k)+'.txt'
    global weightFile
    weightFile = dbPrefix+'_weight.txt'
    global profileFile
    profileFile = dbPrefix+'_profile.txt'
    global kmerDict
    kmerDict = {}
    kmerDict[k] = load_kmer_dict(dbFile)
    global weightDict
    weightDict = load_weight_dict(weightFile)
    global stProfile
    stProfile = load_st_from_file(profileFile)
    load_config(__config__)
def load_st_from_file(profileF):
    """
    Function   : load_st_from_file
    Input      : profile definition file
    Output     : Updates the DB dictionary variables
    Description: Used in loading the DB as set of variables
    """
    with open(profileF, 'r') as definitionFile:
        st = {}
        for line in definitionFile:
            if not line.startswith("marker"):
                cols = line.rstrip().rsplit('\t')
                if cols[0] not in st:
                    st[cols[0]] = {}
                st[cols[0]][cols[1]] = cols[2]
    return st
def load_kmer_dict(dbFile):
    """
    Function   : load_kmer_dict
    Input      : DB prefix
    Output     : Updates the DB dictionary variables
    Description: Used in loading the DB as set of variables
    """
    kmerTableDict = {}
    with open(dbFile, 'r') as kmerTableFile:
        lines = kmerTableFile.readlines()
        for line in lines:
            array = line.rstrip().rsplit('\t')
            kmerTableDict[array[0]] = {}
            kmerTableDict[array[0]][array[1]] = array[2][1:-1].rsplit(',')
    return kmerTableDict
def load_weight_dict(weightFile):
    """
    Function   : load_weight_dict
    Input      : Weight file prefix
    Output     : Updates the DB dictionary variables
    Description: Used in loading the DB as set of variables
    """
    weightDict = {}
    with open(weightFile, 'r') as weightTableFile:
        lines = weightTableFile.readlines()
        for line in lines:
            array = line.rstrip().rsplit('\t')
            try:
                (loc, allele) = array[0].replace('-', '_').rsplit('_', 1)
            except ValueError:
                print("Error : Allele name in locus file should be seperated by '_' or '-'")
                exit(0)
            if loc not in weightDict:
                weightDict[loc] = {}
            weightDict[loc][allele] = float(array[1])
    return weightDict

def load_config(config):
    """
    Function   : load_config
    Input      : config file path from getopts
    Output     : Updates configDict
    Description: Used to find allele fasta files for getCoverage
    """
    global configDict
    configDict = {}
    with open(config) as configFile:
        lines = configFile.readlines()
        head = ''
        for line in lines:
            if line.rstrip() == '':
                continue
            if line.rstrip() == '[loci]':
                head = 'loci'
                configDict[head] = {}
            elif line.rstrip() == '[profile]':
                head = 'profile'
                configDict[head] = {}
            else:
                arr = line.strip().split()
                configDict[head][arr[0]] = arr[1]
    for head in configDict:
        for element in configDict[head]:
            if not os.path.isfile(configDict[head][element]):
                print("ERROR: %s file does not exist at %s" % (element, configDict[head][element]))
                exit(0)
    return configDict

def print_results(results, output_filename, overwrite):
    """
    Function   : print_results
    Input      : results, output file, overwrite?
    Output     : Prints on the screen or in a file
    Description: Prints the results in the format asked by the user
    """
    # pp.pprint(results)
    if output_filename != None:
        if overwrite is False:
            outfile = open(output_filename, "a")
        else:
            outfile = open(output_filename, "w")
    # pp.pprint(stProfile)
    output = {}
    output["sample"] = []
    outString = 'Sample'
    logging.debug("Post-processing: Finding most likely phenotypes and markers")
    for s in results:
        sample = s
        output["sample"].append(s)
        for loc in results[s]:
            if loc == "genericmarkers":
                for mID in results[s][loc]:
                    if stProfile[loc][mID] not in output:
                        output[stProfile[loc][mID]] = {}
                    output[stProfile[loc][mID]][s] = results[s][loc][mID]
            else:
                max_mID = max(results[s][loc].items(), key=operator.itemgetter(1))[0]
                if stProfile[loc][max_mID] not in output:
                    output[stProfile[loc][max_mID]] = {}
                output[stProfile[loc][max_mID]][s] = results[s][loc][max_mID]
    sorted_samples = sorted(output["sample"])
    for sample in sorted_samples:
        outString += (f"\t{sample}")
    outString += "\n"
    sorted_output_keys = sorted(output)
    for key in sorted_output_keys:
        if key != "sample":
            outString += f"{key}"
            for sample in output["sample"]:
                if sample in output[key]:
                    outString += f"\t{output[key][sample]}"
                else:
                    outString += f"\t0"
            outString += "\n"
    if output_filename != None:
        outfile.write(f"{outString}\n")
    else:
        print(f"{outString}\n")
################################################################################
# Predict part ends here
################################################################################
def reverse_complement(seq):
    """
    Build DB part starts
    Returns the reverse complement of the sequence
    """
    seqU = seq.upper()
    seq_dict = {'A':'T', 'T':'A',
                'G':'C', 'C':'G', 'Y':'R',
                'R':'Y', 'S':'S', 'W':'W',
                'K':'M', 'M':'K', 'N':'N'}
    try:
        return "".join([seq_dict[base] for base in reversed(seqU)])
    except Exception:
        strn = "Reverse Complement Error:" + seqU
        logging.debug(strn)

def get_fasta_dict(fullLocusFile):
    """
    Function   : get_fasta_dict
    Input      : locus file name
    Output     : dictionary with all the allele sequences
    Description: Stores each allele sequence in a dictionary
    """
    logging.debug("Create Fasta Dict")
    logging.debug(fullLocusFile)
    fastaFile = open(fullLocusFile, 'r').read()
    entries = [x for x in fastaFile.split('>') if len(x) != 0]
    fastaDict = {}
    for entry in entries:
        key = [x for x in entry.split('\n')[0].split() if len(x) != 0][0]
        sequence = ''.join(entry.split('\n')[1:]).rstrip()
        fastaDict[key] = {'sequence':sequence}
    return fastaDict

def form_kmer_db(configDict, k, output_filename):
    """
    Function   : form_kmer_db
    Input      : configuration file, k value, output prefix
    Output     : abil_URDOcaller DB
    Description: Constructs the k-mer DB in both strand orientation
    """
    dbFileName = output_filename+'_'+str(k)+'.txt'
    weightFileName = output_filename+'_weight.txt'
    kmerDict = {}
    mean = {}
    for locus in configDict['loci']:
        msgs = "formKmerDB :" +locus
        logging.debug(msgs)
        fastaDict = get_fasta_dict(configDict['loci'][locus])
        total_kmer_length = 0
        seq_location = 0
        for allele in list(fastaDict.keys()):
            seq = fastaDict[allele]['sequence'].strip()
            seq_len = len(seq)
            total_kmer_length += seq_len
            seq_location += 1
            try:
                (loc, num) = allele.replace('-', '_').rsplit('_', 1)
            except ValueError:
                print("Error : Allele name in locus file should be seperated by '_' or '-'")
                exit(0)
            splitId = allele.replace('-', '_').rsplit('_', 1)
            i = 0
            while i+k <= seq_len:
                kmer = seq[i:i+k]
                revCompKmer = reverse_complement(kmer)
                if kmer not in kmerDict:
                    kmerDict[kmer] = {}
                    kmerDict[kmer][splitId[0]] = []
                    kmerDict[kmer][splitId[0]].append(int(splitId[1]))
                else:
                    if splitId[0] not in kmerDict[kmer]:
                        kmerDict[kmer][splitId[0]] = []
                        kmerDict[kmer][splitId[0]].append(int(splitId[1]))
                    else:
                        kmerDict[kmer][splitId[0]].append(int(splitId[1]))
                if revCompKmer not in kmerDict:
                    kmerDict[revCompKmer] = {}
                    kmerDict[revCompKmer][splitId[0]] = []
                    kmerDict[revCompKmer][splitId[0]].append(int(splitId[1]))
                else:
                    if splitId[0] not in kmerDict[revCompKmer]:
                        kmerDict[revCompKmer][splitId[0]] = []
                        kmerDict[revCompKmer][splitId[0]].append(int(splitId[1]))
                    else:
                        kmerDict[revCompKmer][splitId[0]].append(int(splitId[1]))
                i += 1
        mean[locus] = total_kmer_length/seq_location*1.0
    with open(dbFileName, 'w') as kfile:
        for key in kmerDict:
            for key1 in kmerDict[key]:
                string = str(key)+'\t'+str(key1)+'\t'+str(kmerDict[key][key1]).replace(" ", "")+'\n'
                kfile.write(string)
    with open(weightFileName, 'w') as wfile:
        for locus in configDict['loci']:
            fastaDict = get_fasta_dict(configDict['loci'][locus])
            for allele in list(fastaDict.keys()):
                splitId = allele.split('_')
                seq = fastaDict[allele]['sequence']
                l = len(seq)
                fac = (l/mean[locus])
                s = allele  + '\t' + str(fac) + '\n'
                if fac > 1.05 or fac < 0.95:
                    wfile.write(s)

def copy_profile(profileDict, output_filename):
    """
    Function   : copy_profile
    Input      : profileDict
    Output     : None
    Description: Duplicated profile file for db
    """
    profileFileName = output_filename+'_profile.txt'
    with open(profileDict['profile']) as f:
        lines = f.readlines()
        with open(profileFileName, "w") as f1:
            f1.writelines(lines)

def make_custom_db(config, k, output_filename):
    """
    Function   : make_custom_db
    Input      : configuration file, k value, output prefix
    Output     : None
    Description: Processes the config file and calls the relevant function
    """
    configDict = {}
    if output_filename == None:
        output_filename = 'kmerDB'
    with open(config, 'r') as configFile:
        lines = configFile.readlines()
        head = ''
        for line in lines:
            if line.rstrip() == '':
                continue
            if line.rstrip() == '[loci]':
                head = 'loci'
                configDict[head] = {}
            elif line.rstrip() == '[profile]':
                head = 'profile'
                configDict[head] = {}
            else:
                print(line.strip().split())
                arr = line.strip().split()
                configDict[head][arr[0]] = arr[1]
    for head in configDict:
        for element in configDict[head]:
            if not os.path.isfile(configDict[head][element]):
                print("ERROR: %s file does not exist at %s" % (element, configDict[head][element]))
                exit(0)
    form_kmer_db(configDict, k, output_filename)
    copy_profile(configDict['profile'], output_filename)

################################################################################
# Build DB part ends
# Check Parameters
################################################################################


def check_params(buildDB, predict, config, k, batch, directory, fastq1, fastq2, dbPrefix):
    """
    Check input parameters
    """
    if predict and buildDB:
        print(HELP_TEXT_SMALL)
        print("Select either predict or buildDB module")
        exit(0)
    if not predict and not buildDB:
        print(HELP_TEXT_SMALL)
        print("Select either predict or buildDB module")
        exit(0)
    if predict:
        if config is None and coverage:
            print(HELP_TEXT_SMALL)
            print("Config parameter is required.")
            exit(0)
        if not os.path.isfile(dbPrefix+'_'+str(k)+'.txt'):
            print(HELP_TEXT_SMALL)
            print("DB file does not exist : ", dbPrefix, '_', str(k), '.txt or change DB prefix.')
            exit(0)
        if not os.path.isfile(dbPrefix+'_weight.txt'):
            print(HELP_TEXT_SMALL)
            print("DB file does not exist : ", dbPrefix, '_weight.txt or change DB prefix.')
            exit(0)
        if not os.path.isfile(dbPrefix+'_profile.txt'):
            print(HELP_TEXT_SMALL)
            print("DB file does not exist : ", dbPrefix, '_profile.txt or change DB prefix.')
            exit(0)
        elif batch:
            if not os.path.isdir(directory):
                print(HELP_TEXT_SMALL)
                print("Error: Directory ("+directory+") does not exist!")
                exit(0)
        elif paired:
            if not os.path.isfile(fastq1):
                print(HELP_TEXT_SMALL)
                print("Error: FASTQ file ("+fastq1+") does not exist!")
                exit(0)
            if not os.path.isfile(fastq2):
                print(HELP_TEXT_SMALL)
                print("Error: FASTQ file ("+fastq2+") does not exist!")
                exit(0)
        elif not paired:
            if not os.path.isfile(fastq1):
                print(HELP_TEXT_SMALL)
                print("Error: FASTQ file ("+fastq1+") does not exist!")
                exit(0)
    if buildDB:
        try:
            if not os.path.isfile(config):
                print(HELP_TEXT_SMALL)
                print("Error: Configuration file ("+config+") does not exist!")
                exit(0)
        except Exception:
            print(HELP_TEXT_SMALL)
            print("Error: Specify Configuration file")
            exit(0)

################################################################################
# The Program Starts Execution Here
# Default Params
################################################################################

__buildDB__ = False
__predict__ = False
__output_filename__ = None
__batch__ = False
__overwrite__ = False
__paired__ = False
__fastq1__ = None
__fastq2__ = None
__user_k__ = False
__config__ = None
__timeDisp__ = False
__reads__ = False
__dbPrefix__ = 'kmer'
__log__ = ''
__k__ = 35
directory = None
__reads__ = False
"""Input arguments"""
options, remainder = getopt.getopt(sys.argv[1:], 'o:x1:2:kbd:phP:c:rva:', [
    'buildDB',
    'predict',
    'output=',
    'config=',
    'prefix=',
    'overwrite',
    'batch',
    'fastq1=',
    'fastq2=',
    'dir=',
    'directory=',
    'help'])
for opt, arg in options:
    if opt in ('-o', '--output'):
        __output_filename__ = arg
    elif opt in ('-x', '--overwrite'):
        __overwrite__ = True
    elif opt in '--buildDB':
        __buildDB__ = True
    elif opt in ('-P', '--prefix'):
        __dbPrefix__ = arg
    elif opt in '--predict':
        __predict__ = True
    elif opt in ('-c', '--config'):
        __config__ = arg
    elif opt in '-k':
        __user_k__ = True
        try:
            __k__ = int(arg)
        except ValueError:
            print("Error: Enter a numerical k value.")
            exit(0)
        # Check to make sure the arg is an int.
    elif opt in ('-1', '--fastq1'):
        __fastq1__ = arg
    elif opt in ('-2', '--fastq2'):
        __fastq2__ = arg
    elif opt in ('-d', '--dir', '--directory'):
        __directory__ = os.path.abspath(arg)
        __batch__ = True
    elif opt in '-a':
        __log__ = arg
    elif opt in '-r':
        __reads__ = True
    elif opt in '-v':
        print(VERSION)
        exit(0)
    elif opt in ('-h', '--help'):
        print(HELP_TEXT)
        exit(0)
check_params(__buildDB__, __predict__, __config__, __k__, __batch__, __directory__, __fastq1__, __fastq2__, __dbPrefix__)
if __buildDB__:
    try:
        if not __log__:
            log = __dbPrefix__+'.log'
    except TypeError:
        log = 'kmer.log'
    logging.basicConfig(filename=log, level=logging.DEBUG,
                        format='%(asctime)s %(message)s', datefmt='%m/%d/%Y %I:%M:%S %p')
    if os.path.isfile(__config__):
        print("Info: Making DB for k = ", __k__)
        print("Info: Making DB with prefix =", __dbPrefix__)
        print("Info: Log file written to ", log)
        make_custom_db(__config__, __k__, _dbPrefix__)
    else:
        print("Error: The input config file "+__config__ +" does not exist.")
elif __predict__:
    try:
        if not __log__:
            log = __dbPrefix__+'.log'
    except TypeError:
        log = 'kmer.log'
    logging.basicConfig(filename=log, level=logging.DEBUG,
                        format='%(asctime)s %(message)s', datefmt='%m/%d/%Y %I:%M:%S %p')
    logging.debug("==============================================================================")
    logging.debug(f"Command: {' '.join(sys.argv)}")
    logging.debug("Starting Marker Prediction")
    logging.debug(f"Temporary directory: {TMPDIR}")
    load_module(__k__, __dbPrefix__)
    rawCounts = {}
    global kCount
    kCount = {}
    if __batch__:
        rawCounts = batch_tool(__directory__, __k__)
    else:
        rawCounts = single_sample_tool(__fastq1__, __fastq2__, __paired__, __k__, rawCounts)
    weightCounts = weight_profile(rawCounts, weightDict)
    print_results(weightCounts, __output_filename__, __overwrite__)
else:
    print("Error: Please select the mode")
    print("--buildDB (for database building) or --predict (for marker discovery)")
