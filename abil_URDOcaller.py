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
commercial or government usage use for any version of this code/algorithm.
If you are a commercial or governmental user, please contact abil@ihrc.com
for permissions.

For additional terms and conditions for government employees, see
"For Government Employees" section

LICENSE TERMS FOR stringMLST
Adopted from: https://creativecommons.org/licenses/by-nc-sa/4.0/

Creative Commons Attribution-NonCommercial-ShareAlike 4.0 International Public
License

By exercising the Licensed Rights (defined below), You accept and agree to be
bound by the terms and conditions of this Creative Commons Attribution-
NonCommercial-ShareAlike 4.0 International Public License ("Public License"). To
the extent this Public License may be interpreted as a contract, You are granted
the Licensed Rights in consideration of Your acceptance of these terms and
conditions, and the Licensor grants You such rights in consideration of benefits
the Licensor receives from making the Licensed Material available under these
terms and conditions.

Section 1 - Definitions.

Adapted Material means material subject to Copyright and Similar Rights that is
derived from or based upon the Licensed Material and in which the Licensed
Material is translated, altered, arranged, transformed, or otherwise modified in
a manner requiring permission under the Copyright and Similar Rights held by the
Licensor. For purposes of this Public License, where the Licensed Material is a
musical work, performance, or sound recording, Adapted Material is always
produced where the Licensed Material is synched in timed relation with a moving
image. Adapter's License means the license You apply to Your Copyright and
Similar Rights in Your contributions to Adapted Material in accordance with the
terms and conditions of this Public License. BY-NC-SA Compatible License means a
license listed at creativecommons.org/compatiblelicenses, approved by Creative
Commons as essentially the equivalent of this Public License. Copyright and
Similar Rights means copyright and/or similar rights closely related to
copyright including, without limitation, performance, broadcast, sound
recording, and Sui Generis Database Rights, without regard to how the rights are
labeled or categorized. For purposes of this Public License, the rights
specified in Section 2(b)(1)-(2) are not Copyright and Similar Rights. Effective
Technological Measures means those measures that, in the absence of proper
authority, may not be circumvented under laws fulfilling obligations under
Article 11 of the WIPO Copyright Treaty adopted on December 20, 1996, and/or
similar international agreements. Exceptions and Limitations means fair use,
fair dealing, and/or any other exception or limitation to Copyright and Similar
Rights that applies to Your use of the Licensed Material. License Elements means
the license attributes listed in the name of a Creative Commons Public License.
The License Elements of this Public License are Attribution, NonCommercial, and
ShareAlike. Licensed Material means the artistic or literary work, database, or
other material to which the Licensor applied this Public License. Licensed
Rights means the rights granted to You subject to the terms and conditions of
this Public License, which are limited to all Copyright and Similar Rights that
apply to Your use of the Licensed Material and that the Licensor has authority
to license. Licensor means the individual(s) or entity(ies) granting rights
under this Public License. NonCommercial means not primarily intended for or
directed towards commercial advantage or monetary compensation. For purposes of
this Public License, the exchange of the Licensed Material for other material
subject to Copyright and Similar Rights by digital file-sharing or similar means
is NonCommercial provided there is no payment of monetary compensation in
connection with the exchange. Share means to provide material to the public by
any means or process that requires permission under the Licensed Rights, such as
reproduction, public display, public performance, distribution, dissemination,
communication, or importation, and to make material available to the public
including in ways that members of the public may access the material from a
place and at a time individually chosen by them. Sui Generis Database Rights
means rights other than copyright resulting from Directive 96/9/EC of the
European Parliament and of the Council of 11 March 1996 on the legal protection
of databases, as amended and/or succeeded, as well as other essentially
equivalent rights anywhere in the world. You means the individual or entity
exercising the Licensed Rights under this Public License. Your has a
corresponding meaning. Section 2 - Scope.

License grant. Subject to the terms and conditions of this Public License, the
Licensor hereby grants You a worldwide, royalty-free, non-sublicensable, non-
exclusive, irrevocable license to exercise the Licensed Rights in the Licensed
Material to: reproduce and Share the Licensed Material, in whole or in part, for
NonCommercial purposes only; and produce, reproduce, and Share Adapted Material
for NonCommercial purposes only. Exceptions and Limitations. For the avoidance
of doubt, where Exceptions and Limitations apply to Your use, this Public
License does not apply, and You do not need to comply with its terms and
conditions. Term. The term of this Public License is specified in Section 6(a).
Media and formats; technical modifications allowed. The Licensor authorizes You
to exercise the Licensed Rights in all media and formats whether now known or
hereafter created, and to make technical modifications necessary to do so. The
Licensor waives and/or agrees not to assert any right or authority to forbid You
from making technical modifications necessary to exercise the Licensed Rights,
including technical modifications necessary to circumvent Effective
Technological Measures. For purposes of this Public License, simply making
modifications authorized by this Section 2(a)(4) never produces Adapted
Material. Downstream recipients. Offer from the Licensor - Licensed Material.
Every recipient of the Licensed Material automatically receives an offer from
the Licensor to exercise the Licensed Rights under the terms and conditions of
this Public License. Additional offer from the Licensor - Adapted Material.
Every recipient of Adapted Material from You automatically receives an offer
from the Licensor to exercise the Licensed Rights in the Adapted Material under
the conditions of the Adapter's License You apply. No downstream restrictions.
You may not offer or impose any additional or different terms or conditions on,
or apply any Effective Technological Measures to, the Licensed Material if doing
so restricts exercise of the Licensed Rights by any recipient of the Licensed
Material. No endorsement. Nothing in this Public License constitutes or may be
construed as permission to assert or imply that You are, or that Your use of the
Licensed Material is, connected with, or sponsored, endorsed, or granted
official status by, the Licensor or others designated to receive attribution as
provided in Section 3(a)(1)(A)(i). Other rights.

Moral rights, such as the right of integrity, are not licensed under this Public
License, nor are publicity, privacy, and/or other similar personality rights;
however, to the extent possible, the Licensor waives and/or agrees not to assert
any such rights held by the Licensor to the limited extent necessary to allow
You to exercise the Licensed Rights, but not otherwise. Patent and trademark
rights are not licensed under this Public License. To the extent possible, the
Licensor waives any right to collect royalties from You for the exercise of the
Licensed Rights, whether directly or through a collecting society under any
voluntary or waivable statutory or compulsory licensing scheme. In all other
cases the Licensor expressly reserves any right to collect such royalties,
including when the Licensed Material is used other than for NonCommercial
purposes. Section 3 - License Conditions.

Your exercise of the Licensed Rights is expressly made subject to the following
conditions.

Attribution.

If You Share the Licensed Material (including in modified form), You must:

retain the following if it is supplied by the Licensor with the Licensed
Material: identification of the creator(s) of the Licensed Material and any
others designated to receive attribution, in any reasonable manner requested by
the Licensor (including by pseudonym if designated); a copyright notice; a
notice that refers to this Public License; a notice that refers to the
disclaimer of warranties; a URI or hyperlink to the Licensed Material to the
extent reasonably practicable; indicate if You modified the Licensed Material
and retain an indication of any previous modifications; and indicate the
Licensed Material is licensed under this Public License, and include the text
of, or the URI or hyperlink to, this Public License. You may satisfy the
conditions in Section 3(a)(1) in any reasonable manner based on the medium,
means, and context in which You Share the Licensed Material. For example, it may
be reasonable to satisfy the conditions by providing a URI or hyperlink to a
resource that includes the required information. If requested by the Licensor,
You must remove any of the information required by Section 3(a)(1)(A) to the
extent reasonably practicable. ShareAlike. In addition to the conditions in
Section 3(a), if You Share Adapted Material You produce, the following
conditions also apply.

The Adapter's License You apply must be a Creative Commons license with the same
License Elements, this version or later, or a BY-NC-SA Compatible License. You
must include the text of, or the URI or hyperlink to, the Adapter's License You
apply. You may satisfy this condition in any reasonable manner based on the
medium, means, and context in which You Share Adapted Material. You may not
offer or impose any additional or different terms or conditions on, or apply any
Effective Technological Measures to, Adapted Material that restrict exercise of
the rights granted under the Adapter's License You apply. Section 4 - Sui
Generis Database Rights.

Where the Licensed Rights include Sui Generis Database Rights that apply to Your
use of the Licensed Material:

for the avoidance of doubt, Section 2(a)(1) grants You the right to extract,
reuse, reproduce, and Share all or a substantial portion of the contents of the
database for NonCommercial purposes only; if You include all or a substantial
portion of the database contents in a database in which You have Sui Generis
Database Rights, then the database in which You have Sui Generis Database Rights
(but not its individual contents) is Adapted Material, including for purposes of
Section 3(b); and You must comply with the conditions in Section 3(a) if You
Share all or a substantial portion of the contents of the database. For the
avoidance of doubt, this Section 4 supplements and does not replace Your
obligations under this Public License where the Licensed Rights include other
Copyright and Similar Rights. Section 5 - Disclaimer of Warranties and
Limitation of Liability.

Unless otherwise separately undertaken by the Licensor, to the extent possible,
the Licensor offers the Licensed Material as-is and as-available, and makes no
representations or warranties of any kind concerning the Licensed Material,
whether express, implied, statutory, or other. This includes, without
limitation, warranties of title, merchantability, fitness for a particular
purpose, non-infringement, absence of latent or other defects, accuracy, or the
presence or absence of errors, whether or not known or discoverable. Where
disclaimers of warranties are not allowed in full or in part, this disclaimer
may not apply to You. To the extent possible, in no event will the Licensor be
liable to You on any legal theory (including, without limitation, negligence) or
otherwise for any direct, special, indirect, incidental, consequential,
punitive, exemplary, or other losses, costs, expenses, or damages arising out of
this Public License or use of the Licensed Material, even if the Licensor has
been advised of the possibility of such losses, costs, expenses, or damages.
Where a limitation of liability is not allowed in full or in part, this
limitation may not apply to You. The disclaimer of warranties and limitation of
liability provided above shall be interpreted in a manner that, to the extent
possible, most closely approximates an absolute disclaimer and waiver of all
liability. Section 6 - Term and Termination.

This Public License applies for the term of the Copyright and Similar Rights
licensed here. However, if You fail to comply with this Public License, then
Your rights under this Public License terminate automatically. Where Your right
to use the Licensed Material has terminated under Section 6(a), it reinstates:

automatically as of the date the violation is cured, provided it is cured within
30 days of Your discovery of the violation; or upon express reinstatement by the
Licensor. For the avoidance of doubt, this Section 6(b) does not affect any
right the Licensor may have to seek remedies for Your violations of this Public
License. For the avoidance of doubt, the Licensor may also offer the Licensed
Material under separate terms or conditions or stop distributing the Licensed
Material at any time; however, doing so will not terminate this Public License.
Sections 1, 5, 6, 7, and 8 survive termination of this Public License. Section 7
- Other Terms and Conditions.

The Licensor shall not be bound by any additional or different terms or
conditions communicated by You unless expressly agreed. Any arrangements,
understandings, or agreements regarding the Licensed Material not stated herein
are separate from and independent of the terms and conditions of this Public
License. Section 8 - Interpretation.

For the avoidance of doubt, this Public License does not, and shall not be
interpreted to, reduce, limit, restrict, or impose conditions on any use of the
Licensed Material that could lawfully be made without permission under this
Public License. To the extent possible, if any provision of this Public License
is deemed unenforceable, it shall be automatically reformed to the minimum
extent necessary to make it enforceable. If the provision cannot be reformed, it
shall be severed from this Public License without affecting the enforceability
of the remaining terms and conditions. No term or condition of this Public
License will be waived and no failure to comply consented to unless expressly
agreed to by the Licensor. Nothing in this Public License constitutes or may be
interpreted as a limitation upon, or waiver of, any privileges and immunities
that apply to the Licensor or You, including from the legal processes of any
jurisdiction or authority.


For Government Employees

Definitions

"Applied Bioformatics Laboratory" and "ABiL" refer to the Applied
Bioinformatics Laboratory, a public-private partnership between IHRC Inc. and
Georgia Institute of Technology, 950 Atlantic Drive / Engineered Biosystems
Building, Room 2200 /Atlanta, GA 30332.  "You" or "Your" refers to the
state, federal or other governmental institutions and their respective
employees or contractors using this program.
"Employees" refers to persons employed directly by You. "Contractors" refers to
persons hired under contract from a thirdparty for the expressed purpose of
performing research for You. "Program" refers to abil_URDOcaller, and the
source code and algorithms contained within this file.  "Research" generally
refers to normal work duties that may require the use of this Program.
"License" refers to this text and the grants given and limitation imposed on You
for the usage of this Program.  "Python" refers to the Python programming
language (www.python.org), and is under separate terms and conditions.

Additional License Rights and Restrictions

ABiL grants You nonexclusive rights for the use of this Program, subject to the
restrictions within this License.  This Program is provided as-is with no
warranty or assumption of responsibility from ABiL. You use this Program at Your
risk and this Program is not suitable for medical diagnosis or other healthcare
application.  ABiL grants You the right to utilize this Program (a) internally
for the expressed purpose of conducting Research and related duties; (b) to make
copies of and distribute this Program internally for the expressed purpose of
conducting Research. You may allow your Contractors to utilize this program given
(a) You agree to enforce this License; (b) Contractors are under contract in
such a way that ABIL's intellectual property rights are protected; and (c) any
redistribution is to computer systems owned and operated by You and utilized
by Contractors for the sole, expressed purpose of conducting Research for You.

You may not:

1)  remove or otherwise obfuscate this License or other marks and
    attribution to ABiL contained within this Program;
2)  redistribute or otherwise transfer this Program to systems or persons not
    employed or contracted by You;
3)  make use of the source code, algorithms, or other intellectual property
    contained within this Program for any purpose other than running this Program;
4)  cause or allow the unauthorized usage, reproduction or transfer of this
    Program.

ABiL reserves all rights not expressly granted within this License.  If You wish
to use this Program for any purpose other than those expressly granted under this
License, You must obtain additional grants from ABiL.  Additional grants may
require additional contractual agreements, fees and licenses.  ABiL claims no
copyright or other license over the Python modules utilized by this Program.
Please refer to the license terms for each module for additional information.


"""

#predict part starts here
############################################################
HELP_TEXT_SMALL = "help"
HELP_TEXT = "HELP"
TMPDIR = tempfile.mkdtemp()
def batch_tool(fdir, kmer):
    """
    Function   : batch_tool
    Input      : Directory name, paired only, k value
    Output     : STs and allelic profiles for each FASTQ file
    Description: Processes all FASTQ files present in the input directory
    """
    if not directory.endswith('/'):
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
        kCount = single_sample_tool(fastq1_processed, fastq2_processed, paired, kmer, rawCounts)
    shutil.rmtree(TMPDIR)
    return kCount
def single_sample_tool(fastq1, fastq2, paired, k, results):
    """
    Function   : single_sample_tool
    Input      : fastq file 1 and 2, paired or single, k value, output dictionary
    Output     : STs and allelic profiles for each FASTQ file
    Description: Processes both FASTQ files passed to the function
    """
    # pp.pprint(kmerDict)
    if reads:
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
    if reads:
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
                if reads: read_file.write('\n'.join('{}'.format(l) for l in lines))
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
    load_config(config)
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

def print_results(results, output_filename, overwrite, timeDisp):
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
        if listMode:
            if not os.path.isfile(fList):
                print(HELP_TEXT_SMALL)
                print("Error: List file ("+fList+") does not exist!")
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

buildDB = False
predict = False
output_filename = None
batch = False
listMode = False
overwrite = False
paired = False
fastq1 = None
fastq2 = None
user_k = False
config = None
timeDisp = False
reads = False
dbPrefix = 'kmer'
log = ''
k = 35
fuzzy = 5
coverage = True
directory = None

"""Input arguments"""
options, remainder = getopt.getopt(sys.argv[1:], 'o:x1:2:kbd:phP:c:trva:', [
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
        output_filename = arg
    elif opt in ('-x', '--overwrite'):
        overwrite = True
    elif opt in '--buildDB':
        buildDB = True
    elif opt in ('-P', '--prefix'):
        dbPrefix = arg
    elif opt in '--predict':
        predict = True
    elif opt in ('-c', '--config'):
        config = arg
    elif opt in '-k':
        user_k = True
        try:
            k = int(arg)
        except ValueError:
            print("Error: Enter a numerical k value.")
            exit(0)
        # Check to make sure the arg is an int.
    elif opt in ('-l', '--list'):
        listMode = True
        fList = arg
    elif opt in ('-1', '--fastq1'):
        fastq1 = arg
    elif opt in ('-2', '--fastq2'):
        fastq2 = arg
    elif opt in ('-d', '--dir', '--directory'):
        directory = os.path.abspath(arg)
        batch = True
    elif opt in '-t':
        timeDisp = True
    elif opt in '-a':
        log = arg
    elif opt in '-r':
        reads = True
    elif opt in '-v':
        print(VERSION)
        exit(0)
    elif opt in ('-h', '--help'):
        print(HELP_TEXT)
        exit(0)
check_params(buildDB, predict, config, k, batch, directory, fastq1, fastq2, dbPrefix)
if buildDB:
    try:
        if not log:
            log = dbPrefix+'.log'
    except TypeError:
        log = 'kmer.log'
    logging.basicConfig(filename=log, level=logging.DEBUG,
                        format='%(asctime)s %(message)s', datefmt='%m/%d/%Y %I:%M:%S %p')
    if os.path.isfile(config):
        print("Info: Making DB for k = ", k)
        print("Info: Making DB with prefix =", dbPrefix)
        print("Info: Log file written to ", log)
        make_custom_db(config, k, dbPrefix)
    else:
        print("Error: The input config file "+config +" does not exist.")
elif predict:
    try:
        if not log:
            log = dbPrefix+'.log'
    except TypeError:
        log = 'kmer.log'
    logging.basicConfig(filename=log, level=logging.DEBUG,
                        format='%(asctime)s %(message)s', datefmt='%m/%d/%Y %I:%M:%S %p')
    logging.debug("==============================================================================")
    logging.debug(f"Command: {' '.join(sys.argv)}")
    logging.debug("Starting Marker Prediction")
    logging.debug(f"Temporary directory: {TMPDIR}")
    load_module(k, dbPrefix)
    rawCounts = {}
    global kCount
    kCount = {}
    if batch:
        rawCounts = batch_tool(directory, k)
    else:
        rawCounts = single_sample_tool(fastq1, fastq2, paired, k, rawCounts)
    weightCounts = weight_profile(rawCounts, weightDict)
    print_results(weightCounts, output_filename, overwrite, timeDisp)
else:
    print("Error: Please select the mode")
    print("--buildDB (for database building) or --predict (for marker discovery)")
