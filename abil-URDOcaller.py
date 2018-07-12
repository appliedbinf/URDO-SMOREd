#!/usr/bin/env python
import getopt
import sys
import logging
import os
import time
import ast
import gzip
import re
import tempfile
import shutil
import xml.etree.ElementTree as ET
try:
    from urllib.request import urlopen, urlretrieve
except ImportError:
        from urllib import urlopen, urlretrieve
import argparse
import math
from collections import defaultdict
from itertools import islice
import operator
import sys
version = """ abil-genecaller v1 (updated : March 29, 2018) """
"""
abil-genecaller free for academic users and requires permission before any 
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
performing research for You. "Program" refers to abil-genecaller, and the
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


The program has 3 basic modes :
    mainTool: for single sample (both single and paired end)
    batchTool: for multiple samples stored at a common location (both single and paired end samples)
    listTool: for multiple samples with location information stored in a list (both single and paired end samples)
predict part starts here
"""

############################################################
# Function   : batchTool
# Input      : Directory name, paired or single, k value
# Output     : STs and allelic profiles for each FASTQ file
# Description: Processes all FASTQ files present in the input
#              directory
#############################################################
helpTextSmall="help"
helpText="HELP"
def batchTool(fdir, paired, k):
    fileList = []
    if not dir.endswith('/'):
        fdir += '/'
    for inputFile in os.listdir(fdir):
        if paired is True:
            if inputFile.endswith('1.fastq') or inputFile.endswith('1.fq') or inputFile.endswith('1.fq.gz') or inputFile.endswith('1.fastq.gz'):
                fastq1 = fdir+inputFile
                fastq2 = fdir+inputFile.replace('1.', '2.')
                fileList.append((fastq1, fastq2))
        else:
            if inputFile.endswith('.fastq') or inputFile.endswith('.fq') or inputFile.endswith('.fq.gz') or inputFile.endswith('.fastq.gz'):
                fastq1 = fdir + inputFile
                fileList.append(fastq1)
    results = multiSampleTool(fileList, paired, k)
    return results
#############################################################
# Function   : singleSampleTool
# Input      : fastq file 1 and 2, paired or single, k value, output dictionary
# Output     : STs and allelic profiles for each FASTQ file
# Description: Processes both FASTQ files passed to the function
#############################################################
def singleSampleTool(fastq1, fastq2, paired, k, results):
    if reads is True:
        readFileName = fastq1.split('/')[-1].split('.')[0][:-1] + '_reads.fq'
        global readFile
        readFile = open(readFileName, 'w+')
    sName = fastq1.split('_')[1]
    msg = f"Preprocessing: working with {fastq1} and {fastq2}"
    logging.debug(msg)
    global alleleCount
    alleleCount = {}
    global kCount
    kCount = {}
    t1 = time.time()
    tmpdir = tempfile.mkdtemp()
    logging.debug(f"Temporary directory: {tmpdir}")
    logging.debug(f"Preprocessing: Merging {fastq1} and {fastq2}")
    vsearch_cmd = f"vsearch --fastq_mergepairs {fastq1} --reverse {fastq2} --fastqout {tmpdir}/reads.fq 2> {fastq1.split('/')[-1].split('.')[0][:-1]}.merge.log 1>/dev/null"
    logging.debug("Preprocessing: VSEARCH command\n\t{}".format(vsearch_cmd))
    os.system(vsearch_cmd)
    singleFileTool(tmpdir, k, sName)
    shutil.rmtree(tmpdir)
    if kCount == {}:
        string = f"No k-mer matches were found for the sample {fastq1} and {fastq2}\n\tProbable cause of the error:  low quality data/too many N's in the data"
        logging.error(f"ERROR: {string}")
        print(string)
#           exit(0)
    profileCount = alleleCount
    logging.debug("singleSampleTool : weightedProfile start")
    weightedProfile = weightedProf(profileCount, weightDict)
    logging.debug("singleSampleTool : weightedProfile finished")
    logging.debug("singleSampleTool : getMaxCount start")
    # finalProfile = getMaxCount(weightedProfile, fileName)
    # print(weightedProfile)
    logging.debug("singleSampleTool : getMaxCount end")
    st = 0
    if profileFile != '':
        logging.debug("singleSampleTool : findST start")
        # st = findST(finalProfile, stProfile)
        logging.debug("singleSampleTool : findST end")
    if reads is True:
        readFile.close()
    return kCount
#############################################################
# Function   : singleFileTool
# Input      : fastq file, k value
# Output     : Edits a global dictionary - results
# Description: Processes the single fastq file
#######################################q######################
def singleFileTool(fastq, k, sName):
    fastq = "/".join([fastq,"reads.fq"])
    msg = f"Analysis: Begin processing merged reads ({fastq})"
    logging.debug(msg)
    if os.path.isfile(fastq):
        logging.debug(f"Analysis: Kmer counting using k={k}")
        finalProfile = {}
        t1 = time.time()
        # fileExplorer(fastq, k)
        if sName not in kCount:
            kCount[sName] = {}
        f = open(fastq)
        for lines in iter(lambda: tuple(islice(f, 4)), ()):
            if len(lines) < 4:
                dbg = "ERROR: Input file is truncated.  Please verify the input FASTQ files are correct"
                logging.debug(dbg)
            try:
                if len(lines[1]) < k:
                    m1 = f"ERROR: Read ID: {[0][1:]}  is  length {len(lines[1])} in {file} smaller than {k}"
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
            # print("{}\t{}\t{}".format(len(firstKmer),len(midKmer),len(lastKmer)))
            if firstKmer in kmerDict[k] or midKmer in kmerDict[k] or lastKmer in kmerDict[k]:
                countKmers(lines[1], k, sName)
                if reads: readFile.write('\n'.join('{}'.format(l) for l in lines))
        # print(kCount)
            t3 = time.time()
    else:
        logging.error(f"ERROR: File does not exist: {fastq}")
#############################################################
# Function   : goodReads
# Input      : sequence read, k, step size
# Output     : Edits the count of global variable alleleCount
# Description: Increment the count for each k-mer match
#############################################################
def countKmers(read, k, sName):
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
    # print(max(alleleCount['allele'].items(), key=operator.itemgetter(1))[0])
    # print(alleleCount['allele'][max(alleleCount['allele'].items(), key=operator.itemgetter(1))[0]])
    # print()
    alleleNumber = max(alleleCount['allele'].items(), key=operator.itemgetter(1))[0]
    alleleKcount = alleleCount['allele'][max(alleleCount['allele'].items(), key=operator.itemgetter(1))[0]]
    if alleleNumber in kCount[sName]:
        # kCount[sName][alleleNumber] += alleleKcount
        kCount[sName][alleleNumber] += 1
    else:
        kCount[sName][alleleNumber] = 1
#############################################################
# Function   : weightedProf
# Input      : allele count global var, weight factors
# Output/Desc: Normalizes alleleCount by weight factor
#############################################################
def weightedProf(alleleCount, weightDict):
    logging.debug("weightedProf")
    weightedDict = {}
    for loc in alleleCount:
        weightedDict[loc] = {}
        for allele in alleleCount[loc]:
            if loc in weightDict:
                if allele in weightDict[loc]:
                    weightedDict[loc][allele] = (alleleCount[loc][allele] / weightDict[loc][allele])
                else:
                    weightedDict[loc][allele] = alleleCount[loc][allele]
            else:
                weightedDict[loc][allele] = alleleCount[loc][allele]
    return weightedDict

#############################################################
# Function   : loadModule
# Input      : k value and prefix of the DB file
# Output     : Updates the DB dictionary variables
# Description: Used in loading the DB as set of variables
#              by calling other functions
#############################################################
def loadModule(k, dbPrefix):
    global dbFile
    dbFile = dbPrefix+'_'+str(k)+'.txt'
    global weightFile
    weightFile = dbPrefix+'_weight.txt'
    global profileFile
    profileFile = dbPrefix+'_profile.txt'
    global kmerDict
    kmerDict = {}
    kmerDict[k] = loadKmerDict(dbFile)
    global weightDict
    weightDict = loadWeightDict(weightFile)
    global stProfile
    stProfile = loadSTfromFile(profileFile)
    loadConfig(config)
#############################################################
# Function   : loadSTfromFile
# Input      : profile definition file
# Output     : Updates the DB dictionary variables
# Description: Used in loading the DB as set of variables
#############################################################
def loadSTfromFile(profileF):
    with open(profileF, 'r') as definitionFile:
        st = {}
        index = {}
        lines = definitionFile.readlines()
        heads = lines[0].rstrip().split('\t')
        for locus in heads:
            index[locus] = heads.index(locus)
        for line in lines:
            pro = line.rstrip().split('\t')
            l = {}
            for locus in heads[1:]:
                try:
                    l[locus] = pro[index[locus]]
                except LookupError:
                    logging.debug("ERROR while loading ST")
                    pass
            st[pro[0]] = l
    return st
#############################################################
# Function   : loadKmerDict
# Input      : DB prefix
# Output     : Updates the DB dictionary variables
# Description: Used in loading the DB as set of variables
#############################################################
def loadKmerDict(dbFile):
    kmerTableDict = {}
    with open(dbFile, 'r') as kmerTableFile:
        lines = kmerTableFile.readlines()
        for line in lines:
            array = line.rstrip().rsplit('\t')
            kmerTableDict[array[0]] = {}
            kmerTableDict[array[0]][array[1]] = array[2][1:-1].rsplit(',')
    return kmerTableDict
#############################################################
# Function   : loadWeightDict
# Input      : Weight file prefix
# Output     : Updates the DB dictionary variables
# Description: Used in loading the DB as set of variables
#############################################################
def loadWeightDict(weightFile):
    weightDict = {}
    with open(weightFile, 'r') as weightTableFile:
        lines = weightTableFile.readlines()
        for line in lines:
            array = line.rstrip().rsplit('\t')
            try:
               (loc, allele) =  array[0].replace('-','_').rsplit('_',1)
            except ValueError:
                print("Error : Allele name in locus file should be seperated by '_' or '-'")
                exit(0) 
            if loc not in weightDict:
                weightDict[loc] = {}
            weightDict[loc][allele] = float(array[1])
    return weightDict
#############################################################
# Function   : loadConfig
# Input      : config file path from getopts
# Output     : Updates configDict
# Description: Used to find allele fasta files for getCoverage
#############################################################
def loadConfig(config):
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

#############################################################
# Function   : printResults
# Input      : results, output file, overwrite?
# Output     : Prints on the screen or in a file
# Description: Prints the results in the format asked by the user
#############################################################
def printResults(results, output_filename, overwrite, timeDisp):
    if output_filename != None:
        if overwrite is False:
            outfile = open(output_filename, "a")
        else:
            outfile = open(output_filename, "w")
    heading = "Sample"
    for gene in sorted(stProfile):
        heading += '\t' + gene
    if output_filename != None:
        outfile.write(heading)
        outfile.write('\n')
    else:
        print(heading)
    # print(results)
    for s in results:
        sample = s
        for gene in sorted(stProfile):
            # print(gene)
            # print(str(stProfile[gene]['allele']))
            if stProfile[gene]['allele'] in results[s]:
                sample += "\t" + str(results[s][stProfile[gene]['allele']])
            else:
                sample += "\t0" 
        
        if output_filename != None:
            outfile.write(sample)
            outfile.write('\n')
        else:
            print(sample)
"""Predict part ends here"""
"""Build DB part starts"""
"""Returns the reverse complement of the sequence"""
def reverseComplement(seq):
    seqU = seq.upper()
    seq_dict = {'A':'T', 'T':'A', 'G':'C', 'C':'G', 'Y':'R', 'R':'Y', 'S':'S', 'W':'W', 'K':'M', 'M':'K', 'N':'N'}
    try:
        return "".join([seq_dict[base] for base in reversed(seqU)])
    except Exception:
        strn = "Reverse Complement Error:" + seqU
        logging.debug(strn)
        pass
#############################################################
# Function   : getFastaDict
# Input      : locus file name
# Output     : dictionary with all the allele sequences
# Description: Stores each allele sequence in a dictionary
#############################################################
def getFastaDict(fullLocusFile):
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
#############################################################
# Function   : formKmerDB
# Input      : configuration file, k value, output prefix
# Output     : abil-genecaller DB
# Description: Constructs the k-mer DB in both strand orientation
#############################################################
def formKmerDB(configDict, k, output_filename):
    dbFileName = output_filename+'_'+str(k)+'.txt'
    weightFileName = output_filename+'_weight.txt'
    kmerDict = {}
    mean = {}
    for locus in configDict['loci']:
        msgs = "formKmerDB :" +locus
        logging.debug(msgs)
        fastaDict = getFastaDict(configDict['loci'][locus])
        sum = 0
        n = 0
        for allele in list(fastaDict.keys()):
            seq = fastaDict[allele]['sequence'].strip()
            l = len(seq)
            sum += l
            n += 1
            try:
                (loc, num) =  allele.replace('-','_').rsplit('_',1)
            except ValueError:
                print("Error : Allele name in locus file should be seperated by '_' or '-'")
                exit(0) 
            splitId =  allele.replace('-','_').rsplit('_',1)
            i = 0
            while i+k <= l:
                kmer = seq[i:i+k]
                revCompKmer = reverseComplement(kmer)
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
        mean[locus] = sum/n*1.0
    with open(dbFileName, 'w') as kfile:
        for key in kmerDict:
            for key1 in kmerDict[key]:
                string = str(key)+'\t'+str(key1)+'\t'+str(kmerDict[key][key1]).replace(" ", "")+'\n'
                kfile.write(string)
    with open(weightFileName, 'w') as wfile:
        for locus in configDict['loci']:
            fastaDict = getFastaDict(configDict['loci'][locus])
            for allele in list(fastaDict.keys()):
                splitId = allele.split('_')
                seq = fastaDict[allele]['sequence']
                l = len(seq)
                fac = (l/mean[locus])
                s = allele  + '\t' + str(fac) + '\n'
                if fac > 1.05 or fac < 0.95:
                    wfile.write(s)
"""Copies the profile definition file as a new file"""
def copyProfileFile(profileDict, output_filename):
    profileFileName = output_filename+'_profile.txt'
    with open(profileDict['profile']) as f:
        lines = f.readlines()
        with open(profileFileName, "w") as f1:
            f1.writelines(lines)
#############################################################
# Function   : makeCustomDB
# Input      : configuration file, k value, output prefix
# Output     : None
# Description: Processes the config file and calls the relevant
#              function
#############################################################
def makeCustomDB(config, k, output_filename):
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
    formKmerDB(configDict, k, output_filename)
    copyProfileFile(configDict['profile'], output_filename)
"""Build DB part ends"""
"""Check Parameters"""
def checkParams(buildDB, predict, config, k, batch, dir, fastq1, fastq2, dbPrefix):
    if predict is True and buildDB is True:
        print(helpTextSmall)
        print("Select either predict or buildDB module")
        exit(0)
    if predict is False and buildDB is False:
        print(helpTextSmall)
        print("Select either predict or buildDB module")
        exit(0)
    if predict is True:
        if config is None and coverage is True:
            print(helpTextSmall)
            print("Config parameter is required.")
            exit(0)
        if not os.path.isfile(dbPrefix+'_'+str(k)+'.txt'):
            print(helpTextSmall)
            print("DB file does not exist : ", dbPrefix, '_', str(k), '.txt or change DB prefix.')
            exit(0)
        if not os.path.isfile(dbPrefix+'_weight.txt'):
            print(helpTextSmall)
            print("DB file does not exist : ", dbPrefix, '_weight.txt or change DB prefix.')
            exit(0)
        if not os.path.isfile(dbPrefix+'_profile.txt'):
            print(helpTextSmall)
            print("DB file does not exist : ", dbPrefix, '_profile.txt or change DB prefix.')
            exit(0)
        if listMode is True:
            if not os.path.isfile(fList):
                print(helpTextSmall)
                print("Error: List file ("+fList+") does not exist!")
                exit(0)
        elif batch is True:
            if not os.path.isdir(dir):
                print(helpTextSmall)
                print("Error: Directory ("+dir+") does not exist!")
                exit(0)
        elif paired is True:
            if not os.path.isfile(fastq1):
                print(helpTextSmall)
                print("Error: FASTQ file ("+fastq1+") does not exist!")
                exit(0)
            if not os.path.isfile(fastq2):
                print(helpTextSmall)
                print("Error: FASTQ file ("+fastq2+") does not exist!")
                exit(0)
        elif paired is False:
            if not os.path.isfile(fastq1):
                print(helpTextSmall)
                print("Error: FASTQ file ("+fastq1+") does not exist!")
                exit(0)
    if buildDB is True:
        try:
            if not os.path.isfile(config):
                print(helpTextSmall)
                print("Error: Configuration file ("+config+") does not exist!")
                exit(0)
        except Exception:
            print(helpTextSmall)
            print("Error: Specify Configuration file")
            exit(0)

"""The Program Starts Execution Here"""
"""Default Params"""

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
        dir = arg
        batch = True
    elif opt in '-t':
        timeDisp = True
    elif opt in '-a':
        log = arg
    elif opt in '-r':
        reads = True
    elif opt in '-v':
        print(version)
        exit(0)
    elif opt in ('-h', '--help'):
        print(helpText)
        exit(0)
checkParams(buildDB, predict, config, k, batch, dir, fastq1, fastq2, dbPrefix)
if buildDB is True:
    try:
        if not log:
            log = dbPrefix+'.log'
    except TypeError:
        log = 'kmer.log'
    logging.basicConfig(filename=log, level=logging.DEBUG, format='%(asctime)s %(message)s', datefmt='%m/%d/%Y %I:%M:%S %p')
    if os.path.isfile(config):
        print("Info: Making DB for k = ", k)
        print("Info: Making DB with prefix =", dbPrefix)
        print("Info: Log file written to ", log)
        makeCustomDB(config, k, dbPrefix)
    else:
        print("Error: The input config file "+config +" does not exist.")
elif predict is True:
    try:
        if not log:
            log = dbPrefix+'.log'
    except TypeError:
        log = 'kmer.log'
    logging.basicConfig(filename=log, level=logging.DEBUG, format='%(asctime)s %(message)s', datefmt='%m/%d/%Y %I:%M:%S %p')
    logging.debug("================================================================================")
    logging.debug("Command: {}".format(" ".join(sys.argv)))
    logging.debug("Starting Marker Prediction")
    loadModule(k, dbPrefix)
    if batch is True:
        results = batchTool(dir, paired, k)
    else:
        results = {}
        results = singleSampleTool(fastq1, fastq2, paired, k, results)
    printResults(results, output_filename, overwrite, timeDisp)
else:
    print("Error: Please select the mode: buildDB (for database building) or predict (for marker discovery)")



