[![CircleCI](https://circleci.com/gh/ar0ch/URDO-predictor.svg?style=svg&circle-token=e5cfa9a42dfa40c1a26051677c8ee81d58bfbc01)](https://circleci.com/gh/ar0ch/URDO-predictor) 
![PyPI - Python Version](https://img.shields.io/pypi/pyversions/Django.svg)

# Readme for abil\_URDOcaller
=============================================================================================
### Usage
```
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
[-v]
[-h][--help]
```
==============================================================================================

### There are two steps to predicting ST using abil\_URDOcaller.
1. Create DB : `abil_URDOcaller.py --buildDB`
2. Predict : `abil_URDOcaller --predict`

#### 1. `abil_URDOcaller.py --buildDB`

**Synopsis:**
`abil_URDOcaller.py --buildDB -c <config file> -k <kmer length(optional)> -P <DB prefix(optional)>`  
config file : is a tab delimited file which has the information for typing scheme ie loci, its multifasta file and profile definition file.
    Format : 
```
[loci]  
locus1    locusFile1
locus2    locusFile2
[profile]
profile   profile_file
```
kmer length : is the kmer length for the db. Note, while processing this should be smaller than the read length.  
    - We suggest kmer lengths of 35, 66 depending on the read length.  
DB prefix(optional) : holds the information for DB files to be created and their location. This module creates 3 files with this prefix.  
    - You can use a folder structure with prefix to store your db at particular location.

**Required arguments**
`--buildDB`  
    Identifier for build db module  
`-c,--config = <configuration file>`  
    Config file in the format described above.   

**Optional arguments**  
`-k = <kmer length>`
    Kmer size for which the db has to be formed(Default k = 35). Note the tool works best with kmer length in between 35 and 66
  for read lengths of 55 to 150 bp. Kmer size can be increased accordingly. It is advised to keep lower kmer sizes 
  if the quality of reads is not very good.  
`-P,--prefix = <prefix>`  
  Prefix for db and log files to be created(Default = kmer). Also you can specify folder where you want the dbb to be created.  
`-a`
    File location to write build log  
`-h,--help`  
    Prints the help manual for this application  

 --------------------------------------------------------------------------------------------
 
#### 2. `abil_URDOcaller.py --predict`
  
`abil_URDOcaller --predict` : can run in two modes
  1) single sample (default mode)
  2) multi-sample: run abil\_URDOcaller for all the samples in a folder (for a particular species)

**Synopsis**
`abil_URDOcaller.py --predict -1 <fastq file> -2 <fastq file> -d <directory location> - -P <DB prefix(optional)> -k <kmer length(optional)> -o <output file> -x`

**Required arguments** 
`--predict`
  Identifier for predict module
`-c,--config = <configuration file>`
  Config file in the format described above. 
  
**Optional arguments**
`-1,--fastq1 = <fastq1_filename>`  
Path to first fastq file for paired end sample and path to the fastq file for single end file.
    - Should have extension fastq or fq.
`-2,--fastq2 = <fastq2_filename>`  
  Path to second fastq file for paired end sample.
    - Should have extension fastq or fq.
`-d,--dir,--directory = <directory>`  
  Directory containing paired end read files for multi-sample prediction  
`-k = <kmer_length>`  
  Kmer length for which the db was created (Default k = 35).  
`-o,--output = <output_filename>`  
  Prints the output to a file instead of stdout.  
`-P,--prefix = <prefix>`  
  Prefix using which the db was created (Defaults = kmer). The location of the db could also be provided.  
`-r`  
  A FASTQ file is generated in the current directory for each sample containing reads with kmer matches.  
`-R, --readsdir = < output directory >`  
  A FASTQ file is generated in the specified directory for each sample containing reads with kmer matches. 
`-v`  
  Prints the version of the software.  
`-x,--overwrite`  
  By default abil\_URDOcaller appends the results to the output\_filename if same name is used.
  This argument overwrites the previously specified output file.
`-h,--help`
  Prints the help manual for this application  

