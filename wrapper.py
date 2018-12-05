#!/usr/bin/python

import gzip, os, numpy, argparse, linecache, sys, multiprocessing
from numpy import median
from Bio import SeqIO

def cpu_threads(max_threads):  
  if multiprocessing.cpu_count() > max_threads:
    return max_threads
  else:
    return multiprocessing.cpu_count()

parser = argparse.ArgumentParser()
parser.add_argument("--pe", help="specify paired-end data, default is single end", \
action="store_true", default=False, required=False)
parser.add_argument("--inp", help="input files in fastq.gz format, if paired-end input both files with a space between them", nargs="+", action="store", required=True)
parser.add_argument("--threads", help="specify number of threads to be used, default is max available threads up to 16 threads", default=cpu_threads(16), required=False, type=int)
parser.add_argument("--kaiju", help="use kaiju for taxonomic identification", action="store_true", default=False)
parser.add_argument("--kraken", help="use kraken2 for taxonomic identification", action="store_true", default=False)
parser.add_argument("--groot", help="use groot for resistome analysis", action="store_true", default=False)
parser.add_argument("--tax_rank", help="set taxonomic rank for output. choose one: phylum, class, order, family, genus, species, default is all ranks.", action="store")
parser.add_argument("--prefix", help="prefix for all output files, default is name of input file(s)", action="store")
parser.add_argument("--trimming", help="turn off trimming with fastp", action="store_false", default=True)
parser.add_argument("--screening", help="turn off host contamination screening with mash and BioBloomCategorizer", action="store_false", default=True)
args = parser.parse_args()

if (args.pe):
  if len(args.inp) != 2:
    print "please use two input files"
    sys.exit()


if not (args.pe):
  if len(args.inp) != 1:
    print "please use one input file or use --pe for paired end files"
    sys.exit()

noPrefix = False
BioBloomCategorizer = False

if args.prefix is None:
  args.prefix = args.inp[0].replace('.fastq.gz', '')  
  noPrefix = True
  
if noPrefix == True:
    
  filtered_R1 = args.inp[0].replace(".fastq.gz", "") + "_noMatch_1.fastq"  
  filtered_R2 = args.inp[0].replace(".fastq.gz", "") + "_noMatch_2.fastq"     
  filtered_SE = args.inp[0].replace(".fastq.gz", "") + "_noMatch.fastq"
  
  trimmed_SE = args.inp[0].replace(".fastq.gz", "") + "_trimmed.fastq.gz"  
  trimmed_R1 = args.inp[0].replace(".fastq.gz", "") + "_trimmed_1.fastq.gz"
  trimmed_R2 = args.inp[1].replace(".fastq.gz", "") + "_trimmed_2.fastq.gz"
 
if noPrefix == False:  
  
  filtered_R1 = args.prefix + "_noMatch_1.fastq"  
  filtered_R2 = args.prefix + "_noMatch_2.fastq"         
  filtered_SE = args.prefix + "_noMatch.fastq"

  trimmed_SE = args.prefix + "_trimmed.fastq.gz"  
  trimmed_R1 = args.prefix + "_trimmed_1.fastq.gz"
  trimmed_R2 = args.prefix + "_trimmed_2.fastq.gz"

if (args.trimming) == False and (args.screening) == True:
  if noPrefix == True:
      
    filtered_R1 = args.inp[0].replace(".fastq.gz", "") + "_noMatch_1.fastq"  
    filtered_R2 = args.inp[0].replace(".fastq.gz", "") + "_noMatch_2.fastq"     
    filtered_SE = args.inp[0].replace(".fastq.gz", "") + "_noMatch.fastq"
    
    trimmed_SE = args.inp[0]
    trimmed_R1 = args.inp[0]
    trimmed_R2 = args.inp[1]
   
  if noPrefix == False:  
    
    filtered_R1 = args.prefix + "_noMatch_1.fastq"  
    filtered_R2 = args.prefix + "_noMatch_2.fastq"         
    filtered_SE = args.prefix + "_noMatch.fastq"
  
    trimmed_SE = args.inp[0]
    trimmed_R1 = args.inp[0]
    trimmed_R2 = args.inp[1]

if (args.trimming) == False and (args.screening) == False:
  trimmed_SE = args.inp[0]
  trimmed_R1 = args.inp[0]
  trimmed_R2 = args.inp[1] 
 
script_location = os.path.dirname(sys.argv[0])
script_dir = os.path.abspath(script_location) 
bin_path = os.path.join(script_dir, "scripts/")
sys.path.append(bin_path)

    
if (args.trimming) == True:
  import QC

if (args.screening) == True:
  import filtering

if (args.kaiju) == True:
  import Kaiju
  
if (args.kraken) == True:
  import Kraken

if (args.groot) == True:
  import groot


