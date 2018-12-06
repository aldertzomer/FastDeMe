#!/usr/bin/python

import gzip, os, numpy, argparse, linecache, sys, re, multiprocessing
from numpy import median
from Bio import SeqIO
from __main__ import *
if (args.trimming) == True:
  from QC import *
if (args.screening) == True:
  from filtering import *

if (args.trimming) == False:
  lines = []
  lengte_reads = []
  def extr_25000_reads(in_file):
    if(in_file.find(".gz") != -1):
      i = 1
      with gzip.open(in_file, "r") as f:  
        while i <= 100000:
          line = f.readline()
          lines.append(line)
          i += 1
    else:
      i = 1
      with open(in_file, "r") as f:  
        while i <= 100000:
          line = f.readline()
          lines.append(line)
          i += 1
      in_file = in_file.replace(".fastq", ".fastq.gz")
      
  def write_reads(in_file):
    with open(in_file, 'w') as f:
      for item in lines:  
          f.write("%s" % item)   

  extr_25000_reads(args.inp[0])
  name_output = os.path.basename(args.inp[0].replace(".fastq.gz", "")) + "_25000_reads" + ".fastq"
  write_reads(name_output)
  
  for seq_record in SeqIO.parse(name_output, "fastq"):
    lengte_reads.append(int((len(seq_record))))
    med = median(lengte_reads)


print "Running Kraken2"

def Kraken2(prefix, input_R1, input_R2, input_SE, input_filtered_R1, input_filtered_R2, input_filtered_SE):
  if (args.pe) and BioBloomCategorizer == False:
    os.system('{} --db {} --threads {} --report {}_kraken_report.txt --output {}_kraken_output.txt --paired {} {}'.format(os.path.join(script_dir, "binaries/kraken2/kraken2"), os.path.join(script_dir, "db/kraken2_bacteria_db"), args.threads, prefix, os.path.basename(args.inp[0].replace(".fastq.gz", "")), input_R1, input_R2))
  
  if not (args.pe) and BioBloomCategorizer == False:
    os.system('{} --db {} --threads {} --report {}_kraken_report.txt --output {}_kraken_output.txt {}'.format(os.path.join(script_dir, "binaries/kraken2/kraken2"), os.path.join(script_dir, "db/kraken2_bacteria_db"), args.threads, prefix, os.path.basename(args.inp[0].replace(".fastq.gz", "")), input_SE))  
  
  
  if (args.pe) and BioBloomCategorizer == True:
    os.system('{} --db {} --threads {} --report {}_kraken_report.txt --output {}_kraken_output.txt --paired {} {}'.format(os.path.join(script_dir, "binaries/kraken2/kraken2"), os.path.join(script_dir, "db/kraken2_bacteria_db"), args.threads, prefix, os.path.basename(args.inp[0].replace(".fastq.gz", "")), input_filtered_R1, input_filtered_R2))
    
  if not (args.pe) and BioBloomCategorizer == True:
    os.system('{} --db {} --threads {} --report {}_kraken_report.txt --output {}_kraken_output.txt {}'.format(os.path.join(script_dir, "binaries/kraken2/kraken2"), os.path.join(script_dir, "db/kraken2_bacteria_db"), args.threads, prefix, os.path.basename(args.inp[0].replace(".fastq.gz", "")), input_filtered_SE))

def Bracken(prefix, tax_rank1, tax_rank2, read):
  os.system("{} -d {} -i {}_kraken_report.txt -o {}_{}.bracken -l {} -r {}".format(os.path.join(script_dir, "binaries/bracken"), os.path.join(script_dir, "db/kraken2_bacteria_db"), prefix, prefix, tax_rank2, tax_rank1, read))


Kraken2(args.prefix, trimmed_R1, trimmed_R2, trimmed_SE, filtered_R1, filtered_R2, filtered_SE)
  
  
if args.tax_rank == "phylum":
  args.tax_rank = "P"
if args.tax_rank == "class":
  args.tax_rank = "C"
if args.tax_rank == "order":
  args.tax_rank = "O"
if args.tax_rank == "family":
  args.tax_rank = "F"
if args.tax_rank == "genus":
  args.tax_rank = "G"
if args.tax_rank == "species":
  args.tax_rank = "S"



rem = med % 10

if rem < 5:
  med = int(med/10)*10
else:
  med = int((med + 10) / 10) * 10
  
  
if (args.tax_rank) is None:
  Bracken(args.prefix, "S", "species", int(med))
  Bracken(args.prefix, "P", "phylum", int(med))
  Bracken(args.prefix, "C", "class", int(med))
  Bracken(args.prefix, "O", "order", int(med))
  Bracken(args.prefix, "F", "family", int(med))
  Bracken(args.prefix, "G", "genus", int(med))

if (args.tax_rank) is not None:
  Bracken(args.prefix, args.tax_rank, args.tax_rank, int(med))
  

if os.path.exists(os.path.basename(args.inp[0].replace(".fastq.gz", "")) + "_kraken_output.txt"):
  os.remove(os.path.basename(args.inp[0].replace(".fastq.gz", "")) + "_kraken_output.txt")
if os.path.exists(os.path.basename(args.inp[0].replace(".fastq.gz", "")) + "_25000_reads.fastq"):
  os.remove(os.path.basename(args.inp[0].replace(".fastq.gz", "")) + "_25000_reads.fastq")
if os.path.exists(os.path.basename(args.inp[1].replace(".fastq.gz", "")) + "_25000_reads.fastq"):
  os.remove(os.path.basename(args.inp[1].replace(".fastq.gz", "")) + "_25000_reads.fastq")  


  


print "Done"

