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
  name_output = args.inp[0].replace(".fastq.gz", "") + "_25000_reads" + ".fastq"
  write_reads(name_output)
  
  lengte_reads= []
  for seq_record in SeqIO.parse(name_output, "fastq"):
    lengte_reads.append(int((len(seq_record))))
    med = median(lengte_reads)

def groot(input_R1, input_R2, prefix, input_SE, input_filtered_R1, input_filtered_R2, input_filtered_SE):
  if (args.pe) and BioBloomCategorizer == False:
    os.system('{} align -p {} -i {} -f {} {} | groot report -c 0.2 > {}_groot_report.txt'.format(os.path.join(script_dir, "binaries/groot"), args.threads, os.path.join(script_dir, "db/groot-db-32/groot_index_{}".format(int(med))), input_R1, input_R2, prefix))
  
  
  if not (args.pe) and BioBloomCategorizer == False:
    os.system('{} align -p {} -i {} -f {} | groot report -c 0.2 > {}_groot_report.txt'.format(os.path.join(script_dir, "binaries/groot"), args.threads, os.path.join(script_dir, "db/groot-db-32/groot_index_{}".format(int(med))), input_SE, prefix))
  
  
  if (args.pe) and BioBloomCategorizer == True:
    os.system('{} align -p {} -i {} -f {} {} | groot report -c 0.2 > {}_groot_report.txt'.format(os.path.join(script_dir, "binaries/groot"), args.threads, os.path.join(script_dir, "db/groot-db-32/groot_index_{}".format(int(med))), input_filtered_R1, input_filtered_R2, prefix))
  
  
  if not (args.pe) and BioBloomCategorizer == True:
    os.system('{} align -p {} -i {} -f {} | groot report -c 0.2 > {}_groot_report.txt'.format(os.path.join(script_dir, "binaries/groot"), args.threads, os.path.join(script_dir, "db/groot-db-32/groot_index_{}".format(int(med))), input_filtered_SE, prefix))
  

rem = med % 10

if rem < 5:
  med = int(med/10)*10
else:
  med = int((med + 10) / 10) * 10
  

print "Running groot"

if noPrefix == True:

  groot(trimmed_R1, trimmed_R2, args.inp[0].replace(".fastq.gz", ""), trimmed_SE, filtered_R1, filtered_R2, filtered_SE)


if noPrefix == False:

  groot(trimmed_R1, trimmed_R2, args.prefix, trimmed_SE, filtered_R1, filtered_R2, filtered_SE)


print "Done"

if os.path.exists(args.inp[0].replace(".fastq.gz", "") + "_25000_reads.fastq"):
  os.remove(args.inp[0].replace(".fastq.gz", "") + "_25000_reads.fastq")
if os.path.exists(args.inp[1].replace(".fastq.gz", "") + "_25000_reads.fastq"):
  os.remove(args.inp[1].replace(".fastq.gz", "") + "_25000_reads.fastq")  