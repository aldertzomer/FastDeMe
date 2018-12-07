#!/usr/bin/python

import gzip, os, numpy, argparse, linecache, sys, multiprocessing
from numpy import median
from Bio import SeqIO
from __main__ import *
if (args.screening) == True:
  from filtering import *
if (args.trimming) == True:
  from QC import *

print "Running Kaiju"

#Running Kaiju definition

def Kaiju(input_R1, input_R2, prefix, input_SE, input_filtered_R1, input_filtered_R2, input_filtered_SE):
  if (args.pe) and BioBloomCategorizer == False:
    os.system('{} -v -t {} -f {} -i {} -j {} -a mem -z {} -o {}_kaiju.txt'.format(os.path.join(script_dir, "binaries/kaiju"), os.path.join(script_dir, "db/DB_RefSeq_Virus/nodes.dmp"), os.path.join(script_dir, "db/DB_RefSeq_Virus/kaiju_db.fmi"), trimmed_R1, trimmed_R2, args.threads, os.path.basename(args.inp[0].replace(".fastq.gz", ""))))
    
  if not (args.pe) and BioBloomCategorizer == False:
    os.system('{} -v -t {} -f {} -i {} -a mem -z {} -o {}_kaiju.txt'.format(os.path.join(script_dir, "binaries/kaiju"), os.path.join(script_dir, "db/DB_RefSeq_Virus/nodes.dmp"), os.path.join(script_dir, "db/DB_RefSeq_Virus/kaiju_db.fmi"), input_SE, args.threads, os.path.basename(args.inp[0].replace(".fastq.gz", ""))))
    
  if (args.pe) and BioBloomCategorizer == True:
    os.system('{} -v -t {} -f {} -i {} -j {} -a mem -z {} -o {}_kaiju.txt'.format(os.path.join(script_dir, "binaries/kaiju"), os.path.join(script_dir, "db/DB_RefSeq_Virus/nodes.dmp"), os.path.join(script_dir, "db/DB_RefSeq_Virus/kaiju_db.fmi"), input_filtered_R1, input_filtered_R2, args.threads, os.path.basename(args.inp[0].replace(".fastq.gz", ""))))
  
  if not (args.pe) and BioBloomCategorizer == True:
    os.system('{} -v -t {} -f {} -i {} -a mem -z {} -o {}_kaiju.txt'.format(os.path.join(script_dir, "binaries/kaiju"), os.path.join(script_dir, "db/DB_RefSeq_Virus/nodes.dmp"), os.path.join(script_dir, "db/DB_RefSeq_Virus/kaiju_db.fmi"), input_filtered_SE, args.threads, os.path.basename(args.inp[0].replace(".fastq.gz", ""))))

def kaijuReport(tax_rank, prefix):
  os.system('{} -t {} -n {} -r {} -i {}_kaiju.txt -o {}_{}_kaiju_summary.txt'.format(os.path.join(script_dir, "binaries/kaijuReport"), os.path.join(script_dir, "db/DB_RefSeq_Virus/nodes.dmp"), os.path.join(script_dir, "db/DB_RefSeq_Virus/names.dmp"), tax_rank, os.path.basename(args.inp[0].replace(".fastq.gz", "")), prefix, tax_rank))

Kaiju(trimmed_R1, trimmed_R2, args.prefix, trimmed_SE, filtered_R1, filtered_R2, filtered_SE)  

if (args.tax_rank) is None:
  kaijuReport('species', args.prefix)
  kaijuReport('phylum', args.prefix)
  kaijuReport('class', args.prefix) 
  kaijuReport('order', args.prefix)
  kaijuReport('family', args.prefix)
  kaijuReport('genus', args.prefix)
if (args.tax_rank) is not None:
  kaijuReport(args.tax_rank, args.prefix)

    
if os.path.exists(os.path.basename(args.inp[0].replace(".fastq.gz", "")) + "_kaiju.txt"):
  os.remove(os.path.basename(args.inp[0].replace(".fastq.gz", "")) + "_kaiju.txt")    
  
print "Done"