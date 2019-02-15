#!/usr/bin/python

import gzip, os, numpy, argparse, linecache, sys, multiprocessing
from numpy import median
from Bio import SeqIO
from __main__ import *
if (args.screening) == True:
  from filtering import *
if (args.trimming) == True:
  from QC import *
  
def KMA(input_R1, input_R2, input_SE, input_filtered_R1, input_filtered_R2, input_filtered_SE, prefix):
  if (args.pe) and BioBloomCategorizer == False:
    os.system("{} -ipe {} {} -t_db {} -o {} -t {}".format(os.path.join(script_dir, "binaries/kma"), os.path.join(output_dir, input_R1), os.path.join(output_dir, input_R2), os.path.join(script_dir, "db/KMA_ResFinder/ResFinder"), os.path.join(output_dir, prefix + "_kma"), args.threads))
  
  if not (args.pe) and BioBloomCategorizer == False:
    os.system("{} -i {} -t_db {} -o {} -t {}".format(os.path.join(script_dir, "binaries/kma"), os.path.join(output_dir, input_SE), os.path.join(script_dir, "db/KMA_ResFinder/ResFinder"), os.path.join(output_dir, prefix + "_kma"), args.threads))
  
  if (args.pe) and BioBloomCategorizer == True:
    os.system("{} -ipe {} {} -t_db {} -o {} -t {}".format(os.path.join(script_dir, "binaries/kma"), os.path.join(output_dir, input_filtered_R1), os.path.join(output_dir, input_filtered_R2), os.path.join(script_dir, "db/KMA_ResFinder/ResFinder"), os.path.join(output_dir, prefix + "_kma"), args.threads))
  
  if not (args.pe) and BioBloomCategorizer == True:
    os.system("{} -i {} -t_db {} -o {} -t {}".format(os.path.join(script_dir, "binaries/kma"), os.path.join(output_dir, input_filtered_SE), os.path.join(script_dir, "db/KMA_ResFinder/ResFinder"), os.path.join(output_dir, prefix + "_kma"), args.threads))

print "Running KMA"

KMA(trimmed_R1, trimmed_R2, trimmed_SE, filtered_R1, filtered_R2, filtered_SE, args.prefix)

print "Done"