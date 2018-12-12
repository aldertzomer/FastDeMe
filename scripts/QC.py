#!/usr/bin/python

import gzip, os, numpy, argparse, linecache, sys, multiprocessing
from numpy import median
from Bio import SeqIO
from __main__ import *

lines = []
lengte_reads = []
quality = []

#Used for calculating a round number from the median 
def iround(x):
  return int(round(x) - .5) + (x > 0)

#Exracting reads from the input file to calculate the statistics for the QC_calc module. Also includes .gz recognition.
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

        
#Writing reads from the extr_25000_reads module to a file
def write_reads(in_file):
  with open(in_file, 'w') as f:
    for item in lines:  
        f.write("%s" % item) 

#Makes lists of read length and PHRED score of individual bases so it can be processed in the QC_calc module
def QC(in_file):
  for seq_record in SeqIO.parse(in_file, "fastq"):
    lengte_reads.append(int((len(seq_record))))
  
  for record in SeqIO.parse(in_file, "fastq"):
    quality.extend(record.letter_annotations["phred_quality"])

#Calculate median, first quartile, third quartile and average PHRED score of the reads. Average PHRED score is one number, so this does not output the average PHRED score per read. 
def QC_calc(reads, phred):
  global med
  med = median(reads)
  global first_quartile
  first_quartile = numpy.percentile(reads, 25)
  global third_quartile 
  third_quartile = numpy.percentile(reads, 75)
  global av_phred 
  av_phred = numpy.mean(phred)

#Writing the calculated statistics (for untrimmed data) to an output file
def write_results_untrimmed(file_name):
  with open(args.prefix + "_QC_results.txt", "a") as f:
    f.write("Statistics for file %s \n" % file_name)
    f.write("Median untrimmed: %s \n" % med)
    f.write("First quartile untrimmed: %s \n" % first_quartile)
    f.write("Third quartile untrimmed: %s \n" % third_quartile)
    f.write("Average PHRED score untrimmed: %s \n" % av_phred)
    f.write("\n")

#Writing the calculated statistics (for trimmed data) to an output file    
def write_results_trimmed(file_name):
  with open(args.prefix + "_QC_results.txt", "a") as f:
    f.write("Statistics for file %s after trimming \n" % file_name)
    f.write("Median trimmed: %s \n" % med)
    f.write("First quartile trimmed: %s \n" % first_quartile)
    f.write("Third quartile trimmed: %s \n" % third_quartile)
    f.write("Average PHRED score trimmed: %s \n" % av_phred)
    f.write("\n")

#Writing certain parts of the trimming report outputted by trim-galore to an output file    
def write_trimming_report_SE(file_name): 
  total_reads_unfiltered = linecache.getline(os.path.basename(file_name.replace(".fastq.gz", "") + "_trimming_report.json"), 4)
  total_reads_filtered = linecache.getline(os.path.basename(file_name.replace(".fastq.gz", "") + "_trimming_report.json"), 14)
  total_bases_unfiltered = linecache.getline(os.path.basename(file_name.replace(".fastq.gz", "") + "_trimming_report.json"), 5)
  total_bases_filtered = linecache.getline(os.path.basename(file_name.replace(".fastq.gz", "") + "_trimming_report.json"), 15)
  adapters = linecache.getline(os.path.basename(file_name.replace(".fastq.gz", "") + "_trimming_report.json"), 40)
  
  with open(args.prefix + "_QC_results.txt", "a") as f:
    f.write("Untrimmed results for file %s \n" % args.inp[0])
    f.write(total_reads_unfiltered.replace('"', "").replace(",", "").replace("_", " ") + "\n")
    f.write(total_bases_unfiltered.replace('"', "").replace(",", "").replace("_", " ") + "\n")
    f.write(adapters.replace('"', "").replace("_", " ") + "\n")
    f.write("\n")
    f.write("Trimmed results for file %s \n" % args.inp[0])
    f.write(total_reads_filtered.replace('"', "").replace(",", "").replace("_", " ") + "\n")
    f.write(total_bases_filtered.replace('"', "").replace(",", "").replace("_", " ") + "\n")  

def write_trimming_report_PE(file_name):
  total_reads_unfiltered = linecache.getline(os.path.basename(file_name.replace(".fastq.gz", "") + "_trimming_report.json"), 52)
  total_reads_filtered = linecache.getline(os.path.basename(file_name.replace(".fastq.gz", "") + "_trimming_report.json"), 142)
  total_bases_unfiltered = linecache.getline(os.path.basename(file_name.replace(".fastq.gz", "") + "_trimming_report.json"), 53)
  total_bases_filtered = linecache.getline(os.path.basename(file_name.replace(".fastq.gz", "") + "_trimming_report.json"), 143)
  adapters = linecache.getline(os.path.basename(file_name.replace(".fastq.gz", "") + "_trimming_report.json"), 48)
  
  total_reads_unfiltered_2 = linecache.getline(os.path.basename(file_name.replace(".fastq.gz", "") + "_trimming_report.json"), 232)
  total_reads_filtered_2 = linecache.getline(os.path.basename(file_name.replace(".fastq.gz", "") + "_trimming_report.json"), 322)
  total_bases_unfiltered_2 = linecache.getline(os.path.basename(file_name.replace(".fastq.gz", "") + "_trimming_report.json"), 233)
  total_bases_filtered_2 = linecache.getline(os.path.basename(file_name.replace(".fastq.gz", "") + "_trimming_report.json"), 323)
  adapters_2 = linecache.getline(os.path.basename(file_name.replace(".fastq.gz", "") + "_trimming_report.json"), 49)
  
  with open(args.prefix + "_QC_results.txt", "a") as f:
    f.write("Untrimmed results for file %s \n" % args.inp[0])
    f.write(total_reads_unfiltered.replace('"', "").replace(",", "").replace("_", " ") + "\n")
    f.write(total_bases_unfiltered.replace('"', "").replace(",", "").replace("_", " ") + "\n")
    f.write(adapters.replace('"', "").replace("_", " ") + "\n")
    f.write("\n")
    f.write("Trimmed results for file %s \n" % args.inp[0])
    f.write(total_reads_filtered.replace('"', "").replace(",", "").replace("_", " ") + "\n")
    f.write(total_bases_filtered.replace('"', "").replace(",", "").replace("_", " ") + "\n")
    f.write("\n")
    f.write("Untrimmed results for file %s \n" % args.inp[1])
    f.write(total_reads_unfiltered_2.replace('"', "").replace(",", "").replace("_", " ") + "\n")
    f.write(total_bases_unfiltered_2.replace('"', "").replace(",", "").replace("_", " ") + "\n")
    f.write(adapters_2.replace('"', "").replace("_", " ") + "\n")
    f.write("\n")
    f.write("Trimmed results for file %s \n" % args.inp[1])
    f.write(total_reads_filtered_2.replace('"', "").replace(",", "").replace("_", " ") + "\n")
    f.write(total_bases_filtered_2.replace('"', "").replace(",", "").replace("_", " ") + "\n")

#Running fastp definition

def fastp(input_SE, output_SE, input_R1, input_R2, output_R1, output_R2):
  if not (args.pe) and (args.threads) <= 16:
    os.system('{} --thread {} -i {} -o {} -l {} -j {} -q 20 >/dev/null 2>/dev/null'.format(os.path.join(script_dir, "binaries/fastp"), args.threads, input_SE, os.path.join(output_dir, output_SE), iround(med*0.8), os.path.basename(args.inp[0].replace(".fastq.gz", "")) + "_trimming_report.json"))
    
  if not (args.pe) and (args.threads) > 16:
    os.system('{} --thread {} -i {} -o {} -l {} -j {} -q 20 >/dev/null 2>/dev/null'.format(os.path.join(script_dir, "binaries/fastp"), int(16), input_SE, os.path.join(output_dir, output_SE), iround(med*0.8), os.path.basename(args.inp[0].replace(".fastq.gz", "")) + "_trimming_report.json"))
  
  if (args.pe) and (args.threads) <= 16:
    os.system('{} --thread {} -i {} -I {} -o {} -O {} -l {} -j {} -q 20 >/dev/null 2>/dev/null'.format(os.path.join(script_dir, "binaries/fastp"), args.threads, input_R1, input_R2, os.path.join(output_dir, output_R1), os.path.join(output_dir, output_R2), iround(med*0.8), os.path.basename(args.inp[0].replace(".fastq.gz", "")) + "_trimming_report.json"))
  
  if (args.pe) and (args.threads) > 16:
    os.system('{} --thread {} -i {} -I {} -o {} -O {} -l {} -j {} -q 20 >/dev/null 2>/dev/null'.format(os.path.join(script_dir, "binaries/fastp"), int(16), input_R1, input_R2, os.path.join(output_dir, output_R1), os.path.join(output_dir, output_R2), iround(med*0.8), os.path.basename(args.inp[0].replace(".fastq.gz", "")) + "_trimming_report.json"))
  


#Extracting 25000 reads

if not (args.pe):
  extr_25000_reads(args.inp[0])
  name_output = os.path.join(output_dir, os.path.basename(args.inp[0].replace(".fastq.gz", "")) + "_25000_reads.fastq")
  write_reads(name_output)
if (args.pe):
  extr_25000_reads(args.inp[0])
  name_output = os.path.join(output_dir, os.path.basename(args.inp[0].replace(".fastq.gz", "")) + "_25000_reads.fastq")
  write_reads(name_output)
  lines = []
  extr_25000_reads(args.inp[1])
  name_output_1 = os.path.join(output_dir, os.path.basename(args.inp[1].replace(".fastq.gz", "")) + "_25000_reads.fastq")
  write_reads(name_output_1)


#Run fastp
print "Running fastp"
for seq_record in SeqIO.parse(name_output, "fastq"):
  lengte_reads.append(int((len(seq_record))))
med = median(lengte_reads)

fastp(args.inp[0], trimmed_SE, args.inp[0], args.inp[1], trimmed_R1, trimmed_R2)

lines = []
lengte_reads = []
quality = []

print "Done"

#Extracting 25000 reads from the trimmed files

if not (args.pe):
  if noPrefix == True:
    extr_25000_reads(os.path.join(output_dir, os.path.basename(args.inp[0].replace(".fastq.gz", ""))) + "_trimmed.fastq.gz")
    name_output_trimmed = os.path.join(output_dir, os.path.basename(args.inp[0].replace(".fastq.gz", "")) + "_25000_reads" + "_trimmed.fastq")
    write_reads(name_output_trimmed)
  if noPrefix == False:
    extr_25000_reads(os.path.join(output_dir, args.prefix) + "_trimmed.fastq.gz")
    name_output_trimmed = os.path.join(output_dir, os.path.basename(args.inp[0].replace(".fastq.gz", "")) + "_25000_reads" + "_trimmed.fastq")
    write_reads(name_output_trimmed)
if (args.pe):
  if noPrefix == True:
    extr_25000_reads(os.path.join(output_dir, os.path.basename(args.inp[0].replace(".fastq.gz", ""))) + "_trimmed_1.fastq.gz")
    name_output_trimmed = os.path.join(output_dir, os.path.basename(args.inp[0].replace(".fastq.gz", "")) + "_25000_reads" + "_trimmed.fastq")
    write_reads(name_output_trimmed)
    lines = []
    extr_25000_reads(os.path.join(output_dir, os.path.basename(args.inp[1].replace(".fastq.gz", ""))) + "_trimmed_2.fastq.gz")
    name_output_trimmed1 = os.path.join(output_dir, os.path.basename(args.inp[1].replace(".fastq.gz", "")) + "_25000_reads" + "_trimmed.fastq")
    write_reads(name_output_trimmed1)
  if noPrefix == False:
    extr_25000_reads(os.path.join(output_dir, args.prefix) + "_trimmed_1.fastq.gz")
    name_output_trimmed = os.path.join(output_dir, os.path.basename(args.inp[0].replace(".fastq.gz", "")) + "_25000_reads" + "_trimmed.fastq")
    write_reads(name_output_trimmed)
    lines = []
    extr_25000_reads(os.path.join(output_dir, args.prefix) + "_trimmed_2.fastq.gz")
    name_output_trimmed1 = os.path.join(output_dir, os.path.basename(args.inp[1].replace(".fastq.gz", "")) + "_25000_reads" + "_trimmed.fastq")
    write_reads(name_output_trimmed1)    
  

#Calculating and writing the statistics to an output file
if not (args.pe):
  QC(name_output)
  QC_calc(lengte_reads, quality)
  write_results_untrimmed(args.inp[0]) 
  QC(name_output_trimmed)
  QC_calc(lengte_reads, quality)
  write_results_trimmed(args.inp[0])
  write_trimming_report_SE(args.inp[0])

if (args.pe):
  QC(name_output)
  QC_calc(lengte_reads, quality)
  write_results_untrimmed(args.inp[0])
  lengte_reads = []
  quality = [] 
  QC(name_output_trimmed)
  QC_calc(lengte_reads, quality)
  write_results_trimmed(args.inp[0])
  lengte_reads = []
  quality = []
  QC(name_output_1)
  QC_calc(lengte_reads, quality)
  write_results_untrimmed(args.inp[1])
  lengte_reads = []
  quality = []
  QC(name_output_trimmed1)
  QC_calc(lengte_reads, quality)
  write_results_trimmed(args.inp[1])
  write_trimming_report_PE(args.inp[0])


#Deleting the intermediate 25000 reads files    
if not (args.pe):
  if os.path.exists(os.path.join(output_dir, os.path.basename(args.inp[0].replace(".fastq.gz", "")) + "_25000_reads.fastq")):
    os.remove(os.path.join(output_dir, os.path.basename(args.inp[0].replace(".fastq.gz", "")) + "_25000_reads.fastq"))
  if os.path.exists(os.path.join(output_dir, os.path.basename(args.inp[0].replace(".fastq.gz", "")) + "_25000_reads_trimmed.fastq")):
    os.remove(os.path.join(output_dir, os.path.basename(args.inp[0].replace(".fastq.gz", "")) + "_25000_reads_trimmed.fastq"))
  if os.path.exists("fastp.html"):
    os.remove("fastp.html")
  if os.path.exists(os.path.basename(args.inp[0].replace(".fastq.gz", "")) + "_trimming_report.json"):
    os.remove(os.path.basename(args.inp[0].replace(".fastq.gz", "")) + "_trimming_report.json")
if (args.pe):
  if os.path.exists(os.path.join(output_dir, os.path.basename(args.inp[0].replace(".fastq.gz", "")) + "_25000_reads.fastq")):
    os.remove(os.path.join(output_dir, os.path.basename(args.inp[0].replace(".fastq.gz", "")) + "_25000_reads.fastq"))
  if os.path.exists(os.path.join(output_dir, os.path.basename(args.inp[0].replace(".fastq.gz", "")) + "_25000_reads_trimmed.fastq")):
    os.remove(os.path.join(output_dir, os.path.basename(args.inp[0].replace(".fastq.gz", "")) + "_25000_reads_trimmed.fastq")) 
  if os.path.exists(os.path.join(output_dir, os.path.basename(args.inp[1].replace(".fastq.gz", "")) + "_25000_reads.fastq")):
    os.remove(os.path.join(output_dir, os.path.basename(args.inp[1].replace(".fastq.gz", "")) + "_25000_reads.fastq"))  
  if os.path.exists(os.path.join(output_dir, os.path.basename(args.inp[1].replace(".fastq.gz", "")) + "_25000_reads_trimmed.fastq")):
    os.remove(os.path.join(output_dir, os.path.basename(args.inp[1].replace(".fastq.gz", "")) + "_25000_reads_trimmed.fastq"))
  if os.path.exists("fastp.html"):
    os.remove("fastp.html")
  if os.path.exists(os.path.basename(args.inp[0].replace(".fastq.gz", "")) + "_trimming_report.json"):
    os.remove(os.path.basename(args.inp[0].replace(".fastq.gz", "")) + "_trimming_report.json")


