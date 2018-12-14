# Project-UU

A fast, easy solution for metagenomic data analysis.

## Installation

The program can be downloaded as an archive or with the following git command:

`git clone <github_link>`

Next, the databases should be downloaded. This can be done by running `getdb.sh` located in the `db/` directory. This will download the four databases that are needed for the program to run and unzip them. Please make sure enough space is available on the drive, since the combined size of the databases is ~165 GB. The combined download size is ~116 GB.

After downloading, the program can be invoked with `wrapper.py`.

## Usage
The program has two mandatory arguments, `--inp` and `--output`. To use the basic version of the program, the following command can be ran:

`./wrapper.py --inp file.fastq.gz --output /path/to/output/folder/`

This will only result in the file getting trimmed and screened for host contamination, which are on by default.

To run perform taxonomic classification using Kaiju and resistome analysis using GROOT with paired end files without host contamination screening, issue the following command:

`./wrapper.py --pe --inp file_R1.fastq.gz file_R2.fastq.gz --output /path/to/output/folder/ --kaiju --groot --screening`

The flag `--pe` is needed when using paired end files. `--kaiju` and `--groot` turn on the respective programs, and `--screening` will skip the screening process.

## Database information
### GROOT

The GROOT database consist of a mixture of the ResFinder, ARG-ANNOT and CARD databases. See the [GROOT documentation](https://groot-documentation.readthedocs.io/en/latest/groot-databases.html) for more details.

### Kaiju

The Kaiju database was made with assembled and annotated bacterial, archaeal and viral reference genomes from the NCBI RefSeq database.

### Kraken

The Kraken database was made with the complete bacterial/archaeal reference genomes from RefSeq.

### Mash/BioBloomCategorizer

The Mash database was made from the complete vertebrate_mammalian and vertebrate_other databases. Since the bloom filters BioBloomCategorizer uses for filtering the host reads are quite large, only filters for common host species are included in the standard database to reduce download size.

Species in standard database|GCF ID
--------------|----------------
*Bos indicus*|GCF_000247795.1
*Bos mutus*|GCF_000298355.1
*Bos taurus*|GCF_002263795.1
*Canis lupus familiaris*|GCF_000002285.3
*Capra hircus*|GCF_001704415.1
*Chinchilla lanigera*|GCF_000276665.1
*Equus caballus*|GCF_002863925.1
*Felis catus*|GCF_000181335.3
*Gorilla gorilla gorilla*|GCF_000151905.2
*Homo sapiens*|GCF_000001405.38
*Mus musculus*|GCF_000001635.26
*Ovis aries*|GCF_000298735.2
*Ovis aries musimon*|GCF_000765115.1
*Pan troglodytes*|GCF_002880755.1
*Rattus norvegicus*|GCF_000001895.5, GCF_000002265.2
*Sus scrofa*|GCF_000003025.6
*Danio rerio*|GCF_000002035.6
*Gallus gallus*|GCF_000002315.5
*Meleagris gallopavo*|GCF_000146605.2

Bloom filters for the remaining species in vertebrate_mammalian and vertebrate_other can be downloaded separately. In cases were the contaminating factor is not known in advance, the program will output the GCF ID of the missing species so it can be downloaded and included by the user. This download consist of a folder containing a `.bf` and `.txt` file. The folder should be put in either `db/mash_db/vertebrate/refseq/vertebrate_mammalian` or `db/mash_db/vertebrate/refseq/vertebrate_other` depending on the species.

Update 6-12-2018:
Fixed a crash when inputting files with their filepath (/path/to/file/example.fastq)
Fixed an issue were some files would be put in program folder instead of input file folder
All files now are put in the folder where the input files are located

Update 10-12-2018:
User can now choose own destination folder (mandatory)

Update 11-12-2018
Fixed an issue were the wrong input file would be chosen for Kaiju, Kraken2 and GROOT
Fixed an issue were GROOT wouldn't start on all PCs

Update 12-12-2018
Fixed a crash of the program when not using a prefix
Fixed an issue were not all unnecessary would be deleted after program finishes

To do:
Program crashes when trimming is turned off

