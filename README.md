# Project-UU
Repo voor interne veranderingen voor main release bioinformatica pipeline project UU

Update 6-12-2018:
Fixed a crash when inputting files with their filepath (/path/to/file/example.fastq)
Fixed an issue were some files would be put in program folder instead of input file folder
All files now are put in the folder where the input files are located

Update 10-12-2018:
User can now choose own destination folder (mandatory)

Update 11-12-2018
Fixed an issue were the wrong input file would be chosen for Kaiju, Kraken2 and GROOT
Fixed an issue were GROOT wouldn't start on all PCs

To do:
Not all uncessary files are correctly deleted at the end when the program finishes
When not using a prefix the program will crash

