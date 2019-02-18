#!/bin/bash

echo "This script needs wget"
echo "Downloading databases from http://klif.uu.nl/download/metagenomics_db/"
echo "Please be patient. This is a large download"
echo
echo "Download Kaiju DB"
wget http://klif.uu.nl/download/metagenomics_db/DB_RefSeq_Virus.tar.gz
echo "Downloading Groot resistome db"
wget http://klif.uu.nl/download/metagenomics_db/groot-db-32.tar.gz
echo "Downloading Kraken2 db"
wget http://klif.uu.nl/download/metagenomics_db/kraken2_bacteria_db.tar.gz
echo "Downloading Mash and BioBloomfilter DBs"
wget http://klif.uu.nl/download/metagenomics_db/mash_db.tar.gz
echo "Downloading KMA DB"
wget http://klif.uu.nl/download/metagenomics_db/KMA_ResFinder.tar.gz

echo "Untarring files"
tar xzvf DB_RefSeq_Virus.tar.gz
tar xzvf groot-db-32.tar.gz
tar xzvf kraken2_bacteria_db.tar.gz
tar xzvf mash_db.tar.gz
tar xzvf KMA_ResFinder.tar.gz

echo "Done"
echo ""
echo "do rm *.tar.gz to remove temporary download files"
exit 1
