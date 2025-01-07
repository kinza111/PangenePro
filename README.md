# PangenePro: an automated pipeline for rapid identification and classification of gene family members
PangenePro is a pipeline designed for the quick and efficient identification of gene family members from multiple related genomes. It further classifies the identified members into core, accessory, and unique gene sets, collectively referred to as pangenes. 

## Features

- **Automation:** It streamlines the previously hectic and error-prone genome-wide gene family identification steps, providing faster and more accurate identification.
- **Refinement of Initial Alignment Results:** Provides the filtered results after the initial identification by removing the redundant accessions and the isoforms.
- **Domain Analysis:** Analyzes the presence of specific domains to provide the confirmed and validated gene family members.
- **Pangenes level Identification:** Identifies gene family members from multiple genomes simultaneously and covers the intraspecies diversity by covering the presence absence variations across multiple genomes, thereby, providing identified members being classified into core, accessory, and unique genes.
- **Results Summary and Graphs:** Efficiently processes and analyzes data, and provides detailed summary and graphs for better visualization of results. 

## Installations

### Prerequisites 
- Linux system due to the software dependencies
- Python 3.8 or higher version
- [Blast+](https://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/ncbi-blast-2.16.0+-x64-linux.tar.gz)
```
$ wget https://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/ncbi-blast-2.16.0+-x64-linux.tar.gz
$ tar zxvf ncbi-blast-2.16.0+-x64-linux.tar.gz
$ ./configure --prefix=/path/to/installation
$ make
$ make install
# Add Blast+ tools to your PATH
$ export PATH=/path/to/installation/:$PATH
```
- [Seqtk](https://github.com/lh3/seqtk)
```
$ git clone https://github.com/lh3/seqtk.git
$ cd seqtk; make
```
- [Diamond](https://github.com/bbuchfink/diamond)
```
$ wget http://github.com/bbuchfink/diamond/releases/download/v2.1.10/diamond-linux64.tar.gz
$ tarxzf diamond-linux64.tar.gz
$ ./configure --prefix=/path/to/installation
$ make
$ make install
# Add Diamond to your PATH
$ export PATH=/path/to/installation/:$PATH
```
- [orthAgogue](https://github.com/samyeaman/orthagogue)
```
$ git clone git://github.com/lh3/orthAgogue.git
$ cd orthAgogue; make
```
- [MCL](https://github.com/micans/mcl)
```
$ conda install bioconda::mcl
$ apt-get install mcl
```
# Usage
## Running the Pipeline 
Clone the repository:
```sh
$ git clone https://github.com/yourusername/PangenePro.git
$ cd PangenePro/Example/
$ $ python get_genes.py AtCRK_ref.faa Ahy_genome.fna Ahy_ proteome.faa Ahy_gff.gff Aip_genome.fna Aip_ proteome.faa Aip_gff.gff Adu_genome.fna Adu_ proteome.faa Adu_gff.gff
```
# Contributing/Contact Us
