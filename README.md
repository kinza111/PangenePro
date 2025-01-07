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
PangenePro supports Linux system due to the software dependencies and requires Python 3.8 or higher version. 
1.	Blast+
$ wget https://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/ncbi-blast-2.16.0+-x64-linux.tar.gz
$ tar zxvf ncbi-blast-2.16.0+-x64-linux.tar.gz
# Add Blast+ tools to your PATH
$ export PATH=/path/to/blast+/bin:$PATH
2.	Seqtk
$ git clone https://github.com/lh3/seqtk.git
$ cd seqtk; make
3.	Diamond 
$ wget http://github.com/bbuchfink/diamond/releases/download/v2.1.10/diamond-linux64.tar.gz
$ tar xzf diamond-linux64.tar.gz
# Add Blast+ tools to your PATH
$ export PATH=/path/to/diamond/bin:$PATH
4.	orthAgogue
$ git clone git://github.com/lh3/seqtk.git
$ cd seqtk
$ make
$ cp seqtk /usr/local/bin
5.	MCL
$ conda install bioconda::mcl
$ apt-get install mcl
Usage 
Running the Application
Clone the repository:
git clone https://github.com/yourusername/PangenePro.git
cd PangenePro
Set up your environment and install dependencies. The query protein sequence, proteome, and genomes files must be in fasta format (.fasta, .faa, .fa) format. The subject genome must be assembled at chromosomal level. The annotation file must be in .gff, .gff3, or .gtf format. 
Run the command to execute pipeline by setting query and subject files
./PangenePro.sh

Example Run
Replace the command in PangenePro.sh script with the example dataset and run it.
$ python get_genes.py AtCRK_ref.faa Ahy_genome.fna Ahy_ proteome.faa Ahy_gff.gff Aip_genome.fna Aip_ proteome.faa Aip_gff.gff Adu_genome.fna Adu_ proteome.faa Adu_gff.gff
