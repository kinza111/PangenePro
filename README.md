# PangenePro: an automated pipeline for rapid identification and classification of gene family members
PangenePro is a pipeline designed for the quick and efficient identification of gene family members from multiple related genomes. It further classifies the identified members into core, accessory, and unique gene sets, collectively referred to as pangenes. 

![PangenePro](https://github.com/kinza111/PangenePro/blob/main/PangenePro.png)

## Features

- **Automation:** It streamlines the previously hectic and error-prone genome-wide gene family identification steps, providing faster and more accurate identification.
- **Refinement of Initial Alignment Results:** Provides the filtered results after the initial identification by removing the redundant accessions and the isoforms.
- **Domain Analysis:** Analyzes the presence of specific domains to provide the confirmed and validated gene family members.
- **Pangenes level Identification:** Identifies gene family members from multiple genomes simultaneously and covers the intraspecies diversity by covering the presence absence variations across multiple genomes, thereby, providing identified members being classified into core, accessory, and unique genes.
- **Results Summary and Graphs:** Efficiently processes and analyzes data, and provides detailed summary and graphs for better visualization of results. 

## Download and Usage
download the PangenePro using wget or through git and and uncompress the package. After downloading, put the bin directory into your PATH.
```
# download the PangenePro
wget https://github.com/kinza111/PangenePro/archive/refs/heads/main.zip
or
git clone git@github.com:kinza111/PangenePro.git
# Add the bin to PATH
$ export PATH=/path/to/PangenePro/bin/:$PATH
```
### Dependencies
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
## Input and Output Files 
### Input files 
PangenePro requires the one protein fasta file containing the sequence of the gene family of interest and the one or more subject genome sequence, their corresponding protein sequence fasta and annotation files. 
The genome and protein sequence file should be a fasta file with following format:
```
>chr1
ATCGATCG...

>ABC
ACKLMN...
```
File extension shoud be, '.fa', '.faa', '.fasta'. The prefix name of the sequence files will be used to indicate the temporary files, so we recommend using the organism's or cultivar names as file name 'cultivar.fa (like Arachis.fa)' while running PangenePro.

The annotattion file extension should be '.gff', '.gff3', or '.gtf'. 
```
ctg123 . gene            1000  9000  .  +  .  ID=gene00001;Name=EDEN
ctg123 . mRNA            1050  9000  .  +  .  ID=mRNA00001;Parent=gene00001;Name=EDEN.1
ctg123 . exon            1300  1500  .  +  .  ID=exon00001;Parent=mRNA00003
ctg123 . CDS             1201  1500  .  +  0  ID=cds00001;Parent=mRNA00001;Name=edenprotein.1
```
### Output files
PangenePro first performs the alignment of queried sequences and the subject proteomes and generates two files: 'Nr_NIs_Accessions.txt' and 'Nr_NIs_PSeqs.faa' for each corresponding genome. These files contain the identified and filtered gene family members' accessions and their fasta sequences. The domain analysis step generates the .xml files which are further filtered to generate final gene family members files, 'protein_domains_all_x.txt' and 'protein_domains_all_x_filtered_proteins.txt'. In the final step, 'Pangenes.txt' file contains the Pangenes clusters and 'gene_sets_summary.csv' contains the classified sets.

## Running the Pipeline 
A small dataset in the 'example' directory can be used to test the PangenePro pipeline. Since it is a small part of original sequence files, the results would not have any biological meaning. Access the dataset and replace the file names in the 'PangenePro.sh' file. 
```
#Specify the query protein sequence, proteome, genome, and anotation files
python get_genes.py AtCRK_ref.faa Ahy_genome.fna Ahy_ proteome.faa Ahy_gff.gff Aip_genome.fna Aip_ proteome.faa Aip_gff.gff Adu_genome.fna Adu_ proteome.faa Adu_gff.gff

$./PangenePro.sh
``` 

## Contributing
We welcome contributions! Please fork the repository and submit pull requests for any enhancements or bug fixes.

## Contact us
- Kinza Fatima; kfati002@ucr.edu
- Haifei Hu; huhaifei@gdaas.cn (E-mail can be in Chinese)
- Muhammad Tahir ul Qamar; tahirulqamar@gcuf.edu.pk
  
We are dedicated to supporting our community and enhancing the project with your valuable feedback.

