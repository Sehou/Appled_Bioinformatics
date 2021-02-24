# RBH_search

RBH_search is a bash script for retrieving Reciprocal Best Hits (RBHs) in two genomes using BLASTP. In other words, this script retrieves the proteins encoded by two genes in two genomes  find each other as the best hit in the other genome. This helps to identify homologs and orthologs while comparing two genomes, common task in comparative genomics.

The script first makes a database of the genomes and performs the forward search. The best hits of the forward search are passed and their sequence are written into FASTA files. These FASTA files are used as query in the reverse search. The RBHs are finally identified and written out into a txt file.

**WARNING:** As this script is based on BLASTP, it can only works on protein genomes. Besides, it raise errors if wrong inputs are provided.



## Installation

This script does not require any installation. However, its repository must be cloned before usage. This can be done as follow:

```bash
git clone https://github.com/Sehou/Applied_Bioinformatics.git
```



## Usage

This script can be used by running the command below in a Linux/Unix command line.

```python
./RBH_search.sh <options>
```

In case it fails , permission rights have to be added using the command below.

```bash
sudo chmod u+x ./RBH_search.sh
```



The options -r for the reference and -s for subject, -p for python BLAST XML output to FASTA are mandatory. All the available options can be see by running: 

```bash
python Admixture_Q_estimate_with_altitude.py --help
```



## Output

Once the program finish executing, the output files are saved  into a folder Blast_RBHs inside the specified output direction if any or in the current working directory. The list of the RBHs are in the file list_RBHS.txt.