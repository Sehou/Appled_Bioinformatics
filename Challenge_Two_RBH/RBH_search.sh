#!/bin/bash

############################################################################
# Author: Kouiho, Sehou Romaric
# Last update: February 2021
# Email: romarickouiho@gmail.com
############################################################################
#									      
# Pipeline script for finding Reciprocal Best Hits in two genomes
#                            using BLAST.
#								
############################################################################

set -e
# set -u
set -o pipefail

# Set default values for optional arguments
evalue=0.000001
wsize=3
numalign=100000
outdir=$PWD

# Usage and help message
usage()
{
cat << EOF

DESCRIPTION:
    This script runs a pipeline for finding Reciprocal Best Hits between
    two genomes or a set of protein sequences using the BLASTP flavour 
            of BLAST. Hence, it only works on protein data.
    
    Although BLAT, LAST and UBLAST are hundredth to a 25th faster than
    BLAST. We choose BLAST to avoid the reduction in the number of 
    homologs and RBHs found by the faster algorithms even if UBLAST
           produced between 0.6 and 0.8 of the RBHs as BLAST.
    
    This script takes the path to a query and target genomes, a path to
    the python parser for the BLAST XML output and optional options to 
                         refine the BLAST search.
    It outputs a folder named Blast_RBHs containing the output into
    the specified output directory or the curent working directory.
    

USAGE: ${0} <options>


OPTIONS:

  Mandatory options:
  -h      Show this help message.
  -r      Absolute path to the target genome (relative path 
                                        will raise an error).
  -s      Absolute path to the subject genome (relative path 
                                        will raise an error).
  -p      Absolute path to the python parser script of BLAST
               XML output to FASTA.
  
  Optional options:
  -e      Threshold for e-value in the BLAST search. Default is 0.00001.
  -k      Word size for the BLAST search. Default is 3.
  -n      Number of alignments for the BLAST search. Default is 100000.
  -o      Output directory. It should be an existing directory otherwise,
          the script will raise an error. Default is the current directory.
 
EOF
}

# Get options from command line
while getopts ":r:s:p:e:k:n:o:h" OPTION;
do
  case ${OPTION} in
    h)  usage; exit 1;;
    r)  target=$OPTARG;;
    s)  query=$OPTARG;;
    p)  parser=$OPTARG;;
    e)  evalue=$OPTARG;;
    k)  wsize=$OPTARG;;
    n)  numalign=$OPTARG;;
    o)  outdir=$OPTARG;;
    \?)  echo "ERROR: unknown option: -$OPTARG" >&2; 
        usage;
        exit 1;;
    :)  echo "ERROR: missing option argument for: -$OPTARG" >&2;
        exit 1;;
    *)  echo "ERROR: unimplemented option: -$OPTARG" >&2;
        exit 1;;
  esac
done

# Check if an options was specified
if ((OPTIND == 1))
then
   echo "ERROR: No options specified"
   usage
   exit 1
fi

shift $((OPTIND -1))


# Check if all mandatory options were passed
if [[ -z $target ]] || [[ -z $query ]]  || [[ -z $parser ]]
then
  printf "\n==================================\n 
  ERROR: missing options\n==================================\n\n"
  usage
  exit 1
fi

# Check if the eventually provided output directory exists
if [ ! -d $outdir ]
then
  printf "\nERROR: provided output directory ${outdir} does not exists.
  Please provided a valid existing directory or remove the -o option.\n"
  usage
  exit 1
fi

# Check if required blast programs are installed.
if ! [ -x "$(command -v blastp)" ]; then
  printf '\nERROR: blast programs to be used are not installed.
   Please install using\n "sudo apt install ncbi-blast+" \n' >&2
  exit 1
fi

# Make output directory and move into it
cd $outdir && mkdir Blast_RBHs && cd Blast_RBHs
outdir=${outdir}/Blast_RBHs
# Make database from the target genome for the forward search
target_database="TARGET_DATABASE"
echo "MAKING SEARCH DATABASE OF THE TARGET GENOME..." 1>&2
[ -d ${outdir}/db ] || mkdir ${outdir}/db
makeblastdb -in ${target} -dbtype prot -out ${outdir}/db/${target_database}

# Perform a forward search of query genome in target genome
echo "RUNNING FORWARD SEARCH IN THE TARGET GENOME..." 1>&2
blastp -query ${query} -db ${outdir}/db/${target_database} -out ${outdir}/forward_blast.xml -outfmt 5 -evalue ${evalue} -word_size ${wsize} -num_alignments ${numalign}
 
# Parse forward search output to fasta format
echo "PARSING FORWARD SEARCH OUTPUT..." 1>&2
python ${parser} -hit ${outdir}/forward_blast.xml -genome ${target} -out ${outdir}/forward_hits

# Make database from the query genome for the reverse search
query_database="QUERY_DATABASE"
echo "MAKING SEARCH DATABASE OF THE QUERY GENOME..." 1>&2
makeblastdb -in ${query} -dbtype prot -out ${outdir}/db/${query_database}

# Perform the reverse search of hits in query genome
# num_alignments is set to 1 as we only want the best hit
echo "RUNNING REVERSE SEARCH IN THE QUERY GENOME..." 1>&2

[ -d ${outdir}/reverse_hits ] || mkdir ${outdir}/reverse_hits
find ${outdir}/forward_hits -name "*_hits.faa" | \
xargs basename -s "_hits.faa" | \
xargs -n1 -I{} blastp -query ${outdir}/forward_hits/{}_hits.faa -db ${outdir}/db/${query_database} -out ${outdir}/reverse_hits/{}_reverse_hits.out -outfmt 6 -word_size ${wsize} -num_alignments 1 -evalue ${evalue}

# Process output to retrieve RBHs
echo "PROCESSING OUTPUT AND RETIEVING RBHs..." 1>&2


for hit in ${outdir}/reverse_hits/*_hits.out
do
  # get the name of the original query
  forward_hit=$(basename $hit _reverse_hits.out)
  printf "\n>The RBHs of %s  are:" "${forward_hit}"
  printf "\nAccession_No\tpercent_identity\tbit_score\tevalue\n"
  
  # check if the accession number is the same as the one of the query in the reverse
  sort  -k3,3 -r -n -t $'\t' $hit | \
  awk -v initial_query=$forward_hit '{ if (initial_query == $2) print $1 "\t" $3 "\t" $11 "\t" $NF;}'
  
done > ${outdir}/list_RBHS.txt

# Print output to stdin
cat ${outdir}/list_RBHS.txt

printf "\nThe output text file list_RBHS.txt containg the RBHs is saved to the path %s" "${outdir}"
echo -e "\nDone successfully.\n"



