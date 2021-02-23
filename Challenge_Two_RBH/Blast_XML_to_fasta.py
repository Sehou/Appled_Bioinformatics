#!/usr/bin/env python3

"""
============================================================================
Author: Kouiho, Sehou Romaric
Last update: Feb 2021
email: romarickouiho@gmail.com
============================================================================

Script for parsing BLAST XML output (format 5) to FASTA file for each hits.

USAGE: Blast_XML_to_fasta.py <arguments>

============================================================================
WARNING!!! This script does not raise any error. Errors are handle if the
            required inputs are provided on the right formats.
============================================================================
"""

# import statements
import argparse
import logging
import os
from Bio.Blast import NCBIXML

logging.basicConfig(format="%(levelname)s : %(message)s", level=logging.INFO)


# function definitions
def get_arguments_from_cmd():
    """ Collect all the given arguments from command line

    :return: a JSON object of the mapping of arguments with their values.
    """

    parser = argparse.ArgumentParser(description="\
    Script for parsing BLAST XML output (format 5) to FASTA \
    file for each hits.\n")

    parser.add_argument("-hit", type=str, required=True,
                        help="path to the hits file to be parsed")
    parser.add_argument("-genome", type=str, required=True,
                        help="path to the input genome")
    parser.add_argument("-out", type=str, required=False,
                        help="output directory")

    return parser.parse_args()


def parse_fasta_file(fasta_file):
    """ Parse  FASTA file to dictionary.

    :param fasta_file: --str, path of the fasta file to be parsed.

    :return: dictionary with seqID as keys and sequences as values
    """

    genome, sequence, seq_id = {}, "", ""

    with open(fasta_file) as records:
        for line in records:
            if not line.strip():
                continue

            line = line.strip()
            if line.startswith(">") and sequence != "":
                genome[seq_id] = sequence
                sequence = ""
            if line.startswith(">"):
                seq_id = line.strip('>').split()[0]
            else:
                sequence = sequence + line
        return genome


def parse_query_hits_in_reference(hits_file, genome_path,
                                  out_dir="forward_hits"):
    """Read the BLAST output XML file and extract the protein
    accession number and sequence from the genome FASTA file.

    :param hits_file: --str, path of the BLAST XML file to be parsed.
    :param genome_path: --str, path of the target genome.
    :param out_dir: --str, path of the output directory.
                    Optional, default is forward_hits.

    :return: This function return None. It writes files to disk.
    """

    genome = parse_fasta_file(genome_path)
    os.makedirs(out_dir, exist_ok=True)

    with open(hits_file) as records:
        hit_records = NCBIXML.parse(records)

        for hits in hit_records:
            file_name = f"{out_dir}/{hits.query.split()[0]}_hits.faa"
            logging.info(f"Writing {file_name}...")

            with open(file_name, 'w') as protein_file:
                for hit in hits.alignments:
                    protein_id = hit.title.split(" ", 2)[1]
                    try:
                        protein_seq = genome[protein_id]
                        protein_file.write(f">{str(protein_id)}\n")
                        protein_file.write(f"{str(protein_seq)}\n")
                    except KeyError:
                        logging.info(f"{protein_id} not found in genome.")
                        pass

        # retrieve and delete empty files
        empty_files = [f"{out_dir}/{f}" for f in \
                       os.listdir(out_dir) \
                       if f.endswith("_hits.faa") \
                       and os.path.getsize(f"{out_dir}/{f}") == 0]

        [os.remove(f) for f in empty_files]

        logging.info("Done parsing XML to FASTA.")


if __name__ == '__main__':
    # get arguments from command line
    arguments = get_arguments_from_cmd()

    # unpack arguments into variables
    hit_file_record = arguments.hit
    genome_file = arguments.genome
    outdir = arguments.out

    # Parse BLAST XML file
    parse_query_hits_in_reference(hit_file_record, genome_file, outdir)

