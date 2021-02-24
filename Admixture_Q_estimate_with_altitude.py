#!/usr/bin/env python3

"""
============================================================================
Author: Kouiho, Sehou Romaric
Last updated: Feb 2021
email: romarickouiho@gmail.com
============================================================================

Script for plotting Q estimates and the altitude at which each sample was
                   collected for ONLY ONE specie.

USAGE: Admixture_Q_estimate_with_altitude.py <arguments>

============================================================================
WARNING!!! Please run "python Admixture_Q_estimate_with_altitude.py --help"
         and carefully read the information about the usage of this script.

         The script raise an error if more invalid options are provided.
         Please keep in mind that this script assume that your are passing
         the data of ONLY ONE specie and it does not check it.
============================================================================
"""

# import statements
import argparse
import logging
import os
import pandas as pd
import matplotlib.pyplot as plt
from sklearn.cluster import KMeans

logging.basicConfig(format="%(levelname)s : %(message)s", level=logging.INFO)


# function definitions
def get_arguments_from_cmd():
    """ Collect all the given arguments from command line

    :return: a JSON object of the mapping of arguments with their values.
    """

    parser = argparse.ArgumentParser(description="\
    Script for plotting Q estimates and the altitude at which each sample \
    was collected for ONLY ONE specie.\n The script raise an error if more \
    invalid options are provided. Please keep in mind that this script \
    assumes that your are passing the data of ONLY ONE specie and"
                                                 "it does not check it.")

    parser.add_argument("-m", type=str, required=True,
                        help="Path to the tab delimited metadata file")
    parser.add_argument("-q", type=str, required=True,
                        help="Path to the Q estimate file")
    parser.add_argument("-o", type=str, required=False,
                        help="Existing directory where to save the plot")

    return parser.parse_args()


def load_q_matrix_and_meta_data(meta_data_file, q_matrix_file):
    """Merge Q estimates file and metadata file into a dataframe.

    :param meta_data_file: --str, path to the metadata file
    :param q_matrix_file:--str, path to the Q estimates file

    :return: pandas dataframe of Q estimates and  metadata
    """

    q_matrix = pd.read_table(q_matrix_file, sep=" ", header=None)

    # Rename columns dynamically for the Q estimate file
    q_matrix.columns = [f"Population{i + 1}" for i in
                        range(len(q_matrix.columns))]
    meta_data = pd.read_table(meta_data_file, header=0)
    # Combine/merge the two dataframes
    combined_data = pd.concat([meta_data, q_matrix], axis=1)

    # Remove row  with NA data and warn user
    if combined_data.isnull().values.any():
        combined_data.dropna(axis=1)
        logging.warning("Some NAN values found in file and rows ignored.")

    return combined_data


def perform_clustering_to_order_sample(comb_dataframe, state=26):
    """Perform K-means clustering on Q estimates, group and order samples.

    :param comb_dataframe: pandas dataframe of Q estimates and  metadata
    :param state: --int, K-mean algorithm random state. Default is 26.

    :return: panda dataframe with the samples ordered according to
                    clusters found by the K-mean.
    """

    populations = [col for col in comb_dataframe.columns \
                   if "Population" in col]
    # Select the data of Q estimates only as a new dataframe
    k_means_data = comb_dataframe[populations]

    # Initialize the K-mean model and fit on selected data
    k_means_model = KMeans(n_clusters=len(populations), random_state=state)
    fitted_k_means = k_means_model.fit(k_means_data)

    # Retrieved samples cluster assignment by the clustering and add it
    comb_dataframe["clusters"] = fitted_k_means.labels_

    # Sort the dataframe based on the clusters assignment
    comb_dataframe.sort_values(by=["clusters"], inplace=True)

    return comb_dataframe


def plot_admixture_q_matrix(comb_dataframe, out_dir=None):
    populations = [col for col in comb_dataframe.columns \
                   if "Population" in col]

    q_estimates = comb_dataframe[populations]
    altitude = comb_dataframe[["altitude"]]

    figure, (axe1, axe2) = plt.subplots(2, sharex=True)
    figure.suptitle('Q estimates  and altitude ')

    plot_alt = altitude.plot(ax=axe1, kind='bar', legend=False, color='red')
    plot_alt.set_ylabel("Altitude", fontsize=12)

    plot_q = q_estimates.plot(stacked=True, kind='bar', ax=axe2)
    # ax.legend()
    plot_q.legend(loc='upper center', bbox_to_anchor=(0.5, 1.05),
          ncol=3, fancybox=True, shadow=True)
    plot_q.set_ylabel("Q estimates", fontsize=12)
    plot_q.set_xlabel("Samples", fontsize=12)
    plot_q.set_xticklabels(comb_dataframe.sample_id, fontsize=8)
    plot_q.set_xticklabels([t if not i % 3 else "" for i, t in enumerate(axe2.get_xticklabels())])

    if out_dir:
        plt.savefig(f"{out_dir}/Admixture_altitude.png", dpi=800, bbox_inches = 'tight')
    else:
        plt.savefig(f"Admixture_altitude.png")
    plt.show()


if __name__ == '__main__':
    # Taking arguments from command line
    args = get_arguments_from_cmd()
    # Unpack arguments into variables
    q_file = args.q
    meta = args.m

    logging.info("Checking provided inputs...")
    # Check if the provided path to Q file exists
    if not os.path.exists(q_file) :
        raise ValueError(f"File {q_file} does not exists.")

    # Check if the provided path to metadata file exists
    if not os.path.exists(meta):
        raise ValueError(f"File {meta} does not exists.")

    # Check if the provided output directory exists
    if args.o and not os.path.exists(args.o):
        raise ValueError(f"Provided output directory does not exists."
                         f"Remove -o option or specify a valid directory.")
    elif not args.o:
        output_dir = os.getcwd()
    else:
        output_dir = args.o

    logging.info("Loading and combining data...")
    # Load and combine data
    data = load_q_matrix_and_meta_data(meta, q_file)

    logging.info("Performing K-mean clustering...")
    # Cluster samples based on Q estimates and order them
    plot_data = perform_clustering_to_order_sample(data)

    logging.info(f"Plotting plot will be saved as Admixture_altitude.png"
                 f" in {output_dir}...")
    # Plot
    plot_admixture_q_matrix(data, output_dir)
