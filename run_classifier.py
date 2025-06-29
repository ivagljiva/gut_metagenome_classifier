#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import os
import sys
import pickle
import argparse
import pandas as pd

from sklearn.linear_model import LogisticRegression

__description__ = "A script to load and run the human gut metagenome classifier described in https://doi.org/10.7554/eLife.89862"
__author__ = "Iva Veseli"
__maintainer__ = "Iva Veseli"
__copyright__ = "Copyleft 2025, Iva Veseli"
__license__ = "GPL 3.0"

IBD_ENRICHED_MODULES_FILE="TRAINING_DATA/IBD_ENRICHED_MODULES.txt" # the modules used to train the classifier

def main(args):
    """Sanity checks the arguments, loads the data and the classifier model, runs the classifier and stores predictions."""
    
    input_sanity_checks(args)
    ppcn = load_data(args)
    predictions = classify(ppcn)
    print(f"Distribution of classes predicted by classifier:\n{predictions.groupby(['class_string']).size()}")
    predictions.to_csv(args.output_file, sep="\t", index_label='sample')
    print(f"Predictions saved to {args.output_file}")

    if args.also_output_ppcn:
        if not args.ppcn_output_file:
            args.ppcn_output_file = "PPCN_matrix.txt"
        ppcn.to_csv(args.ppcn_output_file, sep="\t", index_label='sample')
        print(f"PPCN matrix saved to {args.ppcn_output_file}")


def input_sanity_checks(args):
    """Ensures the provided arguments are sensible."""

    if (not args.ppcn_table and not args.copy_numbers) or (args.ppcn_table and args.copy_numbers):
        raise RuntimeError("Input Error: Please choose one of the input options: either provide --ppcn-table, or provide --copy-numbers "
                           "and --populations.")
    if (args.copy_numbers and not args.populations):
        raise RuntimeError("Input Error: Providing module copy numbers with --copy-numbers requires that you also provide population sizes "
                           "with the --populations flag.")

    if args.ppcn_table and args.also_output_ppcn:
        raise RuntimeError("Input Error: the --also-output-ppcn flag doesn't work with --ppcn-table input.")
    if args.ppcn_output_file and not args.also_output_ppcn:
        raise RuntimeError("Input Error: if you specify the --ppcn-output-file parameter, you should probably also specify --also-output-ppcn.")
    
def load_data(args):
    """Reads the input file(s), computes PPCN if necessary, and returns the PPCN matrix."""

    ppcn_matrix = None
    if args.ppcn_table:
        ppcn_matrix = pd.read_csv(args.ppcn_table, sep="\t", index_col=0)
        
    elif args.copy_numbers:
        copies = pd.read_csv(args.copy_numbers, sep="\t", index_col=0)
        pops = pd.read_csv(args.populations, sep="\t", index_col=0)
        if copies.index.name != 'module':
            raise RuntimeError("Input Error: the first column of your copy number matrix should be named 'module'.")
        if 'num_populations' not in pops.columns:
            raise RuntimeError("Input Error: your populations table should include a column called 'num_populations'.")
        sample_set_cn = set(copies.columns.to_list())
        sample_set_pop = set(pops.index.to_list())
        if not sample_set_cn.issubset(sample_set_pop):
            difference = sample_set_cn.difference(sample_set_pop)
            raise RuntimeError(f"Input Error: some of the samples in your copy number matrix are not present in the populations "
                               f"table. Here are the affected samples: {', '.join(list(difference))}")
        ppcn_matrix = copies / pops['num_populations']
        ppcn_matrix = ppcn_matrix.T

    # load set of IBD-enriched modules and use it to subset the data matrix
    ibd_enriched_modules = pd.read_csv(IBD_ENRICHED_MODULES_FILE, sep='\t')['module'].to_list()
    mods_not_in_table = set(ibd_enriched_modules).difference(set(ppcn_matrix.columns.to_list()))
    if mods_not_in_table:
        raise RuntimeError(f"Input Error: some of the IBD-enriched modules used to train the classifier were not found in "
                           f"the input data. Here they are: {', '.join(list(mods_not_in_table))}")
    ppcn_matrix_ibd_enriched = ppcn_matrix[ibd_enriched_modules]
    assert ppcn_matrix_ibd_enriched.shape[1] == len(ibd_enriched_modules)

    return ppcn_matrix_ibd_enriched

def classify(features):
    """Loads classifier and returns predictions as a data frame."""

    final_model = pickle.load(open("classifier.pickle", 'rb'))
    pred = final_model.predict(features)
    predictions_df = pd.DataFrame(pred, index=features.index, columns=["class"])
    predictions_df.loc[predictions_df['class'] == 1, 'class_string'] = 'IBD'
    predictions_df.loc[predictions_df['class'] == 0, 'class_string'] = 'HEALTHY'

    return predictions_df

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__description__)
    groupA = parser.add_argument_group("INPUT OPTION 1", "Start from a matrix of per-population copy numbers (PPCNs) of KEGG modules for "
                                                         "each metagenome sample.")
    groupA.add_argument("--ppcn-table", required=False, metavar="FILE", help="A tab-delimited matrix where rows are sample names, columns are "
                                                         "modules, and values are per-population copy numbers. The header for the module column "
                                                         "should be 'module'.")

    groupB = parser.add_argument_group("INPUT OPTION 2", "Start from a tab-delimited matrix of KEGG module copy numbers and a table of "
                                                         "population sizes for each metagenome sample. The PPCN calculation will be done "
                                                         "before classification.")
    groupB.add_argument("--copy-numbers", required=False, metavar="FILE", help="A matrix where rows are modules, columns are sample "
                                                         "names, and values are unnormalized module copy numbers.")
    groupB.add_argument("--populations", required=False, metavar="FILE", help="A tab-delimited table where the first column contains "
                                                         "sample names, and there is a column called 'num_populations' that contains "
                                                         "the number of microbial populations estimated to be present in each sample.")                                   

    groupC = parser.add_argument_group("OUTPUT", "What you want to get back from this program.")
    groupC.add_argument("-o", "--output-file", required=False, default="predictions.txt", help="The name of the output file in "
                                                          "which to store the classifier predictions.")
    groupC.add_argument("--also-output-ppcn", required=False, action='store_true', help="If you use input option 2, this "
                                                          "flag will ensure that the computed PPCN values are stored in an output file.")
    groupC.add_argument("--ppcn-output-file", required=False, metavar="FILE", help="The name for the computed PPCN output matrix, "
                                                           "if requested.")
    args = parser.parse_args()

    try:
        main(args)
    except Exception as e:
        print(e)
        sys.exit(-1)
