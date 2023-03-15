#!/usr/bin/env python
# -*- coding: utf-8 -*-

import sys
import os
import numpy as np
import pandas as pd
import logging
import argparse


FORMAT = '%(levelname)s %(asctime)-15s %(name)-20s %(message)s'
logging.basicConfig(level=logging.INFO, format=FORMAT)
logger = logging.getLogger(__name__)


def replace2na(ratio_df, count_df, cutoff):
    new_ratio_df = ratio_df.copy()
    cell_list = new_ratio_df.columns[4:]

    for cell in cell_list:
        new_ratio_df.loc[count_df.loc[:, cell]<cutoff, cell] = np.nan
    return new_ratio_df


def extract_info(window_info_dir, barcode_list, column_name, column_index, out_dir):
    first_path = os.path.join(window_info_dir, barcode_list[0]+".window.info.pkl")
    first_df = pd.read_pickle(first_path)

    all_df = first_df.iloc[:, np.r_[0:4,column_index]].copy()
    all_df.rename({column_name: barcode_list[0]}, axis=1, inplace=True)

    concat_list = []
    n_cell = 0

    for bc in barcode_list[1:]:
        pkl_path = os.path.join(window_info_dir, bc+".window.info.pkl")
        df = pd.read_pickle(pkl_path)
        
        a_series = df.iloc[:,column_index]
        a_series.rename(bc, inplace=True)

        concat_list.append(a_series)
        n_cell += 1

        if n_cell == 100:
            concat_list.insert(0, all_df)
            all_df = pd.concat(concat_list, axis=1)
            concat_list = []
            n_cell = 0
            logger.info("finish 100")

    concat_list.insert(0, all_df)
    all_df = pd.concat(concat_list, axis=1)

    for chrom in map(str, range(1,23)):
        out_pkl = os.path.join(out_dir, chrom+"."+column_name+".100.100.pkl")
        chrom_df = all_df[all_df["chrom"]==chrom].copy()
        chrom_df.reset_index(drop=True).to_pickle(out_pkl)


def parse_args(args):
    description = "merge the window information from all the cells, split the tables by chromosomes, and mask windows that has small number informative reads with NA."
    parser = argparse.ArgumentParser(description=description)
    parser.add_argument("-i","--input", required=True, type=str, help="input window information directory containing pickle files for each cell")
    parser.add_argument("-b","--barcode", required=True, type=str, help="barcode list file, each line is a unique barcode")
    parser.add_argument("-c", "--cutoff", type=int, required=False, default=3,  help="lowest cutoff of summation of two allele counts to replace log2 ratio to NA")
    parser.add_argument("-a", "--out_allele_count", type=str, required=True,  help="output direcotry for two allele counts")
    parser.add_argument("-s", "--out_sum_count", type=str, required=True,  help="output direcotry for summation of two allele counts")
    parser.add_argument("-l", "--out_log2ratio", type=str, required=True,  help="output direcotry for log2 ratio of two allele counts")
    parser.add_argument("-o", "--out_abslog2ratio", type=str, required=True,  help="output direcotry for absolute log2 ratio of two allele counts")
    
    args = parser.parse_args(args)

    if not os.path.isfile(args.barcode):
        raise ValueError("barcode list file does not exist!")
    if not os.path.isdir(args.input):
        raise ValueError(f"input directory does not exists: {args.input}")
    if not os.path.isdir(args.out_allele_count):
        raise ValueError(f"output directory for allele count does not exists: {args.out_allele_count}")
    if not os.path.isdir(args.out_sum_count):
        raise ValueError(f"output directory for sum count does not exists: {args.out_sum_count}")
    if not os.path.isdir(args.out_log2ratio):
        raise ValueError(f"output directory for log2 ratio does not exists: {args.out_log2ratio}")
    if not os.path.isdir(args.out_abslog2ratio):
        raise ValueError(f"output directory for absolute log2 ratio does not exists: {args.out_abslog2ratio}")


    return {
        'window_info_dir' : args.input,
        'bc_file' : args.barcode,
        'cutoff' : args.cutoff,
        'allele_count_dir' : args.out_allele_count,
        'sum_count_dir' : args.out_sum_count,
        'log2ratio_dir' : args.out_log2ratio,
        'abslog2ratio_dir': args.out_abslog2ratio
    }


def main(args=None):
    args = parse_args(args)
    logger.info('\n'.join(['Arguments:'] + [f'{a} : {args[a]}' for a in args]))
    
    logger.info("merge the window information from all the cells, split the tables by chromosomes, and mask windows that has small number informative reads with NA.")

    barcode_list = []
    with open(args["bc_file"]) as fin:
        for line in fin:
            barcode = line.strip()
            barcode_list.append(barcode)

    logger.info("allele count...")
    extract_info(args["window_info_dir"], barcode_list, "allele_count", 4, args["allele_count_dir"])

    logger.info("sum count...")
    extract_info(args["window_info_dir"], barcode_list, "sum_count", 5, args["sum_count_dir"])

    logger.info("log2 ratio...")
    extract_info(args["window_info_dir"], barcode_list, "log2ratio", 6, args["log2ratio_dir"])

    logger.info("absolute log2 ratio and filtering...")
    chrom_list = list(map(str, range(1,23)))
    for chrom in chrom_list:
        ratio_pkl = os.path.join(args["log2ratio_dir"], chrom+".log2ratio.100.100.pkl")
        countSum_pkl = os.path.join(args["sum_count_dir"], chrom+".sum_count.100.100.pkl")
        countSum_df = pd.read_pickle(countSum_pkl)
        ratio_df = pd.read_pickle(ratio_pkl)

        abs_df = pd.concat([ratio_df.iloc[:,:4], np.abs(ratio_df.iloc[:,4:])], axis=1)

        abs_pkl = os.path.join(args["abslog2ratio_dir"], chrom+".abslog2ratio.100.100.pkl")
        abs_df.to_pickle(abs_pkl)

        filtered_ratio_df = replace2na(ratio_df, countSum_df, cutoff=args["cutoff"])
        filtered_ratio_df.to_pickle(os.path.join(args["log2ratio_dir"], chrom+"_log2ratio.filtered.100.100.pkl"))

        filtered_absratio_df = replace2na(abs_df, countSum_df, cutoff=args["cutoff"])
        filtered_absratio_df.to_pickle(os.path.join(args["abslog2ratio_dir"], chrom+"_abslog2ratio.filtered.100.100.pkl"))



if __name__ == "__main__":
    main()