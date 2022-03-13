#!/usr/bin/env python
# -*- coding: utf-8 -*-
### split 10x single-cell bam file by barcode (CB tag)

import sys
import os
import logging
import pysam
import argparse


FORMAT = '%(levelname)s %(asctime)-15s %(name)-20s %(message)s'
logging.basicConfig(level=logging.INFO, format=FORMAT)
logger = logging.getLogger(__name__)


def parse_args(args):
    description = "Split CB tag sorted bam file into single-cell bams."
    parser = argparse.ArgumentParser(description=description)
    parser.add_argument("-i","--input", required=True, type=str, help="CB tag sorted bam file which includes all the cells")
    parser.add_argument("-b","--barcode", required=True, type=str, help="barcode list file, each line is a unique barcode")
    parser.add_argument("-t", "--threads", type=int, required=False, default=1, help="number of threads for multi-threading")
    parser.add_argument("-o", "--out", type=str, required=False, default='./', help="Running directory where to write the single-cell bams (default: current directory)")
    args = parser.parse_args(args)

    if not os.path.isfile(args.input):
        raise ValueError("input bam file does not exist!")
    if not os.path.isfile(args.barcode):
        raise ValueError("barcode list file does not exist!")
    if not os.path.isdir(args.outdir):
        raise ValueError(f"Running directory does not exists: {args.outdir}")

    return {
        'bam' : args.input,
        'bc_file' : args.barcode,
        'n_threads' : args.threads,
        'rundir' : os.path.abspath(args.outdir)
    }


# each line is a barcode in the file
def read_barcode(barcode_list):
    barcodes = set()
    with open(barcode_list) as fin:
        for line in fin:
            barcodes.add(line.strip())


def split_bam(bc_sorted_bam, barcode_set, output_dir, thread):
    samfile = pysam.AlignmentFile(bc_sorted_bam, "rb")

    cur_barcode = ""
    i = 0
    barcode_count = 0
    while True:
        try:
            record = next(samfile)
            if record.has_tag("CB"):
                barcode = record.get_tag("CB")
                if barcode in barcode_set:
                    if barcode != cur_barcode:
                        if i > 0:
                            fout.close()
                            out_sort_bam = os.path.join(output_dir, cur_barcode+"_sorted.bam")
                            out_sort_bam_index = os.path.join(output_dir, cur_barcode+"_sorted.bam.bai")
                            pysam.sort("-o", out_sort_bam, "-@", str(thread), out_bam)
                            pysam.index("-@", str(thread), out_sort_bam, out_sort_bam_index)
                        out_bam = os.path.join(output_dir, barcode+".bam")
                        fout = pysam.AlignmentFile(out_bam, "wb", template=samfile)
                        fout.write(record)
                        i += 1
                        barcode_count += 1
                        cur_barcode = barcode
                    else:
                        fout.write(record)
                        i += 1
            
                    if i > 0 and i%1000000==0:
                        logger.info(str(i)+" reads..."+str(barcode_count)+" barcodes ...")

        except StopIteration:
            fout.close()
            out_sort_bam = os.path.join(output_dir, cur_barcode+"_sorted.bam")
            out_sort_bam_index = os.path.join(output_dir, cur_barcode+"_sorted.bam.bai")
            pysam.sort("-o", out_sort_bam, "-@", str(thread), out_bam)
            pysam.index("-@", str(thread), out_sort_bam, out_sort_bam_index)
            break


def main(args=None):
    args = parse_args(args)
    logger.info('\n'.join(['Arguments:'] + [f'{a} : {args[a]}' for a in args]))
    
    logger.info("Split CB tag sorted bam file into single-cell bams.")

    bc_set = read_barcode(args.bc_file)

    logger.info("begin to split...")
    split_bam(args.bam, bc_set, args.rundir, args.n_threads)


if __name__ == "__main__":
    main()