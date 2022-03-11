#!/usr/bin/env python
# -*- coding: utf-8 -*-
### split 10x single-cell bam file by barcode (CB tag)

import sys
import os
import logging
import pysam



FORMAT = '%(levelname)s %(asctime)-15s %(name)-20s %(message)s'
logging.basicConfig(level=logging.INFO, format=FORMAT)
logger = logging.getLogger(__name__)


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

if __name__ == "__main__":
    input_sorted_bam, barcode_list, thread, output_dir = sys.argv[1:]

    bc_set = read_barcode(barcode_list)
    split_bam(input_sorted_bam, bc_set, output_dir, thread)