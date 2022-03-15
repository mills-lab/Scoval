import sys
import os
import pandas as pd
from random import choices


def is_overlap(perm_row, cnv_df):
    sub_df = cnv_df.query("barcode==\'"+perm_row["barcode"]+"\' & chrom==\'"+perm_row["chrom"]+"\'")
    if sub_df.empty:
        return False
    else:
        for index, row in sub_df.iterrows():
            if perm_row["start"] < row.start:
                if perm_row["end"] > row.start:
                    return True
            else:
                if perm_row["start"] < row.end:
                    return True
    return False


if __name__ == "__main__":
    cnv_callset_file, perm_dir, goodCell_txt, out_bed = sys.argv[1:]
    cnv_df = pd.read_csv(cnv_callset_file, sep="\t", header=0, names=["chrom", "start", "end", "bam", "copy_number", "barcode", "size"])

    good_cell_list = []
    with open(goodCell_txt) as fin:
        for line in fin:
            good_cell_list.append(line.strip())

    out_df = pd.DataFrame(columns=["chrom", "start", "end", "copy_number", "size", "barcode"])
    for i in range(10):
        print("perm"+str(i)+"....")
        input_bed = os.path.join(perm_dir, "perm"+str(i)+".bed")

        perm_df = pd.read_csv(input_bed, sep="\t", header=0, names=["chrom", "start", "end", "bam", "copy_number", "barcode", "size"])
        input_df = perm_df[["chrom", "start", "end", "copy_number", "size"]].copy()

        input_df["barcode"] = choices(good_cell_list, k=input_df.shape[0])

        del_index = []
        for index, row in input_df.iterrows():
            if is_overlap(row, cnv_df):
                del_index.append(index)

        new_df = input_df.drop(index = del_index)
        out_df = pd.concat([out_df, new_df], ignore_index=True)
    
    out_df.to_csv(out_bed, index=False, header=None, sep="\t")
    
