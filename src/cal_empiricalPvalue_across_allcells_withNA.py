import pandas as pd
import numpy as np
import sys
import os


def get_start_index(df, cnv_start):
    if cnv_start <= df.iloc[0].start:
        return 0
    else:
        start_index = df.index[df["start"] <= cnv_start].tolist()[-1]
        if df.iloc[start_index].end <= cnv_start:
            start_index += 1
        return start_index


def get_end_index(df, cnv_end):
    if cnv_end >= df.iloc[-1].end:
        return df.shape[0]-1
    else:
        end_index = df.index[df["end"] >= cnv_end].tolist()[0]
        if df.iloc[end_index].start >= cnv_end:
            end_index -= 1
        return end_index


def validate(cnv_call_file, ratio_pkl_dir):
    cnv_df = pd.read_csv(cnv_call_file, sep="\t")

    ratio_dict = {}
    chr_list = list(map(str, range(1,23)))
    for chrom in chr_list:
        pkl_path = os.path.join(ratio_pkl_dir, chrom+"_abslog2ratio.filtered.100.100.pkl")
        df = pd.read_pickle(pkl_path)
        df.reset_index(drop=True)
        ratio_dict["chr"+chrom] = df.copy()
    
    out_df = cnv_df.copy()
    out_df["median_log2"] = 0.0
    out_df["pvalue"] = 0.0
    for index, row in cnv_df.iterrows():
        if "bam" in row:
            barcode = row["bam"].split(".")[0]
        else:
            barcode = row["barcode"].split("-")[0]
        ratio_df = ratio_dict[row["chrom"]]
        try:
            start_index = get_start_index(ratio_df, row["start"])
        except IndexError:
            print(row)
        end_index = get_end_index(ratio_df, row["end"])

        sub_df = ratio_df.iloc[start_index:end_index+1]
        if sub_df.empty:
            cell_ratio = np.nan
            pvalue = np.nan
        else:
            bg_ratio_values = ratio_df.iloc[start_index:end_index+1, 4:].apply(np.nanmedian, 0).values
                        
            cell_ratio = np.nanmedian(ratio_df.iloc[start_index:end_index+1].loc[:, barcode])
            if np.isnan(cell_ratio):
                pvalue = np.nan
            else:
                pvalue = sum(bg_ratio_values > cell_ratio)/len(bg_ratio_values)

        out_df.at[index, "median_log2"] = cell_ratio
        out_df.at[index, "pvalue"] = pvalue
    
    return out_df


if __name__ == "__main__":
    cnv_call_file, ratio_pkl_dir, out_pkl = sys.argv[1:]
    out_df = validate(cnv_call_file, ratio_pkl_dir)
    out_df.to_pickle(out_pkl)
    # out_df.to_csv(out_csv, sep="\t", index=False)