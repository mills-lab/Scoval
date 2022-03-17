import os, sys
import re
import pandas as pd
import numpy as np
from collections import Counter,defaultdict
import logging
import argparse

FORMAT = '%(levelname)s %(asctime)-15s %(name)-20s %(message)s'
logging.basicConfig(level=logging.DEBUG, format=FORMAT)
logger = logging.getLogger("window")


def parse_args(args):
    description = "count informative reads and make the 100 het-SNP window."
    parser = argparse.ArgumentParser(description=description)
    parser.add_argument("-s","--snp", required=True, type=str, help="pickle file of phased heterozygous SNPs information from the previous step")
    parser.add_argument("-m","--mpileup", required=True, type=str, help="mpileup results from the previous step")
    parser.add_argument("-w", "--window", type=int, required=False, default=100, help="window size by het-SNP")
    parser.add_argument("-p", "--step", type=int, required=False, default=100, help="step size by het-SNP")
    parser.add_argument("--out", required=True, type=str, help="output pickle file path")
    args = parser.parse_args(args)

    if not os.path.isfile(args.snp):
        raise ValueError("SNP file does not exist!")
    if not os.path.isfile(args.mpileup):
        raise ValueError("mpileup file does not exists!")

    return {
        'snp' : args.snp,
        'mpileup' : args.mpileup,
        'window' : args.window,
        'step' : args.step,
        'out' : os.path.abspath(args.out)
    }


def process_count(m):
    m1 = m.replace('$', '')
    m2 = m1.replace('*', '')
    m3 = re.sub(r'\^.', '', m2)

    all_match = re.findall(r'([\+-]\d+)', m3)

    if all_match:
        for t in all_match:
            n_repeat = t[1:]
            if t[0] == '+':
                reg = "\\+"+ n_repeat + "\\w{" + n_repeat + "}"
            elif t[0] == '-':
                reg = "\\-"+ n_repeat + "\\w{" + n_repeat + "}"
            
            m3 = re.sub(r'{}'.format(reg), '', m3)
            
    return m3.upper()


def allele_to_count(allele, ref, alt):
    c = Counter(allele)
    return str(c[ref])+":"+str(c[alt])


def convert(het_df, mpileup_df):
    mpileup_df["new_allele"] = mpileup_df["allele"].apply(process_count)
    mpileup_df = mpileup_df[["chrom", "pos", "new_allele"]].copy()
    mpileup_df["chrom"] = mpileup_df["chrom"].astype(str)
    het_df["chrom"] = het_df["chrom"].astype(str)
    out_df  = het_df.merge(mpileup_df, how="left", on=["chrom", "pos"])
    out_df["new_allele"] = out_df["new_allele"].replace(np.nan, '', regex=True)
    
    out_df["allele_ratio"] = out_df.apply(lambda x: allele_to_count(x["new_allele"], x["ref"], x["alt"]), axis=1)

    return out_df


def get_h0(ratio):
    a, b = map(int, ratio.split(":"))
    return a


def get_h1(ratio):
    a, b = map(int, ratio.split(":"))
    return b


def get_info(df, pseudocount=0.1):
    h0_sum = df.iloc[:,7].apply(get_h0).sum()
    h1_sum = df.iloc[:,7].apply(get_h1).sum()
    return str(h0_sum)+":"+str(h1_sum), h0_sum+h1_sum, np.log2((h0_sum+pseudocount)/(h1_sum+pseudocount))


def make_window(count_df, window_size, step_size):
    out_df = pd.DataFrame(columns=["chrom", "start", "end", "ps", "allele_count", "sum_count", "log2ratio"])

    ps_list = count_df.ps.unique()

    for ps in ps_list:
        ps_df = count_df[count_df["ps"]==ps].copy()
        if ps_df.shape[0] <= window_size:
            allele_count, sum_count, log2ratio = get_info(ps_df)
            tmp_series = pd.Series({"chrom": ps_df.iloc[0,0], "start": ps_df.iloc[0,1], "end": ps_df.iloc[-1, 1], "ps":ps, "allele_count": allele_count, "sum_count": sum_count, "log2ratio": log2ratio})
            out_df = out_df.append(tmp_series, ignore_index = True)
        else:
            start_index = list(range(0, ps_df.shape[0], window_size))
            end_index = list(range(step_size, ps_df.shape[0], window_size)) + [ps_df.shape[0]]

            for i in range(len(start_index)):
                sub_df = ps_df.iloc[start_index[i]:end_index[i]].copy()
                sub_chrom = sub_df.iloc[0,0]
                sub_start = sub_df.iloc[0,1]
                sub_end = sub_df.iloc[-1,1]
                allele_count, sum_count, log2ratio = get_info(sub_df)
                tmp_series = pd.Series({"chrom": sub_chrom, "start": sub_start, "end": sub_end, "ps": ps, "allele_count": allele_count, "sum_count": sum_count, "log2ratio": log2ratio})
                out_df = out_df.append(tmp_series, ignore_index = True)
    return out_df


def main(args=None):
    args = parse_args(args)
    logger.info('\n'.join(['Arguments:'] + [f'{a} : {args[a]}' for a in args]))
    
    logger.info("count informative reads and make the 100 het-SNP window")
    logger.info("reading files...")
    het_df = pd.read_pickle(args["snp"])
    mpileup_df = pd.read_csv(args["mpileup"], sep="\t", header = None, names=["chrom","pos","tmp1","tmp2","allele","tmp3"])

    count_df = convert(het_df, mpileup_df)
    out_df = make_window(count_df, args["window_size"], args["step_size"])
    out_df.to_pickle(args["out"])


if __name__ == '__main__':
    main()
