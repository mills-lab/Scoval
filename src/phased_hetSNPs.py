import os, sys
import argparse
import subprocess as sp
import pandas as pd
from shutil import which
import logging

FORMAT = '%(levelname)s %(asctime)-15s %(name)-20s %(message)s'
logging.basicConfig(level=logging.DEBUG, format=FORMAT)
logger = logging.getLogger("hetSNP")


def parse_args(args):
    description = "Extract phased heterozygous SNPs from VCF file."
    parser = argparse.ArgumentParser(description=description)
    parser.add_argument("-v","--vcf", required=True, type=str, help="VCF file of phased germline SNPs")
    parser.add_argument("-b","--bcftools", required=False, default=None, type=str, help="Path to the directory to \"bcftools\" executable, required in default mode (default: bcftools is directly called as it is in user $PATH)")
    parser.add_argument("-c", "--chromosomes", type=str, required=False, default=','.join([str(i) for i in range(1, 23)]), help="comma-separeted list of,omosomes between apices (default: \"1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22\")")
    parser.add_argument("--outdir", required=False, default='./', type=str, help="Running directory where to write the phased heterozygous SNPs (default: current directory)")
    args = parser.parse_args(args)

    if not os.path.isfile(args.vcf):
        raise ValueError("VCF file does not exist!")
    if not os.path.isdir(args.outdir):
        raise ValueError(f"Running directory does not exists: {args.outdir}")
    
    bcftools = args.bcftools
    if not bcftools:
        bcftools = "bcftools"
    if which(bcftools) is None:
        raise ValueError("bcftools has not been found or is not executable!")

    return {
        'vcf' : args.vcf,
        'bcftools' : bcftools,
        'chrs' : args.chromosomes,
        'rundir' : os.path.abspath(args.outdir)
    }


def get_phased_het_snp(vcf_file, bcftools, chrs, out_tsv):
    logger.info("Running bcftools")
    fmt = "\'%CHROM\t%POS\t%REF\t%ALT[\t%GT\t%PS]\n\'"
    cmd_bcftools = f"{bcftools} view -m2 -M2 -v snps -p -f 'PASS' -g het -r {chrs} {vcf_file} | {bcftools} query -u -f {fmt}"

    with open(out_tsv, "w") as fout:
        proc = sp.run(cmd_bcftools, shell=True, stdout=fout, stderr=sp.PIPE)
    if proc.returncode == 0:
        logger.info("phased het SNP done")
    else:
        logger.error("error in phased het SNP")
        logger.error(proc.stderr)
        sys.exit(101)


def drop_jump_ps(df):
    drop_ps_set = set()
    known_ps_set = set()

    last_ps = 0
    for index, row in df.iterrows():
        if row["ps"] != last_ps:
            if row["ps"] in known_ps_set:
                drop_ps_set.add(row["ps"])          
            else:
                known_ps_set.add(row["ps"])
        last_ps = row["ps"]
    
    df.drop(df[df["ps"].isin(drop_ps_set)].index, inplace=True)
    df.reset_index(drop=True, inplace=True)

    return df


def main(args=None):
    args = parse_args(args)
    logger.info('\n'.join(['Arguments:'] + [f'{a} : {args[a]}' for a in args]))
    
    logger.info("Processing VCF file and extracting phased het SNPs")
    
    raw_path = os.path.join(args["rundir"], "raw_phased_het_snp.tsv")
    get_phased_het_snp(args["vcf"], args["bcftools"], args["chrs"], raw_path)

    logger.info("dropping duplicate SNPs")
    raw_df = pd.read_csv(raw_path, sep="\t", header=None, names=["chrom", "pos", "ref", "alt", "gt", "ps"])
    raw_df.drop_duplicates(subset=["chrom", "pos"], keep=False)

    logger.info("filtering jump PS")
    filtered_path = os.path.join(args["rundir"], "phased_het_snp.pkl")
    pos_tsv_path = os.path.join(args["rundir"], "phased_het_snp_pos.tsv")
    if "." in raw_df["ps"]:
        raw_df.to_pickle(filtered_path, protocol=4)
        raw_df[["chrom", "pos"]].to_csv(pos_tsv_path, index=False, header=False, sep="\t")
    else:
        drop_df = drop_jump_ps(raw_df)
        drop_df.to_pickle(filtered_path, protocol=4)
        drop_df[["chrom", "pos"]].to_csv(pos_tsv_path, index=False, header=False, sep="\t")


if __name__ == '__main__':
    main()