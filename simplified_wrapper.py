#!/usr/bin/env python

import setlog
from argparse import ArgumentParser
import os
import subprocess

import concatmsa
import gbk2sqldb
import getmlstgenes

global log
log=setlog.init(toconsole=True)

def run_minimlst(indir, cpu):
    """Runs a reduced version of autoMLST's de novo workflow.
    Arguments:
        indir: Directory of genbank files from which to create the tree.
        threads: Number of cpu cores to use
    """

    # create SQL database from genbank files
    log.info("1. Creating SQL database from genbank files...")
    out_db = os.path.join(indir, "genbank_to_sql")
    gbk2sqldb.runall(indir, out_db, mcpu=cpu)
    # extract genes from database
    log.info("2. Extracting genes from database...")
    singles_dir = os.path.join(indir, "singles")
    getmlstgenes.findsingles(out_db, outdir=singles_dir)
    # align and trim - currently hardcoded to amino acid seqs
    log.info("3. Aligning and trimming...")
    singles_found = [single for single in os.listdir(singles_dir)
                     if os.path.isfile(os.path.join(singles_dir, single)) and single.endswith("faa")]
    # create alignment directory
    align_dir = os.path.join(indir, "aligned")
    if not os.path.isdir(align_dir):
        os.mkdir(align_dir)
    # create trimmed directory
    trim_dir = os.path.join(indir, "trimmed")
    if not os.path.isdir(trim_dir):
        os.mkdir(trim_dir)
    for single_file in singles_found:
        single_path = os.path.join(singles_dir, single_file)
        aligned_path = os.path.join(align_dir, single_file)
        # align
        log.debug(" - 3.1 Aligning {}".format(single_file))
        mafft_cmd = "mafft --quiet --thread {} {}".format(cpu, single_path)
        mafft_cmd = mafft_cmd.split()
        with open(aligned_path, "w") as outfile:
            subprocess.call(mafft_cmd, stdout=outfile)
        # trim
        log.debug(" - 3.2 Trimming {}".format(single_file))
        trimmed_path = os.path.join(trim_dir, single_file)
        trimal_cmd = "trimal -in {} -out {} -automated1".format(aligned_path, trimmed_path)
        trimal_cmd = trimal_cmd.split()
        subprocess.call(trimal_cmd)
    # concatenate
    log.info("4. Concatenating files...")
    trimmed_files = [os.path.join(trim_dir, trimmed_file) for trimmed_file in os.listdir(trim_dir)]
    concatmsa.concatmsa(trimmed_files, os.path.join(indir, "supermatrix.fa"), os.path.join(indir, "raxmlpart.txt"),
                        " ")  # split on space
    # build tree
    log.info("5. Building tree...")
    iq_tree_cmd = "iqtree -s {} -q {} -bb 1000 -nt {}".format(os.path.join(indir, "supermatrix.fa"),
                                                       os.path.join(indir, "raxmlpart.txt"),
                                                       cpu)
    iq_tree_cmd = iq_tree_cmd.split()
    subprocess.call(iq_tree_cmd)
    log.info("Done!")

if __name__ == "__main__":
    parser = ArgumentParser(description="Simplified workflow for generating a MLSA-based tree from a directory of gbks.")
    parser.add_argument("indir", help="Directory containing gbk files")
    parser.add_argument("threads", help="Number of cpu cores to use")
    args = parser.parse_args()
    input_dir = args.indir
    run_minimlst(input_dir, int(args.threads))