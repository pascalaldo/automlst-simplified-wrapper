#!/usr/bin/env python

from argparse import ArgumentParser
import os
import subprocess

import concatmsa
import gbk2sqldb
import getmlstgenes

def run_minimlst(indir):
    """Runs a reduced version of autoMLST's de novo workflow.
    Arguments:
        indir: Directory of genbank files from which to create the tree.
    """
    # create SQL database from genbank files
    out_db = os.path.join(indir, "genbank_to_sql")
    gbk2sqldb.runall(indir, out_db)
    # extract genes from database
    singles_dir = os.path.join(indir, "singles")
    getmlstgenes.findsingles(out_db, outdir=singles_dir)
    # align and trim - currently hardcoded to amino acid seqs
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
        mafft_cmd = "mafft --quiet --thread 1 {}".format(single_path)
        mafft_cmd = mafft_cmd.split()
        with open(aligned_path, "w") as outfile:
            subprocess.call(mafft_cmd, stdout=outfile)
        # trim
        trimmed_path = os.path.join(trim_dir, single_file)
        trimal_cmd = "trimal -in {} -out {} -automated1".format(aligned_path, trimmed_path)
        trimal_cmd = trimal_cmd.split()
        subprocess.call(trimal_cmd)
    # concatenate
    trimmed_files = [os.path.join(trim_dir, trimmed_file) for trimmed_file in os.listdir(trim_dir)]
    concatmsa.concatmsa(trimmed_files, os.path.join(indir, "supermatrix.fa"), os.path.join(indir, "raxmlpart.txt"),
                        " ")  # split on space
    # build tree
    iq_tree_cmd = "iqtree -s {} -q {} -bb 1000".format(os.path.join(indir, "supermatrix.fa"),
                                                       os.path.join(indir, "raxmlpart.txt"))
    iq_tree_cmd = iq_tree_cmd.split()
    subprocess.call(iq_tree_cmd)

if __name__ == "__main__":
    parser = ArgumentParser(description="Simplified workflow for generating a MLSA-based tree from a directory of gbks.")
    parser.add_argument("indir", help="Directory containing gbk files")
    args = parser.parse_args()
    input_dir = args.indir
    run_minimlst(input_dir)
