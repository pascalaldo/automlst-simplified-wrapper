#!/usr/bin/env python
# Copyright (C) 2015,2016 Mohammad Alanjary
# University of Tuebingen
# Interfaculty Institute of Microbiology and Infection Medicine
# Lab of Nadine Ziemert, Div. of Microbiology/Biotechnology
# Funding by the German Centre for Infection Research (DZIF)
#
# This file is part of ARTS
# ARTS is free software. you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version
#
# License: You should have received a copy of the GNU General Public License v3 with ARTS
# A copy of the GPLv3 can also be found at: <http://www.gnu.org/licenses/>.

import os, argparse, setlog, json, pickle, sqlite3 as sql

log = setlog.init(toconsole=True,level="info")

def writeallgenes(db,glist,ignore,outdir=".",tax="",allgenes=False,genelimit=10000,outgroups=None,pct=0.5,writeaa=True,rename=False):
    outdir = os.path.realpath(outdir)
    if type(glist) is not list and os.path.exists(str(glist)):
        with open(glist,"r") as fil:
            if ".json" in glist:
                temp = json.load(fil)
                glist = [x.get('acc','') for x in temp]
                glist = [x for x in glist if x]
            else:
                glist = [x.strip() for x in fil]
    if type(ignore) is not list and os.path.exists(str(ignore)):
        with open(ignore,"r") as fil:
            ignore = [x.strip() for x in fil]

    glist = glist[:genelimit]

    #Use taxonomy dictionary
    if tax and os.path.exists(tax):
        pass
    else:
        tax = os.path.join(os.path.dirname(os.path.realpath(__file__)),"gcf2names.pkl")
    log.info("Loading taxonomy definitions...")
    with open(tax,"r") as fil:
        taxonomy = pickle.load(fil)

    log.info("Writing genes: %s"%glist)
    if not os.path.exists(outdir):
        os.makedirs(outdir)

    conn = sql.connect(db)

    if allgenes:
        temp = [x[0] for x in conn.execute("Select hmmhit from HMMhits").fetchall()]
        glist = temp + [x[0] for x in conn.execute("Select hmmhit from RNAhits").fetchall()]

    for gene in glist:
        table = "HMMhits"
        if "_rRNA" in gene:
            table = "RNAhits"
        query = "SELECT orgname,seqid,description,naseq,aaseq FROM Seqs WHERE seqid in (SELECT seqid FROM %s WHERE flags>=1 and hmmhit=='%s' and hmmcov >= %.4f)"%(table,gene,pct)
        if len(ignore):
            query += " AND orgname NOT IN ('%s')"%"','".join(ignore)
        results = conn.execute(query)

        #WRITE ALL RESULTS
        with open(os.path.join(outdir,gene+".fna"),"w") as nafil, open(os.path.join(outdir,gene+".faa"),"w") as aafil:
            for row in results:
                #Translate assembly id using taxonomy table if applicable. Defualt to user title and add {}identifiers to title
                orgname = taxonomy.get(row[0],row[0])
                #mark outgroups
                if outgroups and row[0] in outgroups:
                    orgname = "OG--" + orgname

                # if row[0] in taxonomy.keys():
                    # orgname = taxonomy[row[0]]["organism_name"]
                    # if taxonomy[row[0]]["strain"] not in orgname:
                    #     orgname += "[%s]" % taxonomy[row[0]]["strain"]
                    # orgname = orgname.replace(" ","_")
                if rename:
                    nafil.write(">%s\n%s\n" % (orgname, row[3]))
                    if table == "HMMhits":
                        aafil.write(">%s\n%s\n" % (orgname, row[4]))
                else:
                    nafil.write(">%s|%s %s\n%s\n" % (orgname, row[1], row[2], row[3]))
                    if table == "HMMhits":
                        aafil.write(">%s|%s %s\n%s\n" % (orgname, row[1], row[2], row[4]))

        #Cleanup empty files
        for x in os.listdir(outdir):
            fname = os.path.join(outdir,x)
            if os.path.isfile(fname) and os.path.getsize(fname) == 0:
                os.remove(fname)

# Commandline Execution
if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="""Extract genes from sequence database""")
    parser.add_argument("input", help="Input sequence database")
    parser.add_argument("-a", "--allgenes", help="Extract every gene in database (default: FALSE)",action="store_true",default=False)
    parser.add_argument("-waa", "--writeaa", help="Also write .faa file",action="store_true",default=False)
    parser.add_argument("-l","--listfile", help="List file of genes to extract (plaintxt)",default="")
    parser.add_argument("-gl","--genelimit", help="Only print top X genes from list file (default = 10000)",type=int,default=10000)
    parser.add_argument("-i","--ignorefile", help="List file of organisms to ignore (plaintxt)",default="")
    parser.add_argument("-t","--tax", help="Use taxonomy db for organism names",default="")
    parser.add_argument("-od","--outdir", help="List file of genes to extract (default current dir)",default="./")
    args = parser.parse_args()
    writeallgenes(args.input,args.listfile,args.ignorefile,args.outdir,args.tax,args.allgenes,args.genelimit,writeaa=args.writeaa)