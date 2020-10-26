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

import argparse, sys, os, time, setlog, json, sqlite3 as sql

global log
log = setlog.init(toconsole=True)

def copydb(finput,sourcedb,ofil):
    if type(finput) is list:
        flist=[x.strip() for x in finput]
    elif type(finput) is file:
        flist=[x.strip() for x in finput]
    elif os.path.exists(finput) and "json" in finput:
        with open(finput,"r") as ifil:
            obj = json.load(ifil)
            flist = obj["selspecies"]
            flist.extend(obj["seloutgroups"])
            #Filter out query orgs
            flist = [x for x in flist if "query" not in x.lower()]
    elif os.path.exists(finput):
        with open(finput,"r") as ifil:
            flist=[x.strip() for x in ifil]
    elif "," in finput:
        flist=[x for x in finput.split(",") if os.path.exists(x)]
    else:
        log.error("No id list found")
        return False
    log.info("Copying: %s from %s"%(flist,sourcedb))
    numrecs=len(flist)
    if numrecs and os.path.exists(sourcedb):
        conn = sql.connect(sourcedb)
        try:
            conn.execute("ATTACH DATABASE '%s' as selorgs"%os.path.realpath(ofil))
            # Copy all Sequence records with matching ids
            log.info("Copying sequence records form %s organisms: [%s]"%(numrecs,','.join(flist)))
            query = "CREATE TABLE selorgs.Seqs AS SELECT * FROM Seqs WHERE orgname IN (%s)"%', '.join('?'*numrecs)
            conn.execute(query,flist)
            log.info("DONE. Copying HMM hits...")
            #Copy all corresponding HMM and RNA hits
            conn.execute("CREATE TABLE selorgs.HMMHits AS SELECT * FROM HMMHits WHERE seqid IN (SELECT seqid FROM selorgs.Seqs)")
            conn.execute("CREATE TABLE selorgs.RNAHits AS SELECT * FROM RNAHits WHERE seqid IN (SELECT seqid FROM selorgs.Seqs)")

            #Make fast indexing
            log.info("DONE. Creating index...")
            query = "CREATE UNIQUE INDEX selorgs.seqidx ON Seqs (seqid)"
            conn.execute(query)

            log.info("DONE. Commiting and saving db....")
            conn.commit()
            conn.close()
        except conn.OperationalError as e:
            log.error("failed sql operation: %s"%e)
            conn.close()
            return False
        return True
    else:
        log.error("No reccords found")
        return False

# Commandline Execution
if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="""Combine fasta files into single sql db.""")
    parser.add_argument("input", nargs="?", help="Input list of organisms (list, Json, or plain txt)", default=sys.stdin)
    parser.add_argument("-db", "--database", help="Source database to copy from", default="refseq.db")
    parser.add_argument("-o", "--out", help="Output database", default="sql_copy.db")
    args = parser.parse_args()
    copydb(args.input, args.database,args.out)