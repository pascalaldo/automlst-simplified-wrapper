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

import argparse, sys, os, time, setlog, sqlalchemy as sql
from sqlalchemy import exc

global log
log = setlog.init(toconsole=True)

def runlist(finput,ofil,transonly=False,orgname=False,inc=False):
    if type(finput) is list:
        flist=[x.strip() for x in finput if os.path.exists(x.strip())]
    elif type(finput) is file:
        flist=[x.strip() for x in finput if os.path.exists(x.strip())]
    elif os.path.exists(finput):
        with open(finput,"r") as ifil:
            flist=[x.strip() for x in ifil if os.path.exists(x.strip())]
    elif "," in finput:
        flist=[x for x in finput.split(",") if os.path.exists(x)]
    else:
        log.error("No File list found")
        return False
    numrecs=len(flist)
    if numrecs:
        sep = '?'
        pg = False
        if ofil.lower().startswith("postgresql://"):
            pg = True
            sep = '%s'
        elif not ofil.lower().startswith("sqlite://"):
            ofil = "sqlite:///"+os.path.realpath(ofil)

        # #Database specific sql
        # if ofil.lower().startswith("mysql://"):
        #     sep = "%s"
        #     AutoInc = " NOT NULL AUTO_INCREMENT"
        #     # login,serv = ofil[8:].split("@")
        #     # user,passwd = login.split(":")
        #     # host,db = serv.split("/")
        #     # port = 3306
        #     # if ":" in host:
        #     #     host,port = host.split(":")
        #     # conn = mysql.connect(user=user,passwd=passwd,host=host,db=db,port=int(port))
        # else:
        #     sep = "?"
        #     AutoInc = ""
        eng = sql.create_engine(ofil)
        csr = eng.connect()
        addedtime = int(time.time())
        #Make table if doesnt exist:
        try:
            metadata = sql.MetaData()
            Seqs = sql.Table('Seqs',metadata,
                             sql.Column('seqid',sql.Integer,primary_key=True),
                             sql.Column('orgname',sql.Text),
                             sql.Column('gene',sql.Text),
                             sql.Column('description',sql.Text),
                             sql.Column('source',sql.Text),
                             sql.Column('loc_start',sql.Integer),
                             sql.Column('loc_end',sql.Integer),
                             sql.Column('loc_strand',sql.Integer),
                             sql.Column('acc',sql.Text),
                             sql.Column('lastscan',sql.Integer),
                             sql.Column('naseq',sql.Text),
                             sql.Column('aaseq',sql.Text)
                            )
            metadata.create_all(eng)
            # sqlcmd="CREATE TABLE Seqs (seqid INTEGER PRIMARY KEY, orgname text, gene text, description text, source text, loc_start int, loc_end int, loc_strand int, acc text, lastscan int, naseq text, aaseq text)"

#Creating table copy:

#             CREATE TABLE xxx.Seqs (seqid INTEGER PRIMARY KEY, orgname text, gene text, description text, source text, loc_start int, loc_end int, loc_strand int, acc text, lastscan int, naseq text, aaseq text)
#             INSERT INTO xxx.Seqs (orgname, gene, description, source, loc_start, loc_end, loc_strand, acc, lastscan, naseq, aaseq) Select orgname, gene, description, source, loc_start, loc_end, loc_strand, acc, lastscan, naseq, aaseq from Seqs where orgname in (Select assembly_id from taxa.taxonomy where genus_name = 'xxx')

            # csr.execute(sqlcmd)
            log.info("Creating db...")
        except exc.OperationalError as e:
            log.info("Table exists. Adding to database...")

        # ex4istingrecs=set(str(x[0])+"@"+str(x[1]) for x in csr.execute("SELECT orgname,gene FROM Seqs"))
        #read each file
        allorgs = []
        for nr,fname in enumerate(flist):
            recs=[]
            temp=None
            aminod=None
            fp,fn = os.path.split(fname)
            org,ext = os.path.splitext(fn)
            if orgname:
                org = orgname
            org = org.replace(" ","_")
            allorgs.append(org)
            #skip .faa files, test type for gbk or fasta
            if ext.lower()==".faa":
                continue
            elif ext.lower()==".fna":
                #test if .faa exists for .fna pair and get seqs in dictionary
                if len(fp):
                    fp+="/"
                if os.path.isfile(fp+org+".faa"):
                    with open(fp+org+".faa","r") as faafil:
                        aminod=[]
                        reci=-1
                        for line in faafil:
                            if line[0]==">":
                                if " " in line:
                                    idx=line.index(" ")
                                    gn = line[1:idx]
                                else:
                                    gn = line.strip()
                                # if gn not in aminod:
                                #     aminod[gn]=""
                                reci+=1
                                aminod.append([gn,""])
                            else:
                                aminod[reci][1]+=line.strip()

                def addamino(xrow,i,transtab=1):
                    if not transonly and aminod and aminod[i][0] == xrow[1]:
                        xrow.append(aminod[i][1])
                    else:
                        from Bio.Seq import Seq
                        from Bio import Data
                        try:
                            aseq=Seq(xrow[-1]).translate(to_stop=True,table=transtab)
                            xrow.append(str(aseq))
                        except Data.CodonTable.TranslationError, e:
                            log.error("Unable to get amino acid translation for reccord: %s,%s"%e)
                    return xrow
                #Read each file and collect sequences
                with open(fname,"r") as ifil:
                    reci=-1
                    for line in ifil:
                        if line[0]==">":
                            if temp: #add last temp
                                recs.append(addamino(temp,reci))
                            reci+=1
                            gn = ds = line.strip()
                            source = "unknown"
                            line=line.replace('"','').replace("'","")
                            if " " in line and "|" in line:
                                idx=line.index(" ")
                                gn = line[1:idx]
                                source = line[1:idx].split("|")[-1]
                                ds = line[idx+1:].strip()
                            if "_PLASMID_" in ds:
                                source+="_PLASMID"
                            locs = (0,0,0)
                            if "|loc|" in ds:
                                idx = ds.index("|loc|")+5
                                locs = [int(x) for x in ds[idx:].split()[:3]]
                            acc = "none"
                            if "|ACC=" in ds:
                                idx = ds.index("|ACC=")+5
                                acc = ds[idx:].split()[0]
                            temp = [org,gn,ds,source,locs[0],locs[1],locs[2],acc,addedtime,""]
                        elif temp:
                            temp[-1]+=line.strip().upper()
                    if temp:
                        recs.append(addamino(temp,reci))

            # Check records exist
            if not len(recs):
                log.error("No sequence records found in: %s"%fn)
                continue

            if inc:
                # Explictly increment seqid values for insert
                maxid = csr.execute('Select Max(seqid) from "Seqs"').fetchone()[0]
                for i,r in enumerate(recs):
                    maxid += 1
                    recs[i] = [maxid] + r

                csr.execute('INSERT INTO "Seqs" (seqid, orgname, gene, description, source, loc_start, loc_end, loc_strand, acc, lastscan, naseq, aaseq) values (%(x)s,%(x)s,%(x)s,%(x)s,%(x)s,%(x)s,%(x)s,%(x)s,%(x)s,%(x)s,%(x)s,%(x)s)'%{"x":sep},recs)
            else:
                #Assume auto-increment
                csr.execute('INSERT INTO "Seqs" (orgname, gene, description, source, loc_start, loc_end, loc_strand, acc, lastscan, naseq, aaseq) values (%(x)s,%(x)s,%(x)s,%(x)s,%(x)s,%(x)s,%(x)s,%(x)s,%(x)s,%(x)s,%(x)s)'%{"x":sep},recs)
            log.info("added %d of %d files (%s  - %d reccords)" % (nr+1,numrecs,fname,len(recs)))
        try:
            log.info("Removing duplicates...")
            if pg:
                csr.execute('DELETE FROM "Seqs" a WHERE EXISTS (SELECT 1 FROM "Seqs" t WHERE t.orgname = a.orgname and t.gene = a.gene and t.loc_start = a.loc_start and t.loc_end = a.loc_end and t.loc_strand = a.loc_strand and t.seqid > a.seqid)')
            else:
                csr.execute('DELETE FROM "Seqs" WHERE seqid NOT IN (SELECT min(t.seqid) FROM "Seqs" t GROUP BY orgname,gene,description,loc_start,loc_end,loc_strand)')
        except exc.OperationalError as e:
            log.error("failed sql operation: %s"%e)

        try:
            log.info("Rebuilding index if present...")
            csr.execute("REINDEX seqidx")
            log.info("Done indexing")
        except exc.OperationalError:
            log.info("No index present")
        csr.close()
        return allorgs
    else:
        log.error("No reccords found")
        return False

# Commandline Execution
if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="""Combine fasta files into single sql db.""")
    parser.add_argument("input", nargs="?", help="Input filelist", default=sys.stdin)
    parser.add_argument("-o", "--out", help="Output sqlite file (default: seqsql.db) or postgre sql db connection 'postgresql://user:pass@host:5432/dbname'", default="seqsql.db")
    parser.add_argument("-org", "--orgname", help="Organism name (default: filename)", default="")
    parser.add_argument("-t", "--trans", help="Only store translation of DNA for protein seqs (default: False)", action='store_true')
    parser.add_argument("-i", "--inc", help="Explicitly increment the sequence id (necessary for table generated with 'CREATE ... AS' )", action='store_true')
    args = parser.parse_args()
    runlist(args.input,args.out,args.trans,args.orgname,args.inc)