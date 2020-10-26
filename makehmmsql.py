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
import argparse, os, setlog, sqlalchemy as sql
import itertools
from sqlalchemy import exc

global log
log = setlog.init(toconsole=True)

def getcoverage(recs):
    for seqid,rows in recs.items():
        covs = {}
        for row in rows:
            if row[0] not in covs:
                covs[row[0]] = [set(),set(),row[5],row[8]]
            #merge hmm coords into set
            covs[row[0]][0] |= set(xrange(row[3],row[4]))
            #merge gene coords into set
            covs[row[0]][1] |= set(xrange(row[6],row[7]))
        for i,row in enumerate(rows):
            #calculate hmm coverage ratio
            recs[seqid][i][-2] = float(len(covs[row[0]][0]))/float(covs[row[0]][2])
            #calculate gene coverage ratio
            recs[seqid][i][-1] = float(len(covs[row[0]][1]))/float(covs[row[0]][3])
    return recs

def markbest(recs):
    for seqid,rows in recs.items():
        mineval = min([x[-7] for x in rows])
        bestmodel = [x[0] for x in rows if x[-7] == mineval][0]
        for i,row in enumerate(rows):
            if row[0] == bestmodel:
                recs[seqid][i][-3] = 10
    return recs

def run(fname,ofil,ev=1e-4,bs=0.0,rna=False,filt2=None):
    sep = '?'
    pg = False
    if ofil.lower().startswith("postgresql://"):
        pg = True
        sep = '%s'
    elif not ofil.lower().startswith("sqlite://"):
        ofil = "sqlite:///"+os.path.realpath(ofil)

    eng = sql.create_engine(ofil)
    csr = eng.connect()
    tabletitle="HMMhits"
    if rna:
        tabletitle="RNAhits"
    #Create tables
    metadata = sql.MetaData()
    hittable = sql.Table(tabletitle,metadata,
                     sql.Column('hmmhit',sql.Text),
                     sql.Column('orgname',sql.Text),
                     sql.Column('seqid',sql.Integer),
                     sql.Column('hmmstart',sql.Integer),
                     sql.Column('hmmend',sql.Integer),
                     sql.Column('hmmlen',sql.Integer),
                     sql.Column('geneStart',sql.Integer),
                     sql.Column('geneEnd',sql.Integer),
                     sql.Column('genelen',sql.Integer),
                     sql.Column('evalue',sql.Float),
                     sql.Column('score',sql.Float),
                     sql.Column('bias',sql.Float),
                     sql.Column('iscore',sql.Float),
                     sql.Column('flags',sql.Integer),
                     sql.Column('hmmcov',sql.Float),
                     sql.Column('genecov',sql.Float)
                    )
    metadata.create_all(eng)
    if filt2 and os.path.isfile(filt2):
        temp={}
        with open(filt2,"r") as fil:
            for line in fil:
                x=line.strip().split()
                temp[x[0]]=float(x[1])
        filt2=temp
    #Add hmmhits:
    with open(fname,"r") as ifil:
        recs={}
        for line in ifil:
            if line[0]!="#":
                x = line.split()
                org,seqid = x[0].split("|")
                if filt2 and x[3] in filt2.keys() and filt2[x[3]]>float(x[7]):
                    continue
                if float(x[6]) < ev and float(x[7]) > bs:  #global quality filter
                    if seqid not in recs:
                        recs[seqid] = []
                    recs[seqid].append([x[3],org,int(seqid),int(x[15]),int(x[16]),int(x[5]),int(x[19]),int(x[20]),int(x[2]),float(x[6]),float(x[7]),float(x[8]),float(x[13]),0,0,0])

        if not len(recs.keys()):
            log.warning("No %s were found"%tabletitle)
            return False
        #Mark Best Hits in hmm search
        recs = markbest(recs)
        #Calculate gene and hmm coverage
        recs = getcoverage(recs)
        recs = list(itertools.chain.from_iterable(recs.values()))
        csr.execute("INSERT INTO \"%(t)s\" VALUES (%(x)s,%(x)s,%(x)s,%(x)s,%(x)s,%(x)s,%(x)s,%(x)s,%(x)s,%(x)s,%(x)s,%(x)s,%(x)s,%(x)s,%(x)s,%(x)s)"%{"x":sep,"t":tabletitle},recs)
        log.info("Added %d records"%len(recs))
        del recs
    # Remove duplicate rows
    try:
        log.info("Removing duplicates...")
        if pg:
            csr.execute('DELETE FROM "'+tabletitle+'" a WHERE EXISTS (SELECT 1 FROM "'+tabletitle+'" t WHERE t.hmmhit = a.hmmhit and t.seqid = a.seqid and t.hmmstart = a.hmmstart and t.hmmend = a.hmmend and t.ctid > a.ctid)')
        else:
            csr.execute('DELETE FROM "'+tabletitle+'" WHERE rowid NOT IN (SELECT min(t.rowid) FROM "'+tabletitle+'" t GROUP BY hmmhit,seqid,hmmstart,hmmend)')
        # log.info("Marking best hits...")
        # csr.execute('UPDATE "'+tabletitle+'" SET flags=0') # clear all previous
        # csr.execute('UPDATE "'+tabletitle+'" SET flags=10 WHERE seqid || hmmhit IN (SELECT seqid || hmmhit FROM "'+tabletitle+"\" GROUP BY seqid HAVING evalue = MIN(evalue))")
    except exc.OperationalError as ex:
        log.error("exception:"%ex)
    #
    # temptable="temp_%d%d"%(os.getpid(),time.time()) # make temporary table name
    # try:
    #     log.info("Calculating total coverage...")
    #
    #     #get hmmalign coords:
    #     result = csr.execute("SELECT seqid,hmmhit,Group_concat('('||hmmstart||','||hmmend||')','&'),Group_concat('('||geneStart||','||geneEnd||')','&') FROM \""+tabletitle+'" GROUP BY seqid,hmmhit')
    #     seqcovs=list(gethmmcov(result))
    #     csr.execute("CREATE TEMP TABLE \"%s_cov\" (seqid int, hmmhit text, hmmcov real, genecov real)"%temptable)
    #     csr.execute("INSERT INTO \"%(t)s_cov\" VALUES (%(x)s,%(x)s,%(x)s,%(x)s)"%{"x":sep,"t":temptable},seqcovs)
    #
    #     sqlcmd="""
    #         CREATE TABLE "%s" AS SELECT A.hmmhit,A.orgname,A.seqid,A.hmmstart,A.hmmend,A.hmmlen,A.geneStart,A.geneEnd,A.genelen,A.evalue,A.score,A.bias,A.iscore,A.flags,1.*B.hmmcov/A.hmmlen as hmmcov,1.*B.genecov/A.genelen as genecov
    #         FROM "%s" A INNER JOIN %s_cov B ON A.seqid=B.seqid AND A.hmmhit=B.hmmhit""" % (temptable,tabletitle,temptable)
    #     csr.execute(sqlcmd)
    #     csr.execute("DROP TABLE \"%s\""%tabletitle)
    #     # csr.execute("DROP TABLE %s_cov"%temptable)
    #     csr.execute('ALTER TABLE "'+temptable+"\" RENAME TO \"%s\""%tabletitle)
    #     # print "Compacting db size..."
    #     # csr.execute("VACUUM")
    # except exc.OperationalError as ex:
    #     log.error("Aborting coverage calculation, Error:%s"%ex)

    csr.close()
    return True

# Commandline Execution
if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="""Combine hmmresults into single sql db""")
    parser.add_argument("input", help="Input Hmm result file")
    parser.add_argument("out", help="Output database")
    parser.add_argument("-r", "--rna", help="Add to RNA table instead of HK gene table", action='store_true')
    parser.add_argument("-e", "--evalue", help="Remove hits > evalue default(1e-4)", type=float, default=1e-4)
    parser.add_argument("-b", "--bitscore", help="Remove hits < bitscore default(0)", type=float, default=0)
    parser.add_argument("-f2", "--filter", help="Use per-model bitscore cutoffs dictated in file (format: modelname value)", default=None)
    args = parser.parse_args()
    run(args.input,args.out,args.evalue,args.bitscore,args.rna,args.filter)