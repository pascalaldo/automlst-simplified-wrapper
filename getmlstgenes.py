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

import os, argparse, setlog, json, sqlite3 as sql, numpy as np

def get16S(db,cov=800):
    conn=sql.connect(db)
    cur=conn.cursor()
    #Combine tables into view and write full hits, fragment hits, genematrix
    cur.execute("CREATE TEMP VIEW rnagenes AS SELECT A.*,B.* FROM RNAhits A INNER JOIN Seqs B ON A.seqid=B.seqid ORDER BY orgname, LENGTH(naseq) DESC")
    #Get longest 16S rRNA hit from each org
    cur.execute("SELECT orgname,seqid,naseq FROM rnagenes WHERE hmmhit='16S_rRNA' GROUP BY orgname")
    seqdict = {}
    for x in cur:
        if x[0] not in seqdict and len(x[2]) >= cov:
            seqdict[x[0]] = {}
            log.debug("RNA gene: %s - %s"%(x[0],len(x[2])))
            #Add sequences (seqid,naseq,aaseq)
            seqdict[x[0]]["RNA_16S_rRNA"] = [[x[1],x[2],False]]
    conn.close()
    return seqdict

def getcoregenes(db,pct=0.4,pct2=0.4,seqdict=None):
    conn = sql.connect(db)
    cur = conn.cursor()
    hkcore=[str(x[0]) for x in cur.execute("SELECT DISTINCT hmmhit FROM HMMhits")]
    orgs=[str(x[0]) for x in cur.execute("SELECT DISTINCT orgname FROM HMMhits")]
    orgs = list(orgs)
    oarry = np.zeros((len(hkcore), len(orgs)))
    if not seqdict:
        seqdict = {}
    if any("RNA_16S_rRNA" in x.keys() for k,x in seqdict.items()):
        log.info("Adding 16S to matrix results")
        #add 16S genes to matrix if seqdict exists
        newrow = np.zeros(len(orgs))
        for org in seqdict.keys():
            j = orgs.index(org)
            newrow[j] = 1
        hkcore.append("RNA_16S_rRNA")
        oarry = np.vstack((oarry,newrow))

    # Combine tables into view and write full hits, fragment hits, genematrix
    cur.execute("CREATE TEMP VIEW allgenes AS SELECT A.*,B.* FROM HMMhits A INNER JOIN Seqs B ON A.seqid=B.seqid WHERE flags>0")
    for i, hkgn in enumerate(hkcore):
        if "RNA_16S_rRNA" in hkgn:
            #Skip 16S already extracted
            continue
        results = cur.execute(
            "SELECT orgname,seqid,genecov,hmmcov,MAX(iscore),genelen,source,loc_start,loc_end,loc_strand,description,score,naseq,aaseq FROM allgenes WHERE hmmhit=? GROUP BY seqid",
            (hkgn,)).fetchall()
        for x in results:
            # add to genematrix
            j = orgs.index(x[0])
            # separate fragments and full length
            if x[2] >= pct or x[3] >= pct2:
                oarry[i, j] += 1
                if x[0] not in seqdict:
                    seqdict[x[0]] = {}
                if hkgn not in seqdict[x[0]]:
                    seqdict[x[0]][hkgn] = []
                seqdict[x[0]][hkgn].append([x[1],x[-2],x[-1]])

    return orgs,hkcore,oarry,seqdict

def findsingles(db,minnum=7,minorg=0.8,maxgenes=30,dnds="",outdir="",keepgenes="",lf=""):
    global log
    if lf:
        log = setlog.init(logfile=lf,level="info")
    else:
        log = setlog.init(toconsole=True,level="info")

    log.info("Getting 16S genes...")
    seqdict = get16S(db)
    log.info("Done. Getting Core genes...")
    orgs,hks,gmat,seqdict = getcoregenes(db,seqdict=seqdict)

    if len(keepgenes):
        keepgenes = keepgenes.split(",")
    else:
        keepgenes = []

    njvals = {}
    if dnds and os.path.isfile(dnds):
         with open(dnds,"r") as fil:
             njvals = json.load(fil)
    numorgs = len(orgs)
    filtinds = [i for i,x in enumerate(gmat) if float(list(x).count(1))/numorgs >= 0.95]
    filtinds.extend([i for i,x in enumerate(hks) if x in keepgenes])
    filtinds = list(set(filtinds))
    hkinds = []
    #Start with organisms that have all keepgenes
    orginds = [i for i in range(numorgs) if all(gmat[hks.index(x),i]==1 for x in keepgenes)]
    hksfound = []

    #Remove organisms with missing/duplicate genes until minimum single copy genes are found
    while float(len(orginds))/numorgs >= minorg:
        hkinds = [i for i in filtinds if float(list(gmat[i][orginds]).count(1)) == len(orginds)]
        hksfound = [hks[i] for i in hkinds]
        if len(hkinds) >= minnum and all(x in hksfound for x in keepgenes):
            log.info("Found minimum genes")
            break
        tophks = sorted([[x,i,list(x).count(1)] for i,x in enumerate(gmat) if i not in hkinds],key=lambda row: row[2],reverse=True)
        maxcount = tophks[0][2]
        tempmat = np.vstack([x[0] for x in tophks if x[2]==maxcount])
        orgcounts = sorted([[list(tempmat[:,j]).count(1),j] for j in orginds], key=lambda row: row[0])
        orginds = [x[1] for x in orgcounts[numorgs-maxcount:]]

    remorg = [x for i,x in enumerate(orgs) if i not in orginds]
    singhks = [x for i,x in enumerate(hks) if i in hkinds]
    #Prioritize by lowest dNdS values
    singhks = sorted(singhks, key=lambda k: njvals.get(k,9))
    #Move keep genes at the top of list
    for x in keepgenes:
        if x in singhks:
            singhks.insert(0, singhks.pop(singhks.index(x)))

    if len(hkinds) >= minnum and all(x in hksfound for x in keepgenes):
        log.info("# of single copy genes found: %s\t# removed organisms: %s"%(len(singhks),len(remorg)))
        log.info("Removed orgs:%s"%remorg)
        for i,x in enumerate(singhks):
            if i<maxgenes:
                log.info("Top Singles: %s\tdNdS=%s"%(x,njvals.get(x,"NA")))
            else:
                log.info("Other Singles: %s\tdNdS=%s"%(x,njvals.get(x,"NA")))
       #WRITE GENES TO FASTA
        if outdir:
            if not os.path.exists(outdir):
                os.makedirs(outdir)
            for hkgn in singhks[:maxgenes]:
                allnucrecs = [">%s\n%s"%(org,x[hkgn][0][1]) for org,x in seqdict.items() if hkgn in x.keys() and x[hkgn][0][1] and org not in remorg]
                allaarecs = [">%s\n%s"%(org,x[hkgn][0][2]) for org,x in seqdict.items() if hkgn in x.keys() and x[hkgn][0][2] and org not in remorg]
                with open(os.path.join(outdir,hkgn+".fna"),"w") as fnafil, open(os.path.join(outdir,hkgn+".faa"),"w") as faafil:
                    if len(allnucrecs):
                        fnafil.write("\n".join(allnucrecs)+"\n")
                    if len(allnucrecs):
                        faafil.write("\n".join(allaarecs)+"\n")
        return remorg,singhks[:maxgenes]
    else:
        log.error("Failed: Only %s single copy genes found after removing orgs: %s\nKeepgenes found: %s"%(len(singhks),remorg,set(singhks)&set(keepgenes)))
        log.info("Try different set or allow more deletions")
        return False,False

# Commandline Execution
if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="""Find and extract MLST genes""")
    parser.add_argument("input", help="Input sequence database")
    parser.add_argument("-dnds", help="Json file with dN/dS values for further filtering",default="")
    parser.add_argument("-log", help="Log file to save status",default="")
    parser.add_argument("-od", "--outdir", help="Directory to rewrite singles with removed organisms (default: currentdir/singles)",default="singles")
    parser.add_argument("-kg", "--keepgenes", help="list of genes that must be kept. (comma seperated, ex: RNA_16S_rRNA,TIGR02013)",default="")
    parser.add_argument("-mx","--maxgenes", help="maximum number of single copy housekeeping genes, uses lowest dNdS valued genes (default:30)",type=int,default=30)
    parser.add_argument("-mn","--minnum", help="minimum number of single copy housekeeping genes (default:7)",type=int,default=7)
    parser.add_argument("-mo","--minorg", help="minimum percent of organisms to remain (default:0.8)",type=float,default=0.8)
    args = parser.parse_args()
    findsingles(args.input,args.minnum,args.minorg,args.maxgenes,args.dnds,args.outdir,args.keepgenes,args.log)