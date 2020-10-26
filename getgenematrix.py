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
#

import argparse, os, setlog, numpy as np, json, sqlite3 as sql

log = setlog.init(toconsole=True,level="info")

def getsingles(genemat,pct2=1.0):
    """Get rows with only single counts"""
    gmat = np.vstack(genemat.values())
    #Store only values of 1 in matrix
    gmat = ((gmat < 2)*(gmat > 0)).astype(int)
    thresh = len(gmat[0])*pct2
    counts = gmat.sum(axis=1)
    #get index of ubiquitus singles and send gene list
    inds = [i for i,count in enumerate(counts) if count >= thresh]
    genenames = np.array(genemat.keys())
    return list(genenames[inds]),counts

def getorgsets(genemat,orgs):
    """Store sets of singles per organism"""
    orgsets = {}
    # get transposed matrix
    gmat = np.vstack(genemat.values()).T
    genenames = np.array(genemat.keys())
    for i,org in enumerate(orgs):
        inds = [j for j,count in enumerate(gmat[i]) if count == 1]
        if org not in orgsets:
            orgsets[org] = set(genenames[inds])
        else:
            log.warning("Duplicate organism name found")
    return orgsets

def getgenesets(genemat,orgs,missing=False):
    """Store sets of singles per organism"""
    genesets = {}
    # get transposed matrix
    gmat = np.vstack(genemat.values())
    orgnames = np.array(orgs)
    allset = set(orgnames)
    for i,gene in enumerate(genemat.keys()):
        inds = [j for j,count in enumerate(gmat[i]) if count == 1]
        if gene not in genesets:
            #Store missing organisms in dictionary
            if missing:
                genesets[gene] = allset - set(orgnames[inds])
            else:
                genesets[gene] = set(orgnames[inds])
        else:
            log.warning("Duplicate gene name found")
    return genesets

def getTopFuncs(genelist,maxgenes=100,exfuncs=[]):
    #get list of functions of top genes (<maxgenes)
    topfuncs = [x["func"] for x in genelist[:maxgenes] if x["func"] not in exfuncs]
    return {f:float(topfuncs.count(f))/maxgenes for f in set(topfuncs)}

def rebalancefuncs(genelist,OVthresh=1.1,maxgenes=100,maxiter=100, exfuncs=[], EQratio=0, ubiq=False):
    """Increase diversity of functions in group < maxgenes by replacing
        lowest priority over-represented function (also replace excluded functions)
        with highest priority under-represented function"""

    topfuncs = getTopFuncs(genelist,maxgenes,exfuncs)

    if maxiter<1:
        return genelist
    if EQratio <= 0:
        allfuncs = set(x["func"] for x in genelist if x["func"] not in exfuncs)
        EQratio = 1.0/len(allfuncs)
        log.debug(allfuncs)
    #First replace any exclude functions
    OVfunc = [x["func"] for x in genelist[:maxgenes] if x["func"] in exfuncs]
    #otherwise use most over-represented function
    if not OVfunc:
        OVfunc = [f for f,ratio in topfuncs.items() if ratio > OVthresh*EQratio]
        OVfunc = sorted(OVfunc,key=lambda x: topfuncs.get(x,9),reverse=True)

    #Get under represented functions
    UVfunc = [f for f,ratio in topfuncs.items() if ratio < EQratio]
    # UVfunc = sorted(UVfunc,key=lambda x: topfuncs.get(x,9))

    if OVfunc and UVfunc:
        OVfunc = OVfunc[0]
        # find lowest priority OVfunc and replace with any highest priority under represented func
        # loop index backward from maxgenes
        for i in range(maxgenes-1,-1,-1):
            func = genelist[i].get("func","")
            if func and func in OVfunc:
                #Find highest priority under represented function outside of maxgenes also impose that it must also be ubiquitous
                for j in range(maxgenes,len(genelist)):
                    func = genelist[j].get("func","")
                    dcount = 0
                    if ubiq:
                        dcount = genelist[j].get("delcount",0)
                    if func in UVfunc and dcount==0:
                        #swap values
                        genelist[i],genelist[j] = genelist[j],genelist[i]
                        # log.debug("moved: %d and %d. (iteration %d)"%(j,i,maxiter))
                        break
                break
        #Recursively return reshuffled list
        maxiter -= 1
        genelist = rebalancefuncs(genelist,OVthresh=OVthresh,maxgenes=maxgenes,maxiter=maxiter,exfuncs=exfuncs,EQratio=EQratio,ubiq=ubiq)
        return genelist
    else:
        return genelist

def prioritize(genemat,orgs,dndsfile="",jsonfile="",metadata="",pct=0.5,maxgenes=100,ignoreorgs=None):
    if not ignoreorgs:
        ignoreorgs = []
    gsets = getgenesets(genemat,orgs,missing=True)
    if not metadata or not os.path.exists(metadata):
        metadata = os.path.join(os.path.dirname(os.path.realpath(__file__)),"model_metadata.json")
    if not dndsfile or not os.path.exists(dndsfile):
        dndsfile = os.path.join(os.path.dirname(os.path.realpath(__file__)),"dnds.json")

    #Load metadata and dnds values:
    with open(metadata,"r") as mfil, open(dndsfile,"r") as dfil:
        metadata = json.load(mfil)
        dndsfile = json.load(dfil)
    numorgs = len(genemat.values()[0])
    genelist = []
    for gene,notsingle in gsets.items():
        #Skip gene if number of organisms needing removal are over pct or if any ignore orgs are in delete set
        if float(len(gsets.get(gene,[])))/numorgs >= pct or any([x in ignoreorgs for x in gsets.get(gene,[])]):
            continue

        rec = {"acc":metadata.get(gene,{}).get("acc",gene),
               "dnds":dndsfile.get(gene,9),"name":metadata.get(gene,{}).get("name",gene),
               "func":metadata.get(gene,{}).get("func","N/A"),
               "desc":metadata.get(gene,{}).get("desc","N/A"),
               "orgdel":list(gsets.get(gene,[])),
               "delcount":len(gsets.get(gene,[]))}
        genelist.append(rec)

    #Sort by number orgs needing removal, segment excluded functions to bottom of each group, then sort by lowest dn/ds value
    exfuncs = ['Unknown function','Hypothetical proteins','Mobile and extrachromosomal element functions','Unclassified']
    genelist = sorted(genelist,key=lambda x:(x["delcount"],int(x["func"] in exfuncs),x["dnds"]))
    if jsonfile:
        with open(jsonfile,"w") as fil:
            json.dump(genelist,fil,indent=2)

    # log.debug("Topgenes Before re-balancing: %s"%getTopFuncs(genelist,maxgenes,exfuncs))
    # log.debug("Delete Org counts: %s"%set([x["delcount"] for x in genelist[:maxgenes]]))

    #Rebalance function diversity for MLST singles
    #genelist = rebalancefuncs(genelist,exfuncs=exfuncs,maxgenes=maxgenes,maxiter=900,ubiq=True)

    # log.debug("Topgenes After re-balancing: %s"%getTopFuncs(genelist,maxgenes,exfuncs))
    # log.debug("Delete Org counts: %s"%set([x["delcount"] for x in genelist[:maxgenes]]))

    return genelist

def getmat(db,pct=0.5,pct2=1.0,ev=0.1,bs=0,bh=True,rna=False,savefil="",prifile="",maxgenes=100,dndsfile="",ignoreorgs=None):
    conn=sql.connect(db)
    cur=conn.cursor()
    tabletitle="HMMhits"
    orgs=sorted([x[0] for x in cur.execute("SELECT DISTINCT orgname FROM "+tabletitle)])
    counts=cur.execute("SELECT orgname,hmmhit,COUNT(DISTINCT seqid) FROM %s WHERE hmmcov >= %.4f AND evalue<=%s AND score>=%s AND flags>=%d GROUP BY orgname,hmmhit"%(tabletitle,pct,ev,bs,bh))
    countdict={}
    for x in counts:
        if x[1] not in countdict:
            countdict[x[1]]=[0]*len(orgs)
        countdict[x[1]][orgs.index(x[0])]=int(x[2])

    if rna:
        tabletitle="RNAhits"
        counts=cur.execute("SELECT orgname,hmmhit,COUNT(DISTINCT seqid) FROM %s WHERE hmmcov >= %.4f AND evalue<=%s AND score>=%s AND flags>=%d GROUP BY orgname,hmmhit"%(tabletitle,pct,ev,bs,bh))
        for x in counts:
            if x[1] not in countdict:
                countdict[x[1]]=[0]*len(orgs)
            countdict[x[1]][orgs.index(x[0])]=int(x[2])

    #save genematrix string:
    gmstr="#totalorgs:\t"+str(len(orgs))+"\n"
    gmstr+="#hmmhit\tOrgcount\tAvg\tStdev\t"+"\t".join(orgs)+"\n"
    for k,v in countdict.items():
        temp = [c for c in v if c > 0] #Exclude mising genes in statistics
        norg = len(temp)
        avgcount = float(np.mean(temp))
        stdcount = float(np.std(temp))

        gmstr+="%s\t%d\t%.2f\t%.4f\t"%(k,norg,avgcount,stdcount)+"\t".join([str(c) for c in v])+"\n"

    conn.close()

    if savefil:
        with open(savefil,"w") as fil, open(os.path.splitext(savefil)[0]+".txt","w") as tfil:
            json.dump({"counts":countdict,"orgs":orgs},fil,indent=2)
            tfil.write(gmstr+"\n")

    if prifile:
        mlst = prioritize(countdict,orgs,jsonfile=prifile,dndsfile=dndsfile,maxgenes=maxgenes,pct=pct,ignoreorgs=ignoreorgs)
    else:
        mlst = getsingles(genemat=countdict,pct2=pct2)

    return countdict,orgs,mlst

# Commandline Execution
if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="""Determine gene matrix to prioritize MLST genes""")
    parser.add_argument("input", help="Sql database (READS TABLES: Seqs, HMMhits)")
    parser.add_argument("-p", "--pct", help="percent match of hmm model length to consider as hit (default: 0.5)", type=float, default=0.5)
    parser.add_argument("-p2", "--pct2", help="percent of orgs to consider as ubiquitous single copy gene (default: 1.0)", type=float, default=1.0)
    parser.add_argument("-mxg", "--maxgenes", help="Limit MLST genes to top X genes (default: 100)", type=int, default=100)
    parser.add_argument("-e", "--evalue", help="Evalue threshold (default: 0.1)", type=float, default=0.1)
    parser.add_argument("-b", "--bitscore", help="Bitscore threshold (default: 0)", type=float, default=0)
    parser.add_argument("-sf", "--savefil", help="Save file to json file (default: disabled)", default="")
    parser.add_argument("-pf", "--priorityfile", help="Save list of top priority MLST genes (based on ubiquity and lowest dn/ds values)", default="")
    parser.add_argument("-df", "--dndsfile", help="JSON file with precalculated dN/dS values used for prioritization", default="")
    parser.add_argument("-bh", "--besthit", help="Only write best hit hmm for each gene (only lowest evalue hmm allowed per gene)", action='store_true')
    parser.add_argument("-r", "--rna", help="Also copy RNA table hits", action='store_true')
    args = parser.parse_args()
    result=getmat(args.input, args.pct, args.pct2, args.evalue, args.bitscore, args.besthit, args.rna, args.savefil, args.priorityfile, args.maxgenes)