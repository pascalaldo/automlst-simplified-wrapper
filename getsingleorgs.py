#!/usr/bin/env python
import os, numpy as np, argparse, json, setlog
from Bio import SeqIO

def write16S(db,outdir,allcopy):
    conn=sql.connect(db)
    cur=conn.cursor()
    #Combine tables into view and write full hits, fragment hits, genematrix
    cur.execute("CREATE TEMP VIEW allgenes AS SELECT A.*,B.* FROM RNAhits A INNER JOIN Seqs B ON A.seqid=B.seqid")
    #Get longest 16S rRNA hit
    cur.execute("SELECT orgname,seqid,naseq FROM allgenes WHERE hmmhit=16S_rRNA GROUP BY orgname ORDER BY LENGTH(naseq) DESC LIMIT 1")
    orgs = set()
    with open(outdir+"RNA_16S_rRNA.fna","w") as nafil:
        for x in cur:
            orgs.add(x[0])
            nafil.write(">%s|%s\n%s\n"%(x[0],x[1],x[2]))
    log.info("Finished RNA_16S_rRNA.fna")
    conn.close()
    return orgs

def findsingles(infil,minnum=7,minorg=0.8,maxgenes=30,dnds="",outdir="",indir=".",log=None):
    if log == None:
        log = setlog.init(toconsole=True)
    orgs = []
    hks = []
    gmat = []
    njvals = {}
    if dnds and os.path.isfile(dnds):
         with open(dnds,"r") as fil:
             njvals = json.load(fil)
    if os.path.isfile(infil):
        with open(infil,"r") as fil:
            for line in fil:
                if line.startswith("#Gene"):
                    orgs.extend(line.strip().split()[5:])
                elif not line.startswith("#"):
                    x = line.strip().split()
                    hks.append(x[0])
                    gmat.append([int(float(v)) for v in x[5:]])
    gmat = np.vstack(gmat)
    numorgs = len(orgs)
    filtinds = [i for i,x in enumerate(gmat) if float(list(x).count(1))/numorgs >= 0.95]
    hkinds = []
    orginds = range(numorgs)
    #Remove organisms with missing/duplicate genes until minimum single copy genes are found
    while len(hkinds) < minnum and float(len(orginds))/numorgs >= minorg:
        hkinds = [i for i in filtinds if float(list(gmat[i][orginds]).count(1)) == len(orginds)]
        if len(hkinds) >= minnum:
            break
        tophks = sorted([[x,i,list(x).count(1)] for i,x in enumerate(gmat) if i not in hkinds],key=lambda row: row[2],reverse=True)
        maxcount = tophks[0][2]
        tempmat = np.vstack([x[0] for x in tophks if x[2]==maxcount])
        orgcounts = sorted([[list(tempmat[:,j]).count(1),j] for j in orginds], key=lambda row: row[0])
        orginds = [x[1] for x in orgcounts[numorgs-maxcount:]]

    remorg = [x for i,x in enumerate(orgs) if i not in orginds]
    singhks = [x for i,x in enumerate(hks) if i in hkinds]
    if len(njvals):
        singhks = sorted(singhks, key=lambda k: njvals.get(k,9))
        mdnds = njvals.get(singhks[:maxgenes][-1],"N/A")
        log.info("Lowest %s single copy genes found with dNdS < %s"%(min(maxgenes,len(singhks)),mdnds))
    if len(hkinds) >= minnum:
        log.info("# of single copy genes found: %s\t# removed organisms: %s"%(len(singhks),len(remorg)))
        log.info("Removed orgs:%s\tHKs:%s"%(remorg,singhks[:maxgenes]))
        if outdir:
            if not os.path.exists(outdir):
                os.makedirs(outdir)
            for hkgn in singhks[:maxgenes]:
                #Copy all single copy sequences to outdir after removing organisms
                seqrecs = SeqIO.parse(os.path.join(indir,hkgn+".fna"),"fasta")
                # seqrecs = [rec for rec in seqrecs if not any(x in rec.id for x in remorg)]
                with open(os.path.join(outdir,hkgn+".fna"),"w") as fnafil, open(os.path.join(outdir,hkgn+".faa"),"w") as faafil:
                    for rec in seqrecs:
                        if not any(x in rec.id for x in remorg):
                            fnafil.write(">%s\n%s\n"%(rec.id,rec.seq))
                            faafil.write(">%s\n%s\n"%(rec.id,rec.seq.translate(to_stop=True)))

        return remorg,singhks[:maxgenes],mdnds
    else:
        log.error("Failed: Only %s single copy genes found after removing orgs: %s"%(len(singhks),remorg))
        log.info("Try different set or allow more deletions")
        return False,False,False
# Commandline Execution
if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="""Use gene matrix to find MLST genes""")
    parser.add_argument("input", help="genematrix file")
    parser.add_argument("-dnds", help="Json file with dN/dS values for further filtering",default="")
    parser.add_argument("-od", "--outdir", help="Directory to rewrite singles with removed organisms",default="")
    parser.add_argument("-id", "--indir", help="Directory where sequence files are (default = current dir)",default=".")
    parser.add_argument("-mx","--maxgenes", help="maximum number of single copy housekeeping genes, uses lowest dNdS valued genes (default:30)",type=int,default=30)
    parser.add_argument("-mn","--minnum", help="minimum number of single copy housekeeping genes (default:7)",type=int,default=7)
    parser.add_argument("-mo","--minorg", help="minimum percent of organisms to remain (default:0.8)",type=float,default=0.8)
    args = parser.parse_args()
    findsingles(args.input,args.minnum,args.minorg,args.maxgenes,args.dnds,args.outdir,args.indir)