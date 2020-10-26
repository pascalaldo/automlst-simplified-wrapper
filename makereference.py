#!/usr/bin/env python
import argparse, subprocess, glob, os, setlog
from concatmsa import concatmsa

log = setlog.init(toconsole=True)

def trimal(infil,outfile):
    cmd = ["trimal","-automated1","-in",infil,"-out",outfile]
    with open(os.devnull,"w") as devnull:
        try:
            subprocess.call(cmd,stdout=devnull)
            log.debug("TrimAl: finished %s"%outfile)
            return True
        except subprocess.CalledProcessError as e:
            log.error("TrimAl: error, could not process %s - %s"%(outfile,e))
            return False

def codonalign(inalign,nafile,outputfile):
    cmd = ["perl",os.path.join(os.path.dirname(os.path.realpath(__file__)),"pal2nal.pl"),"-output", "fasta", inalign, nafile]
    with open(outputfile,"w") as ofil:
        try:
            subprocess.call(cmd, stdout=ofil)
            log.debug("Pal2Nal: finished %s"%outputfile)
            return True
        except subprocess.CalledProcessError as e:
            log.error("Pal2Nal: error, could not process %s - %s"%(inalign,e))
            return False

def maftalign(seq,outfile,cpu=1):
    #ofil = tempfile.NamedTemporaryFile(dir=os.path.split(outfile)[0])
    cmd = ["mafft-linsi","--quiet",seq]
    if cpu>1:
        cmd[1:1] = ["--thread","%s"%cpu]
    # if os.path.isfile(outfile+".tmafft"):
    #     cmd[1:1] = ["--treein", outfile+".tmafft"]
    try:
        # ofil = tempfile.NamedTemporaryFile(dir=os.path.split(newseqs)[0])
        with open(outfile,"w") as ofil:
            subprocess.call(cmd,stdout=ofil)
        log.debug("MAFFT: finished %s"%seq)
        return True
    except subprocess.CalledProcessError as e:
        log.error("MAFFT: error, failed to align %s"%outfile)
        return False

def raxml(outdir,inalgn,bootstrap=1000,mcpu=1,part=None):
    #Run raxml EPA
    fname = os.path.split(inalgn)[-1]
    cmd="raxmlHPC-SSE3 -f a -m GTRGAMMAI -p 12345 -x 12345 -N %s -w %s -s %s -n %s"%(bootstrap,os.path.realpath(outdir),inalgn,fname)
    if mcpu >1:
        cmd+=" -T %s"%mcpu
    if part and os.path.exists(part):
        cmd+=" -q %s"%part
    try:
        log.debug("Starting: %s"%cmd)
        with open(os.devnull,"w") as devnull:
            subprocess.call(cmd.split(),stdout=devnull)
        log.debug("RAxML: finished %s"%inalgn)
        return True
    except subprocess.CalledProcessError as e:
        log.error("RAxML: error, could not process %s - %s"%(inalgn,e))
        return False


def makeref(indir,outdir,mcpu=1,concat=False):
    singlist = [os.path.split(x)[-1] for x in glob.glob(os.path.join(indir,"*.faa"))]
    if os.path.exists(os.path.join(indir,"RNA_16S_rRNA.fna")):
        singlist.append("RNA_16S_rRNA.fna")

    os.makedirs(os.path.join(outdir,"trees"))

    log.info("Starting alignments...")
    finished=[]
    for fname in singlist:
        outfil = os.path.join(outdir,fname)
        outfil2 = os.path.join(outdir,os.path.splitext(fname)[0]+".fna")
        fnafil = os.path.join(indir,os.path.splitext(fname)[0]+".fna")
        trimfil = outfil2+".trimmed"
        if "RNA_16S_rRNA.fna" in fname and maftalign(os.path.join(indir,fname),outfil,mcpu):
            if trimal(outfil,trimfil):
                finished.append(trimfil)
        elif maftalign(os.path.join(indir,fname),outfil,mcpu):
            if codonalign(outfil,fnafil,outfil2):
                if trimal(outfil2,trimfil):
                    finished.append(trimfil)
    if concat:
        log.info("Making supermatrix...")
        outfil = os.path.join(outdir,"supermatrix.fa")
        partfil = os.path.join(outdir,"nucpart.txt")
        concatmsa(os.path.join(outdir,"*.trimmed"),outfil,partfil)
        log.info("Done. Building tree...")
        raxml(os.path.join(outdir,"trees"),outfil,mcpu=mcpu,part=partfil)
    else:
        log.info("Starting buildtrees...")
        for fname in finished:
            raxml(os.path.join(outdir,"trees"),fname,mcpu=mcpu)

# Commandline Execution
if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="""Use extracted MLST genes and builds reference folder with single copy codon alignments""")
    parser.add_argument("input", help="Directory of extracted MLST genes")
    parser.add_argument("outdir", help="Output directory of alignments and trees")
    parser.add_argument("-cpu", "--multicpu", help="Turn on Multi processing set # Cpus (default: Off, 1)", type=int, default=1)
    parser.add_argument("-cat", "--concatmsa", help="Build species tree using concatenated supermatrix",action="store_true",default=False)
    args = parser.parse_args()
    makeref(args.input,args.outdir,args.multicpu,args.concatmsa)