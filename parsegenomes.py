#!/usr/bin/env python
import subprocess, os, shutil, setlog, argparse
import gbk2fa
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.Data.CodonTable import TranslationError

global log
log=setlog.init(toconsole=True)

def runprodigal(infasta, outfasta):
    try:
        #crude check for format
        # outfasta2 = os.path.splitext(outfasta)[0]+".faa"
        # cmd = ["prodigal","-d",outfasta,"-a",outfasta2,"-i",infasta,"-f","gff"]
        cmd = ["prodigal","-d",outfasta,"-i",infasta,"-f","gff"]
        with open(os.devnull,"w") as devnull:
            log.info("Prodigal: Starting file %s"%infasta)
            subprocess.call(cmd, stdout=devnull, stderr=devnull)
            log.info("Prodigal: Finished file %s"%infasta)
        return outfasta
    except Exception, e:
        log.error("Problem running prodigal for file %s"%infasta)
        return False

def fastatype(fastafil, tol=0.05):
    """Test if fasta is parsed correctly and is a multi-record ORF document by testing if some records contain multiple stop codons (ORF fasta should only have a few related to RNAs)"""
    counts = 0
    total = 0
    with open(fastafil,"r") as fil:
        try:
            fasta = SeqIO.parse(fil,"fasta")
            for rec in fasta:
                if Seq.translate(rec.seq).count("*") > 1:
                    counts += 1
                total += 1
        except (ValueError,TranslationError) as e:
            log.error("Could not translate file %s : %s"%(os.path.split(fastafil)[-1],e))
            return False
    if not total:
        log.error("No reccords found in file: %s"%fastafil)
        return False
    if float(counts)/total <= tol:
        return "orf"
    return "contig"

def parsegbks(flist, outdir):
    for fname in flist:
        try:
            gbk2fa.convertgenes(fname,outdir=outdir,rename=True,genes=True,f="rRNA",prfix="QS--")
        except Exception as e:
            log.error("Could not parse %s : %s"%(os.path.split(fname)[-1],e))

def parseall(indir,outdir):
    #allow explicit list input instead of directory
    if type(indir) is list:
        allfiles = indir
        indir = ""
    else:
        allfiles = os.listdir(indir)
    #extention filtering
    allowed = [".fasta",".fa",".faa",".fna",".fas"]
    fastafiles = [os.path.join(indir,x) for x in allfiles if os.path.splitext(x)[-1].lower() in allowed]
    allowed = [".gbff",".gbf",".gbk",".gb",".genbank"]
    gbkfiles = [os.path.join(indir,x) for x in allfiles if os.path.splitext(x)[-1].lower() in allowed]

    #Get types for fasta files and attempt to parse as fasta
    fastatypes = [(fastatype(x),x) for x in fastafiles]

    contigfiles = [x[1] for x in fastatypes if x[0] == "contig"]
    orffiles = [x[1] for x in fastatypes if x[0] == "orf"]

    log.info("Input files found: %d,%d,%d (fasta contig, fasta multi-record orf, genbank)"%(len(contigfiles),len(orffiles),len(gbkfiles)))

    #Try to parse all files with gbk ext
    log.info("Parsing all gbk and fasta files...")
    parsegbks(gbkfiles,outdir)

    #Find orfs using prodigal for contig files if exists
    for fname in contigfiles:
        #copy input files to job directory and run gene finding
        outfile = os.path.join(outdir,"QS--"+os.path.splitext(os.path.split(fname)[-1])[0])
        log.info("copying %s to %s"%(fname,outfile))
        shutil.copy(fname,outfile+".fa")
        runprodigal(fname,outfile+".fna")

    #Find all orf documents and copy to output directory
    for fname in orffiles:
        outfile = os.path.join(outdir,"QS--"+os.path.splitext(os.path.split(fname)[-1])[0])
        shutil.copy(fname,outfile+".fa") #Normalize input by using orf multi-gene fasta as whole geneome fasta
        shutil.copy(fname,outfile+".fna")
    return True

# Commandline Execution
if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="""Normalize all genome inputs into fasta whole genome record and multi-fasta orf""")
    parser.add_argument("indir", help="Input directory of sequences (genbank and/or fasta sequences)")
    parser.add_argument("outdir", help="Output directory of save fasta full genome and multi-record orfs")
    args = parser.parse_args()
    parseall(args.indir,args.outdir)
    # runprodigal(args.infasta,args.outfasta)