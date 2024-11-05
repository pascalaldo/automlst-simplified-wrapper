#!/usr/bin/env python
import argparse, subprocess, glob, os, setlog
import gbk2fa, makeseqsql, seqsql2fa, makehmmsql #fasta2genes
from werkzeug.utils import secure_filename
from Bio import Entrez
from backports import tempfile

global log
log=setlog.init(toconsole=True)

def getNCBIgbk(acc, folder):
    try:
        #crude check for format
        if acc.upper().startswith("GCA_") or acc.upper().startswith("GCF_"):
            log.error("Error, cannot use assembly accession number, please use Genbank or Refseq sequence accession")
            filename = False
        elif acc.replace("_","").replace(".","").isalnum() and len(acc) <= 20 and len(acc) >= 6:
            if "." in acc:
                acc = acc[:acc.index(".")]
            Entrez.email = "amlstuser@ziemertlab.com"
            handle = Entrez.efetch(db="nucleotide", rettype="gb", id=acc, retmode="text")
            filename = secure_filename(acc+".gbk")
            with open(filename,"w") as outfile:
                outfile.write(handle.read())
        else:
            log.error("Error with Accession number ")
            filename = False
    except Exception, e:
        log.error("Error retrieving gbk from NCBI")
        filename = False
    return filename

def hmmsearch(outname,hmmdb,infile,mcpu=1,cut="tc"):
    log.info("Searching HMMs...")
    cmd=["hmmsearch", "--domtblout", outname, "--noali", "--notextw", hmmdb, infile]
    if mcpu>1:
        cmd[1:1] = ["--cpu", str(mcpu)]
    if cut and cut.lower() in ("ga","tc","nc"):
        cmd[1:1] = ["--cut_%s"%cut.lower()]
    with open(os.devnull,"w") as devnull:
        subprocess.call(cmd, stdout=devnull, stderr=devnull)
    #subprocess.call(cmd)

def runall(indir,outdb,mcpu=1,remotelist=None,rnamodels=None,aamodels=None):
    with tempfile.TemporaryDirectory() as tempfolder:
        convfolder = os.path.join(tempfolder,"converted")
        if remotelist and os.path.isfile(remotelist):
            log.info("Downloading remote list of accessions...")
            with open(remotelist,"r") as fil:
                for line in fil:
                    if getNCBIgbk(line.strip(),indir):
                        log.info("Downloaded %s"%line.strip())
                    else:
                        log.error("Failed to download %s"%line.strip())
        gbk2fa.runall(os.path.realpath(indir),convfolder,True,mcpu,"rRNA",True)
    
        refdb = os.path.realpath(outdb)
        aaseqs = os.path.join(tempfolder,"allseqs.faa")
        naseqs = os.path.join(tempfolder,"allseqs.fna")
        makeseqsql.runlist(glob.glob(os.path.join(convfolder,"*.fna")),refdb,True)
        seqsql2fa.writefasta(refdb,aaseqs,latest=True)
        seqsql2fa.writefasta(refdb,naseqs,nuc=True,latest=True)
    
        #Run Hmm searches on seqs
        if not rnamodels:
            rnamodels = os.path.join(os.path.dirname(os.path.realpath(__file__)),"barnap_bact_rRna.hmm")
        if not aamodels:
            #aamodels = os.path.join(os.path.dirname(os.path.realpath(__file__)),"allTFequivPFcore.hmm")
            aamodels = os.path.join(os.path.dirname(os.path.realpath(__file__)),"reducedcore.hmm")
    
        aaresult = os.path.join(tempfolder,"aaresult.domhr")
        rnaresult = os.path.join(tempfolder,"rnaresult.domhr")
        hmmsearch(aaresult,aamodels,aaseqs,mcpu=mcpu)
        hmmsearch(rnaresult,rnamodels,naseqs,mcpu=mcpu)
        makehmmsql.run(aaresult,refdb)
        makehmmsql.run(rnaresult,refdb,rna=True)

# Commandline Execution
if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="""Takes collection of genbank files as input and builds sql database with hmm result hits""")
    parser.add_argument("input", help="Directory of gbk files")
    parser.add_argument("outdb", help="Output directory of alignments and trees")
    parser.add_argument("-cpu", "--multicpu", help="Turn on Multi processing set # Cpus (default: Off, 1)", type=int, default=1)
    parser.add_argument("-r", "--remote", help="List of accession numbers to download and add to input folder", default=None)
    args = parser.parse_args()
    runall(args.input,args.outdb,args.multicpu,args.remote)
