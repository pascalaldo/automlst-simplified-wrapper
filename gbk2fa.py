#!/usr/bin/python
import sys
import os
import multiprocessing as mp
from Bio import SeqIO
import argparse
import setlog

global log
log=setlog.init(toconsole=True)

# def addremotecontigs(reclist):
#     """If Contig records exists try to download from NCBI and add to records"""
#     newlist = []
#     for rec in reclist:
#         if 'contig' in rec.annotations:
#             log.info("Genbank file contains supercontig references attempting to download")
#
#         else:
#             newlist.append(rec)

def getheader(seq_record,recnum,userecnum=False):
    """Get record header information"""
    if len(seq_record.id) and not userecnum:
        seqtitle = seq_record.id
    elif len(seq_record.name) and not userecnum:
        seqtitle = seq_record.name
    else:
        log.warning("record has no ID using reccord number")
        seqtitle = "scaffold_"+str(recnum)
    if len(seq_record.description) > 0:
        seqdesc = seq_record.description
    else:
        seqdesc = ""
        log.warning("Genbank does not have description")
    if "source" in seq_record.annotations.keys():
        seqdesc += "|source|"+seq_record.annotations["source"].replace(".","").replace(" ","_")+"|"
    if "organism" in seq_record.annotations.keys():
        seqdesc += "|org|"+seq_record.annotations["organism"].replace(".","").replace(" ","_").replace(",","")+"|"
    return seqtitle, seqdesc


def appendheader(seqtitle, desc, lnum, qual, loc):
    """Add details to fasta header"""
    gi = "lcl|" + str(lnum)
    gdesc = ""
    if "db_xref" in qual.keys():
        gi = qual["db_xref"][0].replace("|", "-").lower()
        if ":" in gi:
            gi = gi.replace(":", "|")
    if "product" in qual.keys():
        gdesc = qual["product"][0]+"|"
    if "gene" in qual.keys():
        gdesc += qual["gene"][0]+"|"
    if "locus_tag" in qual.keys():
        gdesc += "locus_tag:"+qual["locus_tag"][0]+"|"
    return gi + "|" + seqtitle + " " + desc + "|" + gdesc + "loc|" + str(loc[0]) + " " + str(loc[1] + 1) + " " + str(loc[2]) + " # ID=1_" + str(lnum) + ";"


def convertgenes(filename, outdir="./", genes=False,f="",rename=False,usetrans=False,plasmid=False,userecnum=False,prfix=""):
    """Parse all gbk records and output nuc and prot sequences for each CDS in multi-fasta format"""
    log.info("Starting %s..."%filename)
    fpath, fname = os.path.split(filename)
    fname, ext = os.path.splitext(fname)
    basename = os.path.join(outdir,prfix+fname)

    #Fix WGS genbank


    SeqIO.convert(filename, "genbank", basename + ".fa", "fasta")
    if genes:
        reclist = SeqIO.parse(filename, "genbank")
        if not os.path.isdir(outdir):
            os.makedirs(outdir)
        cutoff = 10
        lnum = 1
        cdscount=0
        recnum = 0
        features=["CDS"]
        if len(f.split(","))>0:
            features.extend(f.split(","))
        with open(basename + ".genelist", "w") as gn_handle, \
                open(basename + ".faa", "w") as aa_handle, \
                open(basename + ".fna", "w") as nuc_handle:
            for seq_record in reclist:
                recnum += 1
                seqtitle, seqdesc = getheader(seq_record,recnum,userecnum)
                bpcount = 0
                bpend = 1
                accession = seq_record.annotations.get("accession",seq_record.id.replace(" ","").replace("#",""))
                feature_count = 0
                for seq_feature in seq_record.features:
                    plasmidTitle=""
                    if "source" in seq_feature.type.lower() and "plasmid" in seq_feature.qualifiers.keys():
                        plasmidTitle="_PLASMID_"+str(seq_feature.qualifiers["plasmid"])
                    if seq_feature.type in features:
                        outseq = seq_feature.extract(seq_record.seq)
                        transeq=""
                        if "translation" in seq_feature.qualifiers.keys() and not usetrans:
                            transeq = seq_feature.qualifiers["translation"][0]
                        elif seq_feature.type=="CDS":
                            transeq = outseq.translate(to_stop=True)
                        if len(outseq) > cutoff:
                            bpend = int(seq_feature.location.end)
                            bpstart = int(seq_feature.location.start)
                            bpstrand = int(seq_feature.location.strand)
                            bpcount += bpend - bpstart
                            seqdetails = appendheader(seqtitle, seqdesc, lnum, seq_feature.qualifiers,[bpstart, bpend, bpstrand])
                            #seqdetails = "%s_%s # %s # %s # %s # ID=%s_%s;"%(accession,cdscount,bpstart,bpend,bpstrand,recnum,feature_count)
                            seqdetails+=plasmidTitle
                            seqdetails+="|Type="+seq_feature.type
                            seqdetails+="|ACC="+accession
                            cols = seqdetails.split("|")
                            nuc_handle.write(">%s\n%s\n" % (seqdetails, outseq))
                            aa_handle.write(">%s\n%s\n" % (seqdetails, transeq))
                            gn_handle.write("%s | %s\t%s\t%s\t%s\n" % (cols[0],cols[1], str(bpstart + 1), str(bpend), str(bpstrand)))
                            lnum += 1
                        if seq_feature.type == "CDS":
                            cdscount += 1
                            feature_count += 1
                log.info("Record #%s CDS bp coverage: %s%%"%(recnum,(bpcount * 100 / bpend)))
                validchars = "abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ1234567890_-"
                if "organism" in seq_record.annotations:
                    orgname= "".join([c for c in seq_record.annotations["organism"].replace(" ","_") if c in validchars])
                # elif "source" in seq_record.annotations:
                #     orgname= "".join([c for c in seq_record.annotations["source"].replace(" ","_") if c in validchars])
                else:
                    orgname=""
        if rename and len(orgname)>0:
            os.rename(basename + ".genelist", os.path.join(outdir, prfix + orgname + ".genelist"))
            os.rename(basename + ".faa", os.path.join(outdir, prfix + orgname + ".faa"))
            os.rename(basename + ".fna", os.path.join(outdir, prfix + orgname + ".fna"))
            os.rename(basename + ".fa", os.path.join(outdir, prfix + orgname + ".fa"))
            outfile=orgname
        else:
            outfile=fname
        log.info("Records: %s CDS features: %s"%(recnum,cdscount))
        log.info("### Finished %s ###"%filename)
        return outfile

def runall(finput="./", outdir="./", genes=False, cpu=1, features="",rename=False,filelist=False,usetrans=False,plasmid=False):
    """Run list of gbk files, use multiprocessing if applicible"""
    ext = (".gbk", ".gb", ".gbf", ".gbff", ".genbank", ".gpff")
    if filelist:
        temp = [line.strip() for line in sys.stdin]
        if len(temp)>0:
            finput=temp
    flist = []
    if type(finput) is list or type(finput) is tuple:
        for fitem in finput:
            if type(fitem) is str and os.path.exists(fitem) and fitem.lower().endswith(ext):
                flist.append(fitem)
            else:
                log.error(fitem + " is invalid skipping")
    elif type(finput) is str and os.path.isdir(finput):
        if not finput.endswith("/"):
            finput += "/"
        flist = [finput + x for x in os.listdir(finput) if x.endswith(ext)]
    elif type(finput) is str and finput.lower().endswith(ext):
        flist = [finput]

    if flist:
        if outdir and not outdir.endswith("/"):
            outdir += "/"
        if not os.path.isdir(outdir):
            os.makedirs(outdir)
        if cpu > 1:
            pool = mp.Pool(cpu)
            for filename in flist:
                pool.apply_async(convertgenes, args=(filename, outdir, genes, features, rename, usetrans, plasmid))
            pool.close()
            pool.join()
        else:
            for filename in flist:
                fname=convertgenes(filename, outdir, genes, features, rename, usetrans, plasmid)
    else:
        log.error("Nothing to convert (see -h or --help for usage)")


if __name__ == '__main__':
    # Commandline args
    parser = argparse.ArgumentParser(
        description="Convert Genbank to Fasta, produce multi-fasta gene files from CDS records")
    parser.add_argument("input", nargs="?", help="Genbank file or directory for conversion (default: current dir)",
                        default="./")
    parser.add_argument("-g", "--genes", help="Output all CDS records as multifasta (.fna and .faa)",
                        action="store_true", default=False)
    parser.add_argument("-f", "--features", help="Include features ex: rRNA,mRNA,...",default="")
    parser.add_argument("-fl", "--filelist", help="Use list of files as input (separated by line)",action="store_true",default=False)
    parser.add_argument("-tr", "--trans", help="Use translation for aa instead of gbk value",action="store_true",default=False)
    parser.add_argument("-p", "--plasmid", help="Keep plasmid records (default: false)",action="store_true",default=False)
    parser.add_argument("-rn", "--rename", help="Rename output file to Organism name",action="store_true",default=False)
    parser.add_argument("-ur", "--usereccordnum", help="Rename each scaffold by gbk record number",action="store_true",default=False)
    parser.add_argument("-c", "--cpu", help="Cpu cores to use (default: max cores)", type=int, default=mp.cpu_count())
    parser.add_argument("-od", "--outdir", help="Store results in OUTDIR (default: currentdir/converted)",
                        default="converted")
    args = parser.parse_args()

    # START EXECUTION
    runall(args.input, args.outdir, args.genes, args.cpu, args.features, args.rename, args.filelist, args.trans, args.plasmid)
