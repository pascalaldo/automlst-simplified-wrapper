#!/usr/bin/env python
import os, json, argparse, setlog, parsegenomes, mash, shutil, re
import copyseqsql, makeseqsql, glob, seqsql2fa, subprocess
import multiprocessing as mp
import makehmmsql, getgenematrix, getgenes, concatmsa, ete3helper
import numpy as np
#from scipy.cluster.vq import kmeans2


###
##### Workflow 1: Build query + reference phylogeny from scratch
###

def catTrees(flist,outfile):
    with open(outfile,"w") as ofil:
        for fname in flist:
            if os.path.exists(fname):
                with open(fname,"r") as ifil:
                    for line in ifil:
                        ofil.write(line)
            else:
                log.warning("Tree file not found for MLST: %s"%fname)
    return outfile

def runmafft(finput,output,thread=1,maxit=300,localpair=False,options="",rename="",fast=False):
    if localpair:
        options += " --localpair"
    if rename:
        #Rename fasta headers by ignoring everything after "rename" character
        with open(finput,"r") as infil, open(finput+".renamed","w") as ofil:
            for line in infil:
                if line.startswith(">") and rename in line:
                    ofil.write(line.strip().split(rename)[0]+"\n")
                else:
                    ofil.write(line)
        finput = finput+".renamed"
    if fast:
        cmd = "mafft --quiet --thread %d --treeout %s %s"%(thread,options,finput)
    else:
        cmd = "mafft --quiet --thread %d --localpair --treeout --maxiterate %d %s %s"%(thread,maxit,options,finput)
    cmd = cmd.split()
    with open(os.devnull,"w") as devnull, open(output,"w") as outfil:
        try:
            subprocess.call(cmd, stdout=outfil, stderr=devnull)
            log.info("Finished alignment %s"%output)
        except Exception as e:
            log.error("Failed to run mafft %s (%s)"%(finput,e))

def runtrimal(finput,output,method="automated1"):
    cmd = "trimal -in %s -out %s -%s"%(finput,output,method)
    cmd = cmd.split()
    with open(os.devnull,"w") as devnull:
        try:
            subprocess.call(cmd, stdout=devnull, stderr=devnull)
            log.info("Finished trimming %s"%output)
        except Exception as e:
            log.error("Failed to run trimall %s (%s)"%(finput,e))

def processmlst(indir,aligndir,cpu=1,parallel=True,trim=False,fast=False):
    flist = glob.glob(os.path.join(indir,"*.fna"))
    if not os.path.exists(aligndir):
        os.makedirs(aligndir)

    if cpu > 1 and parallel:
        pool = mp.Pool(cpu)
        try:
            for filename in flist:
                fname = os.path.split(filename)[-1]
                if trim:
                    pool.apply_async(runtrimal, args=(filename, os.path.join(aligndir,fname)))
                else:
                    pool.apply_async(runmafft, (filename, os.path.join(aligndir,fname)),dict(fast=fast))
        finally:
            pool.close()
            pool.join()
    else:
        for filename in flist:
            fname = os.path.split(filename)[-1]
            if trim:
                runtrimal(filename, os.path.join(aligndir,fname))
            else:
                runmafft(filename, os.path.join(aligndir,fname),thread=cpu,fast=fast)

def hmmsearch(outname,hmmdb,infile,mcpu=1,cut="tc"):
    log.info("Searching sequences against hmmdb: %s"%hmmdb)
    cmd=["hmmsearch", "--domtblout", outname, "--noali", "--notextw", hmmdb, infile]
    if mcpu>1:
        cmd[1:1] = ["--cpu", str(mcpu)]
    if cut and cut.lower() in ("ga","tc","nc"):
        cmd[1:1] = ["--cut_%s"%cut.lower()]
    with open(os.devnull,"w") as devnull:
        subprocess.call(cmd, stdout=devnull, stderr=devnull)

def buildseqdb(selorgs,queryseqs,resultdir,sourcedb):
    #Get list of all selected and outgroups
    reflist = list(selorgs.get("selspecies",[]))
    reflist.extend(selorgs.get("seloutgroups",[]))

    #Copy from reference database to query database
    seqsdb = os.path.join(resultdir,"seqs.db")
    copyseqsql.copydb(queryseqs,sourcedb,seqsdb)

def getorgs(resultdir,mashresult,skip="",IGlimit=50,OGlimit=1,minorgs=25):
    #Get user selection if exists
    usersel = os.path.join(resultdir,"userlist.json")
    autosel = os.path.join(resultdir,"autoOrglist.json")
    selection = False
    if os.path.exists(usersel):
        with open(usersel,"r") as fil:
            selection = json.load(fil)
    elif os.path.exists(autosel):
        with open(autosel,"r") as fil:
            selection = json.load(fil)
    elif "skip2" in skip:
        # commonrank = mashresult.get("commonrank",["phyl","",""])
        # if commonrank[0]:
        #     selection = {"selspecies":[x.get("id","") for x in mashresult.get("reforgs",[])[:IGlimit] if str(x.get(commonrank[0]+"id","")) in str(commonrank[1])],
        #              "seloutgroups":[x.get("id","") for x in mashresult.get("outgroups",[])[:OGlimit]]}
        # else:
        qorgs = mashresult.get("queryorgs",[])
        selrefs = mashresult.get("reforgs",[])[:IGlimit-len(qorgs)-OGlimit]
        selection = {"selspecies":[x.get("id","") for x in selrefs]}
        orgrecs = {i:v for i,v in enumerate(qorgs+selrefs)}
        # log.debug("Issue_orgrecs: %s"%orgrecs)
        commonrank = mash.getcommongroup(orgrecs)
        # log.debug("Issue_commonrank: %s"%commonrank)
        OGorgs = mash.getoutgrouporgs(commonrank,mashresult.get("reforgs",[]),glimit=OGlimit)
        # log.debug("Issue_OGorgs: %s"%OGorgs)
        selection["seloutgroups"] = [x.get("id","") for x in OGorgs[:OGlimit]]
        # if not len(selection["selspecies"]) > minorgs:
        #     selection = {"selspecies":[x.get("id","") for x in mashresult.get("reforgs",[])[:minorgs]],
        #              "seloutgroups":[x.get("id","") for x in mashresult.get("outgroups",[])[:OGlimit]]}

        with open(autosel,"w") as fil:
            json.dump(selection,fil)

    return selection

def getmlstselection(resultdir,mlstpriority,maxmlst=100,minmlst=10,skip="",ignoreorgs=None,concat=True):
    if not ignoreorgs:
        ignoreorgs = []
    usersel = os.path.join(resultdir,"usergenes.json")
    autosel = os.path.join(resultdir,"autoMlstlist.json")
    selection = False
    delorgs = set()
    if os.path.exists(usersel):
        with open(usersel,"r") as fil:
            temp = json.load(fil)
            selection = temp["selection"]
            delorgs = temp["delorgs"]
            concat = temp.get("mode","concatenated")
            if concat != "concatenated":
                concat = False
            return selection, list(delorgs), concat

    if os.path.exists(autosel):
        with open(autosel,"r") as fil:
            temp = json.load(fil)
            selection = temp["selection"]
            delorgs = temp["delorgs"]
    else:
        minrecs = [x for x in mlstpriority if not len(x["orgdel"])]
        log.info("MLST genes in all orgs = %s"%len(minrecs))
        if len(minrecs) >= minmlst:
            recs = minrecs[:maxmlst]
        else:
            recs = [x for x in mlstpriority if not any([org in x["orgdel"] for org in ignoreorgs])][:minmlst]
        selection = [x["acc"] for x in recs]
        for x in recs:
            delorgs = delorgs | set(x["orgdel"])
        with open(autosel,"w") as fil:
            json.dump({"selection":selection,"delorgs":list(delorgs)},fil)

    return selection, list(delorgs), concat

def runIQtree(outdir,infasta,partfile="",cpu=1,model="GTR",bs=0,fout="speciestree",titlesep="",outgroup=""):
    #If title seperator exists rename titles in file first
    if titlesep:
        tf = infasta+".renamed"
        with open(tf,"w") as ofil, open(infasta,"r") as ifil:
            for line in ifil:
                lineout = line.strip()
                if line.startswith(">") and titlesep in line:
                    lineout = line.strip().split(titlesep)[0]
                ofil.write(lineout+"\n")
        infasta = tf
    log.info("Building tree: model=%s, bootstrap=%s, file=%s"%(model,bs,fout))
    cmd=["iqtree", "-quiet", "-nt", str(cpu), "-s", infasta, "-m", model,"-pre",os.path.join(outdir,fout)]
    if bs:
        cmd.extend(["-bb",str(bs)])
    if outgroup:
        #get full title of first outgroup to root tree
        with open(infasta,"r") as ifil:
            for line in ifil:
                if line.startswith(">") and outgroup in line:
                    cmd.extend(["-o", "%s"%line[1:].strip()])
                    break
    if partfile:
        cmd.extend(["-spp", partfile])
    with open(os.devnull,"w") as devnull:
        subprocess.call(cmd,stdout=devnull,stderr=devnull)

def runAstral(resultdir,intree,outtree,astjar=""):
    log.info("Started running Astral: %s"%(intree))
    if not astjar or not os.path.realpath(astjar):
        astjar = os.path.join(os.path.dirname(__file__),"astral.5.5.9.jar")
    cmd=["java", "-Xmx3000M", "-jar", str(astjar), "-i", intree, "-o", outtree]
    with open(os.path.join(resultdir,"astral.log"),"w") as astlog:
        subprocess.call(cmd,stdout=astlog,stderr=astlog)

def colphylogeny(resultdir,trimdir,cpu=1,model="MFP",bs=0):
    treedir = os.path.join(resultdir,"trees")
    if not os.path.exists(treedir):
        os.makedirs(treedir)
    flist = glob.glob(os.path.join(trimdir,"*.fna"))
    #build all trees
    if cpu > 1:
        pool = mp.Pool(cpu)
        try:
            for filename in flist:
                pool.apply_async(runIQtree, (treedir, filename), dict(model=model,bs=bs,fout=os.path.split(filename)[-1],titlesep="|"))
        finally:
            pool.close()
            pool.join()
    else:
        for filename in flist:
            runIQtree(treedir, filename, fout=os.path.split(filename)[-1], bs=bs, model=model, titlesep="|")

    #Combine all trees using ASTRAL
    flist = glob.glob(os.path.join(treedir,"*.treefile"))
    alltrees = catTrees(flist,os.path.join(treedir,"alltrees.tree"))
    coltree = os.path.join(treedir,"summaryTree.tree")
    runAstral(resultdir,alltrees,coltree)

def concatphylogeny(resultdir,concatfasta,partfile,cpu=1,model="GTR",bs=0,outgroup=""):
    trimdir = os.path.join(resultdir,"mlst_trimmed")
    log.info("Concatenating all MLST genes")
    concatmsa.concatmsa(os.path.join(trimdir,"*.fna"),concatfasta,partfile,fsplit=" ",checktype=False)
    #Build tree from concatenated fasta
    treedir = os.path.join(resultdir,"trees")
    if not os.path.exists(treedir):
        os.makedirs(treedir)
    runIQtree(treedir,concatfasta,partfile=partfile,cpu=cpu,fout="concatTree.tree",bs=bs,model=model)

def screenmlst(mlstdir,aligndir,cpu=1,mingenes=50):
    """Additional screen for outlier genes with highly discordant parsimony trees"""
    guidetrees = glob.glob(os.path.join(mlstdir,"*.tree"))
    #combine all to tree file
    hkgenes = []
    alltrees = os.path.join(aligndir,"guidetrees.tree")
    with open(alltrees,"w") as gtfil:
        for fname in guidetrees:
            with open(fname,"r") as tfil:
                hk = os.path.split(fname)[-1]
                hk = hk[:hk.find(".")]
                hkgenes.append(hk)
                gtfil.write("".join([x.strip() for x in tfil]).replace(";","")+";\n")
    #Run robinson foulds distance calculation
    cmd = ["iqtree","-nt",str(cpu),"-rf_all",alltrees]
    with open(os.devnull,"w") as devnull:
        subprocess.call(cmd,stdout=devnull,stderr=devnull)

    #Read distance matrix and store median values for each gene
    treedists = []
    with open(alltrees+".rfdist","r") as fil:
        line = fil.next()
        for line in fil:
            temp = np.median([int(x) for x in line.strip().split()[1:] if int(x)>0])
            treedists.append(temp)

    #segment tress into 2 partitions using k-means
    #km,kmgroup = kmeans2(np.array(treedists),2)
    #Ensure group with highest distance is marked as group 1, group 0 is highest group
    #if km[0] > km[1]:
    #    kmgroup = np.ones(len(kmgroup))-kmgroup
    #inds = [i for i,x in enumerate(kmgroup) if x>0]
    #get hkgenes with highest group designation
    # return kmgroup,hkgenes,treedists

    #Assume normal distribution and highlight trees outside 1 standard deviation (~ 15% of distant genes)
    threshold = np.mean(treedists)+np.std(treedists)
    inds = np.array([int(x>=threshold) for x in treedists])


    hkhigh = sorted(np.array(hkgenes)[inds>0], key=lambda x: treedists[hkgenes.index(x)])

    #if lowest group is under min genes remove genes from highest group
    mindiff = mingenes - (len(hkgenes)-len(hkhigh))
    if mindiff > 0:
        hkhigh = hkhigh[mindiff+1:]

    #Mask screened MLST genes
    for hk in hkhigh:
        fname = os.path.join(aligndir,hk+".fna")
        os.rename(fname,fname+".removed")

    #write table of genes and median distances
    fname = os.path.join(aligndir,"removedsumamry.tsv")
    with open(fname,"w") as fil:
        for i,hk in enumerate(hkgenes):
            fil.write("%s\t%s\t%s\n"%(hk,treedists[i],inds[i]))

    return hkhigh

def startwf1(indir,resultdir,checkpoint=False,concat=False,mashmxdist=0.5,cpu=1,
             skip="",refdb="",hmmdb="",rnadb="",maxmlst=100,model="GTR",bs=0,
             kf=False,maxorg=50,filtMLST=True,fast=False,minmlst=10):
    """WORKFLOW 1: Build phylogeny from scratch"""
    if not checkpoint:
        checkpoint = "w1-0"
    queryseqs = os.path.join(resultdir,"queryseqs")

    #Parse all inputs
    if checkpoint == "w1-0":
        log.info("JOB_STATUS::Parsing all genomes...")
        if parsegenomes.parseall(indir,queryseqs):
            checkpoint = "w1-1"
            log.info("JOB_CHECKPOINT::%s"%checkpoint)
            log.info("JOB_PROGRESS::5/100")
        else:
            log.error("Problem parsing input genomes")
            return False

    #Run MASH distances
    mashresult = False
    if checkpoint == "w1-1":
        log.info("JOB_STATUS::Running MASH ANI estimation against reference sequences...")
        mashresult = mash.getdistances(queryseqs,resultdir,cpu=cpu,maxdist=mashmxdist)
        if mashresult:
            checkpoint = "w1-2"
            log.info("JOB_CHECKPOINT::%s"%checkpoint)
            log.info("JOB_PROGRESS::10/100")
        else:
            log.error("MASH distance failed")
            return False

    #Get set of organisms to build seq DB
    selorgs = False
    if checkpoint == "w1-2":
        log.info("JOB_STATUS::Loading mash results...")
        if mashresult:
            pass
        elif os.path.exists(os.path.join(resultdir,"reflist.json")):
            with open(os.path.join(resultdir,"reflist.json"),"r") as fil:
                mashresult = json.load(fil)
                log.info("Loading mash results...")
        else:
            log.error("No Mash results to process")
            return False

        selorgs = getorgs(resultdir,mashresult,skip=skip,IGlimit=maxorg,minorgs=25)
        if selorgs:
            checkpoint = "w1-3"
            log.info("JOB_CHECKPOINT::%s"%checkpoint)
            log.info("JOB_PROGRESS::15/100")
        else:
            checkpoint = "w1-STEP2"
            log.info("JOB_CHECKPOINT::%s"%checkpoint)

    if checkpoint == "w1-STEP2":
        log.info("JOB_STATUS:: Waiting for selected organisms")
        return "waiting"

    #Copy reference sequence database and add query organisms
    orgdb = os.path.join(resultdir,"refquery.db")
    # selorgs = {}
    if not selorgs and os.path.exists(os.path.join(resultdir,"userlist.json")):
        with open(os.path.join(resultdir,"userlist.json"),"r") as fil:
            selorgs = json.load(fil)
    elif not selorgs and os.path.exists(os.path.join(resultdir,"autoOrglist.json")):
        with open(os.path.join(resultdir,"autoOrglist.json"),"r") as fil:
            selorgs = json.load(fil)
    if checkpoint == "w1-3":
        log.info("JOB_STATUS:: Collating selected genomes...")
        #Clear old db if exists
        if os.path.exists(orgdb):
            os.remove(orgdb)

        flist = list(selorgs["selspecies"])
        flist.extend(selorgs["seloutgroups"])
        flist = [x for x in flist if "query" not in x.lower()]
        log.info("refdb: %s"%refdb)
        if copyseqsql.copydb(flist,refdb,orgdb):
            #get file list of query orgs and add to db
            seqlist = glob.glob(os.path.join(queryseqs,"*.fna"))
            allorgs = makeseqsql.runlist(seqlist,orgdb,True,"",True)
            if allorgs:
                checkpoint = "w1-4"
                log.info("JOB_CHECKPOINT::%s"%checkpoint)
        else:
            log.error("Failed at copydb")

    #Run HMM searches on query organism
    if checkpoint == "w1-4":
        if not orgdb:
            log.error("No querydb rerun at checkpoint w1-3")
            return False
        #Write newly added seqences to file
        naseqs = os.path.join(resultdir,"addedseqs.fna")
        aaseqs = os.path.join(resultdir,"addedseqs.faa")
        seqsql2fa.writefasta(orgdb,aaseqs,False,"",True)
        seqsql2fa.writefasta(orgdb,naseqs,True,"",True)

        #Run HMM searches
        log.info("JOB_STATUS:: Searching for MLST genes in query sequences...")
        log.info("JOB_PROGRESS::25/100")
        if not hmmdb:
            hmmdb = os.path.join(os.path.dirname(os.path.realpath(__file__)),"reducedcore.hmm")
        if not rnadb:
            rnadb = os.path.join(os.path.dirname(os.path.realpath(__file__)),"barnap_bact_rRna.hmm")
        hmmsearch(aaseqs+".domhr",hmmdb,aaseqs,mcpu=cpu)
        hmmsearch(naseqs+".domhr",rnadb,naseqs,mcpu=cpu)

        #Add HMMresults to refdb and cleanup unused seqs
        log.info("Adding query HMM results to database")
        status = makehmmsql.run(aaseqs+".domhr",orgdb)
        if not status:
            #If no genes were found report error
            log.error("No MLST genes could be found. Stop processing")
        log.info("Adding query RNA results to database")
        status = makehmmsql.run(naseqs+".domhr",orgdb,rna=True)
        os.remove(aaseqs)
        os.remove(naseqs)
        checkpoint = "w1-5"
        log.info("JOB_CHECKPOINT::%s"%checkpoint)

    #Calculate mlst core genes and export
    mlstdir = os.path.join(resultdir,"mlstgenes")
    mlstpriority = os.path.join(resultdir,"mlstpriority.json")
    genematjson = os.path.join(resultdir,"mlstmatrix.json")
    dndsfile = os.path.join(os.path.dirname(os.path.realpath(__file__)),"dnds.json")
    mlstselection = []
    delorgs = []
    if checkpoint == "w1-5":
        allquery = [os.path.splitext(os.path.split(x)[-1])[0].replace(" ","_") for x in glob.glob(os.path.join(resultdir,"queryseqs","*.fna"))]
        ignoreorgs = list(selorgs.get("seloutgroups",[]))
        ignoreorgs.extend(allquery)

        #getmlstgenes.findsingles(orgdb,maxgenes=500,outdir=mlstdir)
        if os.path.exists(mlstpriority) and os.path.exists(genematjson):
            with open(genematjson,"r") as mfil, open(mlstpriority,"r") as pfil:
                mlstpriority = json.load(pfil)
                temp = json.load(mfil)
                genemat = temp["counts"]
                orgs = temp["orgs"]
                del temp
        else:
            genemat,orgs,mlstpriority = getgenematrix.getmat(orgdb,pct=0.5,pct2=1.0,bh=True,rna=True,savefil=genematjson,prifile=mlstpriority,dndsfile=dndsfile,ignoreorgs=allquery)
        mlstselection, delorgs, concat = getmlstselection(resultdir,mlstpriority,maxmlst,ignoreorgs=ignoreorgs,concat=concat,minmlst=minmlst)
        if "skip3" in skip.lower() or os.path.exists(os.path.join(resultdir,"usergenes.json")):
            #Export selected genes to mlst folder
            log.info("JOB_STATUS:: Writing MLST genes...")
            getgenes.writeallgenes(orgdb,mlstselection,delorgs,outdir=mlstdir,outgroups=selorgs.get("seloutgroups",None),pct=0.5,rename=True)
            checkpoint = "w1-6"
            log.info("JOB_CHECKPOINT::%s"%checkpoint)
        else:
            checkpoint = "w1-STEP3"
            log.info("JOB_CHECKPOINT::%s"%checkpoint)

    if checkpoint == "w1-STEP3":
        log.info("JOB_STATUS:: Waiting for selected MLST genes")
        return "waiting"

    ## Align and trim all MLST genes
    aligndir = os.path.join(resultdir,"mlst_aligned")
    trimdir = os.path.join(resultdir,"mlst_trimmed")

    if checkpoint == "w1-6":
        log.info("JOB_STATUS:: Aligning MLST genes")
        log.info("JOB_PROGRESS::30/100")
        #align all
        processmlst(mlstdir,aligndir,cpu=cpu,fast=fast)
        checkpoint = "w1-6b"
        log.info("JOB_CHECKPOINT::%s"%checkpoint)

    if checkpoint == "w1-6b":
        if filtMLST and filtMLST != "False":
            #Extra screen of MLST genes to remove outliers based on starting tree distance
            log.info("JOB_STATUS:: Screening for inconsistent MLST genes")
            log.info("JOB_PROGRESS::55/100")
            excludemlst = screenmlst(mlstdir,aligndir,cpu=cpu)
            log.info("JOB_STATUS:: Excluded MLST genes: %s"%excludemlst)

        #trim all
        log.info("JOB_STATUS:: Trimming alignments")
        log.info("JOB_PROGRESS::65/100")
        processmlst(aligndir,trimdir,cpu=cpu,trim=True)

        checkpoint = "w1-7"
        log.info("JOB_CHECKPOINT::%s"%checkpoint)

    treedir = os.path.join(resultdir,"trees")
    finishedtree = ""
    #Build trees
    if checkpoint == "w1-7":
        log.info("JOB_PROGRESS::75/100")
        if concat and concat != "False":
            log.info("JOB_STATUS:: Running concatenated supermatrix phylogeny")
            concatfasta = os.path.join(resultdir,"concatMLST.fasta")
            partfile = os.path.join(resultdir,"nucpartition.txt")
            concatphylogeny(resultdir, concatfasta, partfile,cpu=cpu,model=model,bs=bs)
            checkpoint = "w1-F"
            log.info("JOB_CHECKPOINT::%s"%checkpoint)
        else:
            log.info("JOB_STATUS:: Running coalescent tree phylogeny")
            colphylogeny(resultdir,trimdir,cpu=cpu,model=model,bs=bs)
            checkpoint = "w1-F"
            log.info("JOB_CHECKPOINT::%s"%checkpoint)

    if checkpoint == "w1-F":
        finaltree = os.path.join(resultdir,"final.tree")
        if os.path.exists(os.path.join(treedir,"concatTree.tree.treefile")) and concat:
            finishedtree = os.path.join(treedir,"concatTree.tree.treefile")
        elif os.path.exists(os.path.join(treedir,"summaryTree.tree")):
            finishedtree = os.path.join(treedir,"summaryTree.tree")
        #Copy final tree to root dir
        if finishedtree:
            fmat = 2 if bs or not concat else 5
            log.debug("Saving final tree... %s"%fmat)
            if not ete3helper.rerootTree(finishedtree,finaltree,fmat=fmat):
                shutil.copy(finishedtree,finaltree)
        else:
            log.error("Could not find final tree. Tree building failed")

        #Job finished? do some cleanup
        if not kf:
            if os.path.exists(orgdb):
                os.remove(orgdb)
            shutil.rmtree(os.path.join(resultdir,"queryseqs"))
            shutil.rmtree(os.path.join(resultdir,"mlst_trimmed"))
            shutil.rmtree(os.path.join(resultdir,"mlstgenes"))
            for oldfil in glob.glob(os.path.join(treedir,"*.model")):
                os.remove(oldfil)

def raxmlEPA(outdir,inalgn,reftree):
    #Run raxml EPA
    cmd="raxmlHPC-SSE3 -f v -m GTRGAMMA -p 12345 -w %s -t %s -s %s -n %s"%(outdir,reftree,inalgn,os.path.split(reftree)[-1])
    log.debug("Starting: %s"%cmd)
    try:
        with open(os.devnull,"w") as devnull:
            subprocess.call(cmd.split(),stdout=devnull)
        log.debug("RAxML-EPA: finished %s"%inalgn)
        return True
    except subprocess.CalledProcessError as e:
        log.error("RAxML-EPA: error, could not process %s - %s"%(inalgn,e))
        return False

def addalltrees(inlist,famrefdir,outdir,cpu=1):
    if cpu==1:
        for x in inlist:
            fname = os.path.split(x)[1]
            ftree = os.path.join(famrefdir,fname+".tree")
            raxmlEPA(outdir,x,ftree)
    else:
        pool = mp.Pool(cpu)
        try:
            for x in inlist:
                fname = os.path.split(x)[1]
                ftree = os.path.join(famrefdir,fname+".tree")
                pool.apply_async(raxmlEPA, (outdir,x,ftree))
        finally:
            pool.close()
            pool.join()


def addallalign(inlist,famrefdir,outdir,cpu=1,fast=False):
    if cpu==1:
        for x in inlist:
            fname = os.path.split(x)[1]
            addopt = "--add %s"%x
            finput = os.path.join(famrefdir,fname)
            output = os.path.join(outdir,fname)
            runmafft(finput,output,options=addopt,fast=fast)
    else:
        pool = mp.Pool(cpu)
        try:
            for x in inlist:
                fname = os.path.split(x)[1]
                addopt = "--add %s"%x
                finput = os.path.join(famrefdir,fname)
                output = os.path.join(outdir,fname)
                pool.apply_async(runmafft, (finput, output),dict(options=addopt,fast=fast))
        finally:
            pool.close()
            pool.join()

def getreference(mashresult,refdir):
    reffams = [str(x).lower() for x in os.listdir(refdir)]

    fams = list(set([q.get("familyname","N/A") for q in mashresult["queryorgs"]]))
    qfam = str(fams[0]).lower()
    if len(fams) == 1 and qfam in reffams:
        #get exact case as folder
        family = [fam for fam in os.listdir(refdir) if qfam == fam.lower()]
        return str(fams[0])
    elif len(fams) > 1:
        famlist = [(x,mashresult["queryorgs"][i].get("familyname","N/A")) for i,x in enumerate(mashresult["orglist"])]
        log.error("Multiple families detected in query geneomes: %s" % "\t".join( ["%s = %s |"%x for x in famlist]) )
    else:
        log.error("No matching family detected")

def startwf2(indir,resultdir,refdir="",checkpoint=False,reference="",model="GTR",bs=0,kf=False,maxmlst=100,
             cpu=1,mashmxdist=0.5,rnadb="",hmmdb="",fast=False):
    """WORKFLOW 2: Get all query genomes and identify reference tree to add sequences to"""
    if not checkpoint:
        checkpoint = "w2-0"
    if reference:
        checkpoint = "w2-3"
    queryseqs = os.path.join(resultdir,"queryseqs")

    #Parse all inputs
    if checkpoint == "w2-0":
        log.info("JOB_STATUS::Parsing all genomes...")
        if parsegenomes.parseall(indir,queryseqs):
            checkpoint = "w2-1"
            log.info("JOB_CHECKPOINT::%s"%checkpoint)
            log.info("JOB_PROGRESS::5/100")
        else:
            log.error("Problem parsing input genomes")
            return False

    #Run MASH distances
    mashresult = False
    if checkpoint == "w2-1":
        log.info("JOB_STATUS::Running MASH ANI estimation against reference sequences...")
        mashresult = mash.getdistances(queryseqs,resultdir,cpu=cpu,maxdist=mashmxdist)
        if mashresult:
            checkpoint = "w2-2"
            log.info("JOB_CHECKPOINT::%s"%checkpoint)
            log.info("JOB_PROGRESS::10/100")
        else:
            log.error("MASH distance failed")
            return False

    #Detect reference group
    if checkpoint == "w2-2":
        log.info("JOB_STATUS::Loading mash results...")
        if mashresult:
            pass
        elif os.path.exists(os.path.join(resultdir,"reflist.json")):
            with open(os.path.join(resultdir,"reflist.json"),"r") as fil:
                mashresult = json.load(fil)
                log.info("Loading mash results...")
        else:
            log.error("No Mash results to process")
            return False

        reference = getreference(mashresult,refdir)
        if reference:
            checkpoint = "w2-3"
            log.info("JOB_STATUS::Detected reference = %s"%reference)
            log.info("JOB_CHECKPOINT::%s"%checkpoint)
            log.info("JOB_PROGRESS::15/100")
        else:
            checkpoint = "w2-F"
            log.info("JOB_CHECKPOINT::%s"%checkpoint)
            log.error("No matching reference found for all query organisms")

    orgdb = os.path.join(resultdir,"allseqs.db")
    seqlist = glob.glob(os.path.join(queryseqs,"*.fna"))
    #Add to db and run HMM searches. Use model in reference folder, or use global models
    if checkpoint == "w2-3":
        #get file list of query orgs and add to db
        allorgs = makeseqsql.runlist(seqlist,orgdb,True,"",False)

        #Write newly added seqences to file
        naseqs = os.path.join(resultdir,"addedseqs.fna")
        aaseqs = os.path.join(resultdir,"addedseqs.faa")
        seqsql2fa.writefasta(orgdb,aaseqs,False,"",True)
        seqsql2fa.writefasta(orgdb,naseqs,True,"",True)

        #Run HMM searches
        log.info("JOB_STATUS:: Searching for MLST genes in query sequences...")
        log.info("JOB_PROGRESS::25/100")

        if not hmmdb:
            if os.path.exists(os.path.join(refdir,reference,"core.hmm")):
                hmmdb = os.path.join(refdir,reference,"core.hmm")
            else:
                hmmdb = os.path.join(os.path.dirname(os.path.realpath(__file__)),"reducedcore.hmm")
        if not rnadb:
            rnadb = os.path.join(os.path.dirname(os.path.realpath(__file__)),"barnap_bact_rRna.hmm")

        hmmsearch(aaseqs+".domhr",hmmdb,aaseqs,mcpu=cpu)
        hmmsearch(naseqs+".domhr",rnadb,naseqs,mcpu=cpu)

        #Add HMMresults
        log.info("Adding query HMM results to database")
        status = makehmmsql.run(aaseqs+".domhr",orgdb)
        if not status:
            #If no genes were found report error
            log.error("No MLST genes could be found. Stop processing")
        log.info("Adding query RNA results to database")
        status = makehmmsql.run(naseqs+".domhr",orgdb,rna=True)
        os.remove(aaseqs)
        os.remove(naseqs)
        checkpoint = "w2-4"
        log.info("JOB_CHECKPOINT::%s"%checkpoint)

    mlstdir = os.path.realpath(os.path.join(resultdir,"mlstgenes"))
    if not os.path.exists(mlstdir):
        os.makedirs(mlstdir)
    aligndir = os.path.realpath(os.path.join(resultdir,"mlst_aligned"))
    if not os.path.exists(aligndir):
        os.makedirs(aligndir)
    trimdir = os.path.realpath(os.path.join(resultdir,"mlst_trimmed"))
    if not os.path.exists(trimdir):
        os.makedirs(trimdir)

    # extract single copy results from genelist
    if checkpoint == "w2-4":
        genelist = [os.path.splitext(os.path.split(x)[1])[0] for x in glob.glob(os.path.join(refdir,reference,"*.fna"))]
        log.info("JOB_STATUS:: Writing MLST genes...")
        singles = [str(x) for x in getgenematrix.getmat(orgdb,rna=True,bh=True)[2][0] if x in genelist]
        getgenes.writeallgenes(orgdb,singles,[],outdir=mlstdir,outgroups=None,pct=0.5,writeaa=False,rename=True)

        log.info("JOB_STATUS:: Adding MLST genes to alignments...")
        addallalign([os.path.join(mlstdir,x+".fna") for x in singles], os.path.join(refdir,reference), aligndir, cpu=cpu,fast=fast)

        log.info("JOB_STATUS:: Trimming all alignments...")
        processmlst(aligndir,trimdir,cpu=cpu,trim=True)
        checkpoint = "w2-5"
        log.info("JOB_CHECKPOINT::%s"%checkpoint)

    treedir = os.path.realpath(os.path.join(resultdir,"trees"))
    if not os.path.exists(treedir):
        os.makedirs(treedir)

    #Add to prebuilt trees
    if checkpoint == "w2-5":
        log.info("JOB_PROGRESS::55/100")
        log.info("JOB_STATUS:: Adding to reference gene trees...")
        inlist = [os.path.join(trimdir,x) for x in os.listdir(trimdir)]
        addalltrees(inlist,os.path.join(refdir,reference),treedir,cpu=cpu)
        checkpoint = "w2-6"
        log.info("JOB_CHECKPOINT::%s" % checkpoint)

    if checkpoint == "w2-6":
        log.info("JOB_PROGRESS::85/100")
        log.info("JOB_STATUS:: Running coalescent tree phylogeny")
        #Combine all trees using ASTRAL
        flist = glob.glob(os.path.join(treedir,"RAxML_labelledTree.*.tree"))
        #reformat newick trees for astral
        for pbt in flist:
            with open(pbt,"r") as ifil, open(pbt+".newick","w") as ofil:
                x = ifil.next()
                ofil.write(re.sub("\[I\d+?\]|\"|'|QUERY___","",x))
        flist = glob.glob(os.path.join(treedir,"*.newick"))

        alltrees = catTrees(flist,os.path.join(treedir,"alltrees.tree"))
        coltree = os.path.join(treedir,"summaryTree.tree")
        runAstral(resultdir,alltrees,coltree)

        checkpoint = "w2-F"
        log.info("JOB_CHECKPOINT::%s" % checkpoint)

    finishedtree=""
    if checkpoint == "w2-F":
        finaltree = os.path.join(resultdir,"final.tree")
        if os.path.exists(os.path.join(treedir,"summaryTree.tree")):
            finishedtree = os.path.join(treedir,"summaryTree.tree")
        #Copy final tree to root dir
        if finishedtree:
            fmat = 2
            log.debug("Saving final tree... %s"%fmat)
            if not ete3helper.rerootTree(finishedtree,finaltree,fmat=fmat):
                shutil.copy(finishedtree,finaltree)
        else:
            log.error("Could not find final tree. Tree building failed")


def startjob(indir,resultdir,skip="",checkpoint=False,workflow=1,refdb="",cpu=1,concat=False,
             model="GTR",bs=0,kf=False,maxmlst=100,maxorg=50,filtMLST=True,refdir="",fast=False,minmlst=10):
    #Setup working directory
    if not os.path.exists(os.path.join(os.path.realpath(resultdir),"queryseqs")):
        os.makedirs(os.path.join(os.path.realpath(resultdir),"queryseqs")) #query sequence folder

    #Read checkpoint from log if exists or use checkpoint parameter
    if not checkpoint and os.path.exists(os.path.join(resultdir,"automlst.log")):
        with open(os.path.join(resultdir,"automlst.log"),"r") as fil:
            for line in fil:
                x = line.strip().split("::")
                if "JOB_CHECKPOINT" in x[0]:
                    checkpoint = x[1]
                elif "WORKFLOW" in x[0]:
                    workflow = int(x[1])
    #Start log
    global log
    log = setlog.init(os.path.join(resultdir,"automlst.log"),toconsole=True)

    #ensure bootsrtap value is valid and no more than 1000
    try:
        bs = int(bs)
        if bs > 1000:
            bs = 1000
        elif bs < 0:
            bs = 0
    except ValueError:
        log.warning("Invalid bootsrap value found, setting to 0")
        bs = 0

    log.debug("START / RESUME JOB last checkpoint: %s"%checkpoint)
    log.info('JOB_PARAMS::{"resultdir":"%s","skip":"%s","workflow":%s,"concat":"%s","model":"%s"}'%(resultdir,skip,workflow,concat,model))
    if workflow == 1:
        log.info("WORKFLOW::1")
        try:
            return startwf1(indir,resultdir,checkpoint=checkpoint,skip=skip,refdb=refdb,cpu=cpu,concat=concat,
                            model=model,bs=bs,kf=kf,maxmlst=maxmlst,maxorg=maxorg,filtMLST=filtMLST,fast=fast,minmlst=minmlst)
        except Exception as e:
            log.error("Unexpected failure: %s"%e)
            raise
    elif workflow == 2:
        log.info("WORKFLOW::2")
        try:
            return startwf2(indir,resultdir,refdir=refdir,checkpoint=checkpoint,model=model,bs=bs,kf=kf,maxmlst=maxmlst,cpu=cpu,fast=fast)
        except Exception as e:
            log.error("Unexpected failure: %s"%e)
            raise
    else:
        log.error("Improper workflow specified: %s"%workflow)
        return False

# Commandline Execution
if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="""Start from input directory of genomes, get all reference distances and build phylogeny""")
    parser.add_argument("indir", help="Directory of input genomes or list of files")
    parser.add_argument("resultdir", help="Directory to store all results")
    parser.add_argument("-skip","--skip", help="Flag to skip manual intervention of organism selection. Ex: 'skip2.skip3' skips organism and mlst selection (default for cmdline execution)",default="skip2.skip3")
    parser.add_argument("-ref","--refdb", help="Reference database of orgs",default="")
    parser.add_argument("-rd","--refdirectory", help="Reference directory of pre-built trees and alignments",default="")
    parser.add_argument("-bs","--bootstrap", help="Bootstrap replicates (default: None)",default="")
    parser.add_argument("-cat","--concat", help="Use concatenated supermatrix to build tree",action="store_true",default=False)
    parser.add_argument("-kf","--keepfiles", help="Do not delete intermediate files after completion",action="store_true",default=False)
    parser.add_argument("-fm","--filtmlst", help="Filter mlst genes by starting tree congruance",action="store_true",default=False)
    parser.add_argument("-fa","--fastalign", help="Align using FFt-2 mode (faster)",action="store_true",default=False)
    parser.add_argument("-m","--model", help="Use specific model for iqtree parameter (default: GTR)",default="GTR")
    parser.add_argument("-ch","--checkpoint", help="Explicitly start at checkpoint",default=False)
    parser.add_argument("-c","--cpu", help="Number of cpu cores to use",type=int, default=1)
    parser.add_argument("-wf","--workflow", help="Workflow to use 1 or 2",choices=(1,2),type=int, default=1)
    parser.add_argument("-mm","--maxmlst", help="Maximum number of MLST genes to use (default:100)",type=int, default=100)
    parser.add_argument("-mn","--minmlst", help="Minimum number of MLST genes to use (default:10)",type=int, default=10)
    parser.add_argument("-mo","--maxorg", help="Maximum number of organisms to use (default:50)",type=int, default=50)
    args = parser.parse_args()
    startjob(args.indir,args.resultdir,args.skip,refdb=args.refdb,cpu=args.cpu,
             concat=args.concat,model=args.model,checkpoint=args.checkpoint,
             bs=args.bootstrap,kf=args.keepfiles,maxmlst=args.maxmlst,
             filtMLST=args.filtmlst,refdir=args.refdirectory,workflow=args.workflow,
             fast=args.fastalign,minmlst=args.minmlst)