import os, tempfile, models, json, datetime, csv, zipfile, time, shutil, re
from flask import render_template, jsonify, request, redirect, abort, make_response, send_from_directory, flash, g
from flask_mail import Message
from redis import Redis
from redis import ConnectionError as redisConnectError
from app import app
from app import mail
from werkzeug.utils import secure_filename
from Bio import Entrez

def getdb():
    rddb = getattr(g,"_redisdb",False)
    if not rddb:
        rddb = g._redisdb = Redis.from_url('redis://localhost:6379/0')
    try:
        rddb.ping()
    except redisConnectError:
        rddb = False
    return rddb

def getserverstats():
    rddb = getdb()
    status={}
    if rddb:
        status["running"] = rddb.llen("AMLSTRQ")
        status["waiting"] = rddb.llen("AMLSTSQ")
        status["finished"] = rddb.llen("AMLSTDQ")
        status["paused"] = rddb.llen("AMLSTPQ")
    return status

def validatefile(fname, asfil=False):
    validgbkext = ['gbk','genbank','gbff','gb','embl']
    validfakext = ['fasta','fa','fna','faa','fas']
    ext = os.path.splitext(fname)[1]
    if not asfil and ext[1:].lower() in validgbkext+validfakext:
        return True
    elif asfil and ext[1:].lower() in validgbkext:
        return True
    return False

def getNCBIgbk(acc):
    try:
        #crude check for format
        if acc.upper().startswith("GCA_") or acc.upper().startswith("GCF_"):
            flash("Error, cannot use assembly accession number, please use Genbank or Refseq sequence accession")
            filename = False
        elif bool(re.search('[A-Z]{4}0+(\.\d){0,}$', acc)):
            flash("Error, WGS accessions not supported, please use Genbank or Refseq sequence accession")
            filename = False
        elif acc.replace("_","").replace(".","").isalnum() and len(acc) <= 20 and len(acc) >= 6:
            if "." in acc:
                acc = acc[:acc.index(".")]
            Entrez.email = "artsuser@ziemertlab.com"
            handle = Entrez.efetch(db="nucleotide", rettype="gbwithparts", id=acc, retmode="text")
            filename = os.path.join(tempfile.mkdtemp(dir=app.config['UPLOAD_FOLDER']), secure_filename(acc+".gbk"))
            with open(filename,"w") as outfile:
                outfile.write(handle.read())
        else:
            flash("Error with Accession number ")
            filename = False
    except Exception, e:
        flash("Error retrieving gbk from NCBI")
        filename = False
    return filename

def getinfile():
    ufile = False
    filename = []
    # USe empty string for FALSE due to storage as string in redis
#    if 'asjobid' in request.form and request.form['asjobid']:
#        if getASstatus(request.form['asjobid']):
#            filename = getASgbk(request.form['asjobid'])
#            if not filename:
#                return False,""
#        else:
#            return False,""
    if 'ncbiacc1' in request.form and request.form['ncbiacc1'] and request.form.get("filesrc") == 'ncbi':
         filename = [getNCBIgbk(request.form['ncbiacc1'])]
         if not filename:
             return []
#    elif 'asseqfile' in request.files and validatefile(request.files['asseqfile'].filename,True):
#        ufile = request.files['asseqfile']
    elif 'seqfile1' in request.files and request.form.get("filesrc") == 'seqfile':
        for seqfile in request.files.getlist('seqfile1'):
            if validatefile(seqfile.filename):
                ufile = seqfile
            if ufile:
                tmp = os.path.join(tempfile.mkdtemp(dir=app.config['UPLOAD_FOLDER']), secure_filename(ufile.filename))
                ufile.save(tmp)
                filename.append(tmp)
    return filename

def isjob(jobid):
    rddb = getdb()
    if rddb and rddb.keys("automlstjob:%s"%jobid):
        return True
    return False

def readdjob(jobid):
    rddb = getdb()
    if rddb:
        if rddb.keys("automlstjob:%s"%jobid) and os.path.exists(os.path.join(app.config['RESULTS_FOLDER'],jobid)):
            rddb.lrem("AMLSTWQ","%s"%jobid,1)
            rddb.lpush("AMLSTSQ","%s"%jobid,1)

def addjob(**kwargs):
    rddb = getdb()
    automlstjob = models.automlstjob(**kwargs)
    if rddb:
        rddb.hmset("automlstjob:%s"%automlstjob.id, automlstjob.getdict())
        rddb.rpush("AMLSTSQ",automlstjob.id)
        # make a results directory?
    return automlstjob

def updatejob(jobid,newref):
    rddb = getdb()
    if rddb:
        redisid = "automlstjob:"+jobid
        rddb.hset(redisid,"reference",newref)

def getjobstatus(jobid):
    jobstatus = "Waiting in queue"
    mashstatus = ""
    checkpoint = ""
    percent = 0
    workflow = 0
    errorlist = []
    paramdict = {}
    if os.path.exists(os.path.join(app.config['RESULTS_FOLDER'],jobid,'automlst.log')):
        with open(os.path.join(app.config['RESULTS_FOLDER'],jobid,'automlst.log'), 'r') as infile:
            #Get align progress
            mlstfound = 0
            mlstaligned = 0

            for line in infile:
                if 'JOB_STATUS' in line:
                    statlist = line.strip().split('::')
                    jobstatus = statlist[1]
                    if "Aligning MLST genes" in jobstatus:
                        mlstaligned = 0
                elif 'JOB_PROGRESS' in line:
                    proglist = line.strip().split('::')
                    fraction = str(proglist[1]).split('/')
                    percent = 100 * (float(fraction[0]) / float(fraction[1]))
                elif 'MASH_STATUS' in line:
                    mashlist = line.strip().split('::')
                    mashstatus = mashlist[1]
                elif 'JOB_CHECKPOINT' in line:
                    checklist = line.strip().split('::')
                    checkpoint = checklist[1]
                elif 'WORKFLOW' in line:
                    workflowline = line.strip().split('::')
                    workflow= workflowline[1]
                elif 'JOB_PARAMS' in line:
                    paramlist = line.strip().split('::')
                    paramdict = json.loads(paramlist[1])
                elif 'Writing genes:' in line:
                    mlstfound = line[line.find("Writing genes:"):].count(",")+1
                elif "Finished alignment" in line:
                    mlstaligned += 1
                elif "ERROR" in line:
                    errormsgline = line.strip().split(' - ')
                    currerror = errormsgline[3]
                    errorlist.append(currerror)
                elif "JOB_REANALYZE" in line:
                    errorlist = []
            if mlstfound > 0 and mlstaligned > 0 and mlstfound != mlstaligned:
                jobstatus += " (%s/%s complete)"%(mlstaligned,mlstfound)
                # total percent is pct * portion extra
                pct = float(mlstaligned)/mlstfound
                percent = percent + (pct * 25)

    jobstatdict = {"progress": percent,"status":jobstatus, "mash":mashstatus, "checkpoint": checkpoint, "workflow": workflow, "params":paramdict, "errors":errorlist}
    return jobstatdict

def reanalyzejob(jobid):
    paramdict={}
    rddb = getdb()
    with open(os.path.join(app.config['RESULTS_FOLDER'],jobid,'automlst.log'),'r') as jobread:
        for line in jobread:
            if 'JOB_PARAMS' in line:
                paramlist = line.strip().split('::')
                paramdict = json.loads(paramlist[1])
                #Clear Skip options
                paramdict["skip"] = ""
                if rddb:
                    rddb.hmset("automlstjob:%s"%jobid,{"skip":""})
    paramdict["skip"]=[]
    with open(os.path.join(app.config['RESULTS_FOLDER'],jobid,'automlst.log'),'a') as joblog:
        joblog.write("\n"+str(datetime.datetime.now())+" - INFO - JOB_REANALYZE::true \n")
        joblog.write(str(datetime.datetime.now())+" - INFO - JOB_CHECKPOINT::W1-STEP2 \n")
        joblog.write(str(datetime.datetime.now())+"- INFO - JOB_STATUS::Reanalyzing\n")
        joblog.write(str(datetime.datetime.now())+"- INFO - JOB_PROGRESS::15/100\n")
        joblog.write(str(datetime.datetime.now())+" - INFO - JOB_PARAMS::"+json.dumps(paramdict)+"\n")
    #Reset to Step2 state, remove all files except keepfiles
    keepfiles = ["automlst.log","mash_distances.txt","mash_distances.json","queryseqs","reflist.json"]
    jobdir = os.path.join(app.config['RESULTS_FOLDER'],jobid)

    # Store the old files instead of deleting them
    oldfiles = tempfile.mkdtemp(prefix="old_",dir=jobdir)
    keepfiles.append(os.path.split(oldfiles)[1])
    for fname in os.listdir(jobdir):
        if fname not in keepfiles:
            os.rename(os.path.join(jobdir,fname),os.path.join(oldfiles,fname))

    # Delete the old files
    shutil.rmtree(oldfiles)

# converts organism list .json to tab-separated text file, filters for only selected species and outgroups (dependent on presence of autoOrglist.json!)
def jsontotsv(jsonpath,jobid):
    resultdict = {}
    treeseqs = []
    with open(os.path.join(app.config['RESULTS_FOLDER'],jobid, 'autoOrglist.json'),'r') as jsonorglist:
        jsonorgdict = json.load(jsonorglist)
        treeseqs = jsonorgdict.get("selspecies")
        treeseqs += jsonorgdict.get("seloutgroups")
    with open(jsonpath, 'r') as jsondict:
        fulldict = json.load(jsondict)
        orglist = fulldict["orglist"]
        refdict = fulldict["reforgs"]
        for i in range(0,len(refdict)):
            if refdict[i]["id"] in treeseqs:
                currorg = refdict[i]
                distvals=currorg["dlist"]
                pvals=currorg["plist"]
                for j in range (0, len(orglist)):
                    currquery = orglist[j]
                    currd=""
                    currp=""
                    if not j >= len(pvals):
                        currp=pvals[j]
                    if not j >= len(distvals):
                        currd=distvals[j]
                    orgcomp = currorg["orgname"] + currquery
                    if orgcomp not in resultdict:
                        resultdict[orgcomp] = ""
                    resultdict[orgcomp] = {"orgname":currorg["orgname"], "strain":currorg["strain"], "id":currorg["id"], "queryname":currquery,"pval":currp, "dist":currd}
    with open(os.path.join(app.config['RESULTS_FOLDER'],jobid,'reftext.txt'),'wb') as csvfile:
        csvwriter = csv.writer(csvfile,delimiter="\t")
        for resultval in resultdict.values():
            csvwriter.writerow([resultval["orgname"],resultval["strain"],resultval["id"],resultval["queryname"],resultval["pval"],resultval["dist"]])

# converts .txt file of Mash distances to JSON
def tsvtojson(tsvpath):
    jsondict = {}
    linelist = []
    with open(tsvpath, 'r') as tsvsource:
        header = tsvsource.next()
        for line in tsvsource:
            linedata = line.strip().split("\t")
            linelist.append(linedata)
    jsondict["data"] = linelist
    return jsondict


def sendnotifymail(msg="",jobid="",to=""):
    try:
        if not msg:
            msg = "Hello, your autoMLST job has been submitted! Your job id is: "
        assert to, jobid
        msgobj = Message("Your autoMLST Job (%s) has been submitted"%jobid,recipients=[to])
        msgobj.html = "%s %s <br> <a href='%sresults/%s'>%sresults/%s</a>"%(msg,jobid,request.url_root,jobid,request.url_root,jobid)
        mail.send(msgobj)
    except Exception as e:
        print "Warning: Email not sent, check email configuration"
        print e

def findjobinfo(jobid):
    jobtitle = [jobid,jobid]
    if os.path.exists(os.path.join(app.config['RESULTS_FOLDER'],jobid,'jobtitle.txt')):
        with open(os.path.join(app.config['RESULTS_FOLDER'],jobid,'jobtitle.txt'),'r') as namefile:
            jobtitle[1] = namefile.next().strip()
    return jobtitle

def getlastresults():
    lastresult = request.cookies.get('automlst.lastresult', False)
    if lastresult:
        lastresult = lastresult.split(";")
        lastresult.reverse()
        lastresult = [findjobinfo(x) for x in lastresult]
    return lastresult

def zipalignments(jobid):
    jobpath = os.path.join(app.config['RESULTS_FOLDER'],jobid)
    alignpath = os.path.join(jobpath,'mlst_aligned')
    with zipfile.ZipFile(os.path.join(jobpath,jobid+'_alignments.zip'),'w') as alignzip:
        for alignfile in os.listdir(alignpath):
            if os.path.splitext(alignfile)[1] == '.fna':
                alignzip.write(os.path.join(alignpath,alignfile), jobid+'_alignments/'+alignfile)
        if os.path.exists(os.path.join(jobpath, 'concatMLST.fasta')):
            alignzip.write(os.path.join(jobpath, 'concatMLST.fasta'), jobid+'_alignments/concatMLST.fasta')

# identifies MLST genes used, extracts information on these genes from mlstpriority.json, outputs tab-separated .txt file
def mlsttsv(jobid):
    genesused = []
    prioritydict = {}
    jobpath = os.path.join(app.config['RESULTS_FOLDER'], jobid)
    alignpath = os.path.join(jobpath, 'mlst_aligned')
    for alignfile in os.listdir(alignpath):
        if os.path.splitext(alignfile)[1] == '.fna':
            genesused.append(os.path.splitext(alignfile)[0])
    with open(os.path.join(jobpath,'mlstpriority.json'),'r') as jsondict:
        fulldict = json.load(jsondict)
        for x in fulldict:
            if x['acc'] in genesused:
                prioritydict[x['acc']] = {"name":x['name'], "function": x['func'], "description": x['desc']}
    with open(os.path.join(jobpath,'mlst_genes.txt'), 'wb') as tsvfile:
        csvwriter = csv.writer(tsvfile, delimiter="\t")
        for key, value in prioritydict.items():
            csvwriter.writerow([key, prioritydict[key]['name'], prioritydict[key]['function'], prioritydict[key]['description']])
