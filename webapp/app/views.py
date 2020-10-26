
import os, routines, json, uuid, datetime, threading
from flask import render_template, jsonify, request, redirect, abort, make_response, send_from_directory, flash, Response, url_for, copy_current_request_context
from Bio import Phylo
from app import app

@app.route('/')
@app.route('/index')
def index():
    return render_template("front.html")

@app.route('/help')
def help():
    return render_template("help.html")

@app.route('/about')
def about():
    return render_template("about.html")

@app.route('/results')
def results():
    lastresult = routines.getlastresults()
    return render_template("results.html",lastresults=lastresult)

@app.route('/download')
def amlstdownl():
    return render_template("download.html")

@app.route('/serverstatus')
def serverstatus():
    status = routines.getserverstats()
    return jsonify(status)

@app.route('/results/getreport', methods=['POST'])
def getreport():
    jobid = request.form.get("jobid")
    return redirect("/results/"+jobid+"/loading")

@app.route('/results/<jobid>')
@app.route('/results/<jobid>/')
def showresults(jobid):
    #return render_template("startjob.html",jobid=jobid)
    return redirect("/results/"+jobid+"/loading")

@app.route('/results/<jobid>/<step>', methods=['GET'])
@app.route('/results/<jobid>/<step>/', methods = ['GET'])
def showstep(jobid,step):
    jobinfo = routines.getjobstatus(jobid)
    paramdict = jobinfo.get("params")
    errmsgs = jobinfo.get("errors", False)
    #skips = paramdict.get("skip",[])
    laststep = request.args.get('laststep','NONE')
    jobname = routines.findjobinfo(jobid)

    #Redirect if job log doesnt exist or not present in redis
    if not routines.isjob(jobid) and jobinfo.get("status","Waiting in queue") == "Waiting in queue":
        flash('Invalid Jobid')
        return redirect('/analyze')

    if step == "loading":
        return render_template("startjob.html", jobid=jobid, jobname=jobname[1], laststep=laststep)
    elif step == "step2":
        if jobinfo["checkpoint"].upper() == "W1-STEP2":
            return render_template("step2.html",jobid=jobid,jobname=jobname[1])
        else:
            return redirect('/results/'+jobid+'/loading')
    elif step == "step3":
        if jobinfo["checkpoint"].upper() == "W1-STEP3":
            return render_template("step3.html",jobid=jobid,jobname=jobname[1])
        else:
            return redirect('/results/'+jobid+'/loading')
    elif step == "report":
        if "-F" in jobinfo["checkpoint"].upper() and not errmsgs:
            if jobinfo["workflow"] == "1":
                return render_template("report.html",jobid=jobid,jobname=jobname[1])
            else:
                return render_template("report.html",jobid=jobid, jobname=jobname[1],workflow=2)
        elif "-F" in jobinfo["checkpoint"].upper() and errmsgs: # currently means any errors at all prevent tree from being shown; separate checkpoint for fatal errors?
            return render_template("report.html",jobid=jobid,jobname=jobname[1], errmsgs = errmsgs)
        else:
            return redirect('/results/'+jobid+'/loading')


@app.route('/results/<jobid>/reanalyze', methods=['GET'])
def reanalyze(jobid):
    if request.args.get("confirm",False) and jobid != 'example':
        routines.reanalyzejob(jobid)
        jobname = routines.findjobinfo(jobid)
        return render_template("step2.html", jobid=jobid,jobname=jobname[1])
    else:
        return redirect('/results/'+jobid+'/loading')

@app.route('/results/<jobid>/tree', methods=['GET'])
def getTree(jobid):
    format = request.args.get('format','newick')
    resultdir = os.path.join(app.config['RESULTS_FOLDER'],jobid)
    resulttree = 'final'
    if os.path.exists(os.path.join(resultdir,resulttree+'.tree')):
        if format == 'newick':
            return send_from_directory(resultdir,resulttree+'.tree',as_attachment=True)
        # possible tree conversions?
        # elif format == 'nexml':
        #     Phylo.convert(os.path.join(resultdir,resulttree+'.tree'),'newick', os.path.join(resultdir,resulttree+'_nexml.xml'),'nexml')
        #     return send_from_directory(resultdir,resulttree+'_nexml.xml',as_attachment=True)
        # elif format == 'phyloxml':
        #     Phylo.convert(os.path.join(resultdir,resulttree +'.tree'), 'newick', os.path.join(resultdir,resulttree +'_phyloxml.xml'), 'phyloxml')
        #     return send_from_directory(resultdir, resulttree + '_phyloxml.xml', as_attachment=True)
    else:
        return "false"


# wrap all the downloads except for the tree into one page?
@app.route('/results/<jobid>/downloadorgs', methods=['GET'])
def downloadorgs(jobid):
    format = request.args.get('format','json')
    resultdir = os.path.join(app.config['RESULTS_FOLDER'],jobid)
    if format == 'json': # narrow down like the txt version?
        return send_from_directory(resultdir, 'reflist.json', as_attachment=True)
    elif format == 'txt':
        if os.path.exists(os.path.join(resultdir, 'reftext.txt')):
            return send_from_directory(resultdir,'reftext.txt',as_attachment=True)
        else:
            jsonpath = os.path.join(resultdir,'reflist.json')
            routines.jsontotsv(jsonpath,jobid) # make a more versatile version of that?
            return send_from_directory(resultdir, 'reftext.txt', as_attachment=True)

@app.route('/results/<jobid>/downloadmash')
def downloadmash(jobid):
    resultdir = os.path.join(app.config['RESULTS_FOLDER'], jobid)
    return send_from_directory(resultdir, 'mash_distances.txt', as_attachment=True)

@app.route('/results/<jobid>/downloadlists', methods=['GET'])
def downloadlists(jobid):
    downl = request.args.get('downl')
    resultdir = os.path.join(app.config['RESULTS_FOLDER'], jobid)
    if downl == 'mlstlist':
        if os.path.exists(os.path.join(resultdir,'mlst_genes.txt')):
            return send_from_directory(resultdir,'mlst_genes.txt', as_attachment=True)
        else:
            routines.mlsttsv(jobid)
            return send_from_directory(resultdir, 'mlst_genes.txt', as_attachment=True)
    elif downl == 'alignment':
        if os.path.exists(os.path.join(resultdir, jobid+'_alignments.zip')):
            return send_from_directory(resultdir, jobid+'_alignments.zip', as_attachment=True)
        else:
            routines.zipalignments(jobid)
            return send_from_directory(resultdir, jobid + '_alignments.zip', as_attachment=True)

@app.route('/results2/<jobid>/refs')
def getrefs(jobid):
    if os.path.exists(os.path.join(app.config['RESULTS_FOLDER'],'genuslist_example2.json')):
        return send_from_directory(app.config['RESULTS_FOLDER'],'genuslist_example2.json')
@app.route('/results2/selectgenus', methods=['POST'])
def selectgenus():
    jobid = request.form.get("jobinfo")
    newref = request.form.get("genusoptions")
    routines.updatejob(jobid,newref)
    genusdict={}
    with open(os.path.join(app.config['RESULTS_FOLDER'],'genuslist_example2.json'),'r') as genusfile:
        genusdict = json.load(genusfile)
        genusdict["genuslist"] = {newref:genusdict["genuslist"][newref]}
        genusdict["maxgenus"] = newref
    with open(os.path.join(app.config['RESULTS_FOLDER'],'genuslist_example2.json'),'w') as fileout:
        json.dump(genusdict,fileout,indent=2)
    return json.dumps({"status":1,"newmax":newref})

@app.route('/results2/refgenus')
def refgenus():
    if os.path.exists(os.path.join(app.config['RESULTS_FOLDER'],'acceleratedrefs.json')):
        return send_from_directory(app.config['RESULTS_FOLDER'],'acceleratedrefs.json')
    else:
        return jsonify([])

@app.route('/analyze')
def analyze():
    return render_template("analyze.html")

@app.route('/analyze2')
def analyze2():
    return render_template("analyze2.html")

@app.route('/upload', methods=['POST'])
def upload():
    #print request.files.getlist("seqfile1")
    filename = routines.getinfile()
    name = [os.path.split(x)[-1] for x in filename]
    #    filename = routines.getNCBIgbk(request.form["ncbiacc1"])
    return json.dumps({"filename": filename,"name": name})
@app.route('/startjob', methods=['POST'])
def startjob():
    jobid = unicode(uuid.uuid4())
    jobdict = {"id": jobid, "workflow": request.form.get("workflow"), "genomes": request.form.getlist('upfiles'),
               "skip":request.form.get('skip2',"")+","+request.form.get('skip3',""), "reference": request.form.get('genusselect','NA'),
               "bootstr":request.form.get('boots',0), "filtmlst": request.form.get('filtmlst',''),
               "mode":request.form.get('optradio',"concatenated"),"modelfind":request.form.get("modelfinder","GTR"), "fastalign":request.form.get("fastalign","")}

    os.mkdir(os.path.join(app.config['RESULTS_FOLDER'],jobid))
    resp = make_response(redirect('/results/' + jobid + '/loading'))

    #Emailer
    email = request.form.get('email', False)
    if email:
        # Run threaded so mail server resp does not block process
        @copy_current_request_context
        def sendmail(x, y, z):
            routines.sendnotifymail(x, y, z)

        mailer = threading.Thread(name='mail_sender', target=sendmail, args=("", jobid, email))
        mailer.start()

    #Jobtitle
    validchars = "abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ1234567890 _-()"
    jobtitle = request.form.get('jobname', False)
    jobtitle = ''.join([c for c in str(jobtitle[:30]) if c in validchars])
    if jobtitle and len(jobtitle):
        with open(os.path.join(app.config['RESULTS_FOLDER'], jobid, "jobtitle.txt"), "w") as namefile:
            namefile.write(jobtitle + "\n")

    #Recent results
    keeplink = request.form.get('keeplink',False)
    if keeplink and keeplink.lower() != "false":
        lastresult = request.cookies.get('automlst.lastresult')
        if lastresult and jobid not in lastresult:
            #Limit number of results kept
            lastresult = lastresult.split(";")
            lastresult.append(jobid)
            lastresult = ";".join(lastresult[-20:])
        else:
            lastresult = jobid
        resp.set_cookie('automlst.lastresult',lastresult)

    #Workflow redirecting
    if request.form.get("workflow") == "1" or request.form.get("workflow") == "2":
        automlstjob = routines.addjob(**jobdict)
        #with open(os.path.join(app.config['RESULTS_FOLDER'],'examplein.json'),'w') as uploadfile:
        #    json.dump(jobdict,uploadfile)
        return resp
    else:
        flash('Invalid workflow')
        return redirect('/analyze')

@app.route('/results/<jobid>/step2/orgs', methods=['GET'])
def getOrgs(jobid):
    orgstart = int(request.args.get('start',0))
    if os.path.exists(os.path.join(app.config['RESULTS_FOLDER'],jobid,'userlist.json')):
        with open(os.path.join(app.config['RESULTS_FOLDER'],jobid,'reflist.json'),'r') as reffile, open(os.path.join(app.config['RESULTS_FOLDER'],jobid,'userlist.json'),'r') as userfile:
            refdict = json.load(reffile)
            userdict = json.load(userfile)
            tempdict = [ref for ref in refdict["reforgs"] if ref["id"] in userdict["selspecies"] and not ref in refdict["reforgs"][0:(orgstart+500)]]
            refdict["reforgs"] = refdict["reforgs"][orgstart:(orgstart + 200)]
            refdict["outgroups"] = [rec for rec in refdict["outgroups"] if
                                    str(rec["phylid"]) != "N/A" and str(rec["familyid"]) != "N/A" and str(
                                        rec["orderid"]) != "N/A" and str(rec["genusid"]) != "N/A"]

            refdict["reforgs"].extend(tempdict)
            refdict.update(userdict)
            return jsonify(refdict)
    elif os.path.exists(os.path.join(app.config['RESULTS_FOLDER'],jobid,'reflist.json')):
        with open(os.path.join(app.config['RESULTS_FOLDER'],jobid,'reflist.json'),'r') as firstref:
            refdict = json.load(firstref)
            refdict["reforgs"] = refdict["reforgs"][orgstart:(orgstart+200)]
            refdict["outgroups"] = [rec for rec in refdict["outgroups"] if str(rec["phylid"]) != "N/A" and str(rec["familyid"]) != "N/A" and str(rec["orderid"]) != "N/A" and str(rec["genusid"]) != "N/A"]
            return jsonify(refdict)
    return jsonify({"error":"No List found."})

@app.route('/results/<jobid>/step2/orgin', methods=['POST'])
@app.route('/results/<jobid>/reanalyze/orgin', methods=['POST'])
@app.route('/results/<jobid>/orgin', methods=['POST'])
def orgin(jobid):
    species = request.form.getlist('specieslist')
    outgroups = request.form.getlist('outgrlist')
    jobid = request.form.get('jobinfo')
    with open(os.path.join(app.config['RESULTS_FOLDER'],jobid,'userlist.json'),'w') as userfile:
        json.dump({"selspecies":species, "seloutgroups":outgroups},userfile)
    with open(os.path.join(app.config['RESULTS_FOLDER'],jobid,'automlst.log'),'a') as joblog:
        joblog.write('\n'+str(datetime.datetime.now())+' - INFO - JOB_CHECKPOINT::w1-3 \n'+str(datetime.datetime.now())+' - INFO - JOB_STATUS::Resuming job...\n') # is this still right?
    routines.readdjob(jobid)
    return redirect('/results/'+jobid+'/loading?laststep=step2') # when is checkpoint set? Might get user stuck in a loop of submitting/getting redirected

@app.route('/results/<jobid>/step2/outgroups', methods=['GET'])
def outgrs(jobid):
    ingroups=[]
    commongr = request.args.get('group', False)
    multigroups = request.args.get('multiple', False)
    if multigroups:
        commonid = (request.args.get('id',False)).split(',')
        ingroups = commonid
    else:
        commonid = request.args.get('id', False)
        ingroups.append(str(commonid))
    outgrlimit = request.args.get('limit',1000)
    startindex = request.args.get('start',0)
    if commongr and commonid:
        with open(os.path.join(app.config['RESULTS_FOLDER'],jobid,'reflist.json'),'r') as outrefs:
            outdict = json.load(outrefs)
            refrec = outdict['reforgs']
            newlist = [rec for rec in refrec if str(rec[commongr+"id"]) not in ingroups and str(rec[commongr+"id"]) != "N/A" and not rec in refrec[0:(startindex+500)]]
            return jsonify(newlist[0:500])
    return send_from_directory(app.config['RESULTS_FOLDER'],'outgroups.json')

@app.route('/results/<jobid>/step3/genes')
def getgenes(jobid):
    if os.path.exists(os.path.join(app.config['RESULTS_FOLDER'],jobid,'mlstpriority.json')):
        return send_from_directory(os.path.join(app.config['RESULTS_FOLDER'],jobid),'mlstpriority.json')

@app.route('/results/<jobid>/step3/genein', methods=['POST'])
def genein(jobid):
    jobid = request.form.get('jobinfo')
    genes = request.form.getlist('mlstlist')
    radioval = request.form.get('optradio')
    rmorgs = request.form.get('removeorgs','')
    with open(os.path.join(app.config['RESULTS_FOLDER'],jobid,'usergenes.json'),'w') as usergenes:
        json.dump({"selection":genes,"mode":radioval,"delorgs":rmorgs.split(",")},usergenes)
    with open(os.path.join(app.config['RESULTS_FOLDER'],jobid,'automlst.log'),'a') as joblog:
        joblog.write('\n'+str(datetime.datetime.now())+' - INFO - JOB_STATUS::Resuming job...\n')
        joblog.write(str(datetime.datetime.now())+' - INFO - JOB_CHECKPOINT::w1-5\n')
    routines.readdjob(jobid)
    return redirect('/results/'+jobid+'/loading?laststep=step3')

@app.route('/results/<jobid>/mash')
def showmash(jobid):
    resultpath =os.path.join(app.config['RESULTS_FOLDER'],jobid)
    tsvpath = os.path.join(resultpath,'mash_distances.txt')
    if os.path.exists(tsvpath):
        jsondata = routines.tsvtojson(tsvpath)
        jsondata["data"] = [rec for rec in jsondata["data"] if float(rec[4])>=0.65]
        return jsonify(jsondata)
    else:
        nodata = {}
        nodata["data"] = []
        return jsonify(nodata)

@app.route('/aniclades')
def aniclades(): # move to static?
    if os.path.exists(os.path.join(app.config['RESULTS_FOLDER'], 'aniclades.json')):
        with open(os.path.join(app.config['RESULTS_FOLDER'], 'aniclades.json'),'r') as anifile:
            anijson = json.load(anifile)
            return jsonify(anijson)

@app.route('/jobstatus/<jobid>')
@app.route('/jobstatus/<jobid>/')
def status(jobid):
    jobstat = routines.getjobstatus(jobid)
    workflow = jobstat["workflow"]
    paramdict = jobstat.get("params")
    if jobstat["checkpoint"].upper() == "W1-STEP2" and workflow == "1":
        #redirdict = {"redirect":"step2"}
        #return jsonify(redirdict)
        jobstat["redirect"] = "step2"
        return jsonify(jobstat)
    #elif jobstat["checkpoint"].upper() == "W1-2" and "skip2" in skips and workflow == "1":
        #jobstat["skip"] = "c1"
        #return jsonify(jobstat)
    elif jobstat["checkpoint"].upper() == "W1-STEP3" and workflow == "1":
        #redirdict = {"redirect":"step3"}
        #return jsonify(redirdict)
        jobstat["redirect"] = "step3"
        return jsonify(jobstat)
    #elif jobstat["checkpoint"].upper() == "W1-3" and "skip3" in skips and workflow == "1":
        #jobstat["skip"] = "c2"
        #return jsonify(jobstat)
    # elif jobstat["checkpoint"].upper() == "W1-F" or jobstat["checkpoint"].upper() == "W2-F":
    elif "-F" in jobstat["checkpoint"].upper():
        #redirdict = {"redirect":"report"}
        #return jsonify(redirdict)
        jobstat["redirect"] = "report"
        return jsonify(jobstat)
    else:
        return jsonify(jobstat)

# probably needs to be updated; simply redirect to completed example job instead? If so, ensure it can't be reanalyzed
# @app.route('/results/example/report')
# def example():
#     return render_template("example.html")


@app.errorhandler(404)
@app.errorhandler(401)
@app.errorhandler(500)
def page_not_found(e):
    return render_template('error.html',title="",errormsg=e)

