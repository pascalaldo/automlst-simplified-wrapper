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

import glob, argparse, json, itertools
from ete3 import Tree

def checkAniClade(treefil,anidict,shrink=False,unrooted=False):
    """Check the percentage of groups that are monophyletic by calculating max members that are monophyletic out of total in each ANI group"""
    T = Tree(treefil)

    groupdict = {}
    treesize = 0
    #lookup and add groupid
    for leaf in T:
        gcfindex = leaf.name.find('GCF_')
        gcfid = leaf.name[gcfindex:gcfindex+13]
        groupid = anidict.get(gcfid,-1)
        leaf.add_features(groupid=groupid)
        treesize += 1

        if groupid not in groupdict:
            groupdict[groupid] = []
        groupdict[groupid].append(leaf.name)

    #reroot by any singleton to avoid artificial breaks in any ani clade
    if unrooted:
        for k,v in groupdict.items():
            if len(v) == 1:
                T.set_outgroup(T&v[0])
                break

    #ignore singletons
    groupdict = {k:v for k,v in groupdict.items() if len(v)>1}
    groupsize = sum([len(g) for g in groupdict.values()])
    scoredict = {k:0 for k in groupdict.keys()}

    #Adds splits until split with no direct monophyletic child is found
    def shrinknodes(node,maxleafs,testset):
        children = node.get_children()
        if len(children) == 2:
            c1 = [x.name for x in children[0].get_leaves()]
            c2 = [x.name for x in children[1].get_leaves()]
            #add to maxleafs
            if len(set(c1)&testset) == len(c1) and len(set(c2)&testset) != len(c2):
                maxleafs = maxleafs|set(c1)
                node = children[1]
                return shrinknodes(node,maxleafs,testset)
            if len(set(c2)&testset) == len(c2) and len(set(c1)&testset) != len(c1):
                maxleafs = maxleafs|set(c2)
                node = children[0]
                return shrinknodes(node,maxleafs,testset)
        return maxleafs

    #Gets group by traversing ancestors from a start leaf
    def grownodes(node,maxleafs,testset):
        node = node.up
        if node:
            temp = set([x.name for x in node.get_leaves()])
            if len(temp&testset) == len(temp):
                maxleafs = temp
                return grownodes(node,maxleafs,testset)
            if shrink:
                return shrinknodes(node,maxleafs,testset)
        return maxleafs

    def getlargestclade(testset):
        maxset = set()
        # remaining = testset
        for org in list(testset):
            # startnode = T.get_leaves_by_name(list(remaining)[0])[0]
            startnode = T.get_leaves_by_name(org)[0]
            temp = grownodes(startnode,set([startnode.name]),testset)
            if len(temp) > len(maxset):
                maxset = temp
            # remaining = remaining - temp
        return maxset

    scoredict.update( {gid:len(getlargestclade( set(members) )) for gid,members in groupdict.items()} )
    #get max mono for each group
    # for groupid,allmembers in groupdict.items():
        # mcount = len(allmembers)
        # while mcount > 1:
        #     if any( [T.check_monophyly(values=x,target_attr="name",unrooted=True)[0] for x in itertools.combinations(allmembers,mcount)] ):
        #         scoredict[groupid] = mcount
        #         break
        #     mcount -= 1

    #return results #total_monophyl/total_groupsize, total_groupsize/treesize, numgroups, numperfect
    if groupsize and treesize:
        return (float(sum(scoredict.values()))/groupsize, float(groupsize)/treesize, len(groupdict.keys()), len([k for k,v in groupdict.items() if len(v)==scoredict[k]]) )
    else:
        return False

def checkAll(treefiles,aniJson,outfile,shrink=False,unrooted=False):
    """Run checkAniClade for every ANI group in json file"""
    treefiles = glob.glob(treefiles)

    #Get Ani group ids and read info columns if exists
    with open(aniJson,"r") as fil:
        aniJson = json.load(fil)
    if "_info" in aniJson:
        info = [str(x) for x in aniJson["_info"]]
        del aniJson["_info"]
    else:
        info = [str(x) for x in range(len(aniJson.values()[0]))]

    with open(outfile,"w") as ofil:
        #write header:  Treefile   (ANI level avg %monophy) (tree coverage "cov") (total groups "grps") (100% mono groups "pfct")
        ofil.write("#TreeFile\t%s\t%s\t%s\t%s\n" % ( "\t".join(info),
                                                     "\t".join([x+"_cov" for x in info]),
                                                     "\t".join([x+"_grps" for x in info]),
                                                     "\t".join([x+"_pfct" for x in info])))
        #loop through trees
        for treefil in treefiles:
            #loop thorugh all levels
            results = []
            for i in range(len(info)):
                print "Starting tree: %s"%treefil
                result = checkAniClade(treefil, {k:v[i] for k,v in aniJson.items()}, shrink=shrink)
                if result:
                    results.append(result)
                else:
                    print "No multi isolate groups found for %s at ani index %s"%(treefil,i)
                    results.append(["N/A"]*4)
            #write results
            ofil.write("%s\t%s\t%s\t%s\t%s\n"% (treefil,
                                                "\t".join([str(x[0]) for x in results]),
                                                "\t".join([str(x[1]) for x in results]),
                                                "\t".join([str(x[2]) for x in results]),
                                                "\t".join([str(x[3]) for x in results])) )

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="Check tree for proper clading of ANI groups")
    parser.add_argument("input",nargs='?',help="Tree files expression ex: directory/*.tree default=*.tree",default="*.tree")
    parser.add_argument("-o","--out",help="Store results in OUT (default: allscores.tsv)",default="allscores.tsv")
    parser.add_argument("-i","--inc", help="include diverging groups into max monophyletic",action="store_true",default=False)
    parser.add_argument("-u","--uroot", help="Set tree as unrooted",action="store_true",default=False)
    parser.add_argument("-aj","--anijson",help="ANI json file, GCF_xxx to ani group ids",default="")
    args = parser.parse_args()
    checkAll(args.input,args.anijson,args.out,args.inc,args.uroot)