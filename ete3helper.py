#!/usr/bin/env python

from ete3 import Tree

def rerootTree(treefile,output,ogterm="OG--",fmat=3):
    intree = Tree(treefile)
    og = [x for x in intree.get_leaves() if x.name.startswith(ogterm)]
    if not len(og):
        return False

    og = og[0]
    intree.set_outgroup(og.name)
    intree.ladderize(direction=1)
    intree.write(outfile=output,format=fmat)
    return True