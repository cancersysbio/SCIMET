#! /usr/bin/python

######################
#Simulation of the partion history of sampled demes in primary tumor
#4 bulk samples (named "X","Y","Z","W") wered sampled at the time when primary tumor size is 10^9 cells
#Also track the deme where the metastatic founder cell was sampled from
######################

import sys,os,math,random
import numpy as np
from collections import Counter
import sets

class node():
    def __init__(self):
        self.name = ""
        self.internal = 1
        self.isroot = 0
        self.blength = 0
        self.index = 0
        self.isancestor = 0

        self.anc=None
        self.left=None
        self.right=None      

def print_newick(tmp_root):
    if tmp_root.internal==1:
        print "(",
        if tmp_root.left != None:
            print_newick(tmp_root.left)
        print ",",
        if tmp_root.right != None:
            print_newick(tmp_root.right)
        print ")",
        print ":%f"%tmp_root.blength,

    else:

        print "%s:%f"%(tmp_root.name,tmp_root.blength),


def print_newick2file(tmp_root,sendout):
    if tmp_root.internal==1:
        sendout.write("(")
        print_newick2file(tmp_root.left,sendout)
        sendout.write(",")
        print_newick2file(tmp_root.right,sendout)
        sendout.write(")")
        sendout.write(":%f"%tmp_root.blength)
    else:
        sendout.write("%s:%f"%(tmp_root.name,tmp_root.blength))

def copy_nodeTree(nodelist):
    nodelist_copy = []
    for z in range(0,len(nodelist)):
        znode = node()
        nodelist_copy.append(znode)
    for x in range(0,len(nodelist)):
        #nodelist_copy.append(nodelist[x])
        nodelist_copy[x].name = nodelist[x].name
        nodelist_copy[x].internal = nodelist[x].internal
        nodelist_copy[x].isroot = nodelist[x].isroot
        nodelist_copy[x].blength = nodelist[x].blength
        #nodelist_copy[x].tlength = nodelist[x].tlength
        nodelist_copy[x].index = nodelist[x].index
        nodelist_copy[x].isancestor = nodelist[x].isancestor
    for y in range(0,len(nodelist)):
        if nodelist[y].anc != None:
            nodelist_copy[y].anc = nodelist_copy[nodelist[y].anc.index]
        else:
            nodelist_copy[y].anc = None
        if nodelist[y].left != None:
            nodelist_copy[y].left = nodelist_copy[nodelist[y].left.index]
        else:
            nodelist_copy[y].left = None
        if nodelist[y].right != None:
            nodelist_copy[y].right = nodelist_copy[nodelist[y].right.index]
        else:
            nodelist_copy[y].right = None
    return nodelist_copy

def collect_leaves(subtree):
    if subtree.internal==1:
        return collect_leaves(subtree.left) + collect_leaves(subtree.right)
    else:
        return [subtree]
    

def label_ancestors(subtree):
    subtree.isancestor +=1
    #print "labeling node",subtree.index,"as ancestors"
    if subtree.anc !=None:
        label_ancestors(subtree.anc)
        

def label_tree(leaf,remove_leaf):
    for w in leaf:
        if w.name in remove_leaf:
            #print w.name
            #print "tracking ancestor for tip node",w.name
            label_ancestors(w)

def find_MRCA(subtree,samplesize):
    if subtree.isancestor==samplesize:
        if subtree.left.isancestor==samplesize:
            return find_MRCA(subtree.left,samplesize)
        
        elif subtree.right.isancestor==samplesize:
            return find_MRCA(subtree.right,samplesize)
        
        else:
            return subtree

def prune_tree(subtree,MRCA):        # assuming already started with MRCA
    #print "let's deal with node",subtree.index

    if subtree.internal==1:
        #print "doing something good"
        if subtree.left != None:
            prune_tree(subtree.left,MRCA)
        if subtree.right != None:
            prune_tree(subtree.right,MRCA)

        if subtree.left != None and subtree.right != None:
            if subtree==MRCA:
                #print subtree.left.index,subtree.right.index
                #pass
                mws12=0
            if subtree.left.isancestor>=1 and subtree.right.isancestor>=1:
                #print "I don't do anything for node",subtree.index
                mws12=1
            else:
                ancestor=subtree.anc
                if subtree.left.isancestor>=1:
                    if ancestor.left==subtree:
                        #print "case 11"
                        left_kid=subtree.left
                        left_kid.anc=ancestor
                        ancestor.left=left_kid
                        left_kid.blength +=subtree.blength 

                    else:# ancestor.right==subtree:
                        #print "case 12"
                        left_kid=subtree.left
                        left_kid.anc=ancestor
                        ancestor.right=left_kid
                        left_kid.blength +=subtree.blength 

                elif subtree.right.isancestor>=1:
                    if ancestor.left==subtree:
                        #print "case 21"
                        right_kid=subtree.right
                        right_kid.anc=ancestor
                        ancestor.left=right_kid
                        right_kid.blength +=subtree.blength 

                    else:     # ancestor.right==subtree:
                        #print "case 22"
                        right_kid=subtree.right
                        right_kid.anc=ancestor
                        ancestor.right=right_kid
                        right_kid.blength +=subtree.blength 
        
        elif subtree.right == None and subtree.left != None:
            if subtree.left.isancestor>=1:
                ancestor=subtree.anc
                if ancestor.left==subtree:
                    left_kid=subtree.left
                    left_kid.anc=ancestor
                    ancestor.left=left_kid
                    left_kid.blength +=subtree.blength 

                else:   # ancestor.right==subtree:
                    left_kid=subtree.left
                    left_kid.anc=ancestor
                    ancestor.right=left_kid
                    left_kid.blength +=subtree.blength 


        elif subtree.left == None and subtree.right != None:
            if subtree.right.isancestor>=1:
                ancestor=subtree.anc
                if ancestor.left==subtree:
                    right_kid=subtree.right
                    right_kid.anc=ancestor
                    ancestor.left=right_kid
                    right_kid.blength +=subtree.blength 

                else:   # ancestor.right==subtree:
                    right_kid=subtree.right
                    right_kid.anc=ancestor
                    ancestor.right=right_kid
                    right_kid.blength +=subtree.blength 

def createLattice(r):
    """
    This is a function to create 3D cubic lattice space with length 2r
    """
    cspace= {}
    for x in range(0,2*r+1):
        for y in range(0,2*r+1):
            for z in range(0,2*r+1):
                cspace[(x,y,z)] = -1
    return cspace

def neighbor26((a,b,c)):
    """
    Collect the 26 neighbour sites of (a,b,c)
    """
    neighbor =[(a-1, b-1, c-1),(a-1, b-1, c),(a-1, b-1, c+1),(a-1, b, c-1),(a-1, b, c),(a-1, b, c+1),(a-1, b+1, c-1),(a-1, b+1, c),(a-1, b+1, c+1),(a, b-1, c-1),(a, b-1, c),(a, b-1, c+1),(a, b, c-1),(a, b, c+1),(a, b+1, c-1),(a, b+1, c),(a, b+1, c+1),(a+1, b-1, c-1),(a+1, b-1, c),(a+1, b-1, c+1),(a+1, b, c-1),(a+1, b, c),(a+1, b, c+1),(a+1, b+1, c-1),(a+1, b+1, c),(a+1, b+1, c+1)]
    return neighbor

def neighbor6((a,b,c)):
    """
    Collect the 6 neighbour sites of (a,b,c)
    """
    neighbor =[(a-1, b, c),(a+1, b, c),(a, b-1, c),(a, b+1, c),(a, b, c-1),(a, b, c+1)]
    return neighbor
    

def localNeighbor((a,b,c),r):
    """
    This is function to collect the local neighbours of position (a,b,c) within an area of radius r in 3D cubic lattice
    """
    neighbour = []
    for x in range(-r,r+1):
        for y in range(-r,r+1):
            for z in range(-r,r+1):
                if pow(x,2)+pow(y,2)+pow(z,2) < pow(r+1,2):
                    neighbour += [(a+x,b+y,c+z)]
    return neighbour

def samplingInEdge(sp,perisite,radius):
    tx,ty,tz = perisite[0],perisite[1],perisite[2]
    localnei = localNeighbor((tx,ty,tz),radius)
    tsample = []
    for x in localnei:
        if sp[x] != -1:
            tsample += [x]
    return tsample

def sparseSampling(region,sample_number,cutoff):
    success = 0
    while success == 0:
        samples = random.sample(region,sample_number)
        repeat = sample_number*(sample_number-1)
        minall = 999
        for x in range(0,repeat):
            rs = random.sample(samples,2)
            min_distance = min([abs(rs[0][0]-rs[1][0]),abs(rs[0][1]-rs[1][1]),abs(rs[0][2]-rs[1][2])])
            if min_distance < minall:
                minall = min_distance
        if min_distance > 2*cutoff:
            success = 1
    return samples

def check_space(loc,sp):
    nei_loc = neighbor26(loc)
    empty_loc = []
    for ekey in nei_loc:
        if sp[ekey] == -1:
            empty_loc += [ekey]
    return empty_loc
    

############main script#################
met_time = int(sys.argv[1])
repl = int(sys.argv[2])

rd = 60                         #radius of the 3D space
pop_size = 200000               #number of demes in final tumor

send_tree = open("tumor3D_primary4samples_Nd"+str(met_time)+".tree"+str(repl),"w")

for w in range(0,100):
    seeding_node = node()
    seeding_node.index = 0
    seeding_node.isroot = 1

    space = createLattice(rd)
    space[(rd,rd,rd)] = 0
    node_tree = [seeding_node]
    current_keys = [(rd,rd,rd)]
    current_size = 0 #current tumor size in deme number
    gindex = 0
    #surface = [(rd,rd,rd)]

    while current_size < pop_size:
        skey = random.choice(current_keys)
        rx,ry,rz = skey[0],skey[1],skey[2]
        neisites = neighbor26((rx,ry,rz))
        empty_neis = []
        for nkey in neisites:
            if space[nkey] == -1:
                empty_neis += [nkey]
        #if len(empty_neis) == 0:
        #    print "something is wrong!!"
        #    break
        if len(empty_neis) > 0:
            rand_prob = random.random()
            if rand_prob < 1-math.exp(-len(empty_neis)*0.25): #deme division
                nextkey = random.choice(empty_neis)
                #current_keys += [nextkey]
                current_size += 1
            
                left = node()
                right = node()
                gindex += 2
                left.index = int(gindex-1)
                right.index = int(gindex)
                left.blength += 1
                right.blength += 1
                left.anc = node_tree[space[skey]]
                right.anc = node_tree[space[skey]]
                node_tree += [left,right]
                node_tree[space[skey]].left = node_tree[gindex-1]
                node_tree[space[skey]].right = node_tree[gindex]
                space[skey] = int(gindex-1)
                space[nextkey] = int(gindex)

                if current_size == met_time:
                    met_deme_index = int(space[nextkey])
                    node_tree[space[nextkey]].internal = -1
                    node_tree[space[nextkey]].name = "met"
                    #met_deme.name = "met"
                    space[nextkey] = -1
                    #met_deme.internal = -1
                else:
                    current_keys += [nextkey]
                    #surface += [nextkey]
                #nei_next = neighbor26(nextkey)
                #for nekey in nei_next:
                #    if nekey in surface:
                #        if len(check_space(nekey,space)) == 0:
                #            surface.remove(nekey)
                
    current_index = []
    for ckey in current_keys:
        current_index += [int(space[ckey])]
    for leave_index in current_index:
        node_tree[leave_index].internal = -1 #current demes are at the leaves
    
    ########sampling bulk tissues########
    periphery = []
    for key in current_keys:
        keynei = neighbor26((key[0],key[1],key[2]))
        for z in keynei:
            if space[z] == -1:
                periphery +=[key]
                break
    quadrant = []
    for pky in periphery:
        if pky[0] > rd and pky[1] > rd and pky[2] > rd:
            quadrant += [pky]
    
    print "replicate=",w
    print "# of demes in this region=",len(quadrant)
    sample_quadrant = sparseSampling(quadrant,4,3)
    #print sample_quadrant
    tissue_sample1 = samplingInEdge(space,sample_quadrant[0],2)
    tissue_sample2 = samplingInEdge(space,sample_quadrant[1],2)
    tissue_sample3 = samplingInEdge(space,sample_quadrant[2],2)
    tissue_sample4 = samplingInEdge(space,sample_quadrant[3],2)
    
    #print "total # of nodes=",len(node_tree)
    #print "# of demes in the periphery=",len(periphery)
    #print tissue_sample1
    #print tissue_sample2
    #print tissue_sample3
    print "# of demes in the four bulk samples=",len(tissue_sample1),len(tissue_sample2),len(tissue_sample3),len(tissue_sample4)
    print

    node_treeX = copy_nodeTree(node_tree)
    #node_treeY = []
    #node_treeY = copy_nodeTree(node_treeX)
    #sample_cell = random.sample(current_index,kcells)
    cell_name = 1
    leaves = []
    toremove = []
    for xkey in tissue_sample1:
        node_treeX[space[xkey]].name = "X"+str(cell_name)
        cell_name += 1
        leaves.append(node_treeX[space[xkey]])
        toremove.append(node_treeX[space[xkey]].name)
    cell_name = 1
    for xkey in tissue_sample2:
        node_treeX[space[xkey]].name = "Y"+str(cell_name)
        cell_name += 1
        leaves.append(node_treeX[space[xkey]])
        toremove.append(node_treeX[space[xkey]].name)
    cell_name = 1
    for xkey in tissue_sample3:
        node_treeX[space[xkey]].name = "Z"+str(cell_name)
        cell_name += 1
        leaves.append(node_treeX[space[xkey]])
        toremove.append(node_treeX[space[xkey]].name)
    cell_name = 1
    for xkey in tissue_sample4:
        node_treeX[space[xkey]].name = "W"+str(cell_name)
        cell_name += 1
        leaves.append(node_treeX[space[xkey]])
        toremove.append(node_treeX[space[xkey]].name)
    
    leaves.append(node_treeX[met_deme_index])
    toremove.append(node_treeX[met_deme_index].name)
    label_tree(leaves,toremove)
    MRCA=find_MRCA(node_treeX[0],len(toremove))
    prune_tree(node_treeX[MRCA.index],MRCA)
    ##T2_mrca += MRCA.left.blength
    ##print_newick(node_tree_copy[MRCA.index])
    ##print ""
    print_newick2file(node_treeX[MRCA.index],send_tree)
    send_tree.write(";\n")
