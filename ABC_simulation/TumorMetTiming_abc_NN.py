#! /usr/bin/python

######################
#Simulation of evolutionary dynamics in paired priamry tumor and metastasis 
#Model: Neutral/Neutral (N/N)
#Output: 9 summary statistics for ABC inference
#copyright- Zheng Hu
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
        
        self.mutlist1=[] ##neutral mutations
        self.mutlist2=[] ##advantageous mutations


def traceLineage(mlineage,label):
    ancestor_str = label.split(',')
    ancestor_list = [int(t) for t in ancestor_str]
    first_anc = ancestor_list[0]
    trace = []
    while first_anc > 0:
        trace += ancestor_list
        ancestor_str = mlineage[first_anc].split(',')
        ancestor_list = [int(t) for t in ancestor_str]
        first_anc = ancestor_list[0]
    return trace


def lowerORupper(value):
    """
    A function to choose the upper or lower integral value given a non-integral number
    """
    lower_int = int(value)
    upper_int = lower_int+1
    if random.random() < value-lower_int:
        return upper_int
    else:
        return lower_int


def initiateFirsDeme_neutral(initial_cells,lineage,cur_index,division_prob):
    """
    The growth of the initial deme from a single transformed tumor cell via a random discrete-time birth-death process
    model - neutral model
    """
    current_deme_size = len(initial_cells)
    neu_list = list(initial_cells)
    while current_deme_size < 5000: #5000 is the number of cells in a tumor deme
        divcells =  int(current_deme_size*division_prob+1) #number of dividing cells in this generation
        neu_list = random.sample(neu_list,divcells)*2 #the mutation distribution in the sampled cells.
        current_deme_size = divcells*2
        new_mut = np.random.poisson(mu*current_deme_size)
        mut_assignment = Counter(np.random.choice(current_deme_size,new_mut))
        for x in mut_assignment.keys():
            mmut = mut_assignment[x]
            mut_strlist = range(cur_index+1,cur_index+1+mmut)
            mut_str = ",".join(map(str,mut_strlist))
            for xi in range(0,mmut):
                cur_index += 1
                lineage += [neu_list[x]]
            neu_list[x] = mut_str
    return neu_list,lineage,cur_index


def initiateFirsDeme_selection(initial_cells,lineage,cur_index,division_prob,sfit):
    """
    The growth of the initial deme from a single transformed tumor cell via a random discrete-time birth-death process
    model - selection model
    """
    neu_list = list(initial_cells)
    adv_list = []
    current_deme_size = len(initial_cells)
    while current_deme_size < 5000:
        n1,n2 = len(neu_list),len(adv_list)                         #n1 and n2 are the current number of neutral founder cells and advantageous cells, respectively
        neu_divcells =  int(n1*division_prob+1)                        #number of dividing cells of neutral lineage in this generation. The other cells will die in the next generation
        neu_list = random.sample(neu_list,neu_divcells)*2
        if n2 > 0:
            adv_divcells = lowerORupper(n2*division_prob*(1+sfit))   #number of dividing cells of advantageous lineage in this generation        
            adv_list = random.sample(adv_list,adv_divcells)*2
        n1,n2 = len(neu_list),len(adv_list)
        current_deme_size = n1+n2
        if n1 > 0:
            new_mut1 = np.random.poisson(mu*n1)               # the total number of mutations occurring in a generation follows Poission distribution with lambda=u*n
            mut_assig1 = Counter(np.random.choice(n1,new_mut1))
            for x1 in mut_assig1.keys():
                nmut = mut_assig1[x1]
                new_mut1 = range(cur_index+1,cur_index+1+nmut)
                mut_str = ",".join(map(str,new_mut1))
                #if nmut > 1:
                #    for t in new_mut1:
                #        multi_events[str(t)] = mut_str
                for xn in range(0,nmut):
                    cur_index += 1
                    lineage += [neu_list[x1]]
                neu_list[x1] = mut_str
        if n2 > 0:
            new_mut2 = np.random.poisson(mu*n2)
            mut_assig2 = Counter(np.random.choice(n2,new_mut2))
            for x2 in mut_assig2.keys():
                nmut = mut_assig2[x2]
                new_mut2 = range(cur_index+1,cur_index+1+nmut)
                mut_str = ",".join(map(str,new_mut2))
                #if nmut > 1:
                #    for t in new_mut2:
                #        multi_events[str(t)] = mut_str
                for xn in range(0,nmut):
                    cur_index += 1
                    lineage += [adv_list[x2]]
                adv_list[x2] = mut_str
        
        if random.random() < adv_prob*n1:                           # occurence of advantageous mutation on the neutral lineage
            cur_index += 1
            current_n1 = len(neu_list)
            lineage += [str(neu_list[current_n1-1])]
            adv_list += [str(cur_index)]
            neu_list = neu_list[0:current_n1-1]
    
    return neu_list,adv_list,lineage,cur_index


def copy_node(xnode):
    ynode = node()
    ynode.isroot = xnode.isroot
    ynode.blength = xnode.blength
    ynode.isancestor = xnode.isancestor
    ynode.anc=xnode.anc
    ynode.left=xnode.left
    ynode.right=xnode.right
    
    return ynode


def vertical_blength(node1): # distence to root
    vlength = node1.blength
    node2 = copy_node(node1)
    if node2.isroot == 1:
        return vlength
    else:
        while node2.isroot == 0:
            node2 = copy_node(node2.anc)
            vlength += node2.blength
        return vlength


def deme_fission_neutral(cnode,cmutlist,lineage,cur_index,division_prob):
    xlist = list(cmutlist)
    total_blength = vertical_blength(cnode)
    #print "blength=",cnode.blength
    #print "total_blength=",total_blength
    if cnode.isroot == 1:
        further = int(cnode.blength)
    else:
        further = int(cnode.blength)-1
    for bl in range(0,further):
        deme_size = len(xlist)
        while deme_size < 10000:
            Ndivcells =  int(deme_size*division_prob) #number of dividing cells in this generation
            xlist = random.sample(xlist,Ndivcells)*2 #the mutation distribution in the sampled cells.
            deme_size = Ndivcells*2
            if total_blength < 80:                   #top deme fissions
                new_mutation = np.random.poisson(mu*deme_size)
                mut_assignment = Counter(np.random.choice(deme_size,new_mutation))
                for x in mut_assignment.keys():
                    nmut = mut_assignment[x]
                    mut_strlist = range(cur_index+1,cur_index+1+nmut)
                    mut_str = ",".join(map(str,mut_strlist))
                    #mut_str = str(range(cur_index+1,cur_index+1+nmut))
                    #mut_str = mut_str[1:len(mut_str)-1]
                    for xi in range(0,nmut):
                        cur_index += 1
                        #lineage[cur_index] = xlist[x]
                        lineage += [xlist[x]]
                    xlist[x] = mut_str
        random.shuffle(xlist)
        xlist=xlist[0:deme_size/2]
    cnode.mutlist1 = list(xlist)
    deme_size = len(xlist)
    while deme_size < 10000:
        Ndivcells =  int(deme_size*division_prob) #number of dividing cells in this generation
        xlist = random.sample(xlist,Ndivcells)*2 #the mutation distribution in the sampled cells.
        deme_size = Ndivcells*2
        if total_blength < 80:                   #top deme fissions
            new_mutation = np.random.poisson(mu*deme_size)
            mut_assignment = Counter(np.random.choice(deme_size,new_mutation))
            for x in mut_assignment.keys():
                nmut = mut_assignment[x]
                mut_strlist = range(cur_index+1,cur_index+1+nmut)
                mut_str = ",".join(map(str,mut_strlist))
                #mut_str = str(range(cur_index+1,cur_index+1+nmut))
                #mut_str = mut_str[1:len(mut_str)-1]
                for xi in range(0,nmut):
                    cur_index += 1
                    #lineage[cur_index] = xlist[x]
                    lineage += [xlist[x]]
                xlist[x] = mut_str
    random.shuffle(xlist)
    xlist1=xlist[0:deme_size/2]
    xlist2=xlist[deme_size/2:deme_size]
    
    return xlist1,xlist2,cur_index,lineage,cnode


def deme_fission_selection(cnode,cmutlist1,cmutlist2,lineage,cur_index,division_prob,sfit):
    xlist = list(cmutlist1)
    ylist = list(cmutlist2)
    total_blength = vertical_blength(cnode)
    #print "blength=",cnode.blength
    #print "total_blength=",total_blength
    if cnode.isroot == 1:
        further = int(cnode.blength)
    else:
        further = int(cnode.blength)-1
    for bl in range(0,further):
        deme_size = len(xlist)+len(ylist)
        while deme_size < 10000:
            n1,n2 = len(xlist),len(ylist)
            Ndivcells =  int(deme_size*division_prob) #number of dividing cells in this generation
            deme_size = Ndivcells*2
            if n2 > 0:
                ydown = int(division_prob*(1+sfit)*n2)
                yup = ydown+1
                plus_prob=division_prob*(1+sfit)*n2-ydown
                if random.random() < plus_prob:
                    ynumber = int(yup)
                else:
                    ynumber = int(ydown)
                if ynumber > 10000:
                    ynumber = int(Ndivcells)
                ylist = random.sample(ylist,ynumber)*2
            else:
                ynumber = 0
                ylist = []
            xnumber = Ndivcells-ynumber
            if xnumber < 0:
                xlist = []
            else:
                xlist = random.sample(xlist,xnumber)*2
            n1,n2 = len(xlist),len(ylist)
            if total_blength < 80:
                if n1 > 0:
                    new_mut1 = np.random.poisson(mu*n1)
                    mut_assig1 = Counter(np.random.choice(n1,new_mut1))
                    for x1 in mut_assig1.keys():
                        nmut = mut_assig1[x1]
                        mut_strlist = range(cur_index+1,cur_index+1+nmut)
                        mut_str = ",".join(map(str,mut_strlist))
                        for xn in range(0,nmut):
                            cur_index += 1
                            lineage += [xlist[x1]]
                        xlist[x1] = mut_str
                if n2 > 0:
                    new_mut2 = np.random.poisson(mu*n2)
                    mut_assig2 = Counter(np.random.choice(n2,new_mut2))
                    for x2 in mut_assig2.keys():
                        nmut = mut_assig2[x2]
                        mut_strlist = range(cur_index+1,cur_index+1+nmut)
                        mut_str = ",".join(map(str,mut_strlist))
                        #mut_str = str(range(cur_index+1,cur_index+1+nmut))
                        #mut_str = mut_str[1:len(mut_str)-1]
                        for xn in range(0,nmut):
                            cur_index += 1
                            lineage += [ylist[x2]]
                        ylist[x2] = mut_str
                if random.random() < adv_prob*n1:
                    cur_index += 1
                    cur_n1 = len(xlist)
                    lineage += [str(xlist[cur_n1-1])]
                    ylist += [str(cur_index)]
                    xlist = xlist[0:cur_n1-1]
            
        random.shuffle(xlist)
        xlist=xlist[0:len(xlist)/2]
        random.shuffle(ylist)
        ylist=ylist[0:len(ylist)/2]
    
    cnode.mutlist1 = list(xlist)
    cnode.mutlist2 = list(ylist)
    deme_size = len(xlist)+len(ylist)
    while deme_size < 10000:
        n1,n2 = len(xlist),len(ylist)
        Ndivcells =  int(deme_size*division_prob) #number of dividing cells in this generation
        deme_size = Ndivcells*2
        if n2 > 0:
            ydown = int(division_prob*(1+sfit)*n2)
            yup = ydown+1
            plus_prob=division_prob*(1+sfit)*n2-ydown
            if random.random() < plus_prob:
                ynumber = int(yup)
            else:
                ynumber = int(ydown)
            if ynumber > 10000:
                ynumber = int(Ndivcells)
            ylist = random.sample(ylist,ynumber)*2
        else:
            ynumber = 0
            ylist = []
        xnumber = Ndivcells-ynumber
        if xnumber < 0:
            xlist = []
        else:
            xlist = random.sample(xlist,xnumber)*2
        n1,n2 = len(xlist),len(ylist)
        if total_blength < 80:
            if n1 > 0:
                new_mut1 = np.random.poisson(mu*n1)
                mut_assig1 = Counter(np.random.choice(n1,new_mut1))
                for x1 in mut_assig1.keys():
                    nmut = mut_assig1[x1]
                    mut_strlist = range(cur_index+1,cur_index+1+nmut)
                    mut_str = ",".join(map(str,mut_strlist))
                    #mut_str = str(range(cur_index+1,cur_index+1+nmut))
                    #mut_str = mut_str[1:len(mut_str)-1]
                    for xn in range(0,nmut):
                        cur_index += 1
                        lineage += [xlist[x1]]
                    xlist[x1] = mut_str
            if n2 > 0:
                new_mut2 = np.random.poisson(mu*n2)
                mut_assig2 = Counter(np.random.choice(n2,new_mut2))
                for x2 in mut_assig2.keys():
                    nmut = mut_assig2[x2]
                    mut_strlist = range(cur_index+1,cur_index+1+nmut)
                    mut_str = ",".join(map(str,mut_strlist))
                    #mut_str = str(range(cur_index+1,cur_index+1+nmut))
                    #mut_str = mut_str[1:len(mut_str)-1]
                    for xn in range(0,nmut):
                        cur_index += 1
                        lineage += [ylist[x2]]
                    ylist[x2] = mut_str
            if random.random() < adv_prob*n1:
                cur_index += 1
                cur_n1 = len(xlist)
                lineage += [str(xlist[cur_n1-1])]
                ylist += [str(cur_index)]
                xlist = xlist[0:cur_n1-1]
    
    random.shuffle(xlist)
    xlist1=xlist[0:len(xlist)/2]
    xlist2=xlist[len(xlist)/2:len(xlist)]
    random.shuffle(ylist)
    ylist1=ylist[0:len(ylist)/2]
    ylist2=ylist[len(ylist)/2:len(ylist)]
    
    return xlist1,xlist2,ylist1,ylist2,cur_index,lineage,cnode


def mutation_counter(mlist,mlineage,size_par,mean_depth,cutoff):
    """
    Find the pervasive mutations with high frequency in the sample (e.g. >50% of cells)
    """
    all_curindex = list(mlist)
    test_size = 10000
    test_index = random.sample(all_curindex,test_size)
    prob_par=size_par*1.0/(size_par+mean_depth)
    index_dist = Counter(test_index)
    all_mutindex = []
    for y in index_dist.keys():
        ylineage = traceLineage(mlineage,y)
        all_mutindex += ylineage*index_dist[y]
    mut_dist = Counter(all_mutindex)
    #testMAF1,testMAF2,testMAF3,testMAF4 = {},{},{},{}
    allMAF = {}
    for x in mut_dist.keys():
        test_alleleFreq = mut_dist[x]*0.5/test_size
        if test_alleleFreq > 0.005:
            site_depth = np.random.negative_binomial(size_par,prob_par)
            if site_depth >= 20:
                var_reads = np.random.binomial(site_depth,test_alleleFreq)
                site_maf = var_reads*1.0/site_depth
                if var_reads >= 3 and site_maf > cutoff:
                    allMAF[x] = round(site_maf,3)
                #if var_reads >= 3 and site_maf >= cutoff1:
                #    testMAF1[x] = site_maf
                #if var_reads >= 3 and site_maf >= cutoff2:
                #    testMAF2[x] = site_maf
                #if var_reads >= 3 and site_maf >= cutoff3:
                #    testMAF3[x] = site_maf
                #if var_reads >= 3 and site_maf >= cutoff4:
                #    testMAF4[x] = site_maf
    return allMAF
    #return testMAF1,testMAF2,testMAF3,testMAF4


def build_tree(instring):               # build the tree, return the root to the caller:
    k=0                                 # indexing on the list of usable nodes
    root=None                           # set root to empty
    index=1

    i=0
    upper=len(instring)
    while i < upper:
        if instring[i]== "(":
            if root == None:
                current=node()
                root=current
                k=k+1
                
                current.isroot=1
                current.index=index
                current.internal=1
                index=index+1
                previous=current

            else:
                if previous.left== None:
                    current=node()
                    k=k+1
                    previous.left=current
                    current.anc=previous
                    
                    current.index=index
                    index=index+1
                    current.isroot=0
                    current.internal=1

                    previous=current
                    
                elif previous.right==None:
                    current=node()
                    k=k+1
                    previous.right=current
                    current.anc=previous
                    
                    current.index=index
                    index=index+1
                    current.isroot=0
                    current.internal=1

                    previous=current
                    
                else:
                    print "both left and right pointer are full, weird tree format!(1)"
                    sys.exit(0)

            i=i+1
                    
        elif instring[i] not in [")", ",", ":", ";"]:
            if previous.left== None:
                current=node()
                k=k+1
                previous.left=current
                current.anc=previous
                
                current.index=index
                index=index+1
                current.isroot=0
                current.internal=0
                temp=instring[i:].find(":")
                current.name=instring[i:i+temp]

                i=i+temp

            elif previous.right==None:
                current=node()
                k=k+1
                previous.right=current
                current.anc=previous

                current.index=index
                index=index+1
                current.isroot=0
                current.internal=0
                temp=instring[i:].find(":")
                current.name=instring[i:i+temp]

                i=i+temp

            else:
                print "both left and right pointer are full, weird tree format!(2)"
                sys.exit(0)

        elif instring[i]==":":
            temp1=temp2=upper
            if ")" in instring[i:]:
                temp1=instring[i:].find(")")
            if "," in instring[i:]:
                temp2=instring[i:].find(",")
            if ";" in instring[i:]:
                temp3=instring[i:].find(";")
                
            temp=min(temp1,temp2,temp3)
            current.blength=float(instring[i+1:i+temp])
            i=i+temp

        elif instring[i]==",":
            i=i+1

        elif instring[i]==")":
            previous=previous.anc
            current=current.anc
            i=i+1
            
        elif instring[i]==";":
            if current == root:
                break;
            else:
                print "error while building tree"
                sys.exit(0)
            i=i+1

        else:
            print "met character(%s) which I don't expect"%instring[i]
            sys.exit(0)
        
    return root
    
def mut_simulatorMet_neutral(tmp_root,cmlist,clineage,current_index,dprob,pop1,pop2,pop3,pop4):
    if tmp_root.internal==1:
        left_list,right_list,current_index,clineage,tmp_root = deme_fission_neutral(tmp_root,cmlist,clineage,current_index,dprob)
        if tmp_root.left != None:
            mut_simulatorMet_neutral(tmp_root.left,left_list,clineage,current_index,dprob,pop1,pop2,pop3,pop4)
        if tmp_root.right != None:
            mut_simulatorMet_neutral(tmp_root.right,right_list,clineage,current_index,dprob,pop1,pop2,pop3,pop4)
    else:
        left_list,right_list,current_index,clineage,tmp_root = deme_fission_neutral(tmp_root,cmlist,clineage,current_index,dprob)
        if "A" in tmp_root.name:
            pop1 += tmp_root.mutlist1
        if "B" in tmp_root.name:
            pop2 += tmp_root.mutlist1
        if "C" in tmp_root.name:
            pop3 += tmp_root.mutlist1
        if "D" in tmp_root.name:
            pop4 += tmp_root.mutlist1

def mut_simulatorMet_selection(tmp_root,neu_list,sel_list,clineage,current_index,dprob,dcoef,pop1,pop2,pop3,pop4):
    if tmp_root.internal==1:
        left_neu_list,right_neu_list,left_sel_list,right_sel_list,current_index,clineage,tmp_root = deme_fission_selection(tmp_root,neu_list,sel_list,clineage,current_index,dprob,dcoef)
        if tmp_root.left != None:
            mut_simulatorMet_selection(tmp_root.left,left_neu_list,left_sel_list,clineage,current_index,dprob,dcoef,pop1,pop2,pop3,pop4)
        if tmp_root.right != None:
            mut_simulatorMet_selection(tmp_root.right,right_neu_list,right_sel_list,clineage,current_index,dprob,dcoef,pop1,pop2,pop3,pop4)
    else:
        left_neu_list,right_neu_list,left_sel_list,right_sel_list,current_index,clineage,tmp_root = deme_fission_selection(tmp_root,neu_list,sel_list,clineage,current_index,dprob,dcoef)
        if "A" in tmp_root.name:
            pop1 += tmp_root.mutlist1+tmp_root.mutlist2
        if "B" in tmp_root.name:
            pop2 += tmp_root.mutlist1+tmp_root.mutlist2
        if "C" in tmp_root.name:
            pop3 += tmp_root.mutlist1+tmp_root.mutlist2
        if "D" in tmp_root.name:
            pop4 += tmp_root.mutlist1+tmp_root.mutlist2


def mut_simulatorPri_neutral(tmp_root,cmlist,clineage,current_index,dprob,pop1,pop2,pop3,pop4,metpop):
    if tmp_root.internal==1:
        left_list,right_list,current_index,clineage,tmp_root = deme_fission_neutral(tmp_root,cmlist,clineage,current_index,dprob)
        if tmp_root.left != None:
            mut_simulatorPri_neutral(tmp_root.left,left_list,clineage,current_index,dprob,pop1,pop2,pop3,pop4,metpop)
        if tmp_root.right != None:
            mut_simulatorPri_neutral(tmp_root.right,right_list,clineage,current_index,dprob,pop1,pop2,pop3,pop4,metpop)
    else:
        left_list,right_list,current_index,clineage,tmp_root = deme_fission_neutral(tmp_root,cmlist,clineage,current_index,dprob)
        if "X" in tmp_root.name:
            pop1 += tmp_root.mutlist1
        if "Y" in tmp_root.name:
            pop2 += tmp_root.mutlist1
        if "Z" in tmp_root.name:
            pop3 += tmp_root.mutlist1
        if "W" in tmp_root.name:
            pop4 += tmp_root.mutlist1
        if "met" in tmp_root.name:
            metpop += tmp_root.mutlist1


def mut_simulatorPri_selection(tmp_root,neu_list,sel_list,clineage,current_index,dprob,dcoef,pop1,pop2,pop3,pop4,metpop):
    if tmp_root.internal==1:
        left_neu_list,right_neu_list,left_sel_list,right_sel_list,current_index,clineage,tmp_root = deme_fission_selection(tmp_root,neu_list,sel_list,clineage,current_index,dprob,dcoef)
        if tmp_root.left != None:
            mut_simulatorPri_selection(tmp_root.left,left_neu_list,left_sel_list,clineage,current_index,dprob,dcoef,pop1,pop2,pop3,pop4,metpop)
        if tmp_root.right != None:
            mut_simulatorPri_selection(tmp_root.right,right_neu_list,right_sel_list,clineage,current_index,dprob,dcoef,pop1,pop2,pop3,pop4,metpop)
    else:
        left_neu_list,right_neu_list,left_sel_list,right_sel_list,current_index,clineage,tmp_root = deme_fission_selection(tmp_root,neu_list,sel_list,clineage,current_index,dprob,dcoef)
        if "X" in tmp_root.name:
            pop1 += tmp_root.mutlist1+tmp_root.mutlist2
        if "Y" in tmp_root.name:
            pop2 += tmp_root.mutlist1+tmp_root.mutlist2
        if "Z" in tmp_root.name:
            pop3 += tmp_root.mutlist1+tmp_root.mutlist2
        if "W" in tmp_root.name:
            pop4 += tmp_root.mutlist1+tmp_root.mutlist2
        if "met" in tmp_root.name:
            metpop += tmp_root.mutlist1+tmp_root.mutlist2

def set_compare(dict1,dict2,cutoff):
    spec1 = sets.Set(dict1.keys())-sets.Set(dict2.keys())
    spec2 = sets.Set(dict2.keys())-sets.Set(dict1.keys())
    spec1_up,spec2_up = [],[]
    for x in spec1:
        if dict1[x] >= cutoff:
            spec1_up += [x]
    for y in spec2:
        if dict2[y] >= cutoff:
            spec2_up += [y]
    return len(spec1_up),len(spec2_up)


def pri_subclone(pridict,metdict,cutoff1,cutoff2):
    pri_subclonal = 0
    for mk in metdict.keys():
        if metdict[mk] >= cutoff1:
            if mk in pridict.keys():
                if pridict[mk] < cutoff2:
                    pri_subclonal += 1
    return pri_subclonal


############main script#################
mu = float(sys.argv[1])         #neutral mutation rate
met_time = float(sys.argv[2])   #dissemination timing - Nd=met_time*5000 (cells)
ptree = str(sys.argv[3])        #deme partition tree in primary tumor
mtree = str(sys.argv[4])        #deme partition tree in metastasis
repl = int(sys.argv[5])         #replicate

p_prob = 0.55                   #initial cell birth probability in primary tumor
m_prob = 0.55                   #initial cell birth probability in metastasis
cluster_size = 1                #seeding cell cluster size

scoef = 0.0                     #selection coefficient
#adv_prob = pow(10,-5)          #rate of advantageous mutations
if met_time > 1:
    met_time = int(met_time)

ss_file = open("SimParamStats_NN_u"+str(mu)+"_metTime"+str(met_time)+"_"+str(repl)+".txt","w")
ss_file.write("mu"+" "+"cluster.size"+" "+"met.time"+" "+"S1"+" "+"S2"+" "+"S3"+" "+"S4"+" "+"S5"+" "+"S6"+" "+"S7"+" "+"S8"+" "+"S9")
ss_file.write("\n")

pri_tree = open(ptree,'r')
met_tree = open(mtree,'r')

pri_tree_array = []
for tree in pri_tree:
    pri_tree_array += [tree]

met_tree_array = []
for tree in met_tree:
    met_tree_array += [tree]

for itree in range(0,len(pri_tree_array)):
    cur_ptree = pri_tree_array[itree]
    cur_mtree = met_tree_array[itree]
    
    mutlineage = ['0']
    mindex = 0
    flist = ['0'] #mutation index in the first deme
    sel_flist = []
    current_pop_size = 1
    #met_cluster = {} #metastasis ancestor clones
    
    flag = 0
    while current_pop_size < 5000: #5000 is the number of cells in a tumor deme
        Ndivcells =  int(current_pop_size*p_prob+1) #number of dividing cells in this generation
        flist = random.sample(flist,Ndivcells)*2 #the mutation distribution in the sampled cells.
        current_pop_size = Ndivcells*2
        new_mutation = np.random.poisson(mu*current_pop_size)
        mut_assignment = Counter(np.random.choice(current_pop_size,new_mutation))
        for x in mut_assignment.keys():
            mmut = mut_assignment[x]
            mut_strlist = range(mindex+1,mindex+1+mmut)
            mut_str = ",".join(map(str,mut_strlist))
            #mut_str = str(range(mindex+1,mindex+1+mmut))
            #mut_str = mut_str[1:len(mut_str)-1]
            for xi in range(0,mmut):
                mindex += 1
                mutlineage += [flist[x]]
            flist[x] = mut_str
        if met_time < 2 and flag==0 and current_pop_size > 1000:
            fmet = random.sample(flist,cluster_size)
            flag = 1
    
    pmut_index = int(mindex)
    
    proot = build_tree(cur_ptree)
    bulk1,bulk2,bulk3,bulk4,met_deme = [],[],[],[],[]

    mut_simulatorPri_neutral(proot,flist,mutlineage,pmut_index,p_prob,bulk1,bulk2,bulk3,bulk4,met_deme)
    #mut_simulatorPri_selection(proot,flist,sel_flist,mutlineage,pmut_index,p_prob,scoef,bulk1,bulk2,bulk3,bulk4,met_deme)

    pmut_index = len(mutlineage)-1 
    #print "final_mut_index=",pmut_index

    p_bulk_merge = bulk1+bulk2+bulk3+bulk4
    #p_bulk_merge = bulk1
    maf_primary = mutation_counter(p_bulk_merge,mutlineage,2,100,0.05)
    
    ####growth of met####
    if met_time >= 2:
        cluster_size = 1
        fmet = random.sample(met_deme,cluster_size)
    
    starting_cells = list(fmet)
    new_mut_index = int(pmut_index)
    #new_lineage = list(mutlineage)
    neu_mlist,mutlineage,new_mut_index = initiateFirsDeme_neutral(starting_cells,mutlineage,new_mut_index,m_prob)
    sel_mlist = []
    mroot = build_tree(cur_mtree)
    tissue1,tissue2,tissue3,tissue4 = [],[],[],[]
    mut_simulatorMet_neutral(mroot,neu_mlist,mutlineage,new_mut_index,m_prob,tissue1,tissue2,tissue3,tissue4)
    #mut_simulatorMet_selection(mroot,neu_mlist,sel_mlist,mutlineage,new_mut_index,m_prob,scoef,tissue1,tissue2,tissue3,tissue4)
   
    m_bulk_merge = tissue1+tissue2+tissue3+tissue4
    #m_bulk_merge = tissue1
    maf_met = mutation_counter(m_bulk_merge,mutlineage,2,100,0.05)
    
    S1,S5 = set_compare(maf_primary,maf_met,0.05)
    S2,S6 = set_compare(maf_primary,maf_met,0.1)
    S3,S7 = set_compare(maf_primary,maf_met,0.2)
    S4,S8 = set_compare(maf_primary,maf_met,0.3)
    S9 = pri_subclone(maf_primary,maf_met,0.3,0.3)
    
    ss_file.write(str(mu)+" "+str(cluster_size)+" "+str(int(math.log10(met_time*5000)))+" "+str(S1)+" "+str(S2)+" "+str(S3)+" "+str(S4)+" "+str(S5)+" "+str(S6)+" "+str(S7)+" "+str(S8)+" "+str(S9))
    ss_file.write("\n")
    
