#! /usr/bin/python

######################################################################################
## An agent-based framework to simulate colon cancer growth and metastasis          ##
## Spatial model: pripheral growth                                                  ##
## Evolution model: N/N, N/S, S/N, S/S                                              ##
## Output: Lp, Lm - number of primary-private and metastasis-private clonal SSNVs   ##
## Author: Zheng Hu                                                                 ##
## Date: 11/9/2017                                                                  ##
######################################################################################

import sys,os,math,random
import numpy as np
from collections import Counter
import sets

class deme():
    def __init__(self):
        self.present= 0     # whether the deme is empty or occupied: 0-empty;1-occupied
        self.neutral = []   # the background cells during tumor growth
        self.advant = []    # the advantageous cells having one driver mutation


def createLattice(d):
    """
    A function to create 3D cubic lattice with length of 2d where each site contains a empty deme
    """
    lattice = {}
    for x in range(0,2*d+1):
        for y in range(0,2*d+1):
            for z in range(0,2*d+1):
                lattice[(x,y,z)] = deme()
    return lattice


def neighbor26((a,b,c)):
    """
    Moore neighbourhood: 26 neighbour sites of (a,b,c)
    """
    neighbor = [(a+i,b+j,c+k)
                for i in [-1,0,1]
                for j in [-1,0,1]
                for k in [-1,0,1]
                if not (i==0 and j==0 and k==0)]

    return neighbor


def neighbor6((a,b,c)):
    """
    von Neumann neighbourhood: 6 neighbour sites of (a,b,c)
    """
    neighbor = [(a-1, b, c),(a+1, b, c),(a, b-1, c),(a, b+1, c),(a, b, c-1),(a, b, c+1)]
    return neighbor


def localNeighbor((a,b,c),r):
    """
    A function to search the local neighbour sites of (a,b,c) within an area of radius r in the 3D cubic lattice
    """
    neighbor = []
    for x in range(-r,r+1):
        for y in range(-r,r+1):
            for z in range(-r,r+1):
                if pow(x,2)+pow(y,2)+pow(z,2) < pow(r+1,2):
                    neighbor += [(a+x,b+y,c+z)]
    return neighbor


def traceLineage(mlineage,mutid):
    """
    A function to obtain the mutational lineage of a cell from the mutation id of the most recently occurred mutation in the cell. 
    For example, the input id (most recently occurred mutation) of target cell is "100" and the output is "1-12-35-56-100", which is the mutation lineage of the cell
    
    mlineage - the list that could help recover the mutational lineage given the most recent mutation id of a lineage
    mutid - the mutation ID of the most recently occurred mutation in the cell
    """
    recent_muts = mutid.split(',') # it is possible that multiple mutations occur during in a cell division. For instance, the mutation id of most recently occurred mutations is "100,101"
    recent_muts = [int(t) for t in recent_muts]
    first_mut = recent_muts[0] # the first mutation id in a multi-mutation event
    trace = []
    while first_mut > 0:
        trace += recent_muts
        recent_muts = mlineage[first_mut].split(',')
        recent_muts = [int(t) for t in recent_muts]
        first_mut = recent_muts[0]
    return trace

    
def lowerORupper(value):
    lower_int = int(value)
    upper_int = lower_int+1
    if random.random() < value-lower_int:
        return upper_int
    else:
        return lower_int


def initiateFirstDeme_neutral(maxsize,lineage,initial_id,current_id,birth_prob):
    """
    The growth of the initial deme via a neutral process

    maxsize - size limit of a deme
    lineage - a list that stores the lineage information of mutations
    """
    neu_list = [str(initial_id)]
    n1 = 1
    current_deme_size = 1
    while current_deme_size < maxsize:
        neu_divcells =  int(n1*birth_prob+1) #number of dividing cells in this generation
        neu_list = random.sample(neu_list,neu_divcells)*2
        n1 = len(neu_list)
        current_deme_size = n1

        new_mut1 = np.random.poisson(mut_rate*n1)
        mut_assig1 = Counter(np.random.choice(n1,new_mut1))
        for x1 in mut_assig1.keys():
            nmut = mut_assig1[x1]
            new_mut1 = range(current_id+1,current_id+1+nmut)
            mut_str = ",".join(map(str,new_mut1))
            #if nmut > 1:
            #    for t in new_mut1:
            #        multi_events[str(t)] = mut_str
            for xn in range(0,nmut):
                current_id += 1
                lineage += [neu_list[x1]]
            neu_list[x1] = mut_str
    
    return neu_list,current_id,lineage


def initiateFirstDeme_selection_v1(maxsize,lineage,initial_id,current_id,sfit,advrate,birth_prob,met_timing):
    """
    The growth of the initial deme via a selection process.

    maxsize - size limit of a deme
    lineage - a list that stores the lineage information of mutations
    """
    neu_list = [str(initial_id)]
    adv_list = []
    current_deme_size = 1
    flag=0
    met_founder = "x"
    while current_deme_size < maxsize:
        n1,n2 = len(neu_list),len(adv_list)
        neu_divcells =  int(n1*birth_prob+1) #number of dividing cells in this generation
        neu_list = random.sample(neu_list,neu_divcells)*2
        if n2 > 0:
            adv_divcells = lowerORupper(n2*birth_prob*(1+sfit)) #number of dividing cells in this generation        
            adv_list = random.sample(adv_list,adv_divcells)*2
        
        n1,n2 = len(neu_list),len(adv_list)
        current_deme_size = n1+n2
        if n1 > 0:
            new_mut1 = np.random.poisson(mut_rate*n1)
            mut_assig1 = Counter(np.random.choice(n1,new_mut1))
            for x1 in mut_assig1.keys():
                nmut = mut_assig1[x1]
                new_mut1 = range(current_id+1,current_id+1+nmut)
                mut_str = ",".join(map(str,new_mut1))
                #if nmut > 1:
                #    for t in new_mut1:
                #        multi_events[str(t)] = mut_str
                for xn in range(0,nmut):
                    current_id += 1
                    lineage += [neu_list[x1]]
                neu_list[x1] = mut_str
        if n2 > 0:
            new_mut2 = np.random.poisson(mut_rate*n2)
            mut_assig2 = Counter(np.random.choice(n2,new_mut2))
            for x2 in mut_assig2.keys():
                nmut = mut_assig2[x2]
                new_mut2 = range(current_id+1,current_id+1+nmut)
                mut_str = ",".join(map(str,new_mut2))
                #if nmut > 1:
                #    for t in new_mut2:
                #        multi_events[str(t)] = mut_str
                for xn in range(0,nmut):
                    current_id += 1
                    lineage += [adv_list[x2]]
                adv_list[x2] = mut_str
        
        if random.random() < advrate*n1:
            current_id += 1
            current_n1 = len(neu_list)
            lineage += [str(neu_list[current_n1-1])]
            adv_list += [str(current_id)]
            neu_list = neu_list[0:current_n1-1]
        
        if flag==0 and met_timing <= current_deme_size:
            ancestor_list = list(neu_list+adv_list)
            met_founder = random.choice(ancestor_list)
            flag=1
    
    return neu_list,adv_list,current_id,lineage,met_founder


def demeGrowthFission_selection_v1(neu_list,adv_list,lineage,current_id,current_size,sfit,advrate,birth_prob):
    """
    A function to simulate deme growth and fission and keep track of the mutational lineages
    
    """
    current_deme_size = len(neu_list)+len(adv_list)
    while current_deme_size < 2*deme_size:
        n1,n2 = len(neu_list),len(adv_list)
        neu_divcells =  lowerORupper(n1*birth_prob) #number of dividing cells in this generation
        neu_list = random.sample(neu_list,neu_divcells)*2
        if n2 > 0:
            adv_divcells =  lowerORupper(n2*birth_prob*(1+sfit)) #number of dividing cells in this generation
            adv_list = random.sample(adv_list,adv_divcells)*2
        
        n1,n2 = len(neu_list),len(adv_list)
        current_deme_size = n1+n2
        if current_size < pow(10,8)/deme_size: # stop mutation occurence when the tumor size is larger than 10^4*5000 = 5*10^7
            if n1 > 0:
                new_mut1 = np.random.poisson(mut_rate*n1)
                mut_assig1 = Counter(np.random.choice(n1,new_mut1))
                for x1 in mut_assig1.keys():
                    nmut = mut_assig1[x1]
                    new_mut1 = range(current_id+1,current_id+1+nmut)
                    mut_str = ",".join(map(str,new_mut1))
                    #if nmut > 1:
                    #    for t in new_mut1:
                    #        multi_events[str(t)] = mut_str
                    for xn in range(0,nmut):
                        current_id += 1
                        lineage += [neu_list[x1]]
                    neu_list[x1] = mut_str
            if n2 > 0:
                new_mut2 = np.random.poisson(mut_rate*n2)
                mut_assig2 = Counter(np.random.choice(n2,new_mut2))
                for x2 in mut_assig2.keys():
                    nmut = mut_assig2[x2]
                    new_mut2 = range(current_id+1,current_id+1+nmut)
                    mut_str = ",".join(map(str,new_mut2))
                    #if nmut > 1:
                    #    for t in new_mut2:
                    #        multi_events[str(t)] = mut_str
                    for xn in range(0,nmut):
                        current_id += 1
                        lineage += [adv_list[x2]]
                    adv_list[x2] = mut_str
            
            if random.random() < advrate*n1:
                current_id += 1
                current_n1 = len(neu_list)
                lineage += [str(neu_list[current_n1-1])]
                adv_list += [str(current_id)]
                neu_list = neu_list[0:current_n1-1]
        
    random.shuffle(neu_list)
    if len(neu_list) > 0:
        offspring_neu = np.random.binomial(len(neu_list),0.5)
    else:
        offspring_neu = 0
    neu_list1=neu_list[0:offspring_neu]
    neu_list2=neu_list[offspring_neu:len(neu_list)]
    
    random.shuffle(adv_list)
    if len(adv_list) > 0:
        offspring_adv = np.random.binomial(len(adv_list),0.5)
    else:
        offspring_adv = 0
    adv_list1=adv_list[0:offspring_adv]
    adv_list2=adv_list[offspring_adv:len(adv_list)]
    
    return neu_list1,neu_list2,adv_list1,adv_list2,current_id,lineage


def demeGrowthFissionSel2_v1(neu_list,adv_list,lineage,current_id,current_size,sfit,advrate,birth_prob):
    """
    A function to simulate deme growth and fission and keep track of the mutational lineages
    
    assumptions:
    (1) the growth of neutral and advantageous clones are independent
    """
    current_deme_size = len(neu_list)+len(adv_list)
    while current_deme_size < 2*deme_size:
        n1,n2 = len(neu_list),len(adv_list)
        neu_divcells =  lowerORupper(n1*birth_prob) #number of dividing cells in this generation
        neu_list = random.sample(neu_list,neu_divcells)*2
        if n2 > 0:
            adv_divcells =  lowerORupper(n2*birth_prob*(1+sfit)) #number of dividing cells in this generation
            adv_list = random.sample(adv_list,adv_divcells)*2
        
        n1,n2 = len(neu_list),len(adv_list)
        current_deme_size = n1+n2
        if current_size < 5*pow(10,7)/deme_size: # stop mutation occurence when the tumor size is larger than 10^4*5000 = 5*10^7
            if n1 > 0:
                new_mut1 = np.random.poisson(mut_rate*n1)
                mut_assig1 = Counter(np.random.choice(n1,new_mut1))
                for x1 in mut_assig1.keys():
                    nmut = mut_assig1[x1]
                    new_mut1 = range(current_id+1,current_id+1+nmut)
                    mut_str = ",".join(map(str,new_mut1))
                    #if nmut > 1:
                    #    for t in new_mut1:
                    #        multi_events[str(t)] = mut_str
                    for xn in range(0,nmut):
                        current_id += 1
                        lineage += [neu_list[x1]]
                    neu_list[x1] = mut_str
            if n2 > 0:
                new_mut2 = np.random.poisson(mut_rate*n2)
                mut_assig2 = Counter(np.random.choice(n2,new_mut2))
                for x2 in mut_assig2.keys():
                    nmut = mut_assig2[x2]
                    new_mut2 = range(current_id+1,current_id+1+nmut)
                    mut_str = ",".join(map(str,new_mut2))
                    #if nmut > 1:
                    #    for t in new_mut2:
                    #        multi_events[str(t)] = mut_str
                    for xn in range(0,nmut):
                        current_id += 1
                        lineage += [adv_list[x2]]
                    adv_list[x2] = mut_str
            
            if random.random() < advrate*n1:
                current_id += 1
                current_n1 = len(neu_list)
                lineage += [str(neu_list[current_n1-1])]
                adv_list += [str(current_id)]
                neu_list = neu_list[0:current_n1-1]
        
    random.shuffle(neu_list)
    if len(neu_list) > 0:
        offspring_neu = np.random.binomial(len(neu_list),0.5)
    else:
        offspring_neu = 0
    neu_list1=neu_list[0:offspring_neu]
    neu_list2=neu_list[offspring_neu:len(neu_list)]
    
    random.shuffle(adv_list)
    if len(adv_list) > 0:
        offspring_adv = np.random.binomial(len(adv_list),0.5)
    else:
        offspring_adv = 0
    adv_list1=adv_list[0:offspring_adv]
    adv_list2=adv_list[offspring_adv:len(adv_list)]
    
    return neu_list1,neu_list2,adv_list1,adv_list2,current_id,lineage


def seqProcessing(sp,sample_keys,mlineage,size_par,mean_depth,purity):
    """
    Model the random sampling process in next-generation sequencing and report the sequencing allele frequencies in a sample of cells
    
    sp- the lattice space constituting the population of demes
    sample_keys- the sampled demes
    size_par- variance parameter in negative-binomial distribution
    mean_depth- the mean depth of the sequencing
    purity- tumor purity
    """
    all_cur_id = []
    all_mut_id = []
    for key in sample_keys:
        smuts = list(sp[key].neutral + sp[key].advant)
        all_cur_id += smuts
    random_sample = 10000 ## the number of cells for sequencing analysis
    sample_id = random.sample(all_cur_id,random_sample)
    id_count = Counter(sample_id)
    for x in id_count.keys():
        xlineage = traceLineage(mlineage,x)
        all_mut_id += xlineage*id_count[x]
    mut_count = Counter(all_mut_id)
    prob_par=size_par*1.0/(size_par+mean_depth)
    sampleAF = {}
    for x in mut_count.keys():
        true_af = mut_count[x]*0.5*purity/random_sample
        if true_af > 0.005:
            site_depth = np.random.negative_binomial(size_par,prob_par)
            if site_depth >= 20:
                var_reads = np.random.binomial(site_depth,true_af)
                seq_af = var_reads*1.0/site_depth
                if var_reads >= 3:
                    #sampleAF[str(x)] = (site_depth,seq_af)
                    sampleAF[str(x)] = seq_af
    return sampleAF


def highMutCCF(sp,sample_keys,mlineage,cutoff):
    """
    The mutations with cff larger than a cutoff
    """
    all_cur_id = []
    all_mut_id = []
    for key in sample_keys:
        smuts = list(sp[key].neutral + sp[key].advant)
        all_cur_id += smuts
    random_sample = 200000
    #random_sample = len(all_cur_id) ## the number of cells for sequencing analysis
    sample_id = random.sample(all_cur_id,random_sample)
    #sample_id = all_cur_id
    id_count = Counter(sample_id)
    for x in id_count.keys():
        xlineage = traceLineage(mlineage,x)
        all_mut_id += xlineage*id_count[x]
    mut_count = Counter(all_mut_id)
    sampleCCF = {}
    for x in mut_count.keys():
        true_ccf = mut_count[x]*1.0/random_sample
        if true_ccf > cutoff:
            sampleCCF[str(x)] = round(true_ccf,4)
    return sampleCCF

def pubMutGenerator(n,size_par,mean_depth,purity):
    """
    A function to generate the public mutations (number is "n") and frequencies
    
    n- number of samples
    size_par- variation parameter in the negative binomial distribution
    mean_death- mean seq depth
    """
    prob_par=size_par*1.0/(size_par+mean_depth)
    mean_af = 0.5*purity
    depth_pub = []
    maf_pub = []
    for k in range(0,n):
        correct = 0
        while correct == 0:
            site_depth = np.random.negative_binomial(size_par,prob_par)
            if site_depth >= 15:
                correct =1
        var_reads = np.random.binomial(site_depth,mean_af)
        site_maf = var_reads*1.0/site_depth
        depth_pub += [site_depth]
        maf_pub += [site_maf]
    return depth_pub,maf_pub


def locationSampling(region,sample_number,cutoff):
    """
    A function to sampling the locations where the bulk tissue will be sampled locally.
    """
    success = 0
    while success == 0:
        locations = random.sample(region,sample_number)
        repeat = sample_number*(sample_number-1)
        minall = 999
        for x in range(0,repeat):
            rs = random.sample(locations,2)
            min_distance = min([abs(rs[0][0]-rs[1][0]),abs(rs[0][1]-rs[1][1]),abs(rs[0][2]-rs[1][2])])
            if min_distance < minall:
                minall = min_distance
        if min_distance > 2*cutoff:
            success = 1
    return locations


def bulkTissueSampling(sp,location,radius):
    """
    Check if the sampled sites have non-empty demes
    """
    local_region = localNeighbor(location,radius)
    bulk_tissue = []
    for x in local_region:
        if sp[x].present == 1:
            bulk_tissue += [x]
    return bulk_tissue


#########################main script to simulate a primary-metastasis pair###############################
repl = int(sys.argv[1])         #simulation replicate    
model = str(sys.argv[2])        #model - NN, NS. SN, or SS


rd = 60                         #radius of the 3D space
deme_size = 5000
final_tumor_size = pow(10,8)    #final tumor size in cell number
final_deme_number = final_tumor_size/deme_size  #number of demes in final tumor
mut_rate = 0.3                  #mutation rate per cell division in exonic regions
p_birth_rate = 0.55             #cell birth probability for primary tumor
m_birth_rate = 0.55             #cell birth probability for metastatic tumor
#p_birth_rate = round(random.uniform(0.54,0.6),4)
#m_birth_rate = round(random.uniform(0.54,0.6),4)
#p_scoef = 0.1
#m_scoef = 0.1
#p_scoef = round(random.uniform(0.01,0.2),4)
#m_scoef = round(random.uniform(0.01,0.2),4)
#mut_rate = round(random.uniform(0.1,1),4)

#p_adv_rate = pow(10,-5)         #rate of advantageous mutations in primary tumor
#m_adv_rate = pow(10,-5)         #rate of advantageous mutations in metastasis

td = pow(10,random.uniform(2,8)) ##primary tumor size at the time of dissemination

if model == "NN" or model == "NS":
    p_adv_rate = 0              #rate of advantageous mutations in primary tumor
    p_scoef = 0                 #selection coefficient in primary tumor
else:
    p_adv_rate = pow(10,-5)
    p_scoef = 0.1

if model == "NN" or model == "SN":
    m_adv_rate = 0              #rate of advantageous mutations in metastasis
    m_scoef = 0                 #selection coefficient in metastasis
else:
    m_adv_rate = pow(10,-5)
    m_scoef = 0.1

mut_lineage = ['0']             # the lineage tracer
mut_index = 0

neu_cells,adv_cells,mut_index,mut_lineage,met_progenitor = initiateFirstDeme_selection_v1(deme_size*2,mut_lineage,0,mut_index,p_scoef,p_adv_rate,p_birth_rate,td) #the growth of the fisrt deme from single transformed cell

random.shuffle(neu_cells)
if len(neu_cells) > 0:
    offspring_neu_number = np.random.binomial(len(neu_cells),0.5)
else:
    offspring_neu_number = 0
neu_cells1=neu_cells[0:offspring_neu_number]
neu_cells2=neu_cells[offspring_neu_number:len(neu_cells)]

random.shuffle(adv_cells)
if len(adv_cells) > 0:
    offspring_adv_number = np.random.binomial(len(adv_cells),0.5)
else:
    offspring_adv_number = 0
adv_cells1=adv_cells[0:offspring_adv_number]
adv_cells2=adv_cells[offspring_adv_number:len(adv_cells)]

first_neighbors = neighbor26((rd,rd,rd))
second_deme_site = random.choice(first_neighbors)

################
space = createLattice(rd)
space[(rd,rd,rd)].present = 1
space[(rd,rd,rd)].neutral = list(neu_cells1)
space[(rd,rd,rd)].advant = list(adv_cells1)
space[second_deme_site].present = 1
space[second_deme_site].neutral = list(neu_cells2)
space[second_deme_site].advant = list(adv_cells2)

current_keys = [(rd,rd,rd),second_deme_site]
surface_keys = [(rd,rd,rd),second_deme_site]
current_deme_number =2 #current tumor size by # of demes
surface_deme_number =2 #current tumor size by # of demes
generation = 1

while current_deme_number < final_deme_number:
    new_keys = []
    for w in range(0,surface_deme_number):
        skey = random.choice(surface_keys)
        if space[skey].present == 1:
            nei_sites = neighbor26(skey)
            empty_neis = [key for key in nei_sites if space[key].present == 0]                    # empty neighbor sites
            
            if len(empty_neis) > 0:
                rand_prob = random.random()
                if rand_prob < 1-math.exp(-len(empty_neis)*0.25):   #when the deme will be chosen for growth and fission
                    pre_neu = list(space[skey].neutral)
                    pre_adv = list(space[skey].advant)
                    post_neu_1,post_neu_2,post_adv_1,post_adv_2,mut_index,mut_lineage = demeGrowthFission_selection_v1(pre_neu,pre_adv,mut_lineage,mut_index,current_deme_number,p_scoef,p_adv_rate,p_birth_rate)
                    newkey = random.choice(empty_neis)
                    new_keys += [newkey]
                    space[skey].neutral = list(post_neu_1)
                    space[newkey].neutral = list(post_neu_2)
                    space[skey].advant = list(post_adv_1)
                    space[newkey].advant = list(post_adv_2)
                    space[newkey].present = 1
                    current_keys += [newkey]
                    current_deme_number += 1

                    if met_progenitor == "x" and current_deme_number == int(td/deme_size+1):
                        anc_list = list(post_neu_2+post_adv_2)
                        met_progenitor = random.choice(anc_list)
    
    ###update surface
    surface_update = list(surface_keys+new_keys)
    surface_keys = []
    for fkey in surface_update:
        neisites = neighbor26(fkey)
        random.shuffle(neisites)
        for key in neisites:
            if space[key].present == 0:
                surface_keys += [fkey]
                break
    surface_deme_number = len(surface_keys)
    generation += 1
    

#####################################
####growth of met###
####################################

#m_neu_cells,mut_index,mut_lineage = initiateFirstDeme_neutral(deme_size,mut_lineage,met_progenitor,mut_index,m_birth_rate)  #the growth of the fisrt deme from single transformed cell
tmp_td=100  #parameter of dissemination timing, not  used in metastasis growth; for for function input
m_neu_cells,m_adv_cells,mut_index,mut_lineage,tmp_progenitor = initiateFirstDeme_selection_v1(deme_size,mut_lineage,met_progenitor,mut_index,m_scoef,m_adv_rate,m_birth_rate,tmp_td)

mspace = createLattice(rd)
mspace[(rd,rd,rd)].present = 1
mspace[(rd,rd,rd)].neutral = list(m_neu_cells)
mspace[(rd,rd,rd)].advant = list(m_adv_cells)

met_current_keys = [(rd,rd,rd)]
met_surface_keys = [(rd,rd,rd)]
met_current_deme_number =1                                                 #current number of demes
met_surface_deme_number =1                                                 #current number of demes

        
while met_current_deme_number < final_deme_number:
    met_new_keys = []
    for w in range(0,met_surface_deme_number):
        skey = random.choice(met_surface_keys)
        if mspace[skey].present == 1:
            nei_sites = neighbor26(skey)
            empty_neis = [key for key in nei_sites if mspace[key].present == 0]                    # empty neighbor sites
            
            if len(empty_neis) > 0:
                rand_prob = random.random()
                if rand_prob < 1-math.exp(-len(empty_neis)*0.25):   #when the deme will be chosen for growth and fission
                    pre_neu = list(mspace[skey].neutral)
                    pre_adv = list(mspace[skey].advant)
                    post_neu_1,post_neu_2,post_adv_1,post_adv_2,mut_index,mut_lineage = demeGrowthFission_selection_v1(pre_neu,pre_adv,mut_lineage,mut_index,met_current_deme_number,m_scoef,m_adv_rate,m_birth_rate)
                    newkey = random.choice(empty_neis)
                    met_new_keys += [newkey]
                    mspace[skey].neutral = list(post_neu_1)
                    mspace[newkey].neutral = list(post_neu_2)
                    mspace[skey].advant = list(post_adv_1)
                    mspace[newkey].advant = list(post_adv_2)
                    mspace[newkey].present = 1
                    met_current_keys += [newkey]
                    met_current_deme_number += 1
    
    ###update surface
    met_surface_update = list(met_surface_keys + met_new_keys)
    met_surface_keys = []
    for fkey in met_surface_update:
        neisites = neighbor26(fkey)
        random.shuffle(neisites)
        for key in neisites:
            if mspace[key].present == 0:
                met_surface_keys += [fkey]
                break
    met_surface_deme_number = len(met_surface_keys)
        

#print round(math.log10(td),3),mut_rate,p_scoef,m_scoef,p_birth_rate,m_birth_rate,len(p_high_muts),len(m_high_muts)

primet = open("TumorSimul3D_LpLm_"+str(model)+"_"+str(repl)+".txt","w")
primet.write("mut_id"+" "+"pccf"+" "+"mccf")
primet.write("\n")

p_samples = random.sample(current_keys,200)
m_samples = random.sample(met_current_keys,200)

p_high_muts = highMutCCF(space,p_samples,mut_lineage,0.01)
m_high_muts = highMutCCF(mspace,m_samples,mut_lineage,0.01)

#print math.log10(td*5000),mut_rate,p_scoef,m_scoef,p_birth_rate,m_birth_rate

all_muts = sets.Set(p_high_muts.keys()) | sets.Set(m_high_muts.keys())
pm_shared = sets.Set(p_high_muts.keys()) & sets.Set(m_high_muts.keys())
p_specific = all_muts - sets.Set(m_high_muts.keys())
m_specific = all_muts - sets.Set(p_high_muts.keys())
p_specific_freq = map(lambda x: p_high_muts[x], p_specific)
m_specific_freq = map(lambda x: m_high_muts[x], m_specific)
lp_count = len(list(filter(lambda x: x > 0.6, p_specific_freq)))
lm_count = len(list(filter(lambda x: x > 0.6, m_specific_freq)))
h_ratio = lm_count*1.0/(lp_count+1)

for m in pm_shared:
    primet.write(str(m)+" "+str(p_high_muts[m])+" "+str(m_high_muts[m]))
    primet.write("\n")

for m in p_specific:
    primet.write(str(m)+" "+str(p_high_muts[m])+" "+"0")
    primet.write("\n")

for m in m_specific:
    primet.write(str(m)+" "+"0"+" "+str(m_high_muts[m]))
    primet.write("\n")

out_para = open("LmLp_output_"+str(model)+"_"+str(repl)+".txt","w")
out_para.write(str(repl)+" "+str(round(math.log10(td),3))+" "+str(mut_rate)+" "+str(p_scoef)+" "+str(m_scoef)+" "+str(p_birth_rate)+" "+str(m_birth_rate)+" "+str(lp_count)+" "+str(lm_count)+" "+str(h_ratio))
out_para.write("\n")

