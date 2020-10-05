# Import all necessary packages

import numpy as np
import csv
import itertools
import re
import pandas as pd
import warnings
from Bio.SeqUtils import MeltingTemp as mt
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.SeqFeature import SeqFeature, FeatureLocation
from Bio.Seq import Seq

#The following is a copyright notice by the authors of iBioCAD, from which
#substantial portions of this code were adopted and modified:

#Copyright (c) 2019 Scott Weisberg

#Permission is hereby granted, free of charge, to any person obtaining a copy
#of this software and associated documentation files (the "Software"), to deal
#in the Software without restriction, including without limitation the rights
#to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
#copies of the Software, and to permit persons to whom the Software is
#furnished to do so, subject to the following conditions:

#The above copyright notice and this permission notice shall be included in all
#copies or substantial portions of the Software.

#THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
#IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
#FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
#AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
#LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
#OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
#SOFTWARE.


#Define infrastructure for the following code
class Part:
    def __init__(self,name,type,sequence):
        self.name = name
        self.sequence = sequence
        self.type = type
    primer_forward = ""
    primer_reverse = ""
    bridge_with_next_part = ""
    bridge_with_previous_part = ""
    primer_forward_tm = 0
    primer_reverse_tm = 0
    
def builds(parts_list):
    import copy
    builds_list = []
    builds = 1
    for i in range(builds):
        build = []
        for part in parts_list:
            build.append(copy.copy(part))
        builds_list.append(build)
    return builds_list


# Define functions to ease the main computation
def reverse_complement(sequence):
    rev_comp = ""
    Watson_Crick = {"A":"T","C":"G","T":"A","G":"C","a":"t","t":"a","c":"g","g":"c"}
    for base in sequence:
        rev_comp = Watson_Crick[base] + rev_comp
    return rev_comp

def complement(sequence):
    comp = ''
    Watson_Crick = {"A":"T","C":"G","T":"A","G":"C","a":"t","t":"a","c":"g","g":"c"}
    for base in sequence:
        comp += Watson_Crick[base]
    return comp

def reverse(sequence):
    rev = ''
    for base in sequence:
        rev = base + rev
    return rev

def flattn(list_2d):
    flattnd = []
    for subl in list_2d:
        if type(subl) is list:
            flattnd += [s for s in subl]
        else:
            flattnd.append(subl)
    return flattnd


#Optimizes the overhangs used in Golden Gate assembly
def golden_gate_optimization(parts_list,free_overhangsets):

    # Write all variable sequences in same order in a new list
    oh_list = []
    for x in parts_list:
        if x.type == 'pegRNA':
            oh_list.append([x.sequence[:20],x.sequence[96:]])
        elif x.type == 'gRNA':
            oh_list.append(x.sequence[:20])
        elif x.type == 'smRNA':
            oh_list.append(x.sequence[:-len(tRNA)])
        else: # if part is tRNA
            oh_list.append('plc')
    
    # Starting in the middle of the variable sequences and moving outwards, find overhang combinations
    for cov in range(2,max([int(np.ceil(np.true_divide(len(prt),2))) for prt in flattn(oh_list)])):
        
        cov_list = []
        for oh in oh_list:
            if type(oh) is list:
                o_list = []
                for o in oh:
                    if 2*cov >= len(o):
                        o_list.append(o)
                    else:
                        o_list.append(o[int(np.ceil(np.true_divide(len(o),2)))-cov:int(np.ceil(np.true_divide(len(o),2)))+cov])
                cov_list.append(o_list)
            else:
                if 2*cov >= len(oh):
                    cov_list.append(oh)
                else:
                    cov_list.append(oh[int(np.ceil(np.true_divide(len(oh),2)))-cov:int(np.ceil(np.true_divide(len(oh),2)))+cov])

        # Find possible overhang combinations
        for s in free_overhangsets:
            golden_gate_overhangs=s
            seq_matches = []
            for x in range(len(parts_list)-1):
                seq_matches.append([])
                for forwoverhang in golden_gate_overhangs:
                    overhangs = [forwoverhang, reverse_complement(forwoverhang)]
                    for overhang in overhangs:
                        if parts_list[x].type == 'pegRNA':
                            if overhang in cov_list[x][1]:
                                seq_matches[x].append(oh_list[x][1][:oh_list[x][1].find(overhang)+4])
                                break
                        elif parts_list[x].type == 'gRNA':
                            if parts_list[x+1].type == 'gRNA':
                                if overhang in cov_list[x+1]:
                                    seq_matches[x].append(oh_list[x+1][:oh_list[x+1].find(overhang)+4])
                                    break
                            elif parts_list[x+1].type == 'pegRNA':
                                if overhang in cov_list[x+1][0]:
                                    seq_matches[x].append(oh_list[x+1][0][:oh_list[x+1][0].find(overhang)+4])
                                    break
                            elif parts_list[x+1].type == 'smRNA':
                                if overhang in cov_list[x+1]:
                                    seq_matches[x].append(oh_list[x+1][:oh_list[x+1].find(overhang)+4])
                                    break
                        elif parts_list[x].type == 'smRNA':
                            if parts_list[x+1].type == 'gRNA':
                                if overhang in cov_list[x+1]:
                                    seq_matches[x].append(oh_list[x+1][:oh_list[x+1].find(overhang)+4])
                                    break
                            elif parts_list[x+1].type == 'pegRNA':
                                if overhang in cov_list[x+1][0]:
                                    seq_matches[x].append(oh_list[x+1][0][:oh_list[x+1][0].find(overhang)+4])
                                    break
                            elif parts_list[x+1].type == 'smRNA':
                                if overhang in cov_list[x+1]:
                                    seq_matches[x].append(oh_list[x+1][:oh_list[x+1].find(overhang)+4])
                                    break
                        else: # if part is tRNA
                            if parts_list[x+1].type == 'gRNA':
                                if overhang in cov_list[x+1]:
                                    seq_matches[x].append(oh_list[x+1][:oh_list[x+1].find(overhang)+4])
                                    break
                            elif parts_list[x+1].type == 'pegRNA':
                                if overhang in cov_list[x+1][0]:
                                    seq_matches[x].append(oh_list[x+1][0][:oh_list[x+1][0].find(overhang)+4])
                                    break
                            elif parts_list[x+1].type == 'smRNA':
                                if overhang in cov_list[x+1]:
                                    seq_matches[x].append(oh_list[x+1][:oh_list[x+1].find(overhang)+4])
                                    break


            combs = []
            for x in itertools.product(*seq_matches):
                combs.append(x)
            for comb in combs:
                if len([c[-4:] for c in comb]) == len(set([c[-4:] for c in comb])):
                    return comb
    #if there are no possible combinations
    return None


# Set template sequences tRNA and gRNA scaffold
tRNA = 'AACAAAGCACCAGTGGTCTAGTGGTAGAATAGTACCCTGCCACGGTACAGACCCGGGTTCGATTCCCGGCTGGTGCA'
scaffld = 'GTTTTAGAGCTAGAAATAGCAAGTTAAAATAAGGCTAGTCCGTTATCAACTTGAAAAAGTGGCACCGAGTCGGTGC'


# Perform scarless Golden Gate assembly computation with provided parts
def scarless_gg(parts_list, primer_tm_range, max_annealing_len, bb_overlaps, additional_overhangs):
    msg = None
    bb_overlaps = [i.lower() for i in bb_overlaps]
    additional_overhangs = [i.lower() for i in additional_overhangs]
    
    # Format parts list
    builds_list = builds(parts_list)
    
    # Go through parts and write all known annotations into list
    new_builds_list = []
    ftrs = []
    for unpacked_list in builds_list:
        mmry = 13
        for part in unpacked_list:
            part.sequence = part.sequence.lower()
            ftrs.append(SeqFeature(FeatureLocation(mmry, mmry+len(part.sequence), strand=1), type=part.name))
            if part.type == 'pegRNA':
                ftrs.append(SeqFeature(FeatureLocation(mmry, mmry+20, strand=1), type='spacer'))
                ftrs.append(SeqFeature(FeatureLocation(mmry+96, mmry+len(part.sequence)-13, strand=1), type='RT template'))
                ftrs.append(SeqFeature(FeatureLocation(mmry+len(part.sequence)-13, mmry+len(part.sequence), strand=1), type='PBS'))
            elif part.type == 'gRNA':
                ftrs.append(SeqFeature(FeatureLocation(mmry, mmry+20, strand=1), type='spacer'))
                ftrs.append(SeqFeature(FeatureLocation(mmry+20+len(scaffld), mmry+len(part.sequence), strand=1), type='tRNA'))
            elif part.type == 'smRNA':
                ftrs.append(SeqFeature(FeatureLocation(mmry, mmry+len(part.sequence)-len(tRNA), strand=1), type='smRNA'))
                ftrs.append(SeqFeature(FeatureLocation(mmry+len(part.sequence)-len(tRNA), mmry+len(part.sequence), strand=1), type='tRNA'))
            mmry += len(part.sequence)
        
        # Iterate through overhang sets with increasing size until fitting one is found
        gg_opt = None
        breakit = False
        for p in range(10,51):
            free_overhangsets = []
            for q in range(5):
                with open('overhangsets/setsof%s.csv'%p,'r') as f:
                    reader = csv.reader(f, delimiter=",")
                    sets = list(reader)[1:]
                
                temp = []
                for s in sets:
                    if len(s) != 0:
                        temp.append(s)

                forwexisting = additional_overhangs+bb_overlaps
                bothexisting = forwexisting+[reverse_complement(i) for i in forwexisting]
                forwset = [x.lower() for x in temp[q]]
                bothset = forwset+[reverse_complement(i) for i in forwset] # Consider also reverse complement of linkers
                # Only grab sets that include all existing overhangs and delete the existing from the set
                if all(i in bothset for i in forwexisting):
                    print(forwset)
                    free_overhangsets.append([i for i in forwset if i not in bothexisting])
                else:
                    continue
            
            if free_overhangsets:
                gg_opt = golden_gate_optimization(unpacked_list,free_overhangsets)
            if gg_opt is not None:
                breakit = True
            if breakit:
                break
        
        # No sets include all existing overhangs
        if gg_opt is None:
            msg = 'comb_warn'
            warnings.warn('The given combination of existing overhangs is not compatible with an optimal overhang set. '
                          'Computing the best overhang set not including the given existing overhangs. There might be '
                          'interference between overhangs.')
            
            for p in range(10,51):
                free_overhangsets = []
                for q in range(5):
                    with open('overhangsets/setsof%s.csv'%p,'r') as f:
                        reader = csv.reader(f, delimiter=",")
                        sets = list(reader)[1:]

                    temp = []
                    for s in sets:
                        if len(s) != 0:
                            temp.append(s)

                gg_opt = golden_gate_optimization(unpacked_list,temp)
                if gg_opt is not None:
                    breakit = True
                if breakit:
                    break
        
        if gg_opt is None:
            msg = 'comb_error'
            return None,None,msg


        #Modify sequences and design primers
        for i in range(len(unpacked_list)):
            
            # If current part is first part, define forward primer with left backbone overlap and reverse primer ordinarily
            if i == 0:
                unpacked_list[i].primer_forward = 'taggtctcc' + reverse_complement(bb_overlaps[0]) + unpacked_list[i].sequence[:max_annealing_len]
                unpacked_list[i].sequence = 'taggtctcc' + reverse_complement(bb_overlaps[0]) + unpacked_list[i].sequence
                if unpacked_list[i].type == 'pegRNA':
                    unpacked_list[i].primer_reverse = reverse_complement(unpacked_list[i].sequence[unpacked_list[i].sequence.find(scaffld.lower())+76-max_annealing_len:unpacked_list[i].sequence.find(scaffld.lower())+76+len(gg_opt[i])] + "tgagacccg")
                    unpacked_list[i+1].primer_forward = "taggtctcc" + unpacked_list[i].sequence[unpacked_list[i].sequence.find(scaffld.lower())+76+len(gg_opt[i])-4:] + unpacked_list[i+1].sequence[:max_annealing_len]
                    unpacked_list[i+1].sequence = "taggtctcc" + unpacked_list[i].sequence[unpacked_list[i].sequence.find(scaffld.lower())+76+len(gg_opt[i])-4:] + unpacked_list[i+1].sequence
                    unpacked_list[i].sequence = unpacked_list[i].sequence[:unpacked_list[i].sequence.find(scaffld.lower())+76+len(gg_opt[i])] + "tgagacccg"
                elif unpacked_list[i].type == 'gRNA':
                    unpacked_list[i].primer_reverse = reverse_complement(unpacked_list[i].sequence[-max_annealing_len:] + unpacked_list[i+1].sequence[:unpacked_list[i+1].sequence.find(gg_opt[i])+len(gg_opt[i])] + "tgagacccg")
                    unpacked_list[i].sequence = unpacked_list[i].sequence + unpacked_list[i+1].sequence[:unpacked_list[i+1].sequence.find(gg_opt[i])+len(gg_opt[i])] + "tgagacccg"
                    if unpacked_list[i+1].type == 'smRNA':
                        unpacked_list[i+1].primer_forward = "taggtctcc" + unpacked_list[i+1].sequence[unpacked_list[i+1].sequence.find(gg_opt[i])+len(gg_opt[i])-4:unpacked_list[i+1].sequence.find(tRNA.lower())+max_annealing_len]
                        unpacked_list[i+1].sequence = "taggtctcc" + unpacked_list[i+1].sequence[unpacked_list[i+1].sequence.find(gg_opt[i])+len(gg_opt[i])-4:]
                    else:
                        unpacked_list[i+1].primer_forward = "taggtctcc" + unpacked_list[i+1].sequence[unpacked_list[i+1].sequence.find(gg_opt[i])+len(gg_opt[i])-4:unpacked_list[i+1].sequence.find(scaffld.lower())+max_annealing_len]
                        unpacked_list[i+1].sequence = "taggtctcc" + unpacked_list[i+1].sequence[unpacked_list[i+1].sequence.find(gg_opt[i])+len(gg_opt[i])-4:]
                elif unpacked_list[i].type == 'smRNA':
                    unpacked_list[i].primer_reverse = reverse_complement(unpacked_list[i].sequence[-max_annealing_len:] + unpacked_list[i+1].sequence[:unpacked_list[i+1].sequence.find(gg_opt[i])+len(gg_opt[i])] + "tgagacccg")
                    unpacked_list[i].sequence = unpacked_list[i].sequence + unpacked_list[i+1].sequence[:unpacked_list[i+1].sequence.find(gg_opt[i])+len(gg_opt[i])] + "tgagacccg"
                    if unpacked_list[i+1].type == 'smRNA':
                        unpacked_list[i+1].primer_forward = "taggtctcc" + unpacked_list[i+1].sequence[unpacked_list[i+1].sequence.find(gg_opt[i])+len(gg_opt[i])-4:unpacked_list[i+1].sequence.find(tRNA.lower())+max_annealing_len]
                        unpacked_list[i+1].sequence = "taggtctcc" + unpacked_list[i+1].sequence[unpacked_list[i+1].sequence.find(gg_opt[i])+len(gg_opt[i])-4:]
                    else:
                        unpacked_list[i+1].primer_forward = "taggtctcc" + unpacked_list[i+1].sequence[unpacked_list[i+1].sequence.find(gg_opt[i])+len(gg_opt[i])-4:unpacked_list[i+1].sequence.find(scaffld.lower())+max_annealing_len]
                        unpacked_list[i+1].sequence = "taggtctcc" + unpacked_list[i+1].sequence[unpacked_list[i+1].sequence.find(gg_opt[i])+len(gg_opt[i])-4:]
                else: #part is tRNA
                    unpacked_list[i].primer_reverse = reverse_complement(unpacked_list[i].sequence[-max_annealing_len:] + unpacked_list[i+1].sequence[:unpacked_list[i+1].sequence.find(gg_opt[i])+len(gg_opt[i])] + "tgagacccg")
                    unpacked_list[i].sequence = unpacked_list[i].sequence + unpacked_list[i+1].sequence[:unpacked_list[i+1].sequence.find(gg_opt[i])+len(gg_opt[i])] + "tgagacccg"
                    if unpacked_list[i+1].type == 'smRNA':
                        unpacked_list[i+1].primer_forward = "taggtctcc" + unpacked_list[i+1].sequence[unpacked_list[i+1].sequence.find(gg_opt[i])+len(gg_opt[i])-4:unpacked_list[i+1].sequence.find(tRNA.lower())+max_annealing_len]
                        unpacked_list[i+1].sequence = "taggtctcc" + unpacked_list[i+1].sequence[unpacked_list[i+1].sequence.find(gg_opt[i])+len(gg_opt[i])-4:]
                    else:
                        unpacked_list[i+1].primer_forward = "taggtctcc" + unpacked_list[i+1].sequence[unpacked_list[i+1].sequence.find(gg_opt[i])+len(gg_opt[i])-4:unpacked_list[i+1].sequence.find(scaffld.lower())+max_annealing_len]
                        unpacked_list[i+1].sequence = "taggtctcc" + unpacked_list[i+1].sequence[unpacked_list[i+1].sequence.find(gg_opt[i])+len(gg_opt[i])-4:]
            
            # If current part is last part, define reverse primer with right backbone overlap and forward primer ordinarily
            elif i == len(unpacked_list)-1:
                unpacked_list[i].primer_reverse = reverse_complement(unpacked_list[i].sequence[-max_annealing_len:] + bb_overlaps[-1] + 'tgagacccg')
                unpacked_list[i].sequence = unpacked_list[i].sequence + bb_overlaps[-1] + 'tgagacccg'
            
            # If current part is not first or last part, do ordinary computation
            else:
                if unpacked_list[i].type == 'pegRNA':
                    unpacked_list[i].primer_reverse = reverse_complement(unpacked_list[i].sequence[unpacked_list[i].sequence.find(scaffld.lower())+76-max_annealing_len:unpacked_list[i].sequence.find(scaffld.lower())+76+len(gg_opt[i])] + "tgagacccg")
                    unpacked_list[i+1].primer_forward = "taggtctcc" + unpacked_list[i].sequence[unpacked_list[i].sequence.find(scaffld.lower())+76+len(gg_opt[i])-4:] + unpacked_list[i+1].sequence[:max_annealing_len]
                    unpacked_list[i+1].sequence = "taggtctcc" + unpacked_list[i].sequence[unpacked_list[i].sequence.find(scaffld.lower())+76+len(gg_opt[i])-4:] + unpacked_list[i+1].sequence
                    unpacked_list[i].sequence = unpacked_list[i].sequence[:unpacked_list[i].sequence.find(scaffld.lower())+76+len(gg_opt[i])] + "tgagacccg"
                elif unpacked_list[i].type == 'gRNA':
                    unpacked_list[i].primer_reverse = reverse_complement(unpacked_list[i].sequence[-max_annealing_len:] + unpacked_list[i+1].sequence[:unpacked_list[i+1].sequence.find(gg_opt[i])+len(gg_opt[i])] + "tgagacccg")
                    unpacked_list[i].sequence = unpacked_list[i].sequence + unpacked_list[i+1].sequence[:unpacked_list[i+1].sequence.find(gg_opt[i])+len(gg_opt[i])] + "tgagacccg"
                    if unpacked_list[i+1].type == 'smRNA':
                        unpacked_list[i+1].primer_forward = "taggtctcc" + unpacked_list[i+1].sequence[unpacked_list[i+1].sequence.find(gg_opt[i])+len(gg_opt[i])-4:unpacked_list[i+1].sequence.find(tRNA.lower())+max_annealing_len]
                        unpacked_list[i+1].sequence = "taggtctcc" + unpacked_list[i+1].sequence[unpacked_list[i+1].sequence.find(gg_opt[i])+len(gg_opt[i])-4:]
                    else:
                        unpacked_list[i+1].primer_forward = "taggtctcc" + unpacked_list[i+1].sequence[unpacked_list[i+1].sequence.find(gg_opt[i])+len(gg_opt[i])-4:unpacked_list[i+1].sequence.find(scaffld.lower())+max_annealing_len]
                        unpacked_list[i+1].sequence = "taggtctcc" + unpacked_list[i+1].sequence[unpacked_list[i+1].sequence.find(gg_opt[i])+len(gg_opt[i])-4:]
                elif unpacked_list[i].type == 'smRNA':
                    unpacked_list[i].primer_reverse = reverse_complement(unpacked_list[i].sequence[-max_annealing_len:] + unpacked_list[i+1].sequence[:unpacked_list[i+1].sequence.find(gg_opt[i])+len(gg_opt[i])] + "tgagacccg")
                    unpacked_list[i].sequence = unpacked_list[i].sequence + unpacked_list[i+1].sequence[:unpacked_list[i+1].sequence.find(gg_opt[i])+len(gg_opt[i])] + "tgagacccg"
                    if unpacked_list[i+1].type == 'smRNA':
                        unpacked_list[i+1].primer_forward = "taggtctcc" + unpacked_list[i+1].sequence[unpacked_list[i+1].sequence.find(gg_opt[i])+len(gg_opt[i])-4:unpacked_list[i+1].sequence.find(tRNA.lower())+max_annealing_len]
                        unpacked_list[i+1].sequence = "taggtctcc" + unpacked_list[i+1].sequence[unpacked_list[i+1].sequence.find(gg_opt[i])+len(gg_opt[i])-4:]
                    else:
                        unpacked_list[i+1].primer_forward = "taggtctcc" + unpacked_list[i+1].sequence[unpacked_list[i+1].sequence.find(gg_opt[i])+len(gg_opt[i])-4:unpacked_list[i+1].sequence.find(scaffld.lower())+max_annealing_len]
                        unpacked_list[i+1].sequence = "taggtctcc" + unpacked_list[i+1].sequence[unpacked_list[i+1].sequence.find(gg_opt[i])+len(gg_opt[i])-4:]
                else: #part is tRNA
                    unpacked_list[i].primer_reverse = reverse_complement(unpacked_list[i].sequence[-max_annealing_len:] + unpacked_list[i+1].sequence[:unpacked_list[i+1].sequence.find(gg_opt[i])+len(gg_opt[i])] + "tgagacccg")
                    unpacked_list[i].sequence = unpacked_list[i].sequence + unpacked_list[i+1].sequence[:unpacked_list[i+1].sequence.find(gg_opt[i])+len(gg_opt[i])] + "tgagacccg"
                    if unpacked_list[i+1].type == 'smRNA':
                        unpacked_list[i+1].primer_forward = "taggtctcc" + unpacked_list[i+1].sequence[unpacked_list[i+1].sequence.find(gg_opt[i])+len(gg_opt[i])-4:unpacked_list[i+1].sequence.find(tRNA.lower())+max_annealing_len]
                        unpacked_list[i+1].sequence = "taggtctcc" + unpacked_list[i+1].sequence[unpacked_list[i+1].sequence.find(gg_opt[i])+len(gg_opt[i])-4:]
                    else:
                        unpacked_list[i+1].primer_forward = "taggtctcc" + unpacked_list[i+1].sequence[unpacked_list[i+1].sequence.find(gg_opt[i])+len(gg_opt[i])-4:unpacked_list[i+1].sequence.find(scaffld.lower())+max_annealing_len]
                        unpacked_list[i+1].sequence = "taggtctcc" + unpacked_list[i+1].sequence[unpacked_list[i+1].sequence.find(gg_opt[i])+len(gg_opt[i])-4:]
        new_builds_list.append(unpacked_list)
    
    #Optimize primer Tm
    for unpacked_list in new_builds_list:
        for part in unpacked_list:
            
            # If both Tms are already below the range, find the nearest Gs
            if mt.Tm_NN(part.primer_forward, nn_table=mt.DNA_NN3, dnac1=125, dnac2=125, Na=50) < primer_tm_range[0] and mt.Tm_NN(part.primer_reverse, nn_table=mt.DNA_NN3, dnac1=125, dnac2=125, Na=50) < primer_tm_range[0]:
                part.primer_forward_tm = mt.Tm_NN(part.primer_forward[:len(part.primer_forward)-re.search('[g,c]', reverse(part.primer_forward)).start()], nn_table=mt.DNA_NN3, dnac1=125, dnac2=125, Na=50)
                part.primer_reverse_tm = mt.Tm_NN(part.primer_reverse[:len(part.primer_reverse)-re.search('[g,c]', reverse(part.primer_reverse)).start()], nn_table=mt.DNA_NN3, dnac1=125, dnac2=125, Na=50)
                
            # If one Tm is below the range, lower the other to make them most similar and find nearest Gs
            elif mt.Tm_NN(part.primer_forward, nn_table=mt.DNA_NN3, dnac1=125, dnac2=125, Na=50) < primer_tm_range[0] and mt.Tm_NN(part.primer_reverse, nn_table=mt.DNA_NN3, dnac1=125, dnac2=125, Na=50) >= primer_tm_range[0]:
                breakit=False
                part.primer_forward_tm = mt.Tm_NN(part.primer_forward[:len(part.primer_forward)-re.search('[g,c]', reverse(part.primer_forward)).start()], nn_table=mt.DNA_NN3, dnac1=125, dnac2=125, Na=50)
                for j in range(len(part.primer_reverse),len(part.primer_reverse)-(max_annealing_len-18),-1):
                    if abs(part.primer_forward_tm-mt.Tm_NN(part.primer_reverse[:j], nn_table=mt.DNA_NN3, dnac1=125, dnac2=125, Na=50)) <=5 and part.primer_reverse[j-1] in ['g','c']:
                        part.primer_reverse_tm = mt.Tm_NN(part.primer_reverse[:j], nn_table=mt.DNA_NN3, dnac1=125, dnac2=125, Na=50)
                        part.primer_reverse = part.primer_reverse[:j]
                        breakit=True
                if not breakit:
                    part.primer_reverse_tm = mt.Tm_NN(part.primer_reverse, nn_table=mt.DNA_NN3, dnac1=125, dnac2=125, Na=50)
            elif mt.Tm_NN(part.primer_reverse, nn_table=mt.DNA_NN3, dnac1=125, dnac2=125, Na=50) < primer_tm_range[0] and mt.Tm_NN(part.primer_forward, nn_table=mt.DNA_NN3, dnac1=125, dnac2=125, Na=50) >= primer_tm_range[0]:
                breakit=False
                part.primer_reverse_tm = mt.Tm_NN(part.primer_reverse[:len(part.primer_reverse)-re.search('[g,c]', reverse(part.primer_reverse)).start()], nn_table=mt.DNA_NN3, dnac1=125, dnac2=125, Na=50)
                for i in range(len(part.primer_forward),len(part.primer_forward)-(max_annealing_len-18),-1):
                    if abs(part.primer_reverse_tm-mt.Tm_NN(part.primer_forward[:i], nn_table=mt.DNA_NN3, dnac1=125, dnac2=125, Na=50)) <=5 and part.primer_forward[i-1] in ['g','c']:
                        part.primer_forward_tm = mt.Tm_NN(part.primer_forward[:i], nn_table=mt.DNA_NN3, dnac1=125, dnac2=125, Na=50)
                        part.primer_forward = part.primer_forward[:i]
                        breakit=True
                if not breakit:
                    part.primer_forward_tm = mt.Tm_NN(part.primer_forward, nn_table=mt.DNA_NN3, dnac1=125, dnac2=125, Na=50)
            
            # If both Tms are within or above the range, lower them into the range and make them most similar and find nearest Gs
            else:
                breakit=False
                # Do all of the above
                for i in range(len(part.primer_forward),len(part.primer_forward)-(max_annealing_len-18),-1):
                    for j in range(len(part.primer_reverse),len(part.primer_reverse)-(max_annealing_len-18),-1):
                        if abs(mt.Tm_NN(part.primer_forward[:i], nn_table=mt.DNA_NN3, dnac1=125, dnac2=125, Na=50)-mt.Tm_NN(part.primer_reverse[:j], nn_table=mt.DNA_NN3, dnac1=125, dnac2=125, Na=50))<=5 and mt.Tm_NN(part.primer_forward[:i], nn_table=mt.DNA_NN3, dnac1=125, dnac2=125, Na=50) >= primer_tm_range[0] and mt.Tm_NN(part.primer_forward[:i], nn_table=mt.DNA_NN3, dnac1=125, dnac2=125, Na=50) <= primer_tm_range[1] and mt.Tm_NN(part.primer_reverse[:j], nn_table=mt.DNA_NN3, dnac1=125, dnac2=125, Na=50) >= primer_tm_range[0] and mt.Tm_NN(part.primer_reverse[:j], nn_table=mt.DNA_NN3, dnac1=125, dnac2=125, Na=50) <= primer_tm_range[1] and part.primer_reverse[j-1] in ['g','c'] and part.primer_forward[i-1] in ['g','c']:
                            part.primer_forward_tm = mt.Tm_NN(part.primer_forward[:i], nn_table=mt.DNA_NN3, dnac1=125, dnac2=125, Na=50)
                            part.primer_forward = part.primer_forward[:i]
                            part.primer_reverse_tm = mt.Tm_NN(part.primer_reverse[:j], nn_table=mt.DNA_NN3, dnac1=125, dnac2=125, Na=50)
                            part.primer_reverse = part.primer_reverse[:j]
                            breakit=True
                        if breakit:
                            break
                    if breakit:
                        break
                # Do range and Gs        
                if not breakit:
                    for i in range(len(part.primer_forward),len(part.primer_forward)-(max_annealing_len-18),-1):
                        for j in range(len(part.primer_reverse),len(part.primer_reverse)-(max_annealing_len-18),-1):
                            if mt.Tm_NN(part.primer_forward[:i], nn_table=mt.DNA_NN3, dnac1=125, dnac2=125, Na=50) >= primer_tm_range[0] and mt.Tm_NN(part.primer_forward[:i], nn_table=mt.DNA_NN3, dnac1=125, dnac2=125, Na=50) <= primer_tm_range[1] and mt.Tm_NN(part.primer_reverse[:j], nn_table=mt.DNA_NN3, dnac1=125, dnac2=125, Na=50) >= primer_tm_range[0] and mt.Tm_NN(part.primer_reverse[:j], nn_table=mt.DNA_NN3, dnac1=125, dnac2=125, Na=50) <= primer_tm_range[1] and part.primer_reverse[j-1] in ['g','c'] and part.primer_forward[i-1] in ['g','c']:
                                part.primer_forward_tm = mt.Tm_NN(part.primer_forward[:i], nn_table=mt.DNA_NN3, dnac1=125, dnac2=125, Na=50)
                                part.primer_forward = part.primer_forward[:i]
                                part.primer_reverse_tm = mt.Tm_NN(part.primer_reverse[:j], nn_table=mt.DNA_NN3, dnac1=125, dnac2=125, Na=50)
                                part.primer_reverse = part.primer_reverse[:j]
                                breakit=True
                            if breakit:
                                break
                        if breakit:
                            break
                # Do only Gs nearest to range
                if not breakit:
                    for i in range(len(part.primer_forward)-1-(max_annealing_len-18), len(part.primer_forward)):
                        for j in range(len(part.primer_reverse)-1-(max_annealing_len-18), len(part.primer_reverse)):
                            if part.primer_reverse[j-1] in ['g','c'] and part.primer_forward[i-1] in ['g','c']:
                                part.primer_forward_tm = mt.Tm_NN(part.primer_forward[:i], nn_table=mt.DNA_NN3, dnac1=125, dnac2=125, Na=50)
                                part.primer_forward = part.primer_forward[:i]
                                part.primer_reverse_tm = mt.Tm_NN(part.primer_reverse[:j], nn_table=mt.DNA_NN3, dnac1=125, dnac2=125, Na=50)
                                part.primer_reverse = part.primer_reverse[:j]
                                breakit=True
                            if breakit:
                                break
                        if breakit:
                            break
                    
                if not breakit:
                    part.primer_forward_tm = mt.Tm_NN(part.primer_forward[:len(part.primer_forward)-reverse(part.primer_forward).find('g')], nn_table=mt.DNA_NN3, dnac1=125, dnac2=125, Na=50)
                    part.primer_reverse_tm = mt.Tm_NN(part.primer_reverse[:len(part.primer_reverse)-reverse(part.primer_reverse).find('g')], nn_table=mt.DNA_NN3, dnac1=125, dnac2=125, Na=50)
    return new_builds_list[0],ftrs,msg


def PTGbldr(inserts):
    '''
    function which takes all desired parts of PTG GG assembly and gives out the respective inserts for PTG. During the 
    process, each part is appended with the same tRNA. The unit of part and tRNA are then treated as one insert.
    
    parts must be of the form [['name', 'type', 'sequence'],[...]]
    '''
    
    PTG_parts = []
    PTG_parts.append(Part('tRNA', 'tRNA', tRNA))
    
    ## Take each coding sequence of the output and append it with a tRNA to the list
    for c,prt in enumerate(inserts):
        if prt[1] == 'pegRNA':
            PTG_parts.append(Part(prt[0], prt[1], str(prt[2])))
            PTG_parts.append(Part('tRNA', 'tRNA', tRNA))
        elif prt[1] == 'gRNA':
            PTG_parts.append(Part(prt[0], prt[1], str(prt[2]) + scaffld + tRNA))
        elif prt[1] == 'smRNA':
            PTG_parts.append(Part(prt[0], prt[1], str(prt[2]) + tRNA))
    
    return PTG_parts


# Execute computation
def runall(arr, tm_range=[52,72], max_ann_len=30, bb_overlaps=['tgcc','gttt'], additional_overhangs=[]):
    msg = None

    all_existing = bb_overlaps+additional_overhangs
    all_existing_forwrev = all_existing+[reverse_complement(i) for i in all_existing]
    if len(all_existing_forwrev) != len(np.unique(all_existing_forwrev)):
        msg = 'comb_error'
        return None,None,msg
    else:

        full_sequence = ''
        PTG = PTGbldr(arr)
        outpt,feat,msg = scarless_gg(PTG, tm_range, max_ann_len, bb_overlaps, additional_overhangs)
    
        if outpt is not None:
            oligos = []
            for c,o in enumerate(outpt):
                oligos.append(o.primer_forward)
                oligos.append(o.primer_reverse)
                if c == 0:
                    full_sequence += o.sequence[13:-13]
                else:
                    full_sequence += o.sequence[9:-13]
        print(full_sequence)

        return outpt,full_sequence,feat,msg

