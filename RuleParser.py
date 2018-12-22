'''

inputs: node list; enzyme + glycan rule set; compartment list (# compartments); 
        compartment composition (which enzyme in which compartment) 

'''

import json
import operator
import numpy as np
import pandas as pd
from GlycanRep import *

def read_node_file(filename):
    '''specifies list of nodes (parent + children) present in the tree'''
    x = pd.read_csv(filename, header=None)
    nodes = []
    for i in range(x.shape[0]):
        nodes = nodes + x.loc[i,:].tolist()
    return list(set(nodes))


def read_enzyme_file(filename, node_list):
    '''specifies the rule set of enzyme-catalysed bonding of 2 sugars, & the attachment position (left, right, any).
    The rules are read from a structured file into a dictionary'''
    x = pd.read_csv(filename, header=None)
    assert x.shape[1] == 4  # check that all 4 columns are present
    rules = {}
    for i in range(x.shape[0]):
        enz = x.loc[i,0]
        n1 = x.loc[i,1]
        n2 = x.loc[i,2]
        r  = x.loc[i,3]
        assert n1 in node_list  ## validate the inputs by cross-checking with other files
        assert n2 in node_list
        assert r  in ('l', 'r', 'a')
        if enz not in rules:
            rules[enz] = {}
        if n1 not in rules[enz]:
            rules[enz][n1] = []
        rules[enz][n1].append( (n2,r) )
    return rules


def read_compartment_file(filename):
    '''list of compartments'''
    x = pd.read_csv(filename, header=None)
    compartments = []
    for i in range(x.shape[0]):
        compartments = compartments + x.loc[i,:].tolist()
    return compartments


def read_compartment_composition_file(filename, compartments, enzymes):
    '''compartment composition: specifies which compartment has which enzyme, thereby activating appropriate rules'''
    f = open(filename)
    lines = f.readlines()
    composition = {}
    for line in lines:
        fields = [ w.strip()  for w in line.split(',')]
        if len(fields) > 1:
            assert fields[0] in compartments
            if fields[0] not in composition:
                composition[fields[0]] = []
            for i in range(1,len(fields)):
                assert fields[i] in enzymes
                composition[fields[0]].append(fields[1])
    return composition


def read_compartment_residence(filename, compartments):
    '''save residence time in compartments as a dictionary; higher time implies higher saturation, i.e. more sites have reacted'''
    x = pd.read_csv(filename, header=None)
    assert x.shape[1] == 2
    assert x.shape[0] == len(compartments)
    residence_times = {}
    for i in range(x.shape[0]):
        c = x.loc[i,0]
        r = x.loc[i,1]
        assert c in compartments
        residence_times[c] = r
    return residence_times


def run_reaction(root, enzymes, compartments, composition, residence ):
    for c in compartments:
        for n in range(residence[c]):
            e = np.random.choice(composition[c])
            n1 = np.random.choice(list(enzymes[e].keys()))
            poss_rxn = enzymes[e][n1]
            np.random.shuffle(poss_rxn)
            for i in range(len(poss_rxn)):
                (n2, r) = poss_rxn[i]
                if r == 'a':
                    r = 'any'
                elif r == 'r':
                    r = 'right'
                else:
                    r = 'left'
                succ = attach_new_node(root,n1,n2,attach_to=r)
                if succ is True:
                    break


def freq_to_distrib( data ):
    assert isinstance(data, dict)
    n = sum(data.values())
    distrib = {}
    for k,v in data.items():
        distrib[k] = v * 1.0/n
    return distrib


def append_freq( base, newfreq):
    assert isinstance(base, dict)
    assert isinstance(newfreq, dict)
    for k,v in newfreq.items():
        if k in base:
            base[k] = base[k] + v
        else:
            base[k] = v
    return base


def composition_distance(src, tgt):
    assert isinstance(src,dict)
    assert isinstance(tgt,dict)
    src_dist = freq_to_distrib(src)
    tgt_dist = freq_to_distrib(tgt)
    all_keys = set(list(src.keys()) + list(tgt.keys()))
    d = 0
    for k in all_keys:
        v1, v2 = 0, 0
        if k in src_dist:
            v1 = src_dist[k]
        if k in tgt_dist:
            v2 = src_dist[k]
        d += abs(v1 - v2)
    return d


def topn( data , n=5):
    assert isinstance(data, dict)
    sorted_x = sorted(data.items(), key=operator.itemgetter(1), reverse=True)
    return sorted_x[:n]


if __name__ == "__main__":
    node_file = 'node_list.csv'
    enz_file = 'enzyme_defn.csv'
    cmp_file = 'compartments.csv'
    comp_file = 'composition.csv'
    resi_file = 'compartment_residence.csv'

    nodes = read_node_file(node_file)
    enz = read_enzyme_file(enz_file, nodes)
    compartments = read_compartment_file(cmp_file)
    compositions = read_compartment_composition_file(comp_file, compartments, enz.keys())
    residence = read_compartment_residence(resi_file, compartments)
    nop = 10000

    basefreq, freq = {}, {}
    chunks = 5000
    repeats = 10000
    runlength, cutofflength, eps = 0, 5, 0.01

    for i in range(repeats):
        print( 'running iteration: %d' % i)
        freq = {}
        for j in range(chunks):
            root = TreeNode('Glc')
            run_reaction(root,enz,compartments,compositions,residence)
            name = '%s' % root
            if name not in freq:
                freq[name] = 0
            freq[name] = freq[name] + 1
        basefreq = append_freq(basefreq, freq)
        if i > 1:
            d = composition_distance(basefreq, freq)
            runlength = runlength + 1 if d < eps else 0
        if runlength > cutofflength:
            break
    print( json.dumps( topn(basefreq, 10), indent=4) )




