#!/usr/bin/env python
import argparse
import copy
import glob
import os
import re
import sys

import math
import numpy as np
import pandas as pd
import statistics as st

from Bio import Phylo



__description__='''

Determine the gap fragmentation score of an alignment, based on the information in the guide tree.

'''


class MSA:
    def __init__(self, filename):
        self.names, self.seqs = self.msa2list(filename)
        self.states = self.seq2state(self.seqs)
        self.mat = np.array(self.states)
        self.df = pd.DataFrame(self.mat, index=self.names)
        cols = []
        for i in range(self.mat.shape[1]):
            col = self.mat[:,i]
            col = ''.join(col)
            cols.append(col)
        self.sorted_columns = sorted(cols)
    
    def msa2list(self, filename):
        ''' Convert the alignment in two lists: names & seqs. '''
        names = []
        seqs = []
        with open(filename) as f:
            for line in f:
                line = line.strip('\n')
                if line[0] == '>':
                    names.append(line[1:])
                    seqs.append('')
                else:
                    seqs[-1] += line
        return(names, seqs)

    def seq2state(self, seqs):
        ''' Determine the 0/1 residue/gap states given a list of sequences. '''
        states = []
        for seq in seqs:
            state = ''
            for i in seq:
                if i == '-':
                    state += '1'
                else:
                    state += '0'
            state = list(state)
            states.append(state)
        return(states)


class Tree:
    
    def __init__(self, filename):
        tree = Phylo.read(filename, 'newick')
        self.tree, self.nodes, self.leaves, self.maxdepth, self.levels = self.preprocess_tree(tree)
        
    def preprocess_tree(self, tree):
        nodes, leaves = tree.get_nonterminals(), tree.get_terminals()
        nodes = self.set_node_names(nodes)
        tree, nodes, leaves, maxdepth, levels = self.tree2level(tree, nodes, leaves)
        return(tree, nodes, leaves, maxdepth, levels)
        
    def set_node_names(self, nodes):
        ''' Rename the nodes if required '''
        for i in range(len(nodes)):
            if not nodes[i].name:
                if i == 0:
                    nodes[i].name = 'root'
                else:
                    nodes[i].name = 'n' + str(i-1)
        return(nodes)
    
    def tree2level(self, tree, nodes, leaves):
        ''' Classify nodes/leaves according to the depth '''
        depths = tree.depths(leaves[0])
        maxdepth = max(depths.values())
        levels = {}
        for element,level in depths.items():
            if level not in levels:
                levels[level] = []
            levels[level].append(element)
            element.comment = {'depth':level}
        return(tree, nodes, leaves, maxdepth, levels)

# ---------------------
# COMPUTE FRAGMENTATION
# ---------------------

def check_gap(nodes):
    ''' Check if any of the given nodes has gapped descendents. '''
    for node in nodes:
        if node.comment['score'] > 0:
            return(True)
    return(False)

def sum_gaps(nodes):
    score = 0.0
    for node in nodes:
        score += node.comment['score']
    return(score)

def rootsum_gaps(nodes):
    score = 0.0
    for node in nodes:
        score += node.comment['score']**2
    score = math.sqrt(score)
    return(score)

def bjorn_rootsum_gaps(nodes, s1):
    s2 = rootsum_gaps(nodes)
    score = 2*s1 - s2
    return(score)

def col2fragmentation(tree, coldict, bjorn):
    ''' Determine the fragmentation score of the tree, for a given column coldict{seq name : 0/1 gap state} '''
    # update leaves' gap state 
    for leaf in tree.leaves:
        leaf.comment['score'] = int(coldict[leaf.name])
    # update nodes' gap state
    for i in range(tree.maxdepth)[::-1]:
        elements = tree.levels[i]
        for element in elements:
            if element.count_terminals() == 1:
                continue
            elif element.is_preterminal():
                childs = element.get_terminals()
                score = sum_gaps(childs)
            else:
                childs = element.get_nonterminals()
                childs2 = element.get_terminals()
                childs = [node for node in childs if node.comment['depth']==(i+1)] 
                sum_ = sum_gaps(childs)
                if (check_gap(childs)) and (sum_ != len(childs2)):
                    if bjorn:
                        score = bjorn_rootsum_gaps(childs, sum_)
                    else:
                        score = rootsum_gaps(childs)
                else:
                    score = sum_
            element.comment['score'] = score
        if i == 0:
            root = element
    fragmentation = root.comment['score']
    return(fragmentation)

def data2fragmentation(tree, msa, post, bjorn):
    fragmentation = 0.0
    pre_col = ""
    for i,col in enumerate(msa.sorted_columns):
        if pre_col == col:
            score = pre_score
        else:
            coldict = dict(zip(msa.names, list(col)))
            score = col2fragmentation(tree, coldict, bjorn)
        fragmentation += score
        pre_col = col
        pre_score = score
    if post == 'average':
        fragmentation = fragmentation / (i+1)
    return(fragmentation)


if __name__ == '__main__':
    
    # parse arguments
    app=argparse.ArgumentParser(description=__description__)
    app.add_argument('-tree', required=True, type=str, help="Guide tree in newick format")
    app.add_argument('-msa', required=True, type=str, help="MSA in FASTA format")
    app.add_argument('-recursion', type=int, default=10**6, help="Recursion limit")
    app.add_argument('-post', type=str, choices=['none', 'average'], default='average', help='Post processing of the fragmentation score summed over columns.')
    app.add_argument('-bjorn', action='store_true', help='Do 2*sum - root of square sum.')
    args=app.parse_args()

    # set recursion limit
    sys.setrecursionlimit(args.recursion)

    # read msa
    msa = MSA(args.msa)
    # print(msa)

    # read tree
    tree = Tree(args.tree)
    # print(tree)

    # determine fragmentation score
    fragmentation = data2fragmentation(tree, msa, args.post, args.bjorn)
    print(fragmentation)



