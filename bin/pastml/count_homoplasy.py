#!/usr/bin/env python

import argparse
import glob
import os
import re
import sys
from Bio import Phylo


__description__='''

Count the number of homoplasic 0|1 state transitions. 

'''


def count_homoplasic_transitions(tree, action='homo', ncols=None):
    
    # root states
    root_p = tree.clade
    root_comment = root_p.comment
    root_states = re.findall('col[0-9]+=([0-9]+)', root_comment)
    if not ncols:
        ncols = [i for i in range(len(root_states))]
    steps, homo = 0, 0

    # for column
    for i in ncols:
        sti = str(i)

        # checked parent-child pairs
        checked = []
        l_res2gap, l_gap2res, l_homoplasic = [], [], []

        # for leaf in tree
        for leaf in tree.get_terminals():

            # get path to leaf
            path = tree.get_path(leaf)

            # continue when path = 1, 
            # since there is no place for one transition to gap, another back, and another to an homoplasic gap
            # this requires a minimal of 
            ### 2 steps if root = gap,
            ### 3 steps, if root = residue
            #if len(path) <= 1:
                #continue

            # pre state
            pre_p = root_p
            root_state = root_states[i]
            pre_state = root_state

            # set states
            if root_state == '0':
                res2gap = 1
            elif root_state == '1':
                res2gap = 0
            else:
                raise ValueError('Check root states')
            gap2res = 0

            # for next path
            for p in path[:]:

                # current states
                comment = p.comment
                match = re.search(rf'col{re.escape(sti)}=[0-9]+', comment)
                state = match.group().split('=')[1]

                # check combination
                comb = (pre_p.name, p.name)
                if comb in checked:
                    if (res2gap == 0) and (comb in l_res2gap):
                        res2gap = 1
                    if (gap2res == 0) and (comb in l_gap2res):
                        gap2res = 1
                    # if checked, update state and continue
                    pre_p = p
                    pre_state = state
                    continue

                # if homoplasic
                if (res2gap > 0) and (gap2res > 0) and (pre_state != state):
                    l_homoplasic.append(comb)
                # if res -> gap
                elif (res2gap == 0) and (pre_state == '1') and (state == '0'):
                    res2gap = 1
                    l_res2gap.append(comb)
                # if gap -> res
                elif (gap2res == 0) and (pre_state == '0') and (state == '1'):
                    gap2res = 1
                    l_gap2res.append(comb)

                # update state
                pre_p = p
                pre_state = state
                checked.append(comb)
        # update numbers
        #print(l_res2gap, l_gap2res)
        steps += len(l_res2gap) + len(l_gap2res) + len(l_homoplasic)
        homo += len(l_homoplasic)
        
    if action == 'homo':
        return(homo)
    elif action == 'steps':
        return(steps)
    elif action == 'both':
        return(homo, steps)
    else:
        raise ValueError('Incorrect action')



if __name__ == '__main__':

    app=argparse.ArgumentParser(description=__description__)
    app.add_argument('-tree', required=True, type=str, help="PastML's output tree")
    app.add_argument('-action', type=str, choices=['homo','steps','both'], default='homo')
    app.add_argument('-ncols', type=int, nargs='+')
    app.add_argument('-recursion',type=int,default=int(1e8), help="Recursion limit")
    args=app.parse_args()

    # set recursion limit
    sys.setrecursionlimit(args.recursion)

    # read tree
    tree = Phylo.read(args.tree, 'newick')

    # count
    if args.action == 'both':
        homo, steps = count_homoplasic_transitions(tree, action=args.action, ncols=args.ncols)
        print(str(homo) + "\t" + str(steps))
    else:
        counts = count_homoplasic_transitions(tree, action=args.action, ncols=args.ncols)
        print(counts)

