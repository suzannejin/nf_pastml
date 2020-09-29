#!/usr/bin/env python

__description__='''

Prune tree (newick), keeping the sequences in MSA.

'''

def msa2list(filename):
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


if __name__ == '__main__':

    import sys
    import argparse
    from ete3 import Tree
    
    app=argparse.ArgumentParser(description=__description__)
    app.add_argument("-msa", required=True, type=str,help="MSA in fasta format.")
    app.add_argument("-tree", required=True, type=str,  help="Tree in newick format.")
    app.add_argument("-out", required=True, type=str, help="Pruned tree.")
    app.add_argument("-format", default=9, type=int, help="Output tree format. Check http://etetoolkit.org/docs/latest/tutorial/tutorial_trees.html?highlight=prune#sec-newick-formats for further information")
    args=app.parse_args()

    # read int: msa and tree
    names, seqs = msa2list(args.msa)
    tree = Tree(args.tree)

    # prune tree and write output
    tree.prune(names)
    tree.write(format=args.format, outfile=args.out)   