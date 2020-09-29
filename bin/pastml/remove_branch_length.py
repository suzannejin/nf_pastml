#!/usr/bin/env python

__description__='''

Remove branch length from tree, 
and keep only the topology and leaf names (by default).

'''

if __name__ == '__main__':

    import sys
    import argparse
    from ete3 import Tree
    
    app=argparse.ArgumentParser(description=__description__)
    app.add_argument("-tree", required=True, type=str,  help="Tree in newick format.")
    app.add_argument("-out", required=True, type=str, help="Pruned tree.")
    app.add_argument("-format", default=9, type=int, help="Output tree format. Check http://etetoolkit.org/docs/latest/tutorial/tutorial_trees.html?highlight=prune#sec-newick-formats for further information")
    args=app.parse_args()

    # read tree
    tree = Tree(args.tree)

    # write tree
    tree.write(format=args.format, outfile=args.out)   
