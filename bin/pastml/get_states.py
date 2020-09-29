#!/usr/bin/env python

__description__='''

Create .csv states file for PastML.

'''

def readAlign(infile,typ):
    ''' Read alignment.
    Input - infile : input alignment filename
          - typ    : alignment format (clustal, nexus, phylip, emboss, fasta, ...)
    '''
    from Bio import AlignIO
    align=AlignIO.read(infile,typ)

    return(align)


def align2dic(msa):
    from Bio import AlignIO
    msadic={}
    for m in msa:
        name=m.name
        seq=m.seq
        msadic[name]=seq
    return(msadic)


def readTree(infile,typ):
    ''' Read phylogenetic tree.
    Input - infile : input tree filename
          - typ    : tree format (newick, phyloxml, nexus, ...)
    '''
    from Bio import Phylo 
    tree=Phylo.read(infile,typ)

    return(tree)


def getTreenames(infile,typ):
    tree=readTree(infile,typ)
    term_names = [term.name for term in tree.get_terminals()]
    return(term_names)


def msa2states_tree(msa, tree, out):
    ''' Check the sequences in tree, and then write the states to output. '''

    # Read MSA
    msa=readAlign(msa,"fasta")
    msadic=align2dic(msa)
    
    # Get tree names
    tree_names=getTreenames(tree,"newick")

    # Write output
    out=open(out,"w")
    # Header
    out.write("ID")
    for name in msadic:
        seq=msadic[name]
        break
    for col in range(len(seq)):
        out.write(",col"+str(col))
    out.write("\n")
    # Gap/residue as binary 0/1
    for name in tree_names:
        seq=msadic[name]
        out.write(name)
        for col in seq:
            if col=="-":
                state=0
            else:
                state=1
            out.write(","+str(state))
        out.write("\n")
    out.close()


def msa2states(msa, out):

    # Read msa
    print('----Reading msa...')
    names, seqs = msa2list(msa)
    
    # Get states
    print('----Converting states...')
    states = seq2state(seqs)

    # Write output
    print('----Writing output csv...')
    states2out(names, states, out)

    print('----Finished')


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


def seq2state(seqs):
    ''' Determine the 0/1 gap/residue states given a list of sequences. '''

    states = []
    for seq in seqs:
        state = ''
        for i in seq:
            if i == '-':
                state += '0'
            else:
                state += '1'
        states.append(state)
    return(states)


def states2out(names, states, out):
    ''' Write the 0/1 states to the output csv file. '''

    out = open(out, 'w')
    
    # header
    out.write('ID')
    for col in range(len(states[0])):
        out.write(',col'+str(col))
    out.write('\n')

    # states
    for i,name in enumerate(names):
        out.write(name)
        state = states[i]
        for j in state:
            out.write(','+j)
        out.write('\n')

    out.close()


if __name__ == '__main__':

    import sys
    import argparse

    app=argparse.ArgumentParser(description=__description__)
    app.add_argument("-msa", required=True, type=str,help="MSA in fasta format.")
    app.add_argument("-out", required=True, type=str, help="Output csv file.")
    app.add_argument("-tree", type=str,  help="Tree in .nwk format. This is required if you want to check the sequences and sort them according to the tree before writing the states.")
    args=app.parse_args()

    if args.tree:
        # This is prefered to facilitate the analysis
        # In this way, the sequences will be sorted according to the tree, consequently they will appear in the same order as the PastML output.
        msa2states_tree(args.msa, args.tree, args.out)
    else:
        msa2states(args.msa, args.out)