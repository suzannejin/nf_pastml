#!/usr/bin/env python

__description__='''

Count state transitions from PastML's output.

'''


def tab2step(filename):
    ''' Record steps (state transition) from tab file '''

    # Get steps
    with open(filename) as f:
        for line in f:
            line = line.strip('\n')
            fields = line.split('\t')
            if fields[0] == 'steps':
                return(int(fields[1]))
    
    # Raise error status if no 'steps' (state transition) line is in file
    raise ValueError('No steps found in file: ' + filename)


def msa2shape(filename):
    ''' Get the shape of the msa '''

    names, seqs = msa2list(filename)
    shape = (len(names), len(seqs[0]))

    return(shape)
    

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

    import argparse
    import glob
    import os
    import sys
    from pathlib import Path

    app=argparse.ArgumentParser(description=__description__)
    app.add_argument('-dir', required=True, type=str, help="PastML's output directory")
    app.add_argument('-msa', type=str, help='MSA')
    app.add_argument('-bycol', action='store_true', help='Print count normalized by the number of columns, too')
    app.add_argument('-bycolrow', action='store_true', help='Print countnormalized by the number of columns and rows, too')
    args=app.parse_args()

    # Count total steps (state transitions) across columns
    step = 0
    # Get tab files
    pattern = os.path.join(args.dir, 'params.*.tab')
    for i,tab in enumerate(glob.glob(pattern)):
        # Read tab file
        step += tab2step(tab)
        o = open(os.path.join(args.dir, 'steps'), 'w')
        o.write(str(step) +'\n')
        o.close()
    # print("Step: " + str(step))

    # normalize
    if args.bycol:
        step2 = step /(i+1)
        o = open(os.path.join(args.dir, 'steps2'), 'w')
        o.write(str(step2) +'\n')
        o.close()
        # print("Step2: " + str(step2))
    if args.bycolrow:
        shape = msa2shape(args.msa)
        step3 = step / shape[0] / shape[1]
        o = open(os.path.join(args.dir, 'steps3'), 'w')
        o.write(str(step3) +'\n')
        o.close()
        # print("Step3: " + str(step3))
