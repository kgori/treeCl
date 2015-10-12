#!/usr/bin/env python
from treeCl.utils import fileIO
from treeCl.wrappers.phylogenetics import Raxml
import os, re, warnings
import numpy as np

def run_raxml_fast(fl, model='GTRGAMMA'):
    fl = os.path.abspath(fl)
    with fileIO.TempDir() as tmpd, fileIO.ChDir(tmpd):
        rax = Raxml()
        rax(s=fl, m=model, p='123', f='F', n='initial_tree', wait=True)
        # print rax.get_stdout()
        rax(s=fl, m=model, p='123', f='e', R='RAxML_binaryModelParameters.initial_tree', t='RAxML_fastTree.initial_tree', n='model_opt', wait=True)
        # print rax.get_stdout()
        resfl = 'RAxML_info.model_opt'
        treefl = 'RAxML_result.model_opt'
        with open(resfl) as resfl_handle, open(treefl) as treefl_handle:
            info = resfl_handle.read()
            tree = treefl_handle.read().rstrip()
    result = parse_raxml_output(info)
    result['tree'] = tree
    return result

def parse_raxml_output(s):
    ac=float(re.search(r'rate A <-> C:\s+([0-9.]+)', s).groups()[0])
    ag=float(re.search(r'rate A <-> G:\s+([0-9.]+)', s).groups()[0])
    at=float(re.search(r'rate A <-> T:\s+([0-9.]+)', s).groups()[0])
    cg=float(re.search(r'rate C <-> G:\s+([0-9.]+)', s).groups()[0])
    ct=float(re.search(r'rate C <-> T:\s+([0-9.]+)', s).groups()[0])
    gt=float(re.search(r'rate G <-> T:\s+([0-9.]+)', s).groups()[0])
    a=float(re.search(r'freq pi\(A\):\s+([0-9.]+)',s).groups()[0])
    c=float(re.search(r'freq pi\(C\):\s+([0-9.]+)',s).groups()[0])
    g=float(re.search(r'freq pi\(G\):\s+([0-9.]+)',s).groups()[0])
    t=float(re.search(r'freq pi\(T\):\s+([0-9.]+)',s).groups()[0])
    alpha = float(re.search(r'alpha:\s+([0-9.]+)' ,s).groups()[0])
    loglk = float(re.search(r'Final GAMMA  likelihood:\s+([0-9-.]+)', s).groups()[0])

    return {
        'loglk': loglk,
        'alpha': alpha,
        'freqs': [a, c, g, t],
        'rates': [ac, ag, at, cg, ct, gt],
    }

def parse_fasttree_output(filename):

    try:
        with open(filename) as handle:
            s = handle.read()
    except IOError:
        warnings.warn('Couldn\'t open {}'.format(filename))
        return None
    
    try:
        loglk, alpha = (float(x) for x in re.search(r'Gamma\(\d+\) LogLk = ([0-9-.]+) alpha = ([0-9.]+)', s).groups())
    except AttributeError:
        warnings.warn('Couldn\'t parse loglk and alpha from {}'.format(filename))
        return None

    try:
        a, c, g, t = (float(x) for x in re.search(r'GTR Frequencies: ([0-9.]+)\s+([0-9.]+)\s+([0-9.]+)\s+([0-9.]+)', s).groups())
    except AttributeError:
        warnings.warn('Couldn\'t parse base frequencies from {}'.format(filename))
        return None

    try:
        ac, ag, at, cg, ct, gt = (float(x) for x in re.search(r'GTR rates\(ac ag at cg ct gt\) ([0-9.]+)\s+([0-9.]+)\s+([0-9.]+)\s+([0-9.]+)\s+([0-9.]+)\s+([0-9.]+)', s).groups())
    except AttributeError:
        warnings.warn('Couldn\'t parse substitution rates from {}'.format(filename))
        return None

    return {
        'loglk': loglk,
        'alpha': alpha,
        'freqs': np.array([a, c, g, t]),
        'rates': np.array([ [0 , ac, ag, at], 
                            [ac, 0 , cg, ct], 
                            [ag, cg, 0 , gt], 
                            [at, ct, gt, 0 ] ])
    }

if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument('-f', '--filename', type=str)
    parser.add_argument('-o', '--output', type=str)
    args = parser.parse_args()

    import json
    result = run_raxml_fast(args.filename)
    with fileIO.fwriter(args.output) as fl:
        json.dump(result, fl)

