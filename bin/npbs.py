#!/usr/bin/env python

import treeCl
from treeCl import Partition
from treeCl.utils import fileIO
from treeCl.wrappers.phylogenetics import Raxml, FastTree
import os
import glob
import re
import logging
import json
from Bio import AlignIO

def get_dirs(path, i):
    working_dir = os.path.join(path, 'npbs{}'.format(i))
    cache_dir = os.path.join(working_dir, 'cache')
    concat_dir = os.path.join(working_dir, 'concat')
    return {'wdir': working_dir, 'cachedir': cache_dir, 'concdir': concat_dir}

def generate_npbs(path, i):
    c = treeCl.Collection(input_dir=path, file_format='phylip')
    working_dir = get_dirs(path, i)['wdir']
    # Check if work already done
    work_done = True
    for rec in c:
        looking_for = '{}.phy'.format(os.path.join(working_dir, rec.name))
        if not (os.path.exists(looking_for) and os.path.getsize(looking_for) > 0):
            if not (os.path.exists(looking_for + '.bz2') and os.path.getsize(looking_for + '.bz2') > 0):
                logger.error("File not found or is empty: {}".format(looking_for))
                work_done = False

    if not work_done:
        npbs = c.permuted_copy()
        if not os.path.exists(working_dir):
            os.mkdir(working_dir)
        for rec in npbs:
            rec.write_alignment('{}.phy'.format(os.path.join(working_dir, rec.name)), 'phylip', True)
            AlignIO.convert('{}.phy'.format(os.path.join(working_dir, rec.name)), 'phylip-relaxed', '{}.phy_'.format(os.path.join(working_dir, rec.name)), 'phylip-relaxed')
            os.system('mv {} {}'.format('{}.phy_'.format(os.path.join(working_dir, rec.name)), '{}.phy'.format(os.path.join(working_dir, rec.name))))

def get_collection(path, i):
    working_dir = get_dirs(path, i)['wdir']
    return treeCl.Collection(input_dir=working_dir, file_format='phylip')

def read_cache(collection, path, i):
    cache_dir = get_dirs(path, i)['cachedir']
    try:
        collection.read_parameters(cache_dir)
        for rec in collection:
            _ = rec.tree
        return True
    except:
        return False

def calc_trees(collection, path, i, logger, method, fastest):
    working_dir = get_dirs(path, i)['wdir']
    if method == 'fasttree':
        collection.fast_calc_trees()
    elif method == 'raxml':
        collection.calc_trees(tree_search=False)
    else:
        raise ValueError('Unrecognised method {}'.format(method))

    collection.write_parameters(os.path.join(working_dir, 'cache'))

def calc_dists(collection, path, i):
    cache_dir = get_dirs(path, i)['cachedir']
    try:
        geo = treeCl.DistanceMatrix.from_csv(os.path.join(cache_dir, 'geo.csv.gz'))
    except:
        geo = collection.get_inter_tree_distances('geo')
        geo.to_csv(os.path.join(cache_dir, 'geo.csv.gz'))
    return geo

def cluster(dm, path, index, nclust):
    cache_dir = get_dirs(path, index)['cachedir']
    if os.path.exists(os.path.join(cache_dir, 'partitions.txt')):
        with treeCl.utils.fileIO.freader(os.path.join(cache_dir, 'partitions.txt')) as fl:
            ps = [eval(line) for line in fl]
    else:
        cl = treeCl.Spectral(dm)
        ps = [cl.cluster(i) for i in range(2, nclust)]
        with treeCl.utils.fileIO.fwriter(os.path.join(cache_dir, 'partitions.txt')) as fl:
            for p in ps:
                fl.write(repr(p)+'\n')
    return ps

def write_concats(collection, ps, path, index, logger, method, fastest):
    concat_dir = get_dirs(path, index)['concdir']
    if not os.path.exists(concat_dir): os.mkdir(concat_dir)

    # First concat is all records, no clustering
    outf = os.path.join(concat_dir, '1cl0.phy')
    outr = os.path.join(concat_dir, '1cl0.json')
    if not os.path.exists(outf):
        conc = collection.concatenate(range(len(collection)))
        al = conc.alignment
        al.write_alignment(outf, 'phylip', True)
        AlignIO.convert(outf, 'phylip-relaxed', '{}_'.format(outf), 'phylip-relaxed')
        os.system('mv {} {}'.format('{}_'.format(outf), outf))
    if not os.path.exists(outr):
        al = treeCl.alignment.Alignment(outf, 'phylip', True)
        if al.is_dna():
            if method == 'raxml':
                logger.info('Calculating fast raxml tree for dna alignment 1cl0.phy')
                result = run_raxml_fast(outf, fastest)
            else:
                logger.info('Calculating FastTree tree for dna alignment 1cl0.phy')
                result = run_fasttree(outf, True)

        else:
            if method == 'raxml':
                logger.info('Calculating fast raxml tree for protein alignment 1cl0.phy')
                result = run_raxml_fast(outf, 'PROTGAMMALG', fastest)
            else:
                logger.info('Calculating FastTree tree for protein alignment 1cl0.phy')
                result = run_fasttree(outf, False)

        with fileIO.fwriter(outr) as fl:
            json.dump(result, fl)

    # Now do rest of concats
    for ngrp, p in enumerate(ps, start=2):
        for i, grp in enumerate(p.get_membership()):
            outf = os.path.join(concat_dir, '{}cl{}.phy'.format(ngrp, i))
            outr = os.path.join(concat_dir, '{}cl{}.json'.format(ngrp, i))
            if not os.path.exists(outf):
                conc = collection.concatenate(grp)
                al = conc.alignment
                al.write_alignment(outf, 'phylip', True)
                AlignIO.convert(outf, 'phylip-relaxed', '{}_'.format(outf), 'phylip-relaxed')
                os.system('mv {} {}'.format('{}_'.format(outf), outf))
            if not os.path.exists(outr):
                al = treeCl.alignment.Alignment(outf, 'phylip', True)
                if al.is_dna():
                    if method == 'raxml':
                        logger.info('Calculating fast raxml tree for dna alignment {}cl{}.phy'.format(ngrp, i))
                        result = run_raxml_fast(outf, fastest)
                    else:
                        logger.info('Calculating FastTree tree for dna alignment {}cl{}.phy'.format(ngrp, i))
                        result = run_fasttree(outf, True)

                else:
                    if method == 'raxml':
                        logger.info('Calculating fast raxml tree for protein alignment {}cl{}.phy'.format(ngrp, i))
                        result = run_raxml_fast(outf, 'PROTGAMMALG', fastest)
                    else:
                        logger.info('Calculating FastTree tree for protein alignment {}cl{}.phy'.format(ngrp, i))
                        result = run_fasttree(outf, False)

                with fileIO.fwriter(outr) as fl:
                    json.dump(result, fl)


def set_logger():
    logger = logging.getLogger('treeCl::npbs')
    logger.setLevel(logging.DEBUG)
    formatter = logging.Formatter('%(asctime)s - %(name)s - %(levelname)s - %(message)s')
    # create console handler with a higher log level
    console = logging.StreamHandler()
    console.setLevel(logging.INFO)
    console.setFormatter(formatter)
    # add the handlers to the logger
    logger.addHandler(console)
    return logger

def run_fasttree(fl, dna=False):
    fl = os.path.abspath(fl)
    with fileIO.TempDir() as tmpd, fileIO.ChDir(tmpd):
        fst = FastTree(verbose=False)
        cmd = '-gtr -gamma -pseudo -out tree.txt {} {}'.format('-nt' if dna else '', fl)
        fst(cmd, wait=True)
        with open('tree.txt') as treefl_handle:
            tree = treefl_handle.read().rstrip()
    result = parse_fasttree_output(fst.get_stderr())
    result['tree'] = treeCl.Tree(tree).as_string('newick',internal_labels=False, suppress_rooting=True).rstrip()
    return result

def run_raxml_fast(fl, model='GTRGAMMA', fastest=True):
    fl = os.path.abspath(fl)
    with fileIO.TempDir() as tmpd, fileIO.ChDir(tmpd):
        rax = Raxml(verbose=False)
        if fastest:
            rax(s=fl, m=model, p='123', n='initial_tree', y=True, wait=True)
            rax(s=fl, m=model, p='123', f='e', n='model_opt', t='RAxML_parsimonyTree.initial_tree', wait=True)
        else:
            rax(s=fl, m=model, p='123', f='F', n='initial_tree', wait=True)
            rax(s=fl, m=model, p='123', f='e', R='RAxML_binaryModelParameters.initial_tree', t='RAxML_fastTree.initial_tree', n='model_opt', wait=True)
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
        'likelihood': loglk,
        'alpha': alpha,
        'frequencies': [a, c, g, t],
        'rates': [ac, ag, at, cg, ct, gt],
    }

def parse_fasttree_output(s):
    try:
        loglk, alpha = (float(x) for x in re.search(r'Gamma\(\d+\) LogLk = ([0-9-.]+) alpha = ([0-9.]+)', s).groups())
    except AttributeError:
        logger.warn('Couldn\'t parse loglk and alpha')
        logger.info(s)
        return None

    return {
        'likelihood': loglk,
        'alpha': alpha,
    }

if __name__ == '__main__':
    import argparse
    logger = treeCl.logging.getLogger()
    logger.info('Handling args...')
    parser = argparse.ArgumentParser()
    parser.add_argument('-p', '--path', type=str, required=True)
    parser.add_argument('-i', '--index', type=str, required=True)
    parser.add_argument('-n', '--nclust', type=int, required=True)
    parser.add_argument('-f', '--fastest', action='store_true')
    parser.add_argument('-m', '--method', type=str, default='fasttree', choices=['fasttree', 'raxml'])
    args = parser.parse_args()
    logger.info('Path = {}, Index = {}, Method = {}'.format(args.path, args.index, args.method))


    logger.info('Writing npbs...')
    generate_npbs(args.path, args.index)


    logger.info('Getting a treeCl.Collection from the new npbs files...')
    c = get_collection(args.path, args.index)


    logger.info('Looking for precomputed trees...')
    trees_found = read_cache(c, args.path, args.index)
    if not trees_found:
        logger.info('Calc trees...')
        calc_trees(c, args.path, args.index, logger, args.method, args.fastest)
    else:
        logger.info('Loaded trees from cache')


    logger.info('Calc distances...')
    dm = calc_dists(c, args.path, args.index)


    logger.info('Cluster...')
    ps = cluster(dm, args.path, args.index, args.nclust)


    logger.info('Write up...')
    write_concats(c, ps, args.path, args.index, logger, args.method, args.fastest)

