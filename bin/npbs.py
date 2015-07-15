#!/usr/bin/env python

import treeCl
from treeCl.utils import fileIO
from treeCl.wrappers.phylogenetics import Raxml
import os
import glob
import re
import logging
import json

def get_dirs(path, i):
    working_dir = os.path.join(path, 'npbs{}'.format(i))
    cache_dir = os.path.join(working_dir, 'cache')
    concat_dir = os.path.join(working_dir, 'concat')
    return {'wdir': working_dir, 'cachedir': cache_dir, 'concdir': concat_dir}

def generate_npbs(path, i):
    c = treeCl.Collection(input_dir=path, file_format='phylip')
    npbs = c.permuted_copy()
    working_dir = get_dirs(path, i)['wdir']
    if not os.path.exists(working_dir):
        os.mkdir(working_dir)
    for rec in npbs:
        rec.write_alignment('{}.phy'.format(os.path.join(working_dir, rec.name)), 'phylip', True)

def get_collection(path, i):
    working_dir = get_dirs(path, i)['wdir']
    return treeCl.Collection(input_dir=working_dir, file_format='phylip')

def calc_trees(collection, path, i, logger):
    working_dir = get_dirs(path, i)['wdir']
    for rec in collection:
        f = rec.parameters.filename
        if rec.is_dna():
            logger.info('Calculating fast raxml tree for dna alignment {}'.format(f))
            result = run_raxml_fast(f)
        else:
            logger.info('Calculating fast raxml tree for protein alignment {}'.format(f))
            result = run_raxml_fast(f, 'PROTGAMMALGF')
        rec.parameters.ml_tree = result['tree']
        rec.parameters.likelihood = result['likelihood']
        rec.parameters.partitions.append(treeCl.parameters.PartitionParameters())
        rec.parameters.partitions.alpha = result['alpha']
        rec.parameters.partitions.rates = result['rates']
        rec.parameters.partitions.frequencies = result['frequencies']

    # collection.calc_trees()
    collection.write_parameters(os.path.join(working_dir, 'cache'))

def calc_dists(collection, path, i):
    geo = collection.get_inter_tree_distances('geo')
    cache_dir = get_dirs(path, i)['cachedir']
    geo.to_csv(os.path.join(cache_dir, 'geo.csv.gz'))
    return geo

def cluster(dm, path, index, nclust):
    cl = treeCl.Clustering(dm)
    decomp = cl.spectral_decomp(None, 'median')
    ps = [cl.spectral_cluster(i, decomp) for i in range(2, nclust)]
    cache_dir = get_dirs(path, index)['cachedir']
    with treeCl.utils.fileIO.fwriter(os.path.join(cache_dir, 'partitions.txt')) as fl:
        for p in ps:
            fl.write(repr(p)+'\n')
    return ps

def write_concats(collection, ps, path, index, logger):
    concat_dir = get_dirs(path, index)['concdir']
    if not os.path.exists(concat_dir): os.mkdir(concat_dir)
    
    # First concat is all records, no clustering
    conc = collection.concatenate(range(len(collection)))
    al = conc.alignment
    outf = os.path.join(concat_dir, '1cl0.phy')
    outr = os.path.join(concat_dir, '1cl0.json')
    al.write_alignment(outf, 'phylip', True)
    if al.is_dna():
        logger.info('Calculating fast raxml tree for dna alignment 1cl0.phy')
        result = run_raxml_fast(outf)
    else:
        logger.info('Calculating fast raxml tree for protein alignment 1cl0.phy')
        result = run_raxml_fast(outf, 'PROTGAMMALG')
    with fileIO.fwriter(outr) as fl:
        json.dump(result, fl)

    # Now do rest of concats
    for ngrp, p in enumerate(ps, start=2):
        for i, grp in enumerate(p.get_membership()):
            conc = collection.concatenate(grp)
            al = conc.alignment
            outf = os.path.join(concat_dir, '{}cl{}.phy'.format(ngrp, i))
            outr = os.path.join(concat_dir, '{}cl{}.json'.format(ngrp, i))
            al.write_alignment(outf, 'phylip', True)
            if al.is_dna():
                logger.info('Calculating fast raxml tree for dna alignment {}cl{}.phy'.format(ngrp, i))
                result = run_raxml_fast(outf)
            else:
                logger.info('Calculating fast raxml tree for protein alignment {}cl{}.phy'.format(ngrp, i))
                result = run_raxml_fast(outf, 'PROTGAMMALG')
            with fileIO.fwriter(outr) as fl:
                json.dump(result, fl)


def set_logger():
    logger = logging.getLogger('treeCl::npbs')
    logger.setLevel(logging.INFO)
    formatter = logging.Formatter('%(asctime)s - %(name)s - %(levelname)s - %(message)s')
    # create console handler with a higher log level
    console = logging.StreamHandler()
    console.setLevel(logging.INFO)
    console.setFormatter(formatter)
    # add the handlers to the logger
    logger.addHandler(console)
    return logger

def run_raxml_fast(fl, model='GTRGAMMA'):
    fl = os.path.abspath(fl)
    with fileIO.TempDir() as tmpd, fileIO.ChDir(tmpd):
        rax = Raxml(verbose=False)
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

if __name__ == '__main__':
    import argparse
    logger = set_logger()
    logger.info('Handling args...')
    parser = argparse.ArgumentParser()
    parser.add_argument('-p', '--path', type=str, required=True)
    parser.add_argument('-i', '--index', type=str, required=True)
    parser.add_argument('-n', '--nclust', type=int, required=True)
    args = parser.parse_args()

    logger.info('Writing npbs...')
    generate_npbs(args.path, args.index)
    logger.info('Getting a treeCl.Collection from the new npbs files...')
    c = get_collection(args.path, args.index)
    logger.info('Calc trees...')
    calc_trees(c, args.path, args.index, logger)
    logger.info('Calc distances...')
    dm = calc_dists(c, args.path, args.index)
    logger.info('Cluster...')
    ps = cluster(dm, args.path, args.index, args.nclust)
    logger.info('Write up...')
    write_concats(c, ps, args.path, args.index, logger)
