#!/usr/bin/env python
import os

def _get_dirs(path):
    working_dir = path
    cache_dir = os.path.join(working_dir, 'ft_cache')
    concat_dir = os.path.join(working_dir, 'ft_concat')
    return {'wdir': working_dir, 'cachedir': cache_dir, 'concdir': concat_dir}

def handle_args():
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument('-p', '--path', type=str, required=True)
    parser.add_argument('-t', '--threads', type=int, default=1)
    parser.add_argument('-n', '--nclust', type=int, required=True)
    args = parser.parse_args()
    return args

def get_collection(path, jobhandler):
    dirs = _get_dirs(path)
    working_dir = dirs['wdir']
    cache_dir = dirs['cachedir']
    try:
        return treeCl.Collection(input_dir=working_dir, file_format='phylip', param_dir=cache_dir)
    except:
        coll = treeCl.Collection(input_dir=working_dir, file_format='phylip')
        _calc_trees(coll, path, jobhandler)
        return coll

def _calc_trees(collection, path, jobhandler):
    cache_dir = _get_dirs(path)['cachedir']
    collection.fast_calc_trees(jobhandler=jobhandler)
    collection.write_parameters(cache_dir)

def calc_dists(collection, path, jobhandler):
    cache_dir = _get_dirs(path)['cachedir']
    try:
        geo = treeCl.DistanceMatrix.from_csv(os.path.join(cache_dir, 'geo.csv.gz'))
    except:
        geo = collection.get_inter_tree_distances('geo', jobhandler=jobhandler)
        geo.to_csv(os.path.join(cache_dir, 'geo.csv.gz'))
    return geo

def get_jobhandler(nthreads):
    if nthreads == 1:
        return treeCl.parutils.SequentialJobHandler()
    elif nthreads > 1:
        return treeCl.parutils.ThreadpoolJobHandler(nthreads)

def cluster(dm, path, nclust):
    from treeCl import Partition
    cache_dir = _get_dirs(path)['cachedir']
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

def write_concats(collection, ps, path):
    concat_dir = _get_dirs(path)['concdir']
    if not os.path.exists(concat_dir): os.mkdir(concat_dir)

    # First concat is all records, no clustering
    outf = os.path.join(concat_dir, '1cl0.phy')
    if not os.path.exists(outf):
        conc = collection.concatenate(range(len(collection)))
        al = conc.alignment
        al.write_alignment(outf, 'phylip', True)

    # Now do rest of concats
    for ngrp, p in enumerate(ps, start=2):
        for i, grp in enumerate(p.get_membership()):
            outf = os.path.join(concat_dir, '{}cl{}.phy'.format(ngrp, i))
            if not os.path.exists(outf):
                conc = collection.concatenate(grp)
                al = conc.alignment
                al.write_alignment(outf, 'phylip', True)


def calc_concat_trees(path, jobhandler):
    concat_dir = _get_dirs(path)['concdir']
    coll = treeCl.Collection(input_dir=concat_dir, file_format='phylip')
    coll.fast_calc_trees(jobhandler=jobhandler)
    coll.write_parameters(concat_dir)


if __name__ == '__main__':
    args = handle_args()

    import treeCl
    logger = treeCl.logging.getLogger(__name__)

    logger.info('Getting jobhandler')
    jh = get_jobhandler(args.threads)
    # logger.debug(str(jh))

    logger.info('Getting collection')
    coll = get_collection(args.path, jh)
    # logger.debug(str(coll))

    logger.info('Calculating distances')
    dm = calc_dists(coll, args.path, jh)
    # logger.debug(str(dm)[:10])

    logger.info('Clustering')
    ps = cluster(dm, args.path, args.nclust)
    # logger.debug(str(ps))

    logger.info('Writing concats')
    write_concats(coll, ps, args.path)

    logger.info('Calculating concat trees')
    calc_concat_trees(args.path, jh)

    logger.info('DONE.')

