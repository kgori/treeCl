import treeCl


def exhaustive_cluster(dm, n):
    hierarchical = treeCl.clustering.Hierarchical(dm)
    try:
        spectral7 = treeCl.Spectral(dm, scale_option=treeCl.clustering.options.LOCAL_SCALE_MANUAL, manual_scale=7)
    except:
        spectral7 = None
    spectral = treeCl.Spectral(dm)
    mds = treeCl.MultidimensionalScaling(dm)
    
    results = {}
    for n in range(2, n+1):
        d = {}
        d['average'] = hierarchical.cluster(n, treeCl.clustering.linkage.AVERAGE)
        d['centroid'] = hierarchical.cluster(n, treeCl.clustering.linkage.AVERAGE)
        d['complete'] = hierarchical.cluster(n, treeCl.clustering.linkage.AVERAGE)
        d['median'] = hierarchical.cluster(n, treeCl.clustering.linkage.AVERAGE)
        d['single'] = hierarchical.cluster(n, treeCl.clustering.linkage.AVERAGE)
        d['ward'] = hierarchical.cluster(n, treeCl.clustering.linkage.AVERAGE)
        d['spectral'] = spectral.cluster(n)
        d['kpca'] = spectral.cluster(n, algo=treeCl.clustering.spectral.KPCA)
        if spectral7:
            d['spectral7'] = spectral7.cluster(n)
            d['kpca7'] = spectral7.cluster(n, algo=treeCl.clustering.spectral.KPCA)
        else:
            d['spectral7'] = d['spectral']
            d['kpca7'] = d['kpca']
        results[n] = d

    return results
