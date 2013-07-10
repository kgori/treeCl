#!/usr/bin/env python

import treeCl
from treeCl.lib.local.datastructs.trcl_tree import TrClTree
from treeCl.lib.local.externals.alf import ALF
from treeCl.lib.remote.utils import fileIO
from treeCl.lib.remote import errors


class Simulator(object):
    """docstring for Simulator"""
    
    def __init__(
        self, class_list, nspecies, master_tree_generator_method='yule', 
        master_tree=None, tmpdir='/tmp', class_tree_permuter='nni',
        num_permutations=0, gene_length_kappa=1, gene_length_theta=1,
        gene_length_min=10,
        ):
        errors.optioncheck(master_tree_generator_method, 
            ['yule', 'coal', 'rtree'])
        errors.optioncheck(class_tree_permuter, 
            ['nni', 'spr', 'genetree'])
        self.num_classes = len(class_list)
        self.num_genes = sum(class_list)
        self.class_list = class_list
        if master_tree is None:
            tree = self.generate_master_tree(master_tree_generator_method, 
                nspecies)
            self.master_tree = tree
            self.num_species = nspecies
        else:
            self.master_tree = master_tree
            if len(master_tree) != nspecies:
                msg = ['Warning: supplied tree has {0} taxa.'.format(
                    len(master_tree)),
                    'Required number is {0}.\n'.format(nspecies),
                    'Resetting number of species to match the supplied tree.']
                print ''.join(msg)
                self.num_species = nspecies
        self.set_gene_lengths(gene_length_kappa, gene_length_theta, 
            gene_length_min)
        self.tmpdir = tmpdir
        self.permuter = class_tree_permuter
        self.num_permutations = num_permutations
        self._class_trees = None
        
    @property
    def master_tree(self):
        return self._master_tree

    @master_tree.setter
    def master_tree(self, tree):
        self._master_tree = tree

    @property
    def class_trees(self):
        if self._class_trees is None:
            self._class_trees = self.generate_class_trees()
        return self._class_trees



    def generate_master_tree(self, method, nspecies):
        errors.optioncheck(method, ['yule', 'coal', 'rtree'])
        if method == 'yule':
            return TrClTree.new_yule(nspecies)
        elif method == 'coal':
            return TrClTree.new_coal(nspecies)
        elif method == 'rtree':
            return TrClTree.new_rtree(nspecies)

    def set_gene_lengths(self, kappa, theta, min_):
        self.gene_length_kappa = kappa
        self.gene_length_theta = theta
        self.gene_length_min = min_

    def generate_class_trees(self):
        class_trees = []
        if self.permuter == 'nni':            
            for k in range(self.num_classes):
                class_tree = self.master_tree.copy()
                for p in range(self.num_permutations):
                    class_tree.rnni(inplace=True)
                class_trees.append(class_tree)
        elif self.permuter == 'spr':
            for k in range(self.num_classes):
                class_tree = self.master_tree.copy()
                for p in range(self.num_permutations):
                    class_tree.rspr(inplace=True,
                        disallow_sibling_sprs=True)
                class_trees.append(class_tree)
        elif self.permuter == 'genetree':
            for k in range(self.num_classes):
                class_trees.append(self.master_tree.sample_gene_tree(
                    nspecies=self.num_species, scale_to=self.num_permutations))
        return class_trees



       




if __name__ == '__main__':

    import argparse
    prog = fileIO.basename(__file__)
    parser = argparse.ArgumentParser(description='{0}'.format(prog))
    parser.add_argument('-t', '--tree', type=str, default='yule')


    args = parser.parse_args()
    tree = args.tree
    print tree

