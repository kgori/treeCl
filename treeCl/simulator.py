#!/usr/bin/env python
from __future__ import print_function, division
from builtins import zip
from builtins import str
from builtins import range
from builtins import object

# standard library
import numbers
import shutil
import sys

# third party
import numpy as np

# treeCl
from .partition import Partition
from .tree import Tree
from . import errors
from .utils import fileIO, print_and_return


class Simulator(object):
    """
    Simulate alignments from several trees.
    Args:
    class_list          = a list with an entry for each class, which is the
                        (integer) number of genes in that class
    permutations_list   = a list with an entry for each class, which is the
                        (integer) number of permutations the class tree has
                        relative to the master tree (see master_tree)
    num_species         = number of leaves on the master tree
    datatype            = 'dna' or 'protein'
                    """

    def __init__(
            self,
            class_list,
            permutations_list,
            nspecies,
            subst_model,
            rate_model,
            master_tree_generator_method='yule',
            master_tree=None,
            class_tree_permuter='nni',
            gene_length_kappa=1.7719,
            gene_length_theta=279.9,
            gene_length_min=10,
            gamma_rate_param=None,
            outdir='./',
            autocorrelated_relaxed_clock=False,
            uncorrelated_relaxed_clock=False,
            scale_rates=False,
            verbosity=0,
    ):
        # default
        errors.optioncheck(master_tree_generator_method, ['yule', 'coal',
                                                          'rtree', 'custom'])
        errors.optioncheck(class_tree_permuter, ['nni', 'spr', 'lgt', 'genetree'])
        if master_tree is None and master_tree_generator_method == 'custom':
            raise Exception('No custom tree was specified')
        self.num_classes = len(class_list)
        self.num_genes = sum(class_list)
        self.class_list = class_list
        self._master_tree = None
        self.verbosity = verbosity
        self.autocorrelated_relaxed_clock = autocorrelated_relaxed_clock
        self.uncorrelated_relaxed_clock = uncorrelated_relaxed_clock
        self.scale_rates = scale_rates
        self.gene_trees = list()
        if master_tree is None:
            tree = self.generate_master_tree(master_tree_generator_method,
                                             nspecies)
            self.master_tree = tree
            self.num_species = nspecies
        else:
            self.master_tree = master_tree
            if len(master_tree) != nspecies:
                msg = [
                    'Warning: supplied tree has {0} taxa.'.format(
                        len(master_tree)),
                    'Required number is {0}.\n'.format(nspecies),
                    'Resetting number of species to match the supplied tree.'
                ]
                print(''.join(msg))
                self.num_species = nspecies
        self.set_gene_lengths(gene_length_kappa, gene_length_theta,
                              gene_length_min)
        self.gamma_rate_param = gamma_rate_param
        self.permuter = class_tree_permuter
        self.permutations_list = permutations_list
        self.datatype = datatype
        self.tmpdir = errors.directorymake(tmpdir)
        self.outdir = outdir
        self.generate_class_trees()  # sets self.class_trees dict
        self.make_alf_dirs()  # sets self.alf_dirs dict
        self.write_alf_params()
        self.get_true_partition()

    @property
    def master_tree(self):
        return self._master_tree

    @master_tree.setter
    def master_tree(self, tree):
        self._master_tree = tree

    def generate_master_tree(self, method, nspecies):
        if method == 'yule':
            tree = Tree.new_yule(nspecies)
            tree.name = '{}_master_tree'.format(method)
            return tree
        elif method == 'coal':
            tree = Tree.new_coal(nspecies)
            tree.name = '{}_master_tree'.format(method)
            return tree
        elif method == 'rtree':
            tree = Tree.new_rtree(nspecies)
            tree.name = '{}_master_tree'.format(method)
            return tree

    def set_gene_lengths(
            self,
            kappa,
            theta,
            min_,
    ):
        self.gene_length_kappa = kappa
        self.gene_length_theta = theta
        self.gene_length_min = min_

    def generate_class_trees(self):
        class_trees = {}
        if self.permuter == 'genetree':
            for k in range(self.num_classes):
                class_trees[k + 1] = self.master_tree.sample_gene_tree(
                    scale_to=self.permutations_list[k])

        else:
            # Base trees for each class
            for k in range(self.num_classes):
                if self.permuter == 'nni':
                    t = self.master_tree.rnni(times=self.permutations_list[k])
                    t.name = 'class{}'.format(k + 1)
                    class_trees[k + 1] = t
                elif self.permuter == 'spr':
                    t = self.master_tree.rspr(times=self.permutations_list[k],
                                              disallow_sibling_sprs=True, keep_entire_edge=True)
                    t.name = 'class{}'.format(k + 1)
                    class_trees[k + 1] = t
                elif self.permuter == 'lgt':
                    t = self.master_tree.rlgt(times=self.permutations_list[k],
                                              disallow_sibling_lgts=True)
                    t.name = 'class{}'.format(k + 1)
                    class_trees[k + 1] = t

            # Expand base class trees into individual trees
            gene_trees = list()
            for k in range(self.num_classes):
                num_genes = self.class_list[k]
                trees = list()

                # populate the trees list
                for _ in range(num_genes):
                    class_tree = class_trees[k + 1]
                    tree = Tree(class_tree.newick)
                    tree.name = class_tree.name
                    trees.append(tree)

                # do per-tree rates/branch length adjustments
                for i, tree in enumerate(trees, start=1):
                    if self.autocorrelated_relaxed_clock:
                        tree.autocorrelated_relaxed_clock(1, 0.01)
                        for node in tree.postorder_node_iter():
                            node.edge_length *= node.rate
                        tree.name += '_{}'.format(i)

                    elif self.uncorrelated_relaxed_clock:
                        tree.uncorrelated_relaxed_clock(1, 0.01 * tree.length())
                        for node in tree.postorder_node_iter():
                            node.edge_length *= node.rate
                        tree.name += '_{}'.format(i)

                    elif self.scale_rates:
                        coeff = np.random.uniform(0.666, 1.333)
                        tree.scale(coeff, inplace=True)
                        tree.name += '_{}'.format(i)

                    else:
                        tree.name += '_{}'.format(i)

                gene_trees.extend(trees)

        self.class_trees = class_trees
        self.gene_trees = gene_trees

    def make_alf_dirs(self):
        alf_dirs = {}
        for i, g in enumerate(self.gene_trees, start=1):
            dirname = fileIO.join_path(self.tmpdir, g.name)
            alf_dirs[i] = errors.directorymake(dirname)
        self.alf_dirs = alf_dirs

    def make_alf_dirs_(self):
        """ DEPRECATED """
        alf_dirs = {}
        for k in range(self.num_classes):
            dirname = fileIO.join_path(self.tmpdir, 'class{0:0>1}'.format(
                k + 1))
            alf_dirs[k + 1] = errors.directorymake(dirname)
        self.alf_dirs = alf_dirs

    def write_alf_params(self):
        if not hasattr(self, 'alf_dirs'):
            self.make_alf_dirs()

        if not hasattr(self, 'class_trees'):
            self.generate_class_trees()

        alf_params = {}
        for i, tree in enumerate(self.gene_trees, start=1):
            alfdir = self.alf_dirs[i]
            datatype = self.datatype
            name = tree.name
            num_genes = 1
            seqlength = self.gene_length_min
            gene_length_kappa = self.gene_length_kappa
            gene_length_theta = self.gene_length_theta
            alf_obj = ALF(tree=tree,
                          datatype=datatype, num_genes=num_genes,
                          seqlength=seqlength, gene_length_kappa=gene_length_kappa,
                          gene_length_theta=gene_length_theta, name=name, tmpdir=alfdir)
            if isinstance(self.gamma_rate_param, numbers.Number):
                alf_obj.params.rate_variation(self.gamma_rate_param)
            if datatype == 'protein':
                alf_obj.params.one_word_model('WAG')
            else:
                alf_obj.params.jc_model()
            alf_params[i] = alf_obj

        self.alf_params = alf_params

    def write_alf_params_(self):
        """ DEPRECATED """
        if not hasattr(self, 'alf_dirs'):
            self.make_alf_dirs()

        if not hasattr(self, 'class_trees'):
            self.generate_class_trees()

        alf_params = {}
        for k in range(self.num_classes):
            alfdir = self.alf_dirs[k + 1]
            tree = self.class_trees[k + 1]
            datatype = self.datatype
            name = 'class{0}'.format(k + 1)
            num_genes = self.class_list[k]
            seqlength = self.gene_length_min
            gene_length_kappa = self.gene_length_kappa
            gene_length_theta = self.gene_length_theta
            alf_obj = ALF(tree=tree,
                          datatype=datatype, num_genes=num_genes,
                          seqlength=seqlength, gene_length_kappa=gene_length_kappa,
                          gene_length_theta=gene_length_theta, name=name, tmpdir=alfdir)
            if datatype == 'protein':
                alf_obj.params.one_word_model('WAG')
            else:
                alf_obj.params.jc_model()
            alf_params[k + 1] = alf_obj

        self.alf_params = alf_params

    def clean(self):
        if not hasattr(self, 'alf_dirs'):
            return
        for directory in self.alf_dirs.values():
            shutil.rmtree(directory)

    def run(self):
        all_records = []
        total_jobs = len(self.gene_trees)
        for i, tree in enumerate(self.gene_trees, start=1):
            if self.verbosity > 0:
                print_and_return('Simulating {} ({:.1f}%)'.format(tree.name,
                                                                  100 * i / total_jobs),
                                 sys.stderr)
            simulated_record = self.alf_params[i].run()[0]
            simulated_record.name = tree.name
            all_records.append(simulated_record)
        self.result = all_records
        self.clean()
        return all_records

    def run_(self):
        """ DEPRECATED """
        all_records = []
        for k in range(self.num_classes):
            simulated_records = self.alf_params[k + 1].run()
            names = ['class{0}_{1:0>{2}}'.format(k + 1, i,
                                                 len(str(self.class_list[k]))) for i in range(1,
                                                                                              len(
                                                                                                  simulated_records) + 1)]
            for (rec, name) in zip(simulated_records, names):
                rec.name = name
            all_records.extend(simulated_records)
        self.result = all_records
        self.clean()
        return all_records

    def write(self):
        if hasattr(self, 'result'):
            errors.directorymake(self.outdir)
            errors.directorymake(fileIO.join_path(self.outdir,
                                                  'base_class_trees'))
            errors.directorymake(fileIO.join_path(self.outdir,
                                                  'gene_trees'))
            for rec in self.result:
                filename = fileIO.join_path(self.outdir, rec.name) + '.phy'
                rec.write_phylip(filename, interleaved=True)
            for i in range(self.num_classes):
                tree = self.class_trees[i + 1]
                name = 'base_tree_class{0:0>{1}}.nwk'.format(i + 1,
                                                             len(str(self.num_classes)))
                filename = fileIO.join_path(self.outdir, 'base_class_trees',
                                            name)
                tree.write_to_file(filename)
            for i, tree in enumerate(self.gene_trees, start=1):
                filename = fileIO.join_path(self.outdir, 'gene_trees',
                                            tree.name + '.nwk')
                tree.write_to_file(filename)
            self.master_tree.write_to_file(fileIO.join_path(self.outdir,
                                                            'master_tree.nwk'))
            filename = fileIO.join_path(self.outdir, 'true_partition.txt')
            with open(filename, 'w') as partition_file:
                partition_file.write(repr(self.true_partition))

    def get_true_partition(self):
        l = []
        for k in range(len(self.class_list)):
            l.extend([k + 1] * self.class_list[k])
        self.true_partition = Partition(l)
        return self.true_partition


if __name__ == '__main__':

    import argparse

    prog = fileIO.basename(__file__)
    parser = argparse.ArgumentParser(description='{0}'.format(prog))
    parser.add_argument('classes', type=int, nargs='+')
    parser.add_argument('-p', '--permutations', type=int, nargs='+')
    parser.add_argument('-s', '--species', type=int, default=12)
    parser.add_argument('-d', '--datatype', type=str, default='protein')
    parser.add_argument('-g', '--tree_generator', type=str, default='yule')
    parser.add_argument('-t', '--tree', type=str)
    parser.add_argument('--permuter', type=str, default='lgt')
    parser.add_argument('-l', '--gamma_params', type=float, nargs=2,
                        default=(1.7719, 279.9))
    parser.add_argument('-m', '--min_length', type=str, default=10)
    parser.add_argument('--tmp', type=str, default='/tmp')
    parser.add_argument('-o', '--output', type=str)
    args = parser.parse_args()

    if args.permutations is None:
        args.permutations = [1 for _ in args.classes]

    sim = Simulator(
        class_list=args.classes,
        permutations_list=args.permutations,
        nspecies=args.species,
        datatype=args.datatype,
        master_tree_generator_method=args.tree_generator,
        master_tree=args.tree,
        class_tree_permuter=args.permuter,
        gene_length_kappa=args.gamma_params[0],
        gene_length_theta=args.gamma_params[1],
        gene_length_min=args.min_length,
        tmpdir=args.tmp,
        outdir=args.output)

    sim.run()
    recs = sim.result
    if args.output is not None:
        sim.write()
