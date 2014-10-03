#!/usr/bin/env python
from __future__ import print_function

# standard library
import os
import random
import re
import tempfile

# third party
import numpy as np

# treeCl
from seq import Seq
from ..tree import Tree
from ..interfacing import pll create_instance, PLLException, pll_to_dict, set_params_from_dict
from ..constants import TMPDIR
from ..errors import directorycheck
from ..utils.lazyprop import lazyprop


def sample_wr(population, k):
    _int = int
    _random = random.random
    n = len(population)
    return [population[_int(_random() * n)] for _ in range(k)]


class TrClSeq(Seq):
    """ A version of the SequenceRecord class with some extra functionality for
    working with tree inference packages, notably TreeCollection """

    def __init__(self, infile=None, file_format='fasta', name=None, datatype=None, headers=None, sequences=None,
                 dv=None, tree=None, tmpdir=None):

        super(TrClSeq, self).__init__(
            infile,
            file_format,
            name,
            datatype,
            headers,
            sequences,
            tmpdir,
        )

        self.TCfiles = {}
        self.dv = dv or []
        if tree and isinstance(tree, Tree):
            self.tree = tree
        else:
            self.tree = None

    def __add__(self, other):
        """ Allows Records to be added together, concatenating sequences """

        self_set = set(self.headers)
        other_set = set(other.headers)
        if not self.datatype == other.datatype:
            print('Trying to add sequences of different datatypes')
            return self
        union = self_set | other_set
        intersection = self_set & other_set
        only_in_self = self_set - other_set
        only_in_other = other_set - self_set

        d = {}
        for k in union:
            if k in intersection:
                d[k] = self.mapping[k] + other.mapping[k]
            elif k in only_in_self:
                d[k] = self.mapping[k] + 'X' * other.seqlength
            elif k in only_in_other:
                d[k] = 'X' * self.seqlength + other.mapping[k]
        dvsum = self.dv + other.dv
        return_object = self.__class__(headers=d.keys(), sequences=d.values(),
                                       datatype=self.datatype).sort_by_name(in_place=False)
        return_object.dv = dvsum
        return return_object

    def bionj(self, tmpdir, verbosity=0):
        """ Uses phyml (via treeCl.externals.tree_builders.Phyml) to build a
        bioNJ tree for the current record """

        if self.tmpdir is not None:
            tmpdir = self.tmpdir
        else:
            directorycheck(tmpdir)

        p = Phyml(self, tmpdir)
        self.tree = p.run('nj', verbosity)
        return self.tree

    def bionj_plus(self, tmpdir, verbosity=0):
        """ Uses phyml (via treeCl.externals.tree_builders.Phyml) to build a
        bioNJ tree for the current record """
        if self.tmpdir is not None:
            tmpdir = self.tmpdir
        else:
            directorycheck(tmpdir)

        p = Phyml(self, tmpdir)
        self.tree = p.run('lr', verbosity)
        return self.tree

    def bootstrap_sample(self, suffix=''):
        """ Samples with replacement from the columns of the alignment """
        columns = self._pivot(self.sequences)
        bootstrap_columns = sample_wr(columns, len(columns))
        bootstrap_sequences = self._pivot(bootstrap_columns)
        if self.name:
            name = '{}.{}'.format(self.name, suffix)
        else:
            name = None
        return self.__class__(headers=self.headers,
                              sequences=bootstrap_sequences,
                              datatype=self.datatype,
                              name=name)

    @lazyprop
    def distinct_sites(self):
        nsites = self._pivot(self.sequences)
        return len(set(nsites))

    def get_dv_matrix(self, verbosity=0):
        """ Uses darwin (via treeCl.externals.DVWrapper) to calculate pairwise
        distances and variances"""
        runDV(self, verbosity)

    def phyml(self, tmpdir, verbosity=0):
        """ Uses phyml (via treeCl.externals.tree_builders.Phyml) to build a
        full ML tree for the current record """
        if self.tmpdir is not None:
            tmpdir = self.tmpdir
        else:
            directorycheck(tmpdir)

        p = Phyml(self, tmpdir)
        self.tree = p.run('ml', verbosity)
        return self.tree

    def pll_get_instance(self, *args):
        delete = False
        try:
            with open(self.infile):
                pass
            alignment = self.infile
        except (IOError, TypeError):
            tmpdir = tempfile.mkdtemp()
            _, alignment = tempfile.mkstemp(dir=tmpdir)
            self.write_phylip(alignment, interleaved=True)
            delete = True

        try:
            instance = create_instance(alignment, *args)
            return instance
        finally:
            if delete:
                os.remove(alignment)
                os.rmdir(tmpdir)

    def pll_optimise(self, partitions, model, nthreads=1, opt_subst=True, seed=None):
        """
        Runs the full raxml search algorithm. Model parameters are set using a combination of partitions and model.
        Optimisation of substitution model parameters is enabled with opt_subst=True.
        :param partitions: Links substitution models to alignment sites. Format is the same as the file RAxML uses
         with the -q flag (e.g. DNA, gene1codon1 = 1-500/3 - see RAxML manual)
        :param model: Dictionary to set model parameters--rates, frequencies and gamma alpha parameters--for each
         partition.
        :param opt_subst: bool: optimise substitution model parameters (T/F).
        :return: Dictionary of optimisation results.
        """

        instance = None
        seed = seed or random.randint(1, 99999)
        tree = self.tree.newick if self.tree else True
        try:
            instance = self.pll_get_instance(partitions, tree, nthreads, seed)
            set_params_from_dict(instance, model)
            instance.optimise(opt_subst)
            return pll_to_dict(instance)
        except ValueError as exc:
            raise exc
        except Exception as exc:
            raise PLLException(exc.message)
        finally:
            if instance:
                del instance

    def pll_optimise_model(self, partitions, model, nthreads=1, seed=None):
        instance = None
        seed = seed or random.randint(1, 99999)
        tree = self.tree.newick if self.tree else True
        try:
            instance = self.pll_get_instance(partitions, tree, nthreads, seed)
            set_params_from_dict(instance, model)
            instance.optimise_model()
            return pll_to_dict(instance)
        except ValueError as exc:
            raise exc
        except Exception as exc:
            raise PLLException(exc.message)
        finally:
            if instance:
                del instance

    def pll_optimise_branches(self, partitions, model, nthreads=1, n_iter=50, seed=None):
        instance = None
        seed = seed or random.randint(1, 99999)
        tree = self.tree.newick if self.tree else True
        try:
            instance = self.pll_get_instance(partitions, tree, nthreads, seed)
            set_params_from_dict(instance, model)
            instance.optimise_branch_lengths(n_iter)
            return pll_to_dict(instance)
        except ValueError as exc:
            raise exc
        except Exception as exc:
            raise PLLException(exc.message)
        finally:
            if instance:
                del instance

    def pll_get_parsimony_tree(self, seed=None):
        instance = None
        seed = seed or random.randint(1, 99999)
        model = 'DNA' if self.guess_datatype() == 'dna' else 'LG'
        partition = '{}, x = 1-{}'.format(model, self.seqlength)
        try:
            instance = self.pll_get_instance(partition, True, 1, seed)
            return Tree(instance.get_tree())
        except ValueError as exc:
            raise exc
        except Exception as exc:
            raise PLLException(exc.message)

        finally:
            if instance:
                del instance

    def _get_tree_collection_strings(self):
        """ Function to get input strings for tree_collection
        tree_collection needs distvar, genome_map and labels -
        these are returned in the order above
        """
        if self.dv <= []:
            return None

        # aliases
        dv = self.dv
        num_matrices = len(dv)
        all_labels = self.headers
        labels_len = len(all_labels)

        # labels string can be built straight away
        labels_string = '{0}\n{1}\n'.format(labels_len, ' '.join(all_labels))

        # distvar and genome_map need to be built up
        distvar_list = [str(num_matrices)]
        genome_map_list = ['{0} {1}'.format(num_matrices, labels_len)]

        # build up lists to turn into strings
        for (i, (matrix, labels)) in enumerate(self.dv, start=1):
            labels = labels.split()
            dim = len(labels)
            if isinstance(matrix, np.ndarray):
                matrix_string = '\n'.join([' '.join(str(x) for x in row)
                                           for row in matrix]) + '\n'
            else:
                matrix_string = matrix
            distvar_list.append('{0} {0} {1}\n{2}'.format(dim, i,
                                                          matrix_string))
            genome_map_entry = ' '.join((str(labels.index(lab) + 1)
                                         if lab in labels else '-1')
                                        for lab in all_labels)
            genome_map_list.append(genome_map_entry)

        distvar_string = '\n'.join(distvar_list)
        genome_map_string = '\n'.join(genome_map_list)

        return distvar_string, genome_map_string, labels_string

    def tree_collection(self,
                        niters=5,
                        keep_topology=False,
                        quiet=True,
                        tmpdir=None,
                        taxon_set=None,
                        guide_tree=None,
                        set_as_record_tree=True):

        import tree_collection

        if tmpdir is not None:
            tmpdir = tmpdir
        elif self.tmpdir is not None:
            tmpdir = self.tmpdir
        else:
            tmpdir = TMPDIR

        gt = None

        if guide_tree is not None:
            gt = guide_tree

        elif self.tree is not None:
            gt = self.tree

        else:
            self.bionj(tmpdir)
            gt = self.tree

        if gt is None:
            raise Exception('Couldn\'t generate a guide tree')

        if self.dv <= []:
            self.dv_matrix(tmpdir)

        if not gt.is_rooted:
            gt.reroot_at_midpoint()
        if not gt.is_rooted:
            raise Exception('Couldn\'t root the guide tree')

        dv, gm, lab = self._get_tree_collection_strings()
        output_tree, score = tree_collection.compute(dv, gm, lab, gt.newick,
                                                     niters, keep_topology,
                                                     quiet)
        if taxon_set is not None:
            result = Tree(output_tree, score, program='tree_collection',
                          name=self.name, output='',
                          taxon_set=taxon_set)
        else:
            result = Tree(output_tree, score, program='tree_collection',
                          name=self.name, output='')

        if set_as_record_tree:
            self.tree = result

        return result

    def distance_fit(self, tree):
        import tree_collection

        if self.dv <= []:
            self.dv_matrix(tmpdir)

        if not tree.is_rooted:
            tree.reroot_at_midpoint()

        if not tree.is_rooted:
            raise Exception('Couldn\'t root the test tree')

        dv, gm, lab = self._get_tree_collection_strings()
        return tree_collection.fit(dv, gm, lab, tree.newick)

    def _pivot(self, lst):
        new_lst = zip(*lst)
        return [''.join(x) for x in new_lst]

    def sanitise(self):
        self.sort_by_name()
        l = []
        for h in self.headers:
            if '/' in h:
                h = h[:h.index('/')]
            h = h.strip().replace(' ', '_')
            l.append(h)
        self.headers = l
        self.sequences = [seq.upper() for seq in self.sequences]
        self.update()

    def shuffle(self):
        """ Modifies in-place """

        columns = self._pivot(self.sequences)
        random.shuffle(columns)
        self.sequences = self._pivot(columns)
        self.update()

    def sort_by_length(self, in_place=True):
        """ Sorts sequences by descending order of length Uses zip as its own
        inverse [ zip(*zip(A,B)) == (A,B) ] Gaps and 'N' characters are not
        counted """

        (h, s) = zip(*sorted(zip(self.headers, self.sequences),
                             key=lambda item: len(item[1].replace('-', '').replace('N',
                                                                                   '')), reverse=True))
        if in_place:
            self.headers = h
            self.sequences = s
        else:
            return self.__class__(name=self.name, headers=h, sequences=s,
                                  datatype=self.datatype)

    def sort_by_name(self, in_place=True):
        """ Sorts sequences by name, treating numbers as integers (i.e. sorting
        like this: 1, 2, 3, 10, 20 not 1, 10, 2, 20, 3). If in_place = False the
        sorting doesn't mutate the underlying object, and the output is returned
        If in_place = True the sorting mutates the self object """

        items = self.mapping.items()
        if not items:
            return self
        sort_key = lambda item: tuple((int(num) if num else alpha) for (num,
                                                                        alpha) in re.findall(r'(\d+)|(\D+)',
                                                                                             item[0]))
        items = sorted(items, key=sort_key)
        (h, s) = zip(*items)
        if in_place:
            self.headers = h
            self.sequences = s
            return self
        else:
            return self.__class__(name=self.name, headers=h, sequences=s,
                                  datatype=self.datatype)

    def split_by_lengths(self, lengths, names=None):
        assert sum(lengths) == self.seqlength
        columns = self._pivot(self.sequences)
        newcols = []
        for l in lengths:
            newcols.append(columns[:l])
            columns = columns[l:]
        newrecs = []
        for col in newcols:
            newseqs = self._pivot(col)
            newrec = self.__class__(headers=self.headers, sequences=newseqs,
                                    datatype=self.datatype)
            newrecs.append(newrec)
        if names:
            for (i, newrec) in enumerate(newrecs):
                newrec.name = names[i]
        else:
            for (i, newrec) in enumerate(newrecs):
                newrec.name = 'record_{0}'.format(i + 1)
        return newrecs
