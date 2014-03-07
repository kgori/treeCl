#!/usr/bin/env python

from ...remote.datastructs.seq import Seq
from ...remote.externals.phyml import Phyml, runPhyml
from ...remote.externals.treecollection import TreeCollection
from ...remote.errors import directorycheck
from ..externals.DVscript import runDV
from ....constants import TMPDIR
from trcl_tree import Tree, TrClTree
import re
import random

def sample_wr(population, k):
    _int = int
    _random = random.random
    n = len(population)
    return [population[_int(_random() * n)] for _ in range(k)]

class TrClSeq(Seq):

    """ A version of the SequenceRecord class with some extra functionality for
    working with tree inference packages, notably TreeCollection """

    def __init__(
        self,
        infile=None,
        file_format='fasta',
        name=None,
        datatype=None,
        headers=[],
        sequences=[],
        dv=None,
        tree=None,
        tmpdir=None,
        ):

        self.TCfiles = {}
        self.dv = dv or []
        if tree and isinstance(tree, Tree):
            self.tree = TrClTree.cast(tree)
        else:
            self.tree = None


        super(TrClSeq, self).__init__(
            infile,
            file_format,
            name,
            datatype,
            headers,
            sequences,
            tmpdir,
            )

    def __add__(self, other):
        """ Allows Records to be added together, concatenating sequences """

        self_set = set(self.headers)
        other_set = set(other.headers)
        if not self.datatype == other.datatype:
            print 'Trying to add sequences of different datatypes'
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
        self.tree = TrClTree.cast(p.run('nj', verbosity))
        return self.tree

    def bionj_plus(self, tmpdir, verbosity=0):
        """ Uses phyml (via treeCl.externals.tree_builders.Phyml) to build a
        bioNJ tree for the current record """
        if self.tmpdir is not None:
            tmpdir = self.tmpdir
        else:
            directorycheck(tmpdir)

        p = Phyml(self, tmpdir)
        self.tree = TrClTree.cast(p.run('lr', verbosity))
        return self.tree

    def bootstrap_sample(self):
        """ Samples with replacement from the columns of the alignment """
        columns = self._pivot(self.sequences)
        bootstrap_columns = sample_wr(columns, len(columns))
        bootstrap_sequences = self._pivot(bootstrap_columns)
        return self.__class__(headers=self.headers, sequences=bootstrap_sequences,
                        datatype=self.datatype)

    def dv_matrix(self, tmpdir):
        """ Uses darwin (via treeCl.externals.DVWrapper) to calculate pairwise
        distances and variances"""
        if self.tmpdir is not None:
            tmpdir = self.tmpdir
        else:
            directorycheck(tmpdir)

        runDV(self, tmpdir)

    def phyml(self, tmpdir, verbosity=0):
        """ Uses phyml (via treeCl.externals.tree_builders.Phyml) to build a
        full ML tree for the current record """
        if self.tmpdir is not None:
            tmpdir = self.tmpdir
        else:
            directorycheck(tmpdir)

        p = Phyml(self, tmpdir)
        self.tree = TrClTree.cast(p.run('ml', verbosity))
        return self.tree

    def likelihood(self, tree, tmpdir, dry_run=False, set_as_record_tree=False):
        if self.tmpdir is not None:
            tmpdir = self.tmpdir
        else:
            directorycheck(tmpdir)

        result = runPhyml(self, tmpdir, 'lk', tree=tree, dry_run=dry_run,
                          set_as_record_tree=set_as_record_tree)
        if dry_run:
            return result
        else:
            return result.score


    # def likelihood(self, tree, tmpdir):
    #     """ Uses phyml (via treeCl.externals.tree_builders.Phyml) to calculate
    #     the likelihood  for the current record """
    #     if self.tmpdir is not None:
    #         tmpdir = self.tmpdir
    #     else:
    #         directorycheck(tmpdir)

    #     p = Phyml(self, tmpdir)

    #     alignment_file = self.write_phylip('{0}/tmp_alignment.phy'.format(
    #         tmpdir), interleaved=True)
    #     newick_file = tree.write_to_file('{0}/tmp_tree.nwk'.format(tmpdir))

    #     p.add_tempfile(alignment_file)
    #     p.add_tempfile(newick_file)
    #     p.add_flag('-i', alignment_file)
    #     p.add_flag('-u', newick_file)
    #     p.add_flag('-b', '0')   # no bootstraps
    #     if self.datatype == 'protein':
    #         p.add_flag('-m', 'WAG') # evolutionary model
    #     else:
    #         p.add_flag('-m', 'GTR')
    #     p.add_flag('-o', 'n')   # no optimisation
    #     if self.datatype == 'protein':
    #         p.add_flag('-d', 'aa')  # datatype
    #     else:
    #         p.add_flag('-d', 'nt')
    #     p.add_flag('--no_memory_check', '')
    #     p.add_flag('--quiet', '')
    #     p.call() # run phyml
    #     (_, stats) = p.read(alignment_file)
    #     p.clean() # cleanup tempfiles
    #     score = float(re.compile('(?<=Log-likelihood: ).+').search(stats).group(0))
    #     return score

    def tree_collection_deprecated(self, tmpdir):
        """ DEPRECATED:   Uses TreeCollection (via
        treeCl.externals.tree_builders.TreeCollection) to build a least squares
        tree for the current record """
        if self.tmpdir is not None:
            tmpdir = self.tmpdir
        else:
            directorycheck(tmpdir)

        if self.dv <= []:
            self.dv_matrix()
        tc = TreeCollection(self, tmpdir)
        self.tree = TrClTree.cast(tc.run())

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
        labels_string   = '{0}\n{1}\n'.format(labels_len, ' '.join(all_labels))

        # distvar and genome_map need to be built up
        distvar_list    = [str(num_matrices)]
        genome_map_list = ['{0} {1}'.format(num_matrices, labels_len)]

        # build up lists to turn into strings
        for (i, (matrix, labels)) in enumerate(self.dv, start=1):
            labels = labels.split()
            dim = len(labels)
            distvar_list.append('{0} {0} {1}\n{2}'.format(dim, i, matrix))
            genome_map_entry = ' '.join((str(labels.index(lab) + 1)
                                        if lab in labels else '-1')
                                        for lab in all_labels)
            genome_map_list.append(genome_map_entry)

        distvar_string = '\n'.join(distvar_list)
        genome_map_string = '\n'.join(genome_map_list)

        return distvar_string, genome_map_string, labels_string

    def tree_collection(self, niters=5, quiet=True, tmpdir=None):
        import tree_collection

        if tmpdir is not None:
            tmpdir = tmpdir
        elif self.tmpdir is not None:
            tmpdir = self.tmpdir
        else:
            tmpdir = TMPDIR

        if not self.tree:
            self.bionj(tmpdir)

        if self.tree is None:
            raise Exception('Couldn\'t generate a BIONJ guide tree')

        if self.dv <= []:
            self.dv_matrix(tmpdir)

        gt = self.tree
        gt.reroot_at_midpoint()
        if not gt.is_rooted:
            raise Exception('Couldn\'t root the guide tree' )
        dv, gm, lab = self._get_tree_collection_strings()
        output_tree, score = tree_collection.compute(dv, gm, lab, gt.newick,
                                                     niters, quiet)
        result = TrClTree(output_tree, score, program='tree_collection',
            name=self.name, output='').scale(0.01)
        self.tree = result
        return result

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
        self._update()

    def shuffle(self):
        """ Modifies in-place """

        columns = self._pivot(self.sequences)
        random.shuffle(columns)
        self.sequences = self._pivot(columns)
        self._update()

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
        if items == []:
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
