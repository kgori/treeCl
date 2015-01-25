import numpy as np
from treeCl.alignment import Alignment
from treeCl.constants import PLL_RANDOM_SEED
from treeCl.tasks import tasks
from treeCl.tree import Tree
from treeCl.utils.decorators import lazyprop

__author__ = 'kgori'


class Concatenation(object):
    """docstring for Concatenation"""

    def __init__(self, collection, indices):
        super(Concatenation, self).__init__()
        if any((x > len(collection)) for x in indices):
            raise ValueError('Index out of bounds in {}'.format(indices))
        if any((x < 0) for x in indices) < 0:
            raise ValueError('Index out of bounds in {}'.format(indices))
        if any((not isinstance(x, int)) for x in indices):
            raise ValueError('Integers only in indices, please: {}'
                             .format(indices))
        self.collection = collection
        self.indices = sorted(indices)

    @lazyprop
    def distances(self):
        return [self.collection.records[i].parameters.partitions.distances for i in self.indices]

    @lazyprop
    def variances(self):
        return [self.collection.records[i].parameters.partitions.variances for i in self.indices]

    @lazyprop
    def frequencies(self):
        return [self.collection.records[i].parameters.partitions.frequencies for i in self.indices]

    @lazyprop
    def datatypes(self):
        return ['dna' if self.collection.records[i].is_dna() else 'protein' for i in self.indices]

    @lazyprop
    def alignment(self):
        al = Alignment([self.collection[i] for i in self.indices])
        al.fast_compute_distances()
        return al

    @lazyprop
    def names(self):
        return [self.collection.records[i].name for i in self.indices]

    @lazyprop
    def lengths(self):
        return [len(self.collection.records[i]) for i in self.indices]

    @lazyprop
    def headers(self):
        return [self.collection.records[i].get_names() for i in self.indices]

    @lazyprop
    def coverage(self):
        total = float(self.collection.num_species())
        return [len(self.collection.records[i]) / total for i in self.indices]

    @lazyprop
    def trees(self):
        return [self.collection.records[i].tree for i in self.indices]

    @lazyprop
    def mrp_tree(self):
        trees = [tree.newick if hasattr('newick', tree) else tree for tree in self.trees]
        return Alignment().get_mrp_supertree(trees)

    def get_tree_collection_strings(self, scale=1):
        """ Function to get input strings for tree_collection
        tree_collection needs distvar, genome_map and labels -
        these are returned in the order above
        """

        # aliases
        num_matrices = len(self.distances)
        label_set = reduce(lambda x, y: x.union(y), (set(l) for l in self.headers))
        labels_len = len(label_set)

        # labels string can be built straight away
        labels_string = '{0}\n{1}\n'.format(labels_len, ' '.join(label_set))

        # distvar and genome_map need to be built up
        distvar_list = [str(num_matrices)]
        genome_map_list = ['{0} {1}'.format(num_matrices, labels_len)]

        # build up lists to turn into strings
        for i in range(num_matrices):
            labels = self.headers[i]
            dim = len(labels)
            dmatrix = np.array(self.distances[i])
            vmatrix = np.array(self.variances[i])
            matrix = np.zeros(dmatrix.shape)
            matrix[np.triu_indices(len(dmatrix), 1)] = dmatrix[np.triu_indices(len(dmatrix), 1)]
            matrix[np.tril_indices(len(vmatrix), -1)] = vmatrix[np.tril_indices(len(vmatrix), -1)]
            if scale:
                matrix[np.triu_indices(dim, 1)] *= scale
                matrix[np.tril_indices(dim, -1)] *= scale * scale

            if isinstance(matrix, np.ndarray):
                matrix_string = '\n'.join([' '.join(str(x) for x in row)
                                           for row in matrix]) + '\n'
            else:
                matrix_string = matrix
            distvar_list.append('{0} {0} {1}\n{2}'.format(dim, i + 1,
                                                          matrix_string))
            genome_map_entry = ' '.join((str(labels.index(lab) + 1)
                                         if lab in labels else '-1')
                                        for lab in label_set)
            genome_map_list.append(genome_map_entry)

        distvar_string = '\n'.join(distvar_list)
        genome_map_string = '\n'.join(genome_map_list)

        guide_tree = Tree(self.mrp_tree)

        for e in guide_tree.postorder_edge_iter():
            if e.length is None:
                if e.head_node == guide_tree.seed_node:
                    e.length = 0.0
                else:
                    e.length = np.random.uniform()

        if not guide_tree.is_rooted:
            guide_tree.reroot_at_midpoint()
        if not guide_tree.is_rooted:
            raise Exception('Couldn\'t root the guide tree')
        tree_string = guide_tree.scale(scale).newick

        return distvar_string, genome_map_string, labels_string, tree_string

    def qfile(self, dna_model='DNA', protein_model='LG', sep_codon_pos=False,
              ml_freqs=False, eq_freqs=False):
        from_ = 1
        to_ = 0
        qs = list()
        if ml_freqs:
            dna_model += 'X'
            protein_model += 'X'
        if eq_freqs and not ml_freqs:
            protein_model += 'F'

        models = dict(dna=dna_model, protein=protein_model)
        for length, name, datatype in zip(self.lengths, self.names,
                                          self.datatypes):
            to_ += length
            if datatype == 'dna' and sep_codon_pos:
                qs.append('{}, {} = {}-{}/3'.format(models[datatype], name, from_,
                                                    to_))
                qs.append('{}, {} = {}-{}/3'.format(models[datatype], name, from_ + 1,
                                                    to_))
                qs.append('{}, {} = {}-{}/3'.format(models[datatype], name, from_ + 2,
                                                    to_))
            else:
                qs.append('{}, {} = {}-{}'.format(models[datatype], name, from_,
                                                  to_))
            from_ += length
        return '\n'.join(qs)

    def paml_partitions(self):
        return 'G {} {}'.format(len(self.lengths),
                                ' '.join(str(x) for x in self.lengths))
