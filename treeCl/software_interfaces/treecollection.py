#!/usr/bin/env python
from __future__ import print_function

# treeCl
from external import TreeSoftware
from phyml import Phyml
from ..datastructs.tree import Tree
from ..errors import TreeBuildingError
from ..utils import fileIO
from ..utils.printing import print_and_return


class TreeCollection(TreeSoftware):
    """ __init__ takes a Seq sequence record as
    first (only) positional argument, and supplied_binary= and
    tmpdir= as keyword arguments """

    default_binary = 'TreeCollection'
    local_dir = fileIO.path_to(__file__)

    def read(self, output):
        output = output.split()
        score = float(output[-1])
        tree = output[-2]
        return (score, tree)

    def run(self, guidetree=None, verbosity=0, **kwargs):

        guidetree_file = self.write_guidetree(guidetree)
        tmpfiles = self.write()
        self.add_flag('-D', tmpfiles['dv'])
        self.add_flag('-M', tmpfiles['map'])
        self.add_flag('-L', tmpfiles['lab'])
        self.add_flag('-T', guidetree_file)
        if verbosity == 1:
            print_and_return('Running TreeCollection on {0}'.format(
                self.record.name))
        elif verbosity > 1:
            print('Running TreeCollection on {0}'.format(self.record.name))
        try:
            (stdout, stderr) = self.call()
            self.clean()
            if verbosity > 1:
                print(stdout, stderr)
            (score, tree) = self.read(stdout)
            tree_object = Tree(tree, score,
                               program=fileIO.basename(self.binary),
                               name=self.record.name, output=stdout,
                               **kwargs)
            self.record.tree = tree_object
            return tree_object
        except:
            raise TreeBuildingError(stderr, fileIO.basename(self.binary))

    def write(self):
        """ Write the distance-variance (dv) file, the labels file and the map
        file all together, because they share information, and this way we only
        have to look it up once"""

        # Look up info

        dv_info = self.record.dv
        num_matrices = len(dv_info)
        if num_matrices == 0:
            print('No distance-variance matrix available')
            return
        all_labels = self.record.headers
        labels_len = len(all_labels)

        base_filename = self.record.get_name(default='tmp_TC')
        f = {}
        f['dv'] = '{0}/{1}_dv.txt'.format(self.tmpdir, base_filename)
        f['lab'] = '{0}/{1}_labels.txt'.format(self.tmpdir, base_filename)
        f['map'] = '{0}/{1}_map.txt'.format(self.tmpdir, base_filename)
        dv_file = open(f['dv'], 'w')
        labels_file = open(f['lab'], 'w')
        map_file = open(f['map'], 'w')

        # Write headers

        dv_file.write('{0}\n'.format(num_matrices))
        map_file.write('{0} {1}\n'.format(num_matrices, labels_len))
        labels_file.write('{0}\n{1}\n'.format(labels_len,
                                              ' '.join(all_labels)))
        labels_file.close()

        # Write info

        for (i, (matrix, labels)) in enumerate(dv_info, start=1):
            labels = labels.split()
            dim = len(labels)
            dv_file.write('{0} {0} {1}\n{2}\n'.format(dim, i, matrix))
            dv_file.flush()
            for lab in all_labels:
                map_file.write(('{0} '.format(labels.index(lab) + 1) if lab
                                                                        in labels else '-1 '))
            map_file.write('\n')
            map_file.flush()
        dv_file.close()
        map_file.close()

        for k in f:
            self.add_tempfile(f[k])  # for cleanup

        return f

    def nj_tree(self):
        p = Phyml(self.record, self.tmpdir)
        tree = p.run('nj')
        tree.reroot_at_midpoint()
        return tree

    def write_guidetree(self, tree=None):
        tree = tree or self.nj_tree()
        base_filename = self.record.get_name(default='tmp_TC')
        filename = '{0}/{1}_tree.nwk'.format(self.tmpdir, base_filename)
        tree.write_to_file(filename)
        self.add_tempfile(filename)
        return filename


def runTC(rec, tmpdir, guidetrees=None, verbosity=0, **kwargs):
    if not isinstance(guidetrees, list):
        guidetrees = [guidetrees]
    tc = TreeCollection(rec, tmpdir)
    trees = [tc.run(guidetree, verbosity, **kwargs) for guidetree in guidetrees]
    return min(trees, key=lambda x: x.score)
