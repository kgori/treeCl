#!/usr/bin/env python

from . import ExternalSoftware
from ..utils import dpy as utils_dpy, fileIO
from ..errors import FileError, filecheck_and_raise
from treeCl.tree import Tree
from ..utils.printing import print_and_return
import re

local_dir = fileIO.path_to(__file__)


class TreeSoftware(ExternalSoftware):

    def __init__(self, record=None, supplied_binary=''):
        super(TreeSoftware, self).__init__(supplied_binary)
        self.record = record
        self.tmpdir = record.tmpdir or '/tmp'

    @property
    def record(self):
        return self._record

    @record.setter
    def record(self, sequence_record):
        self._record = sequence_record


class Phyml(TreeSoftware):

    default_binary = 'phyml'
    score_regex = re.compile('(?<=Log-likelihood: ).+')

    def read(self, filename):
        tree_filename = filename + '_phyml_tree.txt'
        stats_filename = filename + '_phyml_stats.txt'
        self.add_tempfile(tree_filename)
        self.add_tempfile(stats_filename)
        with open(tree_filename) as treefile:
            with open(stats_filename) as statsfile:
                return (treefile.read(), statsfile.read())

    def run(self, analysis=None, verbosity=0):
        if analysis:
            self.set_default_flags(analysis)
        else:
            analysis = fileIO.basename(self.binary)
        if verbosity > 1:
            print self.flags
            print 'Writing tempfiles to', self.tmpdir
        filename = self.write()
        filecheck_and_raise(filename)
        self.add_flag('-i', filename)
        if verbosity == 1:
            print_and_return('Running phyml on {0}'.format(self.record.name))
        elif verbosity > 1:
            print 'Running phyml on {0}'.format(self.record.name)
        (stdout, stderr) = self.call(verbose=(True if verbosity > 1 else False))
        (tree, stats) = self.read(filename)
        score = float(self.score_regex.search(stats).group(0))
        if verbosity > 1:
            print 'Cleaning tempfiles'
        self.clean()
        tree_object = Tree(tree, score, analysis, self.record.name, stats)
        self.record.tree = tree_object
        if verbosity > 1:
            print 'Done.'
        return tree_object

    def write(self):
        filename = self.record.get_name(default='tmp_phyml_input')
        filename = '{0}/{1}.phy'.format(self.tmpdir, filename)
        self.record.write_phylip(filename)
        self.add_tempfile(filename)
        return filename

    def read_datatype(self, datatype=None):
        datatype = datatype or self.record.datatype
        if datatype == 'protein':
            return {'-d': 'aa', '-m': 'WAG'}
        elif datatype == 'dna':
            return {'-d': 'nt', '-m': 'GTR'}

    def set_default_flags(self, analysis='ml', datatype=None):

        defaults = self.read_datatype(datatype=datatype)
        if defaults:
            defaults['-a'] = 'e'
            defaults['-b'] = 0
            defaults['-c'] = 4
            defaults['-q'] = ''
            defaults['--quiet'] = ''
            if analysis == 'ml' or analysis == 'full':
                defaults['-o'] = 'tlr'
            elif analysis == 'nj' or analysis == 'bionj':
                defaults['-o'] = 'n'

            for flag in defaults:
                self.add_flag(flag, defaults[flag])


class TreeCollection(TreeSoftware):

    default_binary = 'TreeCollection'

    def read(self, output):
        output = output.split()
        score = float(output[-1])
        tree = output[-2]
        return (score, tree)

    def run(self, guidetree=None, verbosity=0):

        guidetree_file = self.write_guidetree(guidetree)
        tmpfiles = self.write()
        self.add_flag('-D', tmpfiles['dv'])
        self.add_flag('-M', tmpfiles['map'])
        self.add_flag('-L', tmpfiles['lab'])
        self.add_flag('-T', guidetree_file)
        if verbosity == 1:
            print_and_return('Running TreeCollection on {0}'.format(self.record.name))
        elif verbosity > 1:
            print 'Running TreeCollection on {0}'.format(self.record.name)
        (stdout, stderr) = self.call()
        self.clean()
        if verbosity > 1:
            print stdout, stderr
        (score, tree) = self.read(stdout)
        tree_object = Tree(tree, score, fileIO.basename(self.binary),
                           self.record.name, stdout).scale(0.01)
        self.record.tree = tree_object
        return tree_object

    def write(self):
        """ Write the distance-variance (dv) file, the labels file and the map
        file all together, because they share information, and this way we only
        have to look it up once"""

        # Look up info

        dv_info = self.record.dv
        num_matrices = len(dv_info)
        if num_matrices == 0:
            print 'No distance-variance matrix available'
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
        labels_file.write('''{0}
{1}
'''.format(labels_len,
                          ' '.join(all_labels)))
        labels_file.close()

        # Write info

        for (i, (matrix, labels)) in enumerate(dv_info, start=1):
            labels = labels.split()
            dim = len(labels)
            dv_file.write('''{0} {0} {1}
{2}
'''.format(dim, i, matrix))
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
        p = Phyml(self.record)
        tree = p.run('nj')
        tree.newick = utils_dpy.bifurcate_base(tree.newick)
        return tree

    def write_guidetree(self, tree=None):
        tree = tree or self.nj_tree()
        base_filename = self.record.get_name(default='tmp_TC')
        filename = '{0}/{1}_tree.nwk'.format(self.tmpdir, base_filename)
        tree.write_to_file(filename)
        self.add_tempfile(filename)
        return filename


# RUNNERS

def runPhyml(rec, analysis, verbosity=0):
    p = Phyml(rec)
    return p.run(analysis, verbosity)


def runTC(rec, guidetrees=None, verbosity=0):
    if not isinstance(guidetrees, list):
        guidetrees = [guidetrees]
    tc = TreeCollection(rec)
    trees = [tc.run(guidetree, verbosity) for guidetree in guidetrees]
    return min(trees, key=lambda x: x.score)
