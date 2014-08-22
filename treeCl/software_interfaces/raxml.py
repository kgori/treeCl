#!/usr/bin/env python
from __future__ import print_function
# standard library
import os
import re
import random

# treeCl
from external import TreeSoftware
from ..datastructs.tree import Tree
from ..datastructs.seq import qfile
from ..errors import filecheck
from ..utils import fileIO
from ..utils.printing import print_and_return

PATTERNS_PER_THREAD_DNA = 500
PATTERNS_PER_THREAD_AA = 200


def rstring(length, numOnly=False, letOnly=False):
    """ Generate a random alphanumeric string of defined length.  """

    numbers = '01234567890123456789'  # double up (bc. up- and lo-case letters)
    letters = 'abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ'
    if numOnly:
        alphabet = numbers
    elif letOnly:
        alphabet = letters
    else:
        alphabet = letters + numbers

    return ''.join(random.choice(alphabet) for _ in range(length))


class Raxml(TreeSoftware):
    """ __init__ takes a Seq sequence record as
    first (only) positional argument, and supplied_binary= and
    tmpdir= as keyword arguments """

    default_binary = 'raxml'
    score_regex = re.compile('(?<=Final GAMMA-based Score of best tree).+')
    local_dir = fileIO.path_to(__file__)

    def get_num_threads(self):
        import multiprocessing

        n_patterns = self.record.n_site_patterns()
        if self.record.datatype == 'dna':
            return min((n_patterns // PATTERNS_PER_THREAD_DNA + 1),
                       multiprocessing.cpu_count())
        return min((n_patterns // PATTERNS_PER_THREAD_AA + 1),
                   multiprocessing.cpu_count())

    def read(self, name_suffix):
        tree_filename = ('{0}/RAxML_bestTree.{1}'
                         .format(self.tmpdir, name_suffix))
        info_filename = ('{0}/RAxML_info.{1}'
                         .format(self.tmpdir, name_suffix))
        pars_filename = ('{0}/RAxML_parsimonyTree.{1}'
                         .format(self.tmpdir, name_suffix))
        log_filename = ('{0}/RAxML_log.{1}'
                        .format(self.tmpdir, name_suffix))
        result_filename = ('{0}/RAxML_result.{1}'
                           .format(self.tmpdir, name_suffix))
        for filename in [tree_filename, info_filename, pars_filename,
                         log_filename, result_filename]:
            self.add_tempfile(filename)
        with open(tree_filename) as treefile:
            with open(info_filename) as infofile:
                return (treefile.read(), infofile.read())

    def run(self, name=None, seed=None, threads=None, bootstrap=False,
            datatype=None, ml_freqs=False, verbosity=0, qfile_string=None,
            **kwargs):
        name = (self.record.name[:10] + rstring(4)
                if self.record.name
                else rstring(14))
        seed = seed or rstring(6, numOnly=True)
        threads = threads or self.get_num_threads()
        datatype = datatype or self.record.datatype
        if verbosity > 1:
            print(self.flags)
            print('Writing tempfiles to', self.tmpdir)
        filename, qfilename = self.write(qfile_string, datatype, ml_freqs)
        filecheck(filename)
        filecheck(qfilename)
        self.set_default_flags()
        self.add_flag('-s', filename)
        self.add_flag('-q', qfilename)
        self.add_flag('-n', name)
        self.add_flag('-p', seed)
        self.add_flag('-T', threads)
        if bootstrap:
            raise Exception('Not implemented bootstrapping yet')

        # OUTPUT TO USER
        if verbosity == 1:
            print_and_return('Running raxml on {0}'.format(self.record.name))
        elif verbosity > 1:
            print('Running raxml on {0}'.format(self.record.name))

        # DRY RUN - just get command string
        if kwargs.get('dry_run', False):
            cmd = self.call(verbose=(True if verbosity > 1 else False),
                            dry_run=True)
            return cmd

        # RUN RAXML
        curr_dir = os.getcwd()
        os.chdir(self.tmpdir)
        (stdout, stderr) = self.call(verbose=(True if verbosity > 1 else False))
        (tree, info) = self.read(name)

        try:
            score = float(self.score_regex.search(info).group(0))
        except:
            print(tree)
            print(info)
        if verbosity > 1:
            print('Cleaning tempfiles')
        self.clean()
        os.chdir(curr_dir)
        tree_object = Tree(newick=tree, score=score, program='raxml',
                           name=self.record.name, output=info, **kwargs)
        if kwargs.get('set_as_record_tree', True):
            self.record.tree = tree_object
        if verbosity > 1:
            print('Done.')
        return tree_object

    def write(self, qfile_string=None, datatype=None, ml_freqs=False):
        """ qfile_string should be a raxml partitions file, given as a string """
        record_name = self.record.get_name(default='tmp_raxml_input')
        filename = '{0}/{1}.phy'.format(self.tmpdir, record_name)
        datatype = datatype or self.record.datatype

        if qfile_string is None:
            model = ('GTR' if datatype == 'dna' else 'WAG')
            qfile_string = qfile([self.record], model, ml_freqs)
        qfilename = '{0}/{1}_qfile.txt'.format(self.tmpdir, record_name)
        with open(qfilename, 'w') as file_:
            file_.write(qfile_string)
        self.record.write_phylip(filename)
        self.add_tempfile(filename)
        self.add_tempfile(qfilename)
        return filename, qfilename

    def read_datatype(self, datatype=None):
        datatype = datatype or self.record.datatype
        if datatype == 'protein':
            return {'-m': 'PROTCATWAG'}
        elif datatype == 'dna':
            return {'-m': 'GTRCAT'}

    def set_default_flags(self, datatype=None):
        defaults = self.read_datatype(datatype=datatype)
        for flag in defaults:
            self.add_flag(flag, defaults[flag])


            # def run_raxml(input, threads, )