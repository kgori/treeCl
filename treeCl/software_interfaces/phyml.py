#!/usr/bin/env python
from __future__ import print_function
# standard library
import os
import re
import shutil
import tempfile
import time

# third party
from bsub import bsub

# treeCl
from external import ExternalSoftware, TreeSoftware
from ..constants import PHYML_MEMORY_MULTIPLIER, PHYML_MEMORY_SPARE, \
    PHYML_MEMORY_MIN, PHYML_MEMORY_STANDARD
from ..datastructs.tree import Tree
from ..errors import filecheck, optioncheck
from ..utils import fileIO
from ..utils.printing import print_and_return


ANALYSES = ['tlr', 'lr', 'l', 'r', 'ml', 'full', 'nj', 'bionj', 'bionj+', 'lk']


class LSFPhyml(ExternalSoftware):
    """ Wraps Phyml objects and runs them on LSF """

    default_binary = 'phyml'
    score_regex = re.compile('(?<=Log-likelihood: ).+')
    local_dir = fileIO.path_to(__file__)

    def __init__(self, records, tmpdir, supplied_binary='', debug=False):
        """ Parent class init method sets temp dir and binary
            Debug=True allows lsf output to be written to tmpdir """
        super(LSFPhyml, self).__init__(tmpdir, supplied_binary, debug)
        self.records = records
        self.temp_dirs = self.setup_temp_dirs()
        self.phyml_objects = self.setup_phyml_objects()
        self.job_ids = set()
        self.debug = debug

    @property
    def records(self):
        return self._records

    @records.setter
    def records(self, records):
        self._records = records

    def write(self):
        pass

    def setup_temp_dirs(self):
        temp_dirs = [tempfile.mkdtemp(dir=self.tmpdir) for rec in self.records]
        return temp_dirs

    def setup_phyml_objects(self):
        phyml_objects = [Phyml(rec, td, self.binary)
                         for (rec, td) in zip(self.records, self.temp_dirs)]
        return phyml_objects

    def get_command_strings(self, analysis='ml'):
        return [phyml.run(analysis, dry_run=True)
                for phyml in self.phyml_objects]

    @staticmethod
    def _memcalc(taxa, categories, sites, states, grace=PHYML_MEMORY_MULTIPLIER,
                 affine=PHYML_MEMORY_SPARE, min=PHYML_MEMORY_MIN):
        """ 18, 4, 65235, 20, grace=1.0) = 2284.8621368408203 ~= 2285 """
        branches = 2 * taxa - 3
        minmem = 4 * (4 * branches * categories * sites * states +
                      2 * branches * categories * states ** 2 +
                      2 * branches * sites * states +
                      6 * branches * sites -
                      2 * categories * sites * states * taxa -
                      sites * taxa) / 1024.0 ** 2
        return int(max(round(minmem * grace, 0) + affine, PHYML_MEMORY_MIN))

    def get_memory_requirements(self):
        l = []
        for record in self.records:
            ntaxa = record.length
            ncategories = 4
            nsites = record.distinct_sites
            nstates = 20 if record.datatype == 'protein' else 4
            memory = self._memcalc(ntaxa, ncategories, nsites, nstates)
            l.append(memory)
        return l

    def launch_lsf(self, command_strings, strategy,
                   minmem=PHYML_MEMORY_STANDARD,
                   verbose=False):
        optioncheck(strategy, ['fixed', 'dynamic'])
        if strategy == 'fixed':
            return self._launch_lsf_fixed_memory(command_strings, minmem,
                                                 verbose)
        elif strategy == 'dynamic':
            return self._launch_lsf_dynamic_memory(command_strings, verbose)

    def _launch_lsf_fixed_memory(self, command_strings, minmem=4096,
                                 verbose=False):
        """ Uses bsub package to send phyml jobs to lsf """
        curr_dir = os.getcwd()
        os.chdir(self.tmpdir)

        job_launcher = bsub('treeCl_static_phyml_task',
                            R='rusage[mem={}]'.format(minmem),
                            M=minmem,
                            verbose=verbose)

        # overwrite kwargs pertaining to output log files
        if not self.debug:
            job_launcher.kwargs['o'] = job_launcher.kwargs['e'] = '/dev/null'

        job_ids = [job_launcher(cmd).job_id
                   for cmd in command_strings]
        self.job_ids.update(job_ids)
        bsub.poll(job_ids)
        os.chdir(curr_dir)

    def _launch_lsf_dynamic_memory(self, command_strings, verbose=False):
        curr_dir = os.getcwd()
        os.chdir(self.tmpdir)

        memory = self.get_memory_requirements()
        job_ids = []
        for i, cmd in enumerate(command_strings):
            memory_reqd = memory[i]
            job_launcher = bsub('treeCl_dynamic_phyml_task',
                                R='rusage[mem={}]'.format(memory_reqd),
                                M=memory_reqd,
                                verbose=verbose)
            if not self.debug:
                job_launcher.kwargs['o'] = '/dev/null'
                job_launcher.kwargs['e'] = '/dev/null'

            job_ids.append(job_launcher(cmd).job_id)
        self.job_ids.update(job_ids)
        bsub.poll(job_ids)
        os.chdir(curr_dir)

    def read(self, analysis, **kwargs):
        self.trees = []
        for phyml in self.phyml_objects:
            (tree, stats) = phyml.read()
            try:
                score = float(self.score_regex.search(stats).group(0))
            except:
                score = 0
            tree_object = Tree(newick=tree, score=score, program=analysis,
                               name=phyml.record.name, output=stats, **kwargs)
            self.trees.append(tree_object)
        return self.trees

    def clean(self):
        """ Delete generated files - TODO: check this works properly """
        for phyml in self.phyml_objects:
            phyml.clean()
        for d in self.temp_dirs:
            if fileIO.can_open(d):
                shutil.rmtree(d)

        deleted = set()
        for job_id in self.job_ids:
            output_file = os.path.join(self.tmpdir,
                                       'treeCl_phyml_task.{}.out'.format(job_id))
            errors_file = os.path.join(self.tmpdir,
                                       'treeCl_phyml_task.{}.err'.format(job_id))
            if (fileIO.delete_if_exists(output_file) and
                    fileIO.delete_if_exists(errors_file)):
                deleted.add(job_id)
        self.job_ids.discard(deleted)

    def run(self, analysis, strategy, minmem, verbose=False, **kwargs):
        """ Run the chosen phyml analysis over LSF """
        command_strings = self.get_command_strings(analysis)
        self.launch_lsf(command_strings, strategy, minmem, verbose)
        trees = self.read(analysis, **kwargs)
        if len(trees) == len(self.records):
            self.clean()
        return trees

    def call(self):
        pass


class Phyml(TreeSoftware):
    """ __init__ takes a Seq sequence record as
    first positional argument, tmpdir as second, and supplied_binary=
    as keyword argument """

    default_binary = 'phyml'
    score_regex = re.compile('(?<=Log-likelihood: ).+')
    local_dir = fileIO.path_to(__file__)

    def read(self, filename=None, tries=5):
        filename = filename or self.filename
        try:
            with open(self.tree_filename) as treefile:
                with open(self.stats_filename) as statsfile:
                    return (treefile.read(), statsfile.read())
        except IOError, e:
            if tries > 0:
                time.sleep(1)
                return self.read(filename, tries - 1)
            print('There was an IOError: {0}'.format(e))
            print('Couldn\'t read PhyML output')
            self.clean()
            raise

    def run(self, analysis=None, verbosity=0,
            **kwargs):
        if analysis:
            self.set_default_flags(analysis)
        else:
            analysis = fileIO.basename(self.binary)
        optioncheck(analysis, ANALYSES)
        if verbosity > 1:
            print(self.flags)
            print('Writing tempfiles to', self.tmpdir)
        filename = self.write()
        filecheck(filename)
        self.add_flag('-i', filename)

        self.tree_filename = filename + '_phyml_tree.txt'
        self.stats_filename = filename + '_phyml_stats.txt'
        self.add_tempfile(self.tree_filename)
        self.add_tempfile(self.stats_filename)

        # OUTPUT TO USER
        if verbosity == 1:
            print_and_return('Running phyml on {0}'.format(self.record.name))
        elif verbosity > 1:
            print('Running phyml on {0}'.format(self.record.name))

        # DRY RUN - just get command string
        if kwargs.get('dry_run', False):
            cmd = self.call(verbose=(True if verbosity > 1 else False),
                            dry_run=True)
            return cmd

        # RUN PHYML
        (stdout, stderr) = self.call(verbose=(True if verbosity > 1 else False))
        (tree, stats) = self.read(filename)
        try:
            score = float(self.score_regex.search(stats).group(0))
        except:
            print(tree)
            print(stats)
        if verbosity > 1:
            print('Cleaning tempfiles')
        self.clean()
        tree_object = Tree(newick=tree,
                           score=score,
                           program=('phyml+' + analysis),
                           name=self.record.name,
                           output=stats, **kwargs)
        if kwargs.get('set_as_record_tree', True):
            self.record.tree = tree_object
        if verbosity > 1:
            print('Done.')
        return tree_object

    def write(self):
        record_name = self.record.get_name(default='tmp_phyml_input')
        with tempfile.NamedTemporaryFile(prefix=record_name,
                                         suffix='.phy',
                                         dir=self.tmpdir,
                                         delete=False) as file_:
            filename = file_.name
            file_.write(self.record.write_phylip('pipe'))
        self.filename = filename
        # filename = '{0}/{1}.phy'.format(self.tmpdir, record_name)
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
            defaults['--no_memory_check'] = ''

            if analysis == 'ml' or analysis == 'full' or analysis == 'tlr':
                defaults['-o'] = 'tlr'
            elif analysis == 'nj' or analysis == 'bionj':
                defaults['-o'] = 'n'
            elif analysis == 'lr' or analysis == 'bionj+':
                defaults['-o'] = 'lr'
            elif analysis == 'l':
                defaults['-o'] = 'l'
            elif analysis == 'r':
                defaults['-o'] = 'r'
            elif analysis == 'lk':
                defaults['-o'] = 'n'

            for flag in defaults:
                self.add_flag(flag, defaults[flag])


def runPhyml(rec, tmpdir, analysis, verbosity=0, tree=None, **kwargs):
    optioncheck(analysis, ANALYSES)
    phyml = Phyml(rec, tmpdir)
    if (analysis == 'lk' or analysis == 'r') and tree is not None:
        tree_name = (tree.name if tree.name else 'tmp_tree')
        tmp_treefile = '{0}/{1}.nwk'.format(tmpdir, tree_name)
        tree.write_to_file(tmp_treefile)
        phyml.add_tempfile(filecheck(tmp_treefile))
        phyml.add_flag('-u', tmp_treefile)
    try:
        return phyml.run(analysis, verbosity, **kwargs)
    finally:
        if not phyml.debug:
            phyml.clean()
            if verbosity > 2:
                print('Cleaned up')


def runLSFPhyml(records, tmpdir, analysis, verbosity, strategy,
                minmem=PHYML_MEMORY_STANDARD, debug=False, **kwargs):
    optioncheck(analysis, ANALYSES)
    optioncheck(strategy, ['fixed', 'dynamic'])
    lsfphyml = LSFPhyml(records, tmpdir, debug=debug)
    verbose = output_files = (True if verbosity > 0 else False)
    try:
        return lsfphyml.run(analysis,
                            strategy=strategy,
                            minmem=minmem,
                            verbose=verbose,
                            output_files=output_files,
                            **kwargs)
    finally:
        if not lsfphyml.debug:
            lsfphyml.clean()
            if verbosity > 2:
                print('Cleaned up')


def phyml_memcalc(taxa, categories, sites, states):
    branches = 2 * taxa - 3
    return 4 * (4 * branches * categories * sites * states +
                2 * branches * categories * states ** 2 +
                2 * branches * sites * states +
                6 * branches * sites -
                2 * categories * sites * states * taxa -
                sites * taxa) / 1024.0 ** 2
