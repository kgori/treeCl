#!/usr/bin/env python

from . import ExternalSoftware
from ..errors import FileError, filecheck_and_raise
from ..utils.fileIO import *
from ..utils.printing import print_and_return

local_dir = path_to(__file__)

class DVWrapper(ExternalSoftware):

    default_binary = 'darwin'
    DV_wrapper = locate_file('DV_wrapper.drw', directory=local_dir)

    def __init__(self, record=None, supplied_binary=''):
        super(DVWrapper, self).__init__(supplied_binary)
        self.record = record
        self.tmpdir = record.tmpdir or '/tmp'

    def __str__(self):
        desc = 'Wrapper for darwin method DV_wrapper.drw'

    def call(
        self,
        filename,
        datatype,
        dryrun=False,
        ):
        datatype = datatype or 'AA'
        cmd = \
            'echo "fil := ReadFastaWithNames(\'{0}\'); seqtype := \'{1}\'; fpath := \'{2}/\'; ReadProgram(\'{3}\');" | {4}'.format(filename,
                datatype, self.tmpdir, self.DV_wrapper, self.binary)
        if dryrun:
            return cmd
        return subprocess(cmd)

    def read_datatype(self):
        if self.record.datatype == 'protein':
            return 'AA'
        elif self.record.datatype == 'dna':
            return 'DNA'

    def read(self):
        output = '{0}/temp_distvar.txt'.format(self.tmpdir)
        filecheck_and_raise(output)
        self.add_tempfile(output)
        return open(output).read().rstrip()

    def run(self, verbosity=0):
        """ Verbosity levels: 0 = No output; 1 = Reports that it is running; 2 =
        More details; 3 or more = Display stdout + stderr

        Level 2 includes level 1, level 3 includes levels 2 and 1"""

        if verbosity == 1:
            print_and_return( 'Calculating distances for {0}'.format(self.record.name) )

        filename = self.write()
        if verbosity > 1:
            print 'Calculating distances for {0}'.format(self.record.name)
            print 'Tempfile {0} written to {1}'.format(filename, self.tmpdir)
        datatype = self.read_datatype()
        if verbosity > 1:
            print 'Datatype inferred to be', datatype
        (stdout, stderr) = self.call(filename, datatype)
        if verbosity > 2:
            print stdout, stderr
        dv_string = self.read()
        labels = ' '.join(self.record.headers)
        if verbosity > 1:
            print 'Cleaning tempfile'
        self.clean()
        self.record.dv = [(dv_string, labels)]
        if verbosity > 1:
            print 'Done.'
        return (dv_string, labels)

    def write(self):
        if self.record.name:
            filename = (self.record.name if len(self.record.name)
                        < 40 else self.record.hashname())
        else:
            filename = 'fasta_tmp'
        filename = '{0}/{1}.tmp'.format(self.tmpdir, filename)
        self.record.write_fasta(filename)
        self.add_tempfile(filename)
        return filename


# RUNNER

def runDV(rec, verbosity=0):
    dv = DVWrapper(rec)
    return dv.run(verbosity)
