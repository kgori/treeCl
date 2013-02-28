#!/usr/bin/env python

import re  # for sorting
from tree import Tree
from externals.dv_wrapper import DVWrapper
from externals.tree_builders import Phyml, TreeCollection
import hashlib
from random import shuffle as shf


class SequenceRecord(object):

    """ Class for reading sequence files in fasta or phylip formats Supports
    writing in fasta, phylip, phylip interleaved, nexus formats, sorting
    sequences by length and name, concatenating sequences when sequence names
    are a perfect match, iterating over records. """

    def __init__(
        self,
        infile=None,
        file_format='fasta',
        name=None,
        datatype=None,
        headers=[],
        sequences=[],
        ):

        self.name = name
        self.headers = headers
        self.sequences = sequences
        self.mapping = {}
        self.length = 0
        self.seqlength = 0
        self.datatype = datatype
        self.is_aligned = False
        if infile:
            if file_format == 'fasta':
                self.read_fasta_file(infile, name=name, datatype=datatype)
            elif file_format == 'phylip':
                self.read_phylip_file(infile, name=name, datatype=datatype)
        self.index = -1
        self._update()

    def _update(self):
        """ For updating the length and mapping attributes of the object after
        reading sequences """

        if self.headers and self.sequences:
            self.mapping = dict(zip(self.headers, self.sequences))
            self.length = len(self.mapping)
            first_seq_length = len(self.sequences[0])
            is_aligned = True
            for seq in self.sequences:
                if len(seq) != first_seq_length:
                    is_aligned = False
                    break
                else:
                    continue
            self.is_aligned = is_aligned
            self.seqlength = first_seq_length

    def __iter__(self):  # Should do this with generators / yield
        return self

    def next(self):  # As above
        self.index += 1
        if self.index == self.length:
            self.index = -1
            raise StopIteration
        return {'header': self.headers[self.index],
                'sequence': self.sequences[self.index]}

    def __str__(self):
        output_string = ''
        if self.is_aligned:
            output_string += 'Aligned_Sequence_Record: {0}\n'.format(self.name)
        else:
            output_string += \
                'Unaligned_Sequence_Record: {0}\n'.format(self.name)
        for i in range(len(self.headers)):
            if len(self.sequences[i]) > 50:
                output_string += '>{0}\n'.format(self.headers[i]) \
                    + '{0}\n'.format((self.sequences[i])[:50]) \
                    + '... ({0}) ...\n'.format(len(self.sequences[i]) - 100) \
                    + '{0}\n'.format((self.sequences[i])[-50:]) + '\n'
            else:
                output_string += '>{0}\n{1}'.format(self.headers[i],
                        self.sequences[i]) + '\n' + '\n'
        output_string += '{0} sequences in record'.format(len(self))
        return output_string

    def __repr__(self):
        return '{0}: Name={1}'.format(self.__class__.__name__, self.name)

    def __len__(self):
        return self.length

    def __add__(self, other):
        """ Allows SequenceRecords to be added together, concatenating sequences
        """

        try:
            assert set(self.headers) == set(other.headers)
        except AssertionError:
            print 'Sequence labels do not match between alignments'
            return self
        d = {}
        for k in self.mapping.keys():
            d[k] = self.mapping[k] + other.mapping[k]
        return SequenceRecord(headers=d.keys(), sequences=d.values())

    def __radd__(self, other):
        return other.__add__(self)

    def __mul__(self, n):
        """ Allows SequenceRecord * x to return a concatenation of the record x
        times """

        if not isinstance(n, int):
            print 'Truncating {0} to {1}'.format(n, int(n))
            n = int(n)
        new_seqs = []
        for s in self.sequences:
            new_seqs.append(s * n)
        return SequenceRecord(headers=self.headers, sequences=new_seqs)

    def __rmul__(self, n):
        if not isinstance(n, int):
            print 'Truncating {0} to {1}'.format(n, int(n))
            n = int(n)
        new_seqs = []
        for s in self.sequences:
            new_seqs.append(s * n)
        return SequenceRecord(headers=self.headers, sequences=new_seqs)

    def __eq__(self, other):
        if type(other) is type(self):
            return set(self.headers) == set(other.headers) \
                and set(self.sequences) == set(other.sequences)
        return False

    def __ne__(self, other):
        return not self.__eq__(other)

    def sort_by_length(self, in_place=True):
        """ Sorts sequences by ungapped length If in_place = False the sorting
        doesn't mutate the underlying object, and the output is returned If
        in_place = True the sorting mutates the self object """

        # Sort sequences by descending order of length
        # Uses zip as its own inverse [ zip(*zip(A,B)) == (A,B) ]

        (h, s) = zip(*sorted(zip(self.headers, self.sequences),
                     key=lambda item: len(item[1].replace('-', '')),
                     reverse=True))
        if in_place:
            self.headers = h
            self.sequences = s
        return SequenceRecord(name=self.name, headers=h, sequences=s)

    def sort_by_name(self, in_place=True):
        """ Sorts sequences by name, treating numbers as integers (i.e. sorting
        like this: 1, 2, 3, 10, 20 not 1, 10, 2, 20, 3). If in_place = False the
        sorting doesn't mutate the underlying object, and the output is returned
        If in_place = True the sorting mutates the self object """

        items = self.mapping.items()
        sort_key = lambda item: tuple((int(num) if num else alpha) for (num,
                                      alpha) in re.findall(r'(\d+)|(\D+)',
                                      item[0]))
        items = sorted(items, key=sort_key)
        (h, s) = zip(*items)
        if in_place:
            self.headers = h
            self.sequences = s
        else:
            return SequenceRecord(name=self.name, headers=h, sequences=s)

    def hashname(self):
        H = hashlib.sha1()
        H.update(self.name)
        return H.hexdigest()

    def get_name(self, maxlen=40, default='noname'):
        if self.name:
            return (self.name if len(self.name) < maxlen else self.hashname())
        else:
            return default

    @staticmethod
    def linebreaker(string, length):
        for i in range(0, len(string), length):
            yield string[i:i + length]

    def make_chunks(self, chunksize):
        num_chunks = self.seqlength / chunksize
        if self.seqlength % chunksize:
            num_chunks += 1

        new_records = []
        generators = [self.linebreaker(s, chunksize) for s in self.sequences]
        for _ in range(num_chunks):
            new_record = type(self)(headers=self.headers)
            new_record.sequences = [next(g) for g in generators]
            new_record._update()
            new_records.append(new_record)
        return new_records

    def read_fasta_file(
        self,
        fasta_file,
        name=None,
        datatype=None,
        ):
        """ FASTA format parser: turns fasta file into Alignment_record object
        """

        headers = []
        sequences = []
        openfile = open(fasta_file, 'r')

        while True:                     # skip over file until first header is 
            line = openfile.readline()  # found,
            if not line:                
                return                  
            if line[0] == '>':          # then break the loop and put the first
                break                   # header into the headers list                                        
        headers.append(line[1:].rstrip())  

        while True:
            line = openfile.readline()
            sequence_so_far = []        # build up sequence a line at a time
            while True:
                if not line:
                    break
                elif not line[0] == '>':
                    sequence_so_far.append(line.rstrip())
                    line = openfile.readline()
                else:
                    break
            sequences.append(''.join(sequence_so_far).replace(',', ''))
            if not line:
                break
            headers.append(line[1:].rstrip())

        first_seq_length = len(sequences[0])    # check all sequences are the 
        is_alignment = True                     # same length
        for seq in sequences:
            if len(seq) != first_seq_length:
                is_alignment = False
                break
            else:
                continue

        if len(headers) != len(sequences):  
            print 'Error matching all headers and sequences'

        if is_alignment:
            self.is_aligned = True
        self.name = name
        self.headers = headers
        self.sequences = sequences
        self.datatype = datatype
        self._update()

    def read_phylip_file(
        self,
        phylip_file,
        name=None,
        datatype=None,
        ):
        """ PHYLIP format parser"""

        openfile = open(phylip_file, 'r')
        info = openfile.readline().split()
        num_taxa = int(info[0])
        seq_length = int(info[1])

        # Initialise lists to hold the headers and sequences
        headers = [None] * num_taxa
        sequences = [None] * num_taxa

        i = 0  # counter monitors how many (informative) lines we've seen
        for line in openfile:
            line = line.rstrip()
            if not line:
                continue  # skip empty lines and don't increment counter

            line = line.split()

            # IF header is None, split line into header / sequence pair
            # ELSE header is not None, we've already seen it
            # so this file must be interleaved, so carry on
            # adding sequence fragments

            if not headers[i % num_taxa]:
                header = line[0]
                sequence_fragment = ''.join(line[1:])

                headers[i % num_taxa] = header
                sequences[i % num_taxa] = [sequence_fragment]
            else:

                sequence_fragment = ''.join(line)
                sequences[i % num_taxa].append(sequence_fragment)

            i += 1  # increment counter

        sequences = [''.join(x) for x in sequences]

        # checks

        try:
            assert len(headers) == len(sequences) == num_taxa
            for sequence in sequences:
                assert len(sequence) == seq_length
        except AssertionError:
            print 'Error reading file'
            return

        self.name = name
        self.datatype = datatype
        (self.headers, self.sequences) = (headers, sequences)
        self._update()

    def write_fasta(
        self,
        outfile='stdout',
        print_to_screen=False,
        linebreaks=None,
        ):
        """ Writes sequences to file in fasta format If outfile = 'stdout' the
        sequences are printed to screen, not written to file If print_to_screen
        = True the sequences are printed to screen whether they are written to
        disk or not """

        if linebreaks:
            try:
                int(linebreaks)
            except ValueError:
                print 'Can\'t use {0} as value for linebreaks'.format(linebreaks)
            sequences = ['\n'.join(self.linebreaker(s, linebreaks)) for s in
                         self.sequences]
        else:
            sequences = self.sequences

        lines = ['>{0}\n{1}'.format(h, seq) for (h, seq) in zip(self.headers,
                 sequences)]
        s = '\n'.join(lines)
        s += '\n'
        if outfile == 'stdout':
            print s
            return s
        elif outfile == 'pipe':
            if print_to_screen:
                print s
            return s
        else:
            open(outfile, 'w').write(s)
            if print_to_screen:
                print s
            return outfile

    def write_nexus(self, outfile='stdout', sequence_type='protein'):
        maxlen = len(max(self.sequences, key=len))
        lines = ['{0:<14} {1:-<{2}}'.format(x, y, maxlen) for (x, y) in
                 zip(self.headers, self.sequences)]
        file_header = '#NEXUS\n' + '\n'
        file_header += 'begin data;\n'
        file_header += \
            '    dimensions ntax={0} nchar={1};\n'.format(self.length, maxlen)
        file_header += \
            '    format datatype={0} interleave=no gap=-;\n'.format(sequence_type)
        file_header += '    matrix\n' + '\n' + '\n'

        file_footer = '    ;\n' + '\n' + 'end;\n'

        s = file_header + '\n'.join(lines) + file_footer
        if outfile == 'stdout':
            print s
            return s
        elif outfile == 'pipe':
            return s
        else:
            open(outfile, 'w').write(s)
            return outfile

    def write_phylip(
        self,
        outfile='stdout',
        print_to_screen=False,
        interleaved=False,
        linebreaks=120,
        ):
        """ Writes sequences to file in phylip format, interleaving optional If
        outfile = 'stdout' the sequences are printed to screen, not written to
        disk If print_to_screen = True the sequences are printed to screen
        whether they are written to disk or not """

        maxlen = len(max(self.sequences, key=len))
        file_header = ' {0} {1}'.format(self.length, maxlen)
        s = [file_header]
        maxheader = len(max(self.headers, key=len))
        label_length = max(maxheader + 1, 10)
        if interleaved:
            seq_length = linebreaks - label_length
            num_lines = maxlen / seq_length
            if maxlen % seq_length:
                num_lines += 1

            for i in range(num_lines):
                for seq_header in self.headers:
                    if i == 0:
                        s.append('{0:<{1}} {2}'.format(seq_header,
                                 label_length, (self.mapping[seq_header])[i
                                 * seq_length:(i + 1) * seq_length]))
                    else:
                        s.append('{0} {1}'.format(' ' * label_length,
                                 (self.mapping[seq_header])[i * seq_length:(i
                                 + 1) * seq_length]))
                s.append('')
        else:
            lines = ['{0:<{1}} {2:-<{3}}'.format(x, label_length, y, maxlen)
                     for (x, y) in zip(self.headers, self.sequences)]
            s.extend(lines)
            s.append('')
        s = '\n'.join(s)
        if outfile == 'stdout':
            print s
            return s
        elif outfile == 'pipe':
            if print_to_screen:
                print s
            return s
        else:
            open(outfile, 'w').write(s)
            if print_to_screen:
                print s
            return outfile


class TCSeqRec(SequenceRecord):

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
        dv=[],
        tree=None,
        tmpdir='/tmp',
        ):

        self.TCfiles = {}
        self.dv = dv
        if tree and isinstance(tree, Tree):
            self.tree = tree
        else:
            self.tree = Tree()
        self.tmpdir = tmpdir

        super(TCSeqRec, self).__init__(
            infile,
            file_format,
            name,
            datatype,
            headers,
            sequences,
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
                d[k] = self.mapping[k] + 'N' * other.seqlength
            elif k in only_in_other:
                d[k] = 'N' * self.seqlength + other.mapping[k]
        dvsum = self.dv + other.dv
        return_object = TCSeqRec(headers=d.keys(), sequences=d.values(),
                                 datatype=self.datatype).sort_by_name(in_place=False)
        return_object.dv = dvsum
        return return_object

    def bionj(self):
        """ Uses phyml (via treeCl.externals.tree_builders.Phyml) to build a
        bioNJ tree for the current record """

        p = Phyml(self)
        self.tree = p.run('nj')

    def dv_matrix(self):
        """ Uses darwin (via treeCl.externals.DVWrapper) to calculate pairwise
        distances and variances"""

        dv = DVWrapper(self)
        self.dv.append(dv.run())

    def phyml(self):
        """ Uses phyml (via treeCl.externals.tree_builders.Phyml) to build a
        full ML tree for the current record """

        p = Phyml(self)
        self.tree = p.run('ml')

    def tree_collection(self):
        """ Uses TreeCollection (via
        treeCl.externals.tree_builders.TreeCollection) to build a least squares
        tree for the current record """

        if self.dv <= []:
            self.get_dv_matrix()
        tc = TreeCollection(self)
        self.tree = tc.run()

    def _pivot(self, lst):
        new_lst = zip(*lst)
        return [''.join(x) for x in new_lst]

    def sanitise(self):
        self.sort_by_name()
        l = []
        for h in self.headers:
            if '/' in h:
                h = h[:h.index('/')]
            while h.startswith(' '):
                h = h[1:]
            h = h.replace(' ', '_')
            l.append(h)
        self.headers = l
        self.sequences = [seq.upper() for seq in self.sequences]
        self._update()

    def shuffle(self):
        """ Modifies in-place """

        columns = self._pivot(self.sequences)
        shf(columns)
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
            return TCSeqRec(name=self.name, headers=h, sequences=s,
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
            return TCSeqRec(name=self.name, headers=h, sequences=s,
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
            newrec = TCSeqRec(headers=self.headers, sequences=newseqs,
                              datatype=self.datatype)
            newrecs.append(newrec)
        if names:
            for (i, newrec) in enumerate(newrecs):
                newrec.name = names[i]
        else:
            for (i, newrec) in enumerate(newrecs):
                newrec.name = 'record_{0}'.format(i + 1)
        return newrecs
