#!/usr/bin/env python
from __future__ import print_function

# standard library
from copy import deepcopy
import hashlib
import itertools
import re

# treeCl
from ..constants import DNA_ACGT_THRESHOLD
from ..errors import optioncheck, directorycheck
from ..utils import fileIO


class Seq(object):
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
            tmpdir=None,
    ):

        self.name = name
        self.headers = headers
        self.sequences = sequences
        self.mapping = {}
        self.length = 0
        self.seqlength = 0
        self.is_aligned = False
        if infile:
            if file_format == 'fasta':
                self.read_fasta_file(infile, name=name)
            elif file_format == 'phylip':
                self.read_phylip_file(infile, name=name)
        self.datatype = (datatype or self.guess_datatype())
        self.index = -1
        self._update()
        self.tmpdir = tmpdir
        if tmpdir is not None:
            directorycheck(tmpdir)

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

    def guess_datatype(self):
        if not self.sequences:
            return
        aa_only = {'E', 'F', 'I', 'L', 'P', 'Q', 'X', 'Z',
                   'e', 'f', 'i', 'l', 'p', 'q', 'x', 'z'}
        all_chars = [x for x in itertools.chain(*self.sequences)
                     if not x in {'-', 'x', 'X'}]
        if aa_only.intersection(all_chars):
            return 'protein'
        acgt_frequency = (sum(all_chars.count(x) for x in 'acgtACGT') /
                          float(len(all_chars)))

        return 'dna' if acgt_frequency > DNA_ACGT_THRESHOLD else 'protein'

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

        return self.__class__(headers=d.keys(), sequences=d.values(),
                              datatype=self.datatype).sort_by_name(in_place=False)

    def __radd__(self, other):
        return other.__add__(self)

    def __mul__(self, n):
        """ Allows SequenceRecord * x to return a concatenation of the record x
        times """

        if not isinstance(n, int):
            print('Truncating {0} to {1}'.format(n, int(n)))
            n = int(n)
        new_seqs = []
        for s in self.sequences:
            new_seqs.append(s * n)
        return self.__class__(headers=self.headers, sequences=new_seqs)

    def __rmul__(self, n):
        if not isinstance(n, int):
            print('Truncating {0} to {1}'.format(n, int(n)))
            n = int(n)
        new_seqs = []
        for s in self.sequences:
            new_seqs.append(s * n)
        return self.__class__(headers=self.headers, sequences=new_seqs)

    def __eq__(self, other):
        if type(other) is type(self):
            return set(self.headers) == set(other.headers) \
                   and set(self.sequences) == set(other.sequences)
        return False

    def __ne__(self, other):
        return not self.__eq__(other)

    def sort_by_length(self, in_place=True, reverse=True):
        """ Sorts sequences by ungapped length If in_place = False the sorting
        doesn't mutate the underlying object, and the output is returned If
        in_place = True the sorting mutates the self object """

        # Sort sequences by descending order of length
        # Uses zip as its own inverse [ zip(*zip(A,B)) == (A,B) ]

        (h, s) = zip(*sorted(zip(self.headers, self.sequences),
                             key=lambda item: len(item[1].replace('-', '')),
                             reverse=reverse))
        if in_place:
            self.headers = h
            self.sequences = s
        return self.__class__(name=self.name, headers=h, sequences=s)

    def sort_by_name(self, in_place=True, reverse=False):
        """ Sorts sequences by name, treating numbers as integers (i.e. sorting
        like this: 1, 2, 3, 10, 20 not 1, 10, 2, 20, 3). If in_place = False the
        sorting doesn't mutate the underlying object, and the output is returned
        If in_place = True the sorting mutates the self object """

        items = self.mapping.items()
        sort_key = lambda item: tuple((int(num) if num else alpha) for (num,
                                                                        alpha) in re.findall(r'(\d+)|(\D+)',
                                                                                             item[0]))
        items = sorted(items, key=sort_key, reverse=reverse)
        (h, s) = zip(*items)
        if in_place:
            self.headers = h
            self.sequences = s
        else:
            return self.__class__(name=self.name, headers=h, sequences=s)

    def hashname(self):
        H = hashlib.sha1()
        H.update(self.name)
        return H.hexdigest()

    def get_name(self, maxlen=40, default='noname'):
        if self.name:
            return (self.name if len(self.name) < maxlen else self.hashname())
        else:
            return default


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
        with fileIO.freader(fasta_file) as openfile:

            while True:  # skip over file until first header
                line = openfile.readline()  # is found,
                if not line:
                    return
                if line[0] == '>':  # then break the loop and put the
                    break  # first header into the headers list
            headers.append(line[1:].rstrip())

            while True:
                line = openfile.readline()
                sequence_so_far = []  # build up sequence a line at a time
                while True:
                    if not line:
                        break
                    elif not line[0] == '>':
                        sequence_so_far.append(line.rstrip())
                        line = openfile.readline()
                    else:
                        break
                sequences.append(''.join(sequence_so_far).translate(None,
                                                                    '\n\r\t ,'))
                if not line:
                    break
                headers.append(line[1:].rstrip())

            first_seq_length = len(sequences[0])  # check if all sequences are
            is_alignment = True  # the same length
            for seq in sequences:
                if len(seq) != first_seq_length:
                    is_alignment = False
                    break
                else:
                    continue

            if len(headers) != len(sequences):
                print('Error matching all headers and sequences')

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

        with fileIO.freader(phylip_file) as openfile:
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
                if headers[i % num_taxa] is None:
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
                print('Error parsing file')
                return

            self.name = name
            self.datatype = datatype
            (self.headers, self.sequences) = (headers, sequences)
            self._update()


    def n_site_patterns(self):
        return len(set(self._pivot(self.sequences)))


    def change_case(self, case):
        optioncheck(case, ['lower', 'upper'])
        if case == 'upper':
            self.sequences = [x.upper() for x in self.sequences]
        else:
            self.sequences = [x.lower() for x in self.sequences]
        self._update()


    @staticmethod
    def linebreaker(string, length):
        for i in range(0, len(string), length):
            yield string[i:i + length]


    @staticmethod
    def _grouper(n, iterable, fillvalue=''):
        "grouper(3, 'ABCDEFG', 'x') --> ABC DEF Gxx"
        args = [iter(iterable)] * n
        return itertools.izip_longest(fillvalue=fillvalue, *args)


    @staticmethod
    def _join(it, sep=''):
        for el in it:
            yield sep.join(el)

    def seqs_to_units(self, units_per_line, unit_size, unit_sep=' '):
        """ Uses nested grouper pattern list comprehension to break sequences
        into `unit_size` length blocks separated by `unit_sep`, with
        `units_per_line` units per line"""
        #
        return ['\n'.join(
            self._join(self._grouper(units_per_line, self._join(self._grouper(unit_size, s))), unit_sep)).strip(
            unit_sep) for s in self.sequences]

    def write_fasta(
            self,
            outfile='stdout',
            print_to_screen=False,
            linebreaks=None,
            units_per_line=10,
            unit_size=10,
            unit_sep='',
    ):
        """ Writes sequences to file in fasta format If outfile = 'stdout' the
        sequences are printed to screen, not written to file If print_to_screen
        = True the sequences are printed to screen whether they are written to
        disk or not """

        if linebreaks is not None:
            unit_size = linebreaks

        if unit_size is None:
            unit_size = max(len(s) for s in self.sequences)
            units_per_line = 1

        sequences = self.seqs_to_units(units_per_line, unit_size, unit_sep)

        lines = ['>{0}\n{1}\n'.format(h, seq) for (h, seq) in zip(self.headers,
                                                                  sequences)]

        if outfile == 'stdout':
            s = ''.join(lines)
            print(s)
            return s
        elif outfile == 'pipe':
            s = ''.join(lines)
            if print_to_screen:
                print(s)
            return s
        else:
            file_writer = fileIO.fwriter(outfile)
            with file_writer as fwr:
                fwr.writelines(lines)
            if print_to_screen:
                print(''.join(lines))
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
            print(s)
            return s
        elif outfile == 'pipe':
            return s
        else:
            with fileIO.fwriter(outfile) as fwr:
                fwr.write(s)
            return outfile

    def write_phylip(
            self,
            outfile='stdout',
            print_to_screen=False,
            interleaved=False,
            linebreaks=None,
            units_per_line=10,
            unit_size=10,
            unit_sep='',
    ):
        """ Writes sequences to file in phylip format, interleaving optional If
        outfile = 'stdout' the sequences are printed to screen, not written to
        disk If print_to_screen = True the sequences are printed to screen
        whether they are written to disk or not """

        linebreaks = linebreaks or 120
        maxlen = max(len(seq) for seq in self.sequences)
        file_header = ' {0} {1}'.format(self.length, maxlen)
        s = [file_header]
        maxheader = max(len(h) for h in self.headers)
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
                                                                                                * seq_length:(
                                                                                                             i + 1) * seq_length]))
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
            print(s)
            return s
        elif outfile == 'pipe':
            if print_to_screen:
                print(s)
            return s
        else:
            with fileIO.fwriter(outfile) as fwr:
                fwr.write(s)
            if print_to_screen:
                print(s)
            return outfile
