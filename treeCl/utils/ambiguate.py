#!/usr/bin/env python
from __future__ import print_function
from builtins import zip
from ..alignment import Alignment

docstring = '''
Combine two unphased sequences into a single sequence with ambiguity codes
Usage: ambiguate.py <FILENAME=eg:phylip.phy>
Output: <phylip_ambig.phy>

'''

ambiguities = {
    'a': frozenset(['a']),
    'c': frozenset(['c']),
    'g': frozenset(['g']),
    't': frozenset(['t']),
    'y': frozenset(['c', 't']),
    'r': frozenset(['a', 'g']),
    'w': frozenset(['a', 't']),
    's': frozenset(['g', 'c']),
    'k': frozenset(['t', 'g']),
    'm': frozenset(['c', 'a']),
    'd': frozenset(['a', 'g', 't']),
    'v': frozenset(['a', 'c', 'g']),
    'h': frozenset(['a', 'c', 't']),
    'b': frozenset(['c', 'g', 't']),
    'x': frozenset(['a', 'c', 'g', 't']),
    'n': frozenset(['a', 'c', 'g', 't']),
}

ambiguities_rev = {v: k for (k, v) in ambiguities.items()}
ambiguities_rev[frozenset(['a', 'c', 'g', 't'])] = 'n'


def get_prefixes(r):
    prefix_list = list()
    for x in r.headers:
        prefix = '.'.join(x.split('.')[:-1])
        if not prefix in prefix_list:
            prefix_list.append(prefix)
    return prefix_list


def get_ambiguity(a, b):
    upper = False
    if a.isupper():
        upper = True
        a = a.lower()
    if b.isupper():
        upper = True
        b = b.lower()
    s1 = ambiguities[a]
    s2 = ambiguities[b]
    union = s1 | s2
    ambig = ambiguities_rev[union]

    return ambig.upper() if upper else ambig


def ambiguate(sequence1, sequence2, delete_ambiguous=False):
    """ delete_ambiguous: Marks sequences for deletion by replacing all
    chars with 'X'. These seqs are deleted later with remove_empty """
    delete = False
    combination = list()
    z = list(zip(sequence1, sequence2))
    for (a, b) in z:
        if a == b:
            combination.append(a)
        else:
            if a == '-' or b == '-':
                combination.append('-')
            else:
                if delete_ambiguous:
                    delete = True
                ambig = get_ambiguity(a, b)
                combination.append(ambig)
    if delete:
        return 'X' * len(combination)
    return ''.join(combination)


def remove_empty(rec):
    """ Deletes sequences that were marked for deletion by convert_to_IUPAC """
    for header, sequence in rec.mapping.items():
        if all(char == 'X' for char in sequence):
            rec.headers.remove(header)
            rec.sequences.remove(sequence)
    rec.update()
    return rec


def get_seqs(rec, pref):
    return rec.mapping[pref + '.1'], rec.mapping[pref + '.2']


if __name__ == '__main__':
    import argparse

    parser = argparse.ArgumentParser()
    parser.add_argument('infile', help='Input file in phylip format')
    parser.add_argument('-d', '--delete_ambiguous', action='store_true',
                        help=('Deletes sequences with any'
                              'ambiguity - could result in empty alignment'))
    parser.add_argument('-c', '--cutoff', type=int, default=2,
                        help=('If number of sequences is less than this, '
                              'no output is written'))
    parser.add_argument('-o', '--outfile', type=str)
    args = parser.parse_args()

    f = args.infile
    rec = Alignment(f, 'phylip')
    prefixes = get_prefixes(rec)
    headers = list()
    sequences = list()
    for pref in prefixes:
        seq1, seq2 = get_seqs(rec, pref)
        combin = ambiguate(seq1, seq2, args.delete_ambiguous)
        headers.append(pref)
        sequences.append(combin)
    newrec = Alignment(headers=headers, sequences=sequences)
    remove_empty(newrec)
    out = (args.outfile if args.outfile
           else f[:f.rindex('.phy')] + '_ambig.phy')
    if len(newrec) >= args.cutoff:
        newrec.write_phylip(out, interleaved=True)
    else:
        print('Empty')
