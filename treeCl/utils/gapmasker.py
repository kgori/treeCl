from __future__ import print_function
from builtins import str
from builtins import object
from collections import defaultdict


class GapMasker(object):
    def __init__(self, template):
        self.template = template
        self.gap_positions = self.get_gap_positions()

    def get_gap_positions(self):
        gap_positions = defaultdict(list)
        names = self.template.headers
        for name in names:
            seq = self.template.mapping[name]
            for pos, char in enumerate(seq):
                if char == '-':
                    gap_positions[name].append(pos)
        return gap_positions

    def mask(self, target):
        try:
            self.check_seqs(target)
            return self.write_gap_positions(target)
        except Exception as e:
            print(e)
            return

    def check_seqs(self, target):
        if len(self.template) != len(target):
            raise Exception('Alignments have different numbers of sequences')

        if set(self.template.headers) != set(target.headers):
            raise Exception('Sequence names don\'t match')

        if self.template.seqlength != target.seqlength:
            raise Exception('Alignments are different lengths')


    def write_gap_positions(self, target):
        for name in self.gap_positions:
            if not name in target.headers:
                raise Exception('Trying to write gaps to non-existent sequence: ',
                                name)
            listseq = list(target.mapping[name])
            for pos in self.gap_positions[name]:
                listseq[pos] = '-'
            seq = ''.join(str(x) for x in listseq)
            target.mapping[name] = seq
        seqs = []
        for name in target.headers:
            seqs.append(target.mapping[name])
        target.sequences = seqs
        target.update()
        return target
