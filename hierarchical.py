#!/usr/bin/env python

from optimiser import Partition, Optimiser, Scorer


class ml_hierarchy(Optimiser):
    def __init__(self, collection, tmpdir='/tmp/', init=None):
        self.collection = collection

        if not self.collection.records[0].tree:
            print 'Calculating NJ trees for collection...'
            self.collection.calc_NJ_trees()

        self.datatype = collection.datatype
        self.Scorer = Scorer(self.collection.records, analysis='nj',
                             datatype=self.datatype,
                             tmpdir=tmpdir)

        self.levels = {}

        if init is None:
            self.initial_assignment = Partition(tuple(range(len(self.collection.records))))
        else:
            self.initial_assignment = init

    def cluster(self, assignment):
        self.levels[max(assignment.partition_vector)] = assignment
        new_assignment = self.merge_closest(assignment)
        if max(new_assignment.partition_vector) > 1:
            self.cluster(new_assignment)
