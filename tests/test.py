#!/usr/bin/env python
import unittest
import treeCl
import os
from treeCl import Partition

thisdir = os.path.dirname(os.path.realpath(__file__))


class PartitionTests(unittest.TestCase):
    def setUp(self):
        self.partition = Partition(['a', 'd', 'a', 'a', 'b', 'a', 'b', 'c', 'c', 'd', 'd', 'd', 'd', 'd'])

    def test_restricted_growth_notation(self):
        """
        Test restricted growth notation
        :return:
        """
        expected = (0, 1, 0, 0, 2, 0, 2, 3, 3, 1, 1, 1, 1, 1)
        self.assertEqual(self.partition.partition_vector, expected)

    def test_membership_to_rgn(self):
        membership = [(0, 2, 3, 5), (1, 9, 10, 11, 12, 13), (4, 6), (7, 8)]
        inferred = Partition.from_membership(membership)
        self.assertEqual(inferred.partition_vector, self.partition.partition_vector)

    def test_numbers(self):
        self.assertEqual(self.partition.num_elements(), 14)
        self.assertEqual(self.partition.num_groups(), 4)

    def test_get_membership(self):
        expected = [(0, 2, 3, 5), (1, 9, 10, 11, 12, 13), (4, 6), (7, 8)]
        self.assertEqual(self.partition.get_membership(), expected)

    def test_maximal(self):
        p1 = Partition([1, 1, 1, 1, 1, 1, 1, 1, 1])
        p2 = Partition([0, 1, 2, 3, 4, 5, 6, 7, 8])
        self.assertFalse(p1.is_maximal())
        self.assertTrue(p2.is_maximal())
        self.assertFalse(self.partition.is_maximal())

    def test_minimal(self):
        p1 = Partition([1, 1, 1, 1, 1, 1, 1, 1, 1])
        p2 = Partition([0, 1, 2, 3, 4, 5, 6, 7, 8])
        self.assertTrue(p1.is_minimal())
        self.assertFalse(p2.is_minimal())
        self.assertFalse(self.partition.is_minimal())

    def test_random(self):
        # Be aware, this has been known to spontaneously fail - problem with testing random things
        p = Partition.random([1, 1, 1], 12)
        self.assertEqual(p.num_groups(), 3)


class CollectionTests(unittest.TestCase):
    def setUp(self):
        self.c = treeCl.Collection(input_dir=os.path.join(thisdir, 'data'), file_format='phylip',
                                   show_progressbars=False)

    def test_len(self):
        self.assertEqual(len(self.c), 15)

    def test_names(self):
        expected = ['class{}_{}'.format(a,b) for a in [1,2,3] for b in [1,2,3,4,5]]
        self.assertEqual(self.c.names, expected)

    def test_lengths(self):
        expected = [511, 420, 139, 635, 229, 593, 228, 218, 404, 167, 227, 650, 104, 1172, 93]
        self.assertEqual(self.c.lengths, expected)


class TreeTests(unittest.TestCase):
    def test_RandomTree_defaultnames(self):
        t = treeCl.tree.RandomTree.new(10)
        expected = ['l1', 'l10', 'l2', 'l3', 'l4', 'l5', 'l6', 'l7', 'l8', 'l9']
        self.assertEqual(sorted(t.labels), expected)

    def test_RandomTree_given_names(self):
        names = ['Jools', 'Jops', 'Stoo', 'Rj', 'Ubik', 'Cj', 'Chris', 'Pete', 'Tadger', 'Hector']
        t = treeCl.tree.RandomTree.new(10, names)
        self.assertEqual(sorted(t.labels), sorted(names))


class DistanceMatrixTests(unittest.TestCase):
    def test_from_csv(self):
        dm = treeCl.DistanceMatrix.from_csv(os.path.join(thisdir, 'data', 'cache', 'geo_dm.csv'))
        self.assertEqual(dm.df.values.sum(), 412.70677069540181)


class EMTests(unittest.TestCase):
    def setUp(self):
        pass

def main():
    unittest.main()


if __name__ == '__main__':
    main()
