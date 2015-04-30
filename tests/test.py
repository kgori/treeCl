#!/usr/bin/env python
import unittest
import treeCl

class PartitionTests(unittest.TestCase):

    def setup(self):
        self.p1 = treeCl.Partition(['a','a','a','b','a','b','c','c'])

    def test_rgn(self):
        """
        Test restricted growth notation
        :return:
        """
        self.setup()
        expected = (0,0,0,1,0,1,2,2)
        self.failUnless(self.p1.partition_vector == expected)

    def test_numbers(self):
        self.setup()
        self.failUnless(self.p1.num_elements() == 8)
        self.failUnless(self.p1.num_groups() == 3)

def main():
    unittest.main()

if __name__ == '__main__':
    main()