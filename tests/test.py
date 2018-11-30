#!/usr/bin/env python
import unittest
import treeCl
import os, shutil
from treeCl import Partition, Alignment
from treeCl.utils.misc import binom_coeff

thisdir = os.path.dirname(os.path.realpath(__file__))

class AlignmentTests(unittest.TestCase):
    def test_read_phylip_file(self):
        filename = os.path.join(thisdir, 'data', 'mini', 'class1_1.phy')
        al = Alignment(filename, 'phylip')
        expected = ['Sp1', 'Sp2', 'Sp3', 'Sp4', 'Sp5']
        self.assertListEqual(expected, al.get_names())

    def test_read_gzip_phylip_file(self):
        filename = os.path.join(thisdir, 'data', 'mini', 'class1_1.phy.gz')
        al = Alignment(filename, 'phylip')
        expected = ['Sp1', 'Sp2', 'Sp3', 'Sp4', 'Sp5']
        self.assertListEqual(expected, al.get_names())

    def test_read_bzip2_phylip_file(self):
        filename = os.path.join(thisdir, 'data', 'mini', 'class1_1.phy.bz2')
        al = Alignment(filename, 'phylip')
        expected = ['Sp1', 'Sp2', 'Sp3', 'Sp4', 'Sp5']
        self.assertListEqual(expected, al.get_names())

    def test_read_fasta_file(self):
        filename = os.path.join(thisdir, 'data', 'mini', 'class1_1.fas')
        al = Alignment(filename, 'fasta')
        expected = ['Sp1', 'Sp2', 'Sp3', 'Sp4', 'Sp5']
        self.assertListEqual(expected, al.get_names())

    def test_read_gzip_fasta_file(self):
        filename = os.path.join(thisdir, 'data', 'mini', 'class1_1.fas.gz')
        al = Alignment(filename, 'fasta')
        expected = ['Sp1', 'Sp2', 'Sp3', 'Sp4', 'Sp5']
        self.assertListEqual(expected, al.get_names())

    def test_read_bzip2_fasta_file(self):
        filename = os.path.join(thisdir, 'data', 'mini', 'class1_1.fas.bz2')
        al = Alignment(filename, 'fasta')
        expected = ['Sp1', 'Sp2', 'Sp3', 'Sp4', 'Sp5']
        self.assertListEqual(expected, al.get_names())

    def compute_distance_correct_result(self):
        pass

    def compute_distance_throws_datatype_mismatch(self):
        self.assertRaises(ValueError)



class PartitionTests(unittest.TestCase):
    def setUp(self):
        self.partition = Partition(['a', 'd', 'a', 'a', 'b', 'a', 'b', 'c', 'c', 'd', 'd', 'd', 'd', 'd'])

    def test_restricted_growth_notation(self):
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
        p = Partition.random([12, 12, 12], 12)
        self.assertEqual(p.num_groups(), 3)

    def test_equal(self):
        p1 = Partition([0,0,0,0,0,2,2,2,2,2,1,1,1,1,1])
        p2 = Partition([2,2,2,2,2,1,1,1,1,1,0,0,0,0,0])
        self.assertEqual(p1, p2)

    def test_not_equal(self):
        p1 = Partition([0,0,0,0,0,0,0,0,1,1,1,1,1,1,1,1])
        p2 = Partition([0,0,0,0,0,0,0,0,1,1,1,1,1,1,1,0])
        self.assertNotEqual(p1, p2)


class CollectionTests(unittest.TestCase):
    def setUp(self):
        self.c = treeCl.Collection(input_dir=os.path.join(thisdir, 'data'), file_format='phylip',
                                   show_progress=False)

    def test_len(self):
        self.assertEqual(len(self.c), 15)

    def test_names(self):
        expected = ['class{}_{}'.format(a,b) for a in [1,2,3] for b in [1,2,3,4,5]]
        self.assertEqual(self.c.names, expected)

    def test_lengths(self):
        expected = [511, 420, 139, 635, 229, 593, 228, 218, 404, 167, 227, 650, 104, 1172, 93]
        self.assertEqual(self.c.lengths, expected)

    def test_read_parameters(self):
        self.c = treeCl.Collection(input_dir=os.path.join(thisdir, 'data'),
                                   param_dir=os.path.join(thisdir, 'data', 'cache'),
                                   file_format='phylip',
                                   show_progress=False)
        rec = self.c[0]
        self.assertEqual(rec.parameters.nj_tree[:72],
                         '((((Sp1:1.47856,(Sp4:1.20999,((Sp8:0.00595845,Sp9:0.00469589):0.27853,Sp')

    def test_read_trees(self):
        self.c = treeCl.Collection(input_dir=os.path.join(thisdir, 'data'),
                                   trees_dir=os.path.join(thisdir, 'data', 'trees'),
                                   file_format='phylip',
                                   show_progress=False)
        rec = self.c[0]
        self.assertEqual(rec.parameters.ml_tree[:72],
                         '((((Sp1:1.48316688535948748573,(Sp4:1.16694627918414717271,((Sp8:0.00749')


class ScorerTests(unittest.TestCase):

    def setUp(self):
        self.workingdir = 'wdir'

    def tearDown(self):
        if os.path.isdir(self.workingdir):
            shutil.rmtree(self.workingdir)

    def test_scorer_can_write(self):
        c = treeCl.Collection(input_dir=os.path.join(thisdir, 'data'),
                              param_dir=os.path.join(thisdir, 'data', 'cache'),
                              file_format='phylip',
                              show_progress=False)

        raxml = treeCl.tasks.RaxmlTaskInterface()
        sc = treeCl.Scorer(c, cache_dir=self.workingdir, task_interface=raxml)
        p = treeCl.Partition([0,0,0,0,0,1,1,1,1,1,2,2,2,2,2])
        sc.write_partition(p)

        # check files were written
        import glob
        files = glob.glob(os.path.join(self.workingdir, '*.phy'))

        self.assertTrue(len(files)>0)


class TreeTests(unittest.TestCase):
    def test_random_tree_default_names(self):
        t = treeCl.tree.RandomTree.new(10)
        expected = ['l1', 'l10', 'l2', 'l3', 'l4', 'l5', 'l6', 'l7', 'l8', 'l9']
        self.assertEqual(sorted(t.labels), expected)

    def test_random_tree_given_names(self):
        names = ['Jools', 'Jops', 'Stoo', 'Rj', 'Ubik', 'Cj', 'Chris', 'Pete', 'Tadger', 'Hector']
        t = treeCl.tree.RandomTree.new(10, names)
        self.assertEqual(sorted(t.labels), sorted(names))

    def test_newick_read(self):
        n = ('((((Sp1:0.0523947839547,Sp2:1.37159604411):2.36974538201,((Sp3:0.179093762783,'
             '(Sp4:0.615505083102,Sp5:0.344065892719):0.0724725996144):0.307962158157,'
             '(Sp6:1.48158479406,Sp7:3.13329090451):1.62357461752):0.62792640958):2.64647302212,'
             '(Sp8:0.145879857199,Sp9:4.33463301328):0.785221836876):0.0653625005117,'
             '((((Sp10:0.0327158596802,Sp11:0.346629825105):0.513499606131,Sp12:0.0931894502388):1.75462968872,'
             '(Sp13:0.0508398281971,Sp14:0.902409030743):0.248348229186):2.66397475192,(Sp15:0.623334704667,'
             'Sp16:0.727987265987):2.45688940891):1.00011564391):0.0;')
        self.assertEqual(treeCl.tree.Tree(n).newick, n)

    def test_rnni(self):
        # Issue 15: Tree.rnni raises a type error because it calls an unimported function
        # (and also shadows the fn name with a bool variable). This test asserts that rnni
        # doesn't throw any error on valid input.

        t = treeCl.tree.Tree(newick="((((T4:42.9474018906,T10:42.9474018906):112.903906732,(T6:14.3433500048,"
                                    "(T2:1.53929863217,T5:1.53929863217):12.8040513726):141.507958618):22.1730692315,"
                                    "T9:178.024377854):34.9689886128,(T3:190.011180702,((T1:0,T8:0):147.182729024,"
                                    "T7:147.182729024):42.8284516785):22.9821857647):2.44503258424;")
        try:
            t.rnni(times=3)
        except:
            self.fail('Tree.rnni() raised an exception unexpectedly')

    def test_rspr(self):
        # Issue 15: Tree.rnni raises a type error because it calls an unimported function
        # (and also shadows the fn name with a bool variable). This test asserts that rspr
        # also doesn't throw any error on valid input.

        t = treeCl.tree.Tree(newick="((((T4:42.9474018906,T10:42.9474018906):112.903906732,(T6:14.3433500048,"
                                    "(T2:1.53929863217,T5:1.53929863217):12.8040513726):141.507958618):22.1730692315,"
                                    "T9:178.024377854):34.9689886128,(T3:190.011180702,((T1:0,T8:0):147.182729024,"
                                    "T7:147.182729024):42.8284516785):22.9821857647):2.44503258424;")
        try:
            t.rspr(times=3)
        except:
            self.fail('Tree.rspr() raised an exception unexpectedly')

    def test_rils(self):
        # Issue 15: Tree.rnni raises a type error because it calls an unimported function
        # (and also shadows the fn name with a bool variable). This test asserts that ILS.rils
        # also doesn't throw any error on valid input.

        t = treeCl.tree.Tree(
            newick="((((T4:42.9474018906,T10:42.9474018906):112.903906732,(T6:14.3433500048,(T2:1.53929863217,"
                   "T5:1.53929863217):12.8040513726):141.507958618):22.1730692315,T9:178.024377854):34.9689886128,"
                   "(T3:190.011180702,((T1:0,T8:0):147.182729024,"
                   "T7:147.182729024):42.8284516785):22.9821857647):2.44503258424;")
        ils = treeCl.tree.ILS(t)
        try:
            ils.rils()
        except:
            self.fail('ILS.rils() raised an exception unexpectedly')


class TreeDistanceTests(unittest.TestCase):
    def setUp(self):
        self.c = treeCl.Collection(input_dir=os.path.join(thisdir, 'data'),
                                   trees_dir=os.path.join(thisdir, 'data', 'trees'),
                                   file_format='phylip',
                                   show_progress=False)
        self.tree1 = treeCl.Tree(self.c[0].tree)
        self.tree2 = treeCl.Tree(self.c[1].tree)

    def test_geo(self):
        distance = treeCl.treedist.geodist(self.tree1, self.tree2, False)
        self.assertAlmostEqual(distance, 0.7051762453884542)

    def test_geo_normalised(self):
        distance = treeCl.treedist.geodist(self.tree1, self.tree2, True)
        self.assertAlmostEqual(distance, 0.11511918913951155)

    def test_euc(self):
        distance = treeCl.treedist.eucdist(self.tree1, self.tree2, False)
        self.assertAlmostEqual(distance, 0.6463384798834909)

    def test_euc_normalised(self):
        distance = treeCl.treedist.eucdist(self.tree1, self.tree2, True)
        self.assertAlmostEqual(distance, 0.10551399341715575)

    def test_wrf(self):
        distance = treeCl.treedist.wrfdist(self.tree1, self.tree2, False)
        self.assertAlmostEqual(distance, 2.1666480728499717)

    def test_wrf_normalised(self):
        distance = treeCl.treedist.wrfdist(self.tree1, self.tree2, True)
        self.assertAlmostEqual(distance, 0.1063259354618597)

    def test_rf(self):
        distance = treeCl.treedist.rfdist(self.tree1, self.tree2, False)
        self.assertAlmostEqual(distance, 2.0)

    def test_rf_normalised(self):
        distance = treeCl.treedist.rfdist(self.tree1, self.tree2, True)
        self.assertAlmostEqual(distance, 0.14285714285714285)

    def test_intertree(self):
        dm = self.c.get_inter_tree_distances('rf', show_progress=False)
        self.assertAlmostEquals(dm.values[0, 1],
                          treeCl.treedist.rfdist(self.tree1, self.tree2, False))

    def test_intertree_partial_overlap(self):
        self.c[0].parameters.ml_tree = treeCl.Tree(self.c[0].parameters.ml_tree).prune_to_subset(
            set(['Sp1', 'Sp2', 'Sp3', 'Sp4', 'Sp5', 'Sp6', 'Sp7'])).newick
        self.c[1].parameters.ml_tree = treeCl.Tree(self.c[1].parameters.ml_tree).prune_to_subset(
            set(['Sp4', 'Sp5', 'Sp6', 'Sp7', 'Sp8', 'Sp9', 'Sp10'])).newick
        dm = self.c.get_inter_tree_distances('geo', show_progress=False)
        self.assertAlmostEquals(dm.values[0, 1], 0.35314280333144532)

    def test_intertree_partial_overlap_value_replacement(self):
        self.c[0].parameters.ml_tree = treeCl.Tree(self.c[0].parameters.ml_tree).prune_to_subset(
            set(['Sp1', 'Sp2', 'Sp3', 'Sp4', 'Sp5', 'Sp6', 'Sp7'])).newick
        self.c[1].parameters.ml_tree = treeCl.Tree(self.c[1].parameters.ml_tree).prune_to_subset(
            set(['Sp4', 'Sp5', 'Sp6', 'Sp7', 'Sp8', 'Sp9', 'Sp10'])).newick
        dm = self.c.get_inter_tree_distances('geo', min_overlap=5, overlap_fail_value=-1, show_progress=False)
        self.assertEquals(dm.values[0, 1], -1)


class DistanceMatrixTests(unittest.TestCase):
    def test_from_csv(self):
        dm = treeCl.DistanceMatrix.from_csv(os.path.join(thisdir, 'data', 'cache', 'geo_dm.csv'))
        self.assertAlmostEqual(dm.df.values.sum(), 412.70677069540181)

    def test_get_names(self):
        dm = treeCl.DistanceMatrix.from_csv(os.path.join(thisdir, 'data', 'cache', 'geo_dm.csv'))
        names = dm.get_names()
        self.assertListEqual(names, ['class1_1', 'class1_2',
                                     'class1_3', 'class1_4',
                                     'class1_5', 'class2_1',
                                     'class2_2', 'class2_3',
                                     'class2_4', 'class2_5',
                                     'class3_1', 'class3_2',
                                     'class3_3', 'class3_4',
                                     'class3_5'])

    def test_set_names(self):
        dm = treeCl.DistanceMatrix.from_csv(os.path.join(thisdir, 'data', 'cache', 'geo_dm.csv'))
        new_names = list('abcdefghijklmno')
        dm.set_names(new_names)
        self.assertListEqual(dm.get_names(), new_names)

    def test_calculation(self):
        c = treeCl.Collection(input_dir=os.path.join(thisdir, 'data'),
                              param_dir=os.path.join(thisdir, 'data', 'cache'),
                              file_format='phylip',
                              show_progress=False)
        dm = c.get_inter_tree_distances('geo', show_progress=False)
        self.assertAlmostEqual(dm.df.values.sum(), 412.70677069540181)

    def test_embed_cmds(self):
        dm = treeCl.DistanceMatrix.from_csv(os.path.join(thisdir, 'data', 'cache', 'geo_dm.csv'))
        embed = dm.embedding(3, 'cmds')
        self.assertEqual(embed.shape, (15, 3))

    def test_embed_kpca(self):
        dm = treeCl.DistanceMatrix.from_csv(os.path.join(thisdir, 'data', 'cache', 'geo_dm.csv'))
        embed = dm.embedding(3, 'kpca')
        self.assertEqual(embed.shape, (15, 3))

    def test_embed_mmds(self):
        dm = treeCl.DistanceMatrix.from_csv(os.path.join(thisdir, 'data', 'cache', 'geo_dm.csv'))
        embed = dm.embedding(3, 'mmds')
        self.assertEqual(embed.shape, (15, 3))

    def test_embed_nmmds(self):
        dm = treeCl.DistanceMatrix.from_csv(os.path.join(thisdir, 'data', 'cache', 'geo_dm.csv'))
        embed = dm.embedding(3, 'nmmds')
        self.assertEqual(embed.shape, (15, 3))

    def test_embed_spectral(self):
        dm = treeCl.DistanceMatrix.from_csv(os.path.join(thisdir, 'data', 'cache', 'geo_dm.csv'))
        embed = dm.embedding(3, 'spectral')
        self.assertEqual(embed.shape, (15, 3))

    def test_embed_tsne(self):
        dm = treeCl.DistanceMatrix.from_csv(os.path.join(thisdir, 'data', 'cache', 'geo_dm.csv'))
        tsne = dm.embedding(3, 'tsne')
        self.assertEqual(tsne.shape, (15, 3))


class ClusteringTests(unittest.TestCase):
    def setUp(self):
        self.dm = treeCl.DistanceMatrix.from_csv(os.path.join(thisdir, 'data', 'cache', 'geo_dm.csv'))

    def test_spectral_kmeans(self):
        cl = treeCl.Spectral(self.dm)
        p = cl.cluster(3, method=treeCl.clustering.methods.KMEANS)
        self.assertEqual(len(p), 3)

    def test_spectral_gmm(self):
        cl = treeCl.Spectral(self.dm)
        p = cl.cluster(3, method=treeCl.clustering.methods.GMM)
        self.assertEqual(len(p), 3)

class RaxmlParserTests(unittest.TestCase):
    def setUp(self):
        self.parser = treeCl.parsers.RaxmlParser()
        self.info = os.path.join(thisdir, 'data', 'parsing', 'RAxML_info.modopt')
        self.result = os.path.join(thisdir, 'data', 'parsing', 'RAxML_result.modopt')
        self.infoq = os.path.join(thisdir, 'data', 'parsing', 'RAxML_info.modoptq')
        self.resultq = os.path.join(thisdir, 'data', 'parsing', 'RAxML_result.modoptq')

    def test_can_parse_default_name(self):
        parse_result = self.parser.to_dict(self.info, self.result, True)
        self.assertEqual(parse_result['partitions'][0]['name'], 'No Name Provided')

    def test_can_parse_provided_name(self):
        parse_result = self.parser.to_dict(self.infoq, self.resultq, True)
        self.assertEqual(parse_result['partitions'][0]['name'], 'class1_1')

    def test_can_parse_dna_multiple_starts(self):
        info = os.path.join(thisdir, 'data', 'parsing', 'RAxML_info.dnamulti')
        tree = os.path.join(thisdir, 'data', 'parsing', 'RAxML_bestTree.dnamulti')
        parse_result = self.parser.to_dict(info, tree, False)
        self.assertEqual(parse_result['partitions'][0]['rates'],
                         [0.713025, 1.907700, 1.010704, 0.900610, 1.628844, 1.000000])

    def test_can_parse_aa_multiple_starts(self):
        info = os.path.join(thisdir, 'data', 'parsing', 'RAxML_info.aamulti')
        tree = os.path.join(thisdir, 'data', 'parsing', 'RAxML_bestTree.aamulti')
        parse_result = self.parser.to_dict(info, tree, False)
        self.assertEqual(parse_result['partitions'][0]['frequencies'],
                         [0.058725, 0.051571, 0.046151, 0.061140, 0.021702,
                          0.041016, 0.072847, 0.073271, 0.029989, 0.065176,
                          0.094756, 0.063071, 0.025892, 0.028991, 0.049385,
                          0.052040, 0.053462, 0.012507, 0.019707, 0.078602])


class RaxmlRunnerTests(unittest.TestCase):
    def setUp(self):
        self.c = treeCl.Collection(input_dir=os.path.join(thisdir, 'data', 'mini'), file_format='phylip',
                                   show_progress=False)

    def test_can_run_gamma(self):
        self.c.calc_trees(model='PROTGAMMAWAG', show_progress=False)
        self.assertFalse(self.c[0].parameters.ml_tree is None)

    def test_can_run_cat(self):
        self.c.calc_trees(model='PROTCATWAG', show_progress=False)
        self.assertFalse(self.c[0].parameters.ml_tree is None)

    def test_can_run_on_dna(self):
        self.c = treeCl.Collection(input_dir=os.path.join(thisdir, 'data', 'dna_alignments'), file_format='phylip',
                                   show_progress=False)
        self.c.calc_trees(indices=[0], model='GTRGAMMA', show_progress=False)
        self.assertFalse(self.c[0].parameters.ml_tree is None)

    def test_can_run_fast_tree(self):
        self.c.calc_trees(indices=[0], fast_tree=True, model='PROTGAMMALGF', show_progress=False)
        self.assertFalse(self.c[0].parameters.ml_tree is None)

    def test_can_run_multiple_starts(self):
        self.c = treeCl.Collection(input_dir=os.path.join(thisdir, 'data', 'dna_alignments'), file_format='phylip',
                                   show_progress=False)
        self.c.calc_trees(indices=[0], model='GTRGAMMA', fast_tree=True, show_progress=False, n_starts=2)
        self.assertFalse(self.c[0].parameters.ml_tree is None)


class ParallelTests(unittest.TestCase):
    def setUp(self):
        self.c = treeCl.Collection(input_dir=os.path.join(thisdir, 'data'),
                                   param_dir=os.path.join(thisdir, 'data', 'cache'),
                                   file_format='phylip',
                                   show_progress=False)

    def test_sequential(self):
        handler = treeCl.parutils.SequentialJobHandler()
        dm = self.c.get_inter_tree_distances('geo', jobhandler=handler, show_progress=False)
        self.assertAlmostEqual(dm.df.values.sum(), 412.70677069540181)

    def test_processpool(self):
        handler = treeCl.parutils.ProcesspoolJobHandler(2)
        dm = self.c.get_inter_tree_distances('geo', jobhandler=handler, show_progress=False)
        self.assertAlmostEqual(dm.df.values.sum(), 412.70677069540181)

    def test_threadpool(self):
        handler = treeCl.parutils.ThreadpoolJobHandler(2)
        dm = self.c.get_inter_tree_distances('geo', jobhandler=handler, show_progress=False)
        self.assertAlmostEqual(dm.df.values.sum(), 412.70677069540181)


class MiscTests(unittest.TestCase):
    def test_binom_coeff(self):
        inputs =  [1, 2,  5, 10,   53,    3289]
        outputs = [0, 1, 10, 45, 1378, 5407116]
        testvals = [binom_coeff(n) for n in inputs]
        self.assertListEqual(testvals, outputs)


def main():
    unittest.main()


if __name__ == '__main__':
    main()
