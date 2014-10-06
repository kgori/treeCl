import bpp
import collections
from numpy import log
import random
import tempfile
from tree import Tree
from interfacing import pll
from utils import fileIO


class Alignment(bpp.Alignment):
    def __init__(self, *args):
        super(Alignment, self).__init__(*args)
        self.infile = None
        self.name = None
        if len(args) > 0 and isinstance(args[0], basestring) and fileIO.can_locate(args[0]):
            self.infile = args[0]

    def __add__(self, other):
        return self.__class__([self, other])

    def __str__(self):
        contents = self.get_sequences()
        output_string = 'Alignment: {}\n'.format(self.name)

        return output_string + ''.join(
            ['>{}\n{}\n... ({}) ...\n{}\n'.format(header, seq[:50], len(seq) - 100, seq[-50:]) for header, seq in
             contents])

    @property
    def tree(self):
        try:
            return Tree(self.get_tree())
        except:
            return Tree(self.get_bionj_tree())

    def read_alignment(self, *args, **kwargs):
        super(Alignment, self).read_alignment(*args, **kwargs)
        self.infile = args[0]

    def pll_get_instance(self, *args):
        tmpdir = tmpfile = None
        try:
            with open(self.infile):
                pass
            alignment = self.infile
            instance = pll.create_instance(alignment, *args)
            return instance

        except (IOError, TypeError):
            tmpdir = tempfile.mkdtemp()
            _, tmpfile = tempfile.mkstemp(dir=tmpdir)
            self.write_alignment(tmpfile, "phylip", interleaved=True)
            instance = pll.create_instance(tmpfile, *args)
            return instance

        finally:
            if tmpfile:
                os.remove(tmpfile)
            if tmpdir:
                os.rmdir(tmpdir)

    def pll_optimise(self, partitions, model=None, nthreads=1, opt_subst=True, seed=None):
        """
        Runs the full raxml search algorithm. Model parameters are set using a combination of partitions and model.
        Optimisation of substitution model parameters is enabled with opt_subst=True.
        :param partitions: Links substitution models to alignment sites. Format is the same as the file RAxML uses
         with the -q flag (e.g. DNA, gene1codon1 = 1-500/3 - see RAxML manual)
        :param model: Dictionary to set model parameters--rates, frequencies and gamma alpha parameters--for each
         partition.
        :param opt_subst: bool: optimise substitution model parameters (T/F).
        :return: Dictionary of optimisation results.
        """

        instance = None
        seed = seed or random.randint(1, 99999)
        tree = self.tree.newick if self.tree else True
        try:
            instance = self.pll_get_instance(partitions, tree, nthreads, seed)
            if model is not None:
                pll.set_params_from_dict(instance, model)
            instance.optimise_tree_search(opt_subst)
            return pll.pll_to_dict(instance)
        except ValueError as exc:
            raise exc
        except Exception as exc:
            raise pll.PLLException(exc.message)
        finally:
            del instance

    def get_unconstrained_likelihood(self):
        weights = collections.Counter()
        sites = self.get_sites()
        n = len(sites)
        ucl = 0

        for site in sites:
            weights[site] += 1

        for v in weights.itervalues():
            if v > 1:
                ucl += v * log(v)
        return ucl - n * log(n)
