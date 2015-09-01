from pyparsing import Suppress, SkipTo, Word, Regex, Literal, OneOrMore, Group, LineEnd, CharsNotIn, nums, alphanums, ParseException
import logging
logger = logging.getLogger(__name__)

FLOAT = Word(nums + '.-').setParseAction(lambda x: float(x[0]))
INT = Word(nums).setParseAction(lambda x: float(x[0]))
WORD = Word(alphanums)

class PhymlParser(object):
    """
    Simple phyml result parser. Assumes GTR model for nucleotide analyses.
    """

    def __init__(self):
        self.MODEL_LABEL = Regex(r'Model of.*substitution:\s+')
        self.ALPHA_LABEL = Regex(r'Gamma shape parameter:\s+')
        self.LNL_LABEL = Regex(r'Log-likelihood:\s+')
        self.F_LABEL = Regex(r'f\(([ACGT])\)=\s+')
        self.R_LABEL = Regex(r'[ACGT]\s+<->\s+[ACGT]\s+')
        self.model = Suppress(SkipTo(self.MODEL_LABEL)) + Suppress(self.MODEL_LABEL) + WORD
        self.lnl = Suppress(SkipTo(self.LNL_LABEL)) + Suppress(self.LNL_LABEL) + FLOAT
        self.alpha = Suppress(SkipTo(self.ALPHA_LABEL)) + Suppress(self.ALPHA_LABEL) + FLOAT
        self.common = self.model + self.lnl + self.alpha
        self.freq = OneOrMore(Group(Suppress(SkipTo(self.F_LABEL)) + Suppress(self.F_LABEL) + FLOAT))
        self.rates = OneOrMore(Group(Suppress(SkipTo(self.R_LABEL)) + Suppress(self.R_LABEL) + FLOAT))
        self.gtr_specific = Group(self.freq) + Group(self.rates)

    def parse(self, filename):
        model = None
        alpha = None
        lnl = None
        freq = None
        rates = None

        with open(filename) as fl:
            s = fl.read()

        try:
            model, lnl, alpha = self.common.parseString(s).asList()

        except ParseException as err:
            logger.error(err)

        if model == 'GTR':
            try:
                freq, rates = self.gtr_specific.parseString(s).asList()
            except ParseException as err:
                logger.error(err)

        return model, alpha, lnl, freq, rates

    def to_dict(self, stats_filename, tree_filename):
        model, alpha, lnl, freq, rates = self.parse(stats_filename)
        try:
            with open(tree_filename) as treefl:
                tree = treefl.read().rstrip()
        except IOError as err:
            logger.error(err)
            return

        result = {'likelihood': lnl[0],
                  'partitions': {0: {'alpha': alpha,
                                     'frequencies': freq},
                                     'rates': rates,
                                     'model': model},
                  'ml_tree': tree}
        return result

pp = PhymlParser()
res_aa = pp.parse('class1_1.phy_phyml_stats')
res_nt = pp.parse('/nfs/research/goldman/kevin/testing_incomplete_taxon_occupancy/scripts/1/concats/415ce5a8740af3a4d5e33dfe5e73f7df00c0d889.phy_phyml_stats')

class RaxmlParser(object):

    def __init__(self):
        self.ALPHA_LABEL = Regex(r'alpha\[\d+\]:')
        self.LNL_LABEL = Literal('GAMMA-based Likelihood:')
        self.FRQ_LABEL = Regex(r'Base frequencies: (?=\d+)') ^ Regex(r'ML estimate base freqs\[\d+\]:')
        self.NAMES_LABEL = Regex(r'Partition: \d+ with name:\s+')
        self.RATES_LABEL = Regex(r'rates\[\d+\].+?:')
        self.alpha = OneOrMore(Suppress(SkipTo(self.ALPHA_LABEL)) + Suppress(self.ALPHA_LABEL) + FLOAT)
        self.lnl = Suppress(SkipTo(self.LNL_LABEL)) + Suppress(self.LNL_LABEL) + FLOAT
        self.frq = OneOrMore(Group(Suppress(SkipTo(self.FRQ_LABEL)) + Suppress(self.FRQ_LABEL) + OneOrMore(FLOAT)))
        self.names = OneOrMore(Suppress(SkipTo(self.NAMES_LABEL)) + Suppress(self.NAMES_LABEL) + CharsNotIn('\n') + Suppress(LineEnd()))
        self.rates = OneOrMore(Group(Suppress(SkipTo(self.RATES_LABEL)) + Suppress(self.RATES_LABEL) + OneOrMore(FLOAT)))

    def parse(self, filename):
        with open(filename) as fl:
            s = fl.read()

        try:
            alphas = self.alpha.parseString(s).asList()
        except ParseException as err:
            logger.error(err)
            alphas = [None]
        try:
            freqs = self.frq.parseString(s).asList()
        except ParseException as err:
            logger.error(err)
            freqs = [None]
        try:
            names = self.names.parseString(s).asList()
        except ParseException as err:
            logger.error(err)
            names = [None]
        try:
            rates = self.rates.parseString(s).asList()
        except ParseException:
            rates = None
        try:
            lnl = self.lnl.parseString(s).asList()
        except ParseException as err:
            logger.error(err)
            lnl = [0]

        return alphas, freqs, names, rates, lnl

    def to_dict(self, info_filename, tree_filename):
        alpha, freqs, names, rates, lnl = self.parse(info_filename)
        try:
            with open(tree_filename) as fl:
                tree = fl.read().rstrip()
        except IOError as err:
            logger.error('No tree file - raxml analysis failed')
            return
        n_parts = len(alpha)
        assert len(freqs) == n_parts
        assert len(names) == n_parts
        if rates is not None:
            assert len(rates) == n_parts

        result = {'likelihood': lnl[0],
                  'partitions': {},
                  'ml_tree': tree}

        for i in range(n_parts):
            subdict = {}
            subdict['alpha'] = alpha[i]
            subdict['frequencies'] = freqs[i]
            subdict['name'] = names[i]
            if rates is not None:
                subdict['rates'] = rates[i]
            result['partitions'][i] = subdict

        return result

rp = RaxmlParser()
