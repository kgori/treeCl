from __future__ import absolute_import
from builtins import zip
from builtins import range
from builtins import object
from pyparsing import Suppress, SkipTo, Word, Regex, Literal, OneOrMore, Group, LineEnd, CharsNotIn, nums, alphanums, ParseException
import logging
logger = logging.getLogger(__name__)

FLOAT = Word(nums + '.-').setParseAction(lambda x: float(x[0]))
INT = Word(nums).setParseAction(lambda x: int(x[0]))
WORD = Word(alphanums+'_')
SPACEDWORD = Word(alphanums+' _')

class PhymlParser(object):
    """
    Simple phyml result parser. Assumes one of the standard models for nucleotide analyses.
    """

    def __init__(self):
        self.MODEL_LABEL = Regex(r'Model of.*substitution:\s+')
        self.ALPHA_LABEL = Regex(r'Gamma shape parameter:\s+')
        self.LNL_LABEL = Regex(r'Log-likelihood:\s+')
        self.F_LABEL = Regex(r'f\(([ACGT])\)=\s+')
        self.R_LABEL = Regex(r'[ACGT]\s+<->\s+[ACGT]\s+')
        self.TSTV_LABEL = Regex(r'Transition/transversion ratio.*:\s+')
        self.model = Suppress(SkipTo(self.MODEL_LABEL)) + Suppress(self.MODEL_LABEL) + WORD
        self.lnl = Suppress(SkipTo(self.LNL_LABEL)) + Suppress(self.LNL_LABEL) + FLOAT
        self.alpha = Suppress(SkipTo(self.ALPHA_LABEL)) + Suppress(self.ALPHA_LABEL) + FLOAT
        self.common = self.model + self.lnl + self.alpha
        self.tstv = OneOrMore(Suppress(SkipTo(self.TSTV_LABEL)) + Suppress(self.TSTV_LABEL) + FLOAT)
        self.freq = OneOrMore(Suppress(SkipTo(self.F_LABEL)) + Suppress(self.F_LABEL) + FLOAT)
        self.rates = OneOrMore(Suppress(SkipTo(self.R_LABEL)) + Suppress(self.R_LABEL) + FLOAT)
        self.gtr_specific = Group(self.freq) + Group(self.rates)
        self.hky_specific = Group(self.tstv) + Group(self.freq)

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

        if model == 'JC69':
            freq = [0.25, 0.25, 0.25, 0.25]
            rates = [1.0, 1.0, 1.0, 1.0, 1.0, 1.0]

        elif model == 'K80':
            freq = [0.25, 0.25, 0.25, 0.25]
            try:
                tstv = self.tstv.parseString(s).asList()
            except ParseException as err:
                logger.error(err)
            
            rates = [1.0, tstv[0], 1.0, 1.0, tstv[0], 1.0]

        elif model == 'F81':
            try:
                freq = self.freq.parseString(s).asList()
            except ParseException as err:
                logger.error(err)
            rates = [1.0, 1.0, 1.0, 1.0, 1.0, 1.0]

        elif model == 'F84' or model == 'HKY85' or model == 'TN93':
            parser = Group(self.tstv) + Group(self.freq)
            try:
                tstv, freq = parser.parseString(s).asList()
            except ParseException as err:
                logger.error(err)
            if model == 'TN93':
                rates = [1.0, tstv[0], 1.0, 1.0, tstv[1], 1.0]
            else:
                rates = [1.0, tstv[0], 1.0, 1.0, tstv[0], 1.0]

        elif model == 'GTR':
            parser = Group(self.freq) + Group(self.rates)
            try:
                freq, rates = parser.parseString(s).asList()
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

        result = {'likelihood': lnl,
                  'partitions': {0: {'alpha': alpha,
                                     'frequencies': freq,
                                     'rates': rates,
                                     'model': model}},
                  'ml_tree': tree}
        return result


class RaxmlParser(object):

    def __init__(self):
        self.ALPHA_LABEL = Regex(r'alpha\[\d+\]:')
        self.LNL_LABEL = Literal('Final GAMMA-based Score of best tree')
        self.FRQ_LABEL = Regex(r'Base frequencies: (?=\d+)') ^ Regex(r'ML estimate base freqs\[\d+\]:')
        self.NAMES_LABEL = Regex(r'Partition: \d+ with name:\s+')
        self.RATES_LABEL = Regex(r'rates\[\d+\].+?:')
        self.MODEL_LABEL = Literal('Substitution Matrix:')
        self.alpha = OneOrMore(Suppress(SkipTo(self.ALPHA_LABEL)) + Suppress(self.ALPHA_LABEL) + FLOAT)
        self.lnl = Suppress(SkipTo(self.LNL_LABEL)) + Suppress(self.LNL_LABEL) + FLOAT
        self.frq = OneOrMore(Group(Suppress(SkipTo(self.FRQ_LABEL)) + Suppress(self.FRQ_LABEL) + OneOrMore(FLOAT)))
        self.names = OneOrMore(Suppress(SkipTo(self.NAMES_LABEL)) + Suppress(self.NAMES_LABEL) + CharsNotIn('\n') + Suppress(LineEnd()))
        self.rates = OneOrMore(Group(Suppress(SkipTo(self.RATES_LABEL)) + Suppress(self.RATES_LABEL) + OneOrMore(FLOAT)))
        self.model = Suppress(SkipTo(self.MODEL_LABEL)) + Suppress(self.MODEL_LABEL) + WORD

        MODEL_LABEL = Literal('Substitution Matrix:')
        SCORE_LABEL = Literal('Final GAMMA  likelihood:')
        DESC_LABEL = Literal('Model Parameters of Partition')
        NAME_LEADIN = Literal(', Name:')
        DATATYPE_LEADIN = Literal(', Type of Data:')
        ALPHA_LEADIN = Literal('alpha:')
        TREELENGTH_LEADIN = Literal('Tree-Length:')
        RATES_LABEL = Regex(r'rate \w <-> \w:')
        FREQS_LABEL = Regex(r'freq pi\(\w\):')

        model = Suppress(SkipTo(MODEL_LABEL)) + Suppress(MODEL_LABEL) + WORD
        likelihood = Suppress(SkipTo(SCORE_LABEL)) + Suppress(SCORE_LABEL) + FLOAT
        description = (Suppress(SkipTo(DESC_LABEL)) +
                       Suppress(DESC_LABEL) + INT +
                       Suppress(NAME_LEADIN) +
                       SPACEDWORD +
                       Suppress(DATATYPE_LEADIN) +
                       WORD)
        alpha = Suppress(ALPHA_LEADIN) + FLOAT
        rates = Suppress(RATES_LABEL) + FLOAT
        freqs = Suppress(FREQS_LABEL) + FLOAT

        self._dash_f_e_parser = (Group(OneOrMore(model)) +
                                 likelihood +
                                 Group(OneOrMore(Group(description +
                                                       alpha +
                                                       Suppress(TREELENGTH_LEADIN) +
                                                       Suppress(FLOAT) +
                                                       Group(OneOrMore(rates)) +
                                                       Group(OneOrMore(freqs))
                                                       ))))

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

    def _dash_f_e_to_dict(self, info_filename, tree_filename):
        """
        Raxml provides an option to fit model params to a tree,
        selected with -f e.
        The output is different and needs a different parser.
        """
        with open(info_filename) as fl:
            models, likelihood, partition_params = self._dash_f_e_parser.parseFile(fl).asList()

        with open(tree_filename) as fl:
            tree = fl.read()

        d = {'likelihood': likelihood, 'ml_tree': tree, 'partitions': {}}

        for model, params in zip(models, partition_params):
            subdict = {}
            index, name, _, alpha, rates, freqs = params
            subdict['alpha'] = alpha
            subdict['name'] = name
            subdict['rates'] = rates
            subdict['frequencies'] = freqs
            subdict['model'] = model
            d['partitions'][index] = subdict

        return d

    def to_dict(self, info_filename, tree_filename, dash_f_e=False):
        """
        Parse raxml output and return a dict
        Option dash_f_e=True will parse the output of a raxml -f e run,
        which has different output
        """
        if dash_f_e:
            return self._dash_f_e_to_dict(info_filename, tree_filename)
        else:
            return self._to_dict(info_filename, tree_filename)

    def _to_dict(self, info_filename, tree_filename):
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
