from pyparsing import Suppress, SkipTo, Word, Regex, Literal, OneOrMore, Group, LineEnd, CharsNotIn, nums, ParseException

FLOAT = Word(nums + '.-').setParseAction(lambda x: float(x[0]))
INT = Word(nums).setParseAction(lambda x: float(x[0]))

class RaxmlParser(object):

    def __init__(self):
        self.ALPHA_LABEL = Regex('alpha\[\d+\]:')
        self.LNL_LABEL = Literal('GAMMA-based Likelihood:')
        self.FRQ_LABEL = Literal('Base frequencies:') ^ Regex('ML estimate base freqs\[\d+\]:')
        self.NAMES_LABEL = Regex('Partition: \d+ with name:\s+')
        self.RATES_LABEL = Regex('rates\[\d+\].+?:')
        self.alpha = OneOrMore(Suppress(SkipTo(self.ALPHA_LABEL)) + Suppress(self.ALPHA_LABEL) + FLOAT)
        self.lnl = Suppress(SkipTo(self.LNL_LABEL)) + Suppress(self.LNL_LABEL) + FLOAT
        self.frq = OneOrMore(Group(Suppress(SkipTo(self.FRQ_LABEL)) + Suppress(self.FRQ_LABEL) + OneOrMore(FLOAT)))
        self.names = OneOrMore(Suppress(SkipTo(self.NAMES_LABEL)) + Suppress(self.NAMES_LABEL) + CharsNotIn('\n') + Suppress(LineEnd()))
        self.rates = OneOrMore(Group(Suppress(SkipTo(self.RATES_LABEL)) + Suppress(self.RATES_LABEL) + OneOrMore(FLOAT)))

    def parse(self, filename):
        with open(filename) as fl:
            s = fl.read()
        alphas = self.alpha.parseString(s).asList()
        freqs = self.frq.parseString(s).asList()
        names = self.names.parseString(s).asList()
        try:
            rates = self.rates.parseString(s).asList()
        except ParseException:
            rates = None
        lnl = self.lnl.parseString(s).asList()
        
        return alphas, freqs, names, rates, lnl

    def to_dict(self, info_filename, tree_filename):
        alpha, freqs, names, rates, lnl = self.parse(info_filename)
        with open(tree_filename) as fl:
            tree = fl.read().rstrip()
        n_parts = len(alpha)
        assert len(freqs) == n_parts
        assert len(names) == n_parts
        if rates is not None:
            assert len(rates) == n_parts

        result = {'likelihood': lnl[0],
                  'partitions': {},
                  'mltree': tree}

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
