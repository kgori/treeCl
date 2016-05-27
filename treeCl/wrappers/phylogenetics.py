from builtins import range
from .abstract_wrapper import AbstractWrapper
from ..utils import smooth_freqs
import re

class Raxml(AbstractWrapper):
    @property
    def _default_exe(self):
        return 'raxmlHPC'

    def _set_help(self):
        self('-h', wait=True)
        self._help = self.get_stdout()


class Phyml(AbstractWrapper):
    @property
    def _default_exe(self):
        return 'phyml'

    def _set_help(self):
        self('-h', wait=True)
        _help = self.get_stdout()

        # Phyml's help string contains extra shell-escape characters to
        # add underlining and boldness to the text. Let's get rid of these
        import re
        rgx = re.compile(r'\x1b\[00;0[0-9]m')
        self._help = rgx.sub('', _help)


class FastTree(AbstractWrapper):
    @property
    def _default_exe(self):
        return 'FastTree'

    def _set_help(self):
        self(wait=True)
        self._help = self.get_stderr()

def parse_fasttree_output(s):
    try:
        loglk, alpha = (float(x) for x in re.search(r'Gamma\(\d+\) LogLk = ([0-9-.]+) alpha = ([0-9.]+)', s).groups())
    except AttributeError:
        raise AttributeError('Couldn\'t parse loglk and alpha')
        return None

    try:
        freqs_match = re.search(r'GTR Frequencies:\s+([0-9.-]+)\s+([0-9.-]+)\s+([0-9.-]+)\s+([0-9.-]+)', s)

    except AttributeError:
        raise AttributeError('Couldn\'t parse GTR frequencies')
        return None

    if freqs_match:
        freqs = []
        for i in range(1, 5):
            freqs.append(float(freqs_match.group(i)))
    else:
        freqs = []

    try:
        rates_match = re.search(
            r'GTR rates\(ac ag at cg ct gt\)\s+([0-9.-]+)\s+([0-9.-]+)\s+([0-9.-]+)\s+([0-9.-]+)\s+([0-9.-]+)\s+([0-9.-]+)\s+',s)

    except AttributeError:
        raise AttributeError('Couldn\'t parse GTR rates')
        return None

    if rates_match:
        rates = []
        for i in range(1, 7):
            rates.append(float(rates_match.group(i)))
    else:
        rates = []


    result = {'likelihood': loglk, 'partitions': {0: {'alpha': alpha}}}
    if freqs:
        result['partitions'][0]['frequencies'] = smooth_freqs(freqs)
    if rates:
        result['partitions'][0]['rates'] = rates
    return result

def parse_raxml_output(s):
    ac=float(re.search(r'rate A <-> C:\s+([0-9.]+)', s).groups()[0])
    ag=float(re.search(r'rate A <-> G:\s+([0-9.]+)', s).groups()[0])
    at=float(re.search(r'rate A <-> T:\s+([0-9.]+)', s).groups()[0])
    cg=float(re.search(r'rate C <-> G:\s+([0-9.]+)', s).groups()[0])
    ct=float(re.search(r'rate C <-> T:\s+([0-9.]+)', s).groups()[0])
    gt=float(re.search(r'rate G <-> T:\s+([0-9.]+)', s).groups()[0])
    a=float(re.search(r'freq pi\(A\):\s+([0-9.]+)',s).groups()[0])
    c=float(re.search(r'freq pi\(C\):\s+([0-9.]+)',s).groups()[0])
    g=float(re.search(r'freq pi\(G\):\s+([0-9.]+)',s).groups()[0])
    t=float(re.search(r'freq pi\(T\):\s+([0-9.]+)',s).groups()[0])
    alpha = float(re.search(r'alpha:\s+([0-9.]+)' ,s).groups()[0])
    loglk = float(re.search(r'Final GAMMA  likelihood:\s+([0-9-.]+)', s).groups()[0])

    return {
        'likelihood': loglk,
        'alpha': alpha,
        'frequencies': [a, c, g, t],
        'rates': [ac, ag, at, cg, ct, gt],
    }
