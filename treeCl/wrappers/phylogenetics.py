from .abstract_wrapper import AbstractWrapper


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
        return 'fasttree'

    def _set_help(self):
        self(wait=True)
        self._help = self.get_stderr()
