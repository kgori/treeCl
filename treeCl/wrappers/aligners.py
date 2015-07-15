from .abstract_wrapper import AbstractWrapper

__all__ = ['Muscle', 'Prank', 'FSA']


class Muscle(AbstractWrapper):
    @property
    def _default_exe(self):
        return 'muscle'

    def _set_help(self):
        self(wait=True)
        self._help = self.get_stderr()


class Prank(AbstractWrapper):
    @property
    def _default_exe(self):
        return 'prank'

    @property
    def _hyphen_policy(self):
        return 1

    def _set_help(self):
        self(help=True, wait=True)
        self._help = self.get_stdout()


class FSA(AbstractWrapper):
    @property
    def _default_exe(self):
        return 'fsa'

    def _set_help(self):
        self(help=True, wait=True)
        self._help = self.get_stdout()


class Mafft(AbstractWrapper):
    @property
    def _default_exe(self):
        return 'mafft'

    def _set_help(self):
        self(help=True, wait=True)
        self._help = self.get_stdout()
