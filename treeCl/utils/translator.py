from builtins import map
#!/usr/bin/env python
from UserDict import UserDict
import re


class Translator(UserDict):
    """ From http://code.activestate.com/recipes/81330-single-pass-multiple-
    replace/ """

    def _make_regex(self):
        """ Build a regular expression object based on the keys of
        the current dictionary """

        return re.compile("(%s)" % "|".join(map(re.escape, list(self.keys()))))

    def __call__(self, mo):
        """ This handler will be invoked for each regex match """

        # Count substitutions
        self.count += 1  # Look-up string

        return self[mo.string[mo.start():mo.end()]]

    def translate(self, text):
        """ Translate text, returns the modified text. """

        # Reset substitution counter
        self.count = 0

        # Process text
        return self._make_regex().sub(self, text)
