"""
If something inherits from Talker, then we can print
text to the terminal in a relatively standard way.
"""

import numpy as np, pprint
import sys

if sys.version_info[0] < 3:
    input = raw_input


class Talker:
    """
    Objects the inherit from Talker have "mute" and "pithy" attributes,
    a report('uh-oh!') method that prints when unmuted,
    a speak('yo!') method that prints only when unmuted and unpithy,
    and an input("what's up?") method that takes input from the prompt.
    """

    _mute = False

    @property
    def _prefix(self):
        s = f"[{self.__class__.__name__.lower()}] "
        return f"{s:>16}"

    def _speak(
        self,
        string="",
        progress=False,
    ):
        """
        Print to the terminal If verbose=True, this will print to terminal. Otherwise, it won't.
        """
        if self._mute == False:
            equalspaces = " " * len(self._prefix)
            toprint = string + ""

            if progress:
                print(
                    "\r" + self._prefix + toprint.replace("\n", "\n" + equalspaces),
                )
            else:
                print(self._prefix + toprint.replace("\n", "\n" + equalspaces))
