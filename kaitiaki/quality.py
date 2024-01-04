import kaitiaki
import numpy as np


class QualityControl:
    def __init__(self, dirname,
                 binary_mode=False,
                 termination_reason='finished'):

        self.binary_mode = binary_mode
        self.dirname = dirname

        self._out = None
        self._out2 = None

        self._plot = kaitiaki.file.plot(f"{dirname}/plot")
        self._plot2 = None

        if self.binary_mode:
            self._plot2 = kaitiaki.file.plot(f"{dirname}/plot2")

        self._comments = []
        self._reason = termination_reason

    def load_outfile(self):
        self._out = kaitiaki.file.out(f"{dirname}/out")

        if self.binary_mode:
            self._out2 = kaitiaki.file.out(f"{dirname}/out2")

    def __score_reason(self):
        return {
            'finished': 0.0,
            'timeout': 0.5,
            'other': 1.0,
        }[self._reason]

    def __score_properties(self):

        # Format: 'property': (min, max, weight)
        return {
            'converged_iterations': (1000, np.infty, 0.5),
            'CO_core': (1e-5, np.infty, 0.5),
            'He_core': (1e-5, np.infty, 0.5)
        }

    def get_property(self, prop):
        if prop in kaitiaki.constants.PLOT_FILE_COLUMNS:
            return self._plot.get(prop).iloc[-1]

        if prop == 'converged_iterations':
            return len(self._plot)

    def score(self):
        return 1.0

        result = 0.0

        result += self.__score_reason()

        self._comments.append(self._reason)

        props = self.__score_properties()

        for prop, (minimum, maximum, weight) in props.items():

            val = self.get_property(prop)

            if not (minimum <= val <= maximum):
                result += weight
                self._comments.append((prop, val))

        return result / sum(prop[2] for prop in props.values())

    def fetch_comments(self):
        return self._comments
