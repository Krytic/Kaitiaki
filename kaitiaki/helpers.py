import numpy as np

import kaitiaki


# Note: following three functions will soon move to binary.py
def RL(q, a, lobe='L1'):
    """Computes the roche lobe radius of a star in a binary.

    Uses the approximation by P. Eggleton for L1, and P. Marchant for L2.
    - Eggleton (1983): https://ui.adsabs.harvard.edu/abs/1983ApJ...268..368E/abstract
    - Marchant et al (2016): https://ui.adsabs.harvard.edu/abs/2016A%26A...588A..50M/abstract

    Note that Eggleton takes q=m2/m1, whereas Marchant takes m1/m2. For consistency, we have elected to take q=m2/m1.


    Args:
        q (float): The mass ratio to compute for. Note that q=m2/m1.
        a (float): The binary separation
        lobe (str): Whether to compute L1 or L2 (default: `'L1'`)

    Returns:
        float: The roche lobe radius (in real units)
    """
    assert lobe in ['L1', 'L2'], f'lobe must be one of L1/L2, not {lobe}'

    def _eggleton(q):
        numerator = 0.49 * (1/q)**(2/3)
        denominator = 0.6*(1/q)**(2/3) + np.log(1+(1/q)**(1/3))
        roche_lobe = numerator / denominator

        return roche_lobe

    match lobe:
        case 'L1':
            roche_lobe = _eggleton(q)
        case 'L2':
            post_factor = (1+0.299*np.arctan(1.84*q**0.397)*q**0.520)
            roche_lobe = _eggleton(q) * post_factor

    return a * roche_lobe


def roche_lobes(donor_plotfile, accretor_plotfile):
    # R1 = 10**donor_plotfile.get('log(R)')
    # R2 = 10**accretor_plotfile.get('log(R)')

    q = accretor_plotfile.get('M') / donor_plotfile.get('M')
    a = accretor_plotfile.get('a')

    get = ['L1', 'L2']

    lobes = {
        lobe: RL(q, a, lobe=lobe) for lobe in get
    }

    return lobes


def contact_phases(donor_plotfile, accretor_plotfile):
    donor = roche_lobes(donor_plotfile, accretor_plotfile)
    accretor = roche_lobes(accretor_plotfile, donor_plotfile)

    donor_fills_L1 = donor['L1'] <= 10**donor_plotfile.get('log(R)')
    accretor_fills_L1 = accretor['L1'] <= 10**accretor_plotfile.get('log(R)')

    in_contact = np.logical_and(donor_fills_L1, accretor_fills_L1)

    return in_contact


def find_runs(x):
    """Find runs of consecutive items in an array.

    https://gist.github.com/alimanfoo/c5977e87111abe8127453b21204c1065"""

    # ensure array
    x = np.asanyarray(x)
    if x.ndim != 1:
        raise ValueError('only 1D array supported')
    n = x.shape[0]

    # handle empty array
    if n == 0:
        return np.array([]), np.array([]), np.array([])

    else:
        # find run starts
        loc_run_start = np.empty(n, dtype=bool)
        loc_run_start[0] = True
        np.not_equal(x[:-1], x[1:], out=loc_run_start[1:])
        run_starts = np.nonzero(loc_run_start)[0]

        # find run values
        run_values = x[loc_run_start]

        # find run lengths
        run_lengths = np.diff(np.append(run_starts, n))

        return run_values, run_starts, run_lengths


class Range:
    def __init__(self, bound_low, bound_high, inclusivity="[]", mutable=False):
        self._bound_low = bound_low
        self._bound_high = bound_high

        self._right_inclusive = inclusivity[1] == "]"

        self._inclusives = [None, None]

        if inclusivity[0] == "[":
            self._left_inclusive = True
            self._inclusives[0] = "["
        else:
            self._left_inclusive = False
            self._inclusives[0] = "("

        if inclusivity[1] == "]":
            self._right_inclusive = True
            self._inclusives[1] = "]"
        else:
            self._right_inclusive = False
            self._inclusives[1] = ")"

        self._mutable = mutable

    def __contains__(self, item):
        part1 = False
        part2 = False

        if self._left_inclusive:
            part1 = self._bound_low <= item
        else:
            part1 = self._bound_low < item

        if self._right_inclusive:
            part2 = self._bound_high >= item
        else:
            part2 = self._bound_high > item

        return (part1 and part2)

    def __getitem__(self, idx):
        if idx == 0:
            return self._bound_low
        elif idx == 1:
            return self._bound_high

        raise ValueError(("idx not allowed -- use 0 to access lower bound, "
                          "or 1 to access upper bound"))

    def inclusives(self):
        return (self._left_inclusive, self._right_inclusive)

    def __setitem__(self, key, val):
        if self._mutable:
            if key == 0:
                self._bound_low = val
            elif key == 1:
                self._bound_high = val
            else:
                raise ValueError(("Key was set improperly. Use 0 to modify "
                                  "the lower bound, or 1 to modify the "
                                  "upper bound."))
        else:
            raise TypeError(("This object was not initialized with "
                             "mutable methods."))

    def __str__(self):
        obj = f"{self._inclusives[0]}{self[0]}, {self[1]}{self._inclusives[1]}"
        return obj
