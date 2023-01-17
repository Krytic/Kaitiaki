import numpy as np

import kaitiaki
import tabulate

@np.vectorize
def compute_separation(P, M, m):
    """Computes the separation of a BSS

    Uses Kepler's third law to compute the separation, in days.

    Decorators:
        np.vectorize

    Arguments:
        P {float} -- The period of the binary star (in days)
        M {float} -- The mass of the primary star (in solar masses)
        m {float} -- The mass of the secondary star (in solar masses)

    Returns:
        {float} -- The orbital separation in Solar Radii
    """

    G = 6.67430e-11
    SOLAR_RADIUS      = 696340000
    SOLAR_MASS        = 1.989e30

    a = ((P * 60 * 60 * 24)**2 * (G * (M+m) * SOLAR_MASS) / (4 * np.pi **2))**(1/3) / SOLAR_RADIUS

    return a

## TODO: This is slow as hell!

def go(outfile_loc: str='out',
       plotfile_loc: str='plot',
       as_string: bool=False,
       explain: bool=False,
       has_he_flash: bool=False):
    """
    Classifies a model according to the modelcheck criteria used by Jan.

    Interpretation:
        result[0] == category
        result[1] == confidence level bounded between 0 and 1.
                     LOWER is better.
    """

    STARS = kaitiaki.STARS.STARSController()

    if has_he_flash and not outfile_loc.endswith('.posthf'):
        outfile_loc = outfile_loc + ".posthf"

    out = STARS.get_last_converged_model(outfile_loc, as_obj=True)

    plot = kaitiaki.plotfile.Plotfile(plotfile_loc)

    """
    Fetch the necessary variables.
    TBD_* means I need to implement this (probably from plotfile)
    """

    MAXIMUM_TEMP = out.get('temp', 'Tmax')
    MASS = out.get('mass')
    AGE = out.get('age')
    RADIUS = 10e0**plot.access()['logR'].to_numpy()[-1]
    LHE = plot.access()['LHe'].to_numpy()

    MODEL_NUM = out.get('modelnum')

    CENTRAL_DENSITY = out.get('dens', 'cntr')
    CENTRAL_HYDROGEN = out.get('H1', 'cntr')
    CENTRAL_HELIUM = out.get('He4', 'cntr')
    CENTRAL_CARBON = out.get('C12', 'cntr')
    CENTRAL_TEMP = out.get('temp', 'cntr')

    LUMINOSITY_HE = out.get('LHe')

    H_CORE_MASS = out.get('MH')
    CO_CORE_MASS = out.get('MC')
    HE_CORE_MASS = out.get('MHe')
    ONE_CORE_MASS = MASS - CO_CORE_MASS - HE_CORE_MASS - H_CORE_MASS
    # H_CORE_MASS = ...

    BINARY_MASS = out.get('Mb')
    BINARY_PERIOD = out.get('P')

    sep = compute_separation(BINARY_PERIOD, MASS, BINARY_MASS - MASS)

    exclude = ['as_string', 'explain', 'STARS', 'out', 'plot', 'exclude', 'LHE']
    explainer = {k: v for k, v in locals().items() if k not in exclude}

    plausible_outcomes = [0.0]

    if (CENTRAL_DENSITY > 8.0):
        # White Dwarf - Central density is high
        plausible_outcomes.append(4.1)

    if ((MASS < 2.0) and (MAXIMUM_TEMP > 8.1)):
        # White Dwarf - total mass is less than 2Msun and max temp is highish
        plausible_outcomes.append(6.6)

    if ((CENTRAL_HYDROGEN < 1e-5) and (CENTRAL_HELIUM < 1e-5)
                                  and (MASS < 1.4)):
        # WHITE DWARF - mass less than Chandrasekhar mass and no central hydrogen and helium
        plausible_outcomes.append(4.2)

    if ((CENTRAL_HYDROGEN < 1e-5) and (CENTRAL_HELIUM < 1e-5)
                                  and (CO_CORE_MASS > 0.3) # actually MCO.
                                  and (HE_CORE_MASS < CO_CORE_MASS+0.1)
                                  and (MASS - HE_CORE_MASS < 2.0)
                                  and (HE_CORE_MASS < 1.3)):
        # CO WHITE DWARF similar to above but when envelope hasn't been removed, also requires that the hydrogen and helium burning shells are close together
        plausible_outcomes.append(4.3)

    if ((CENTRAL_HYDROGEN < 1e-5) and (CENTRAL_HELIUM < 1e-5)
                                  and (CO_CORE_MASS > 0.3)
                                  and (HE_CORE_MASS < CO_CORE_MASS+0.1)):
        # WHITE DWARF - similar to above but without requiring low enevelope mass and low CO core mass
        plausible_outcomes.append(4.4)

    if ((CENTRAL_HYDROGEN < 1e-5) and (CENTRAL_HELIUM < 1e-5)
                                  and (CO_CORE_MASS > 0.5)
                                  and (CO_CORE_MASS < 1.3)
                                  and (HE_CORE_MASS < 1.45)):
        # WHITE DWARF
        plausible_outcomes.append(4.5)

    if ((CENTRAL_HYDROGEN < 1e-5) and (CENTRAL_HELIUM < 1e-5)
                                  and (ONE_CORE_MASS > 1.3)
                                  and (CO_CORE_MASS > 1.3)
                                  and (HE_CORE_MASS) > 1.3):
        # SAGB STAR - ONe WD
        plausible_outcomes.append(5.1)

    if ((CENTRAL_HYDROGEN < 1e-5) and (CENTRAL_HELIUM < 1e-5)
                                  and (HE_CORE_MASS > 1.3)
                                  and (CO_CORE_MASS > 1.3)
                                  and (HE_CORE_MASS < 1.4)
                                  and (CO_CORE_MASS < 1.4)):
        # ONe WHITE DWARF
        plausible_outcomes.append(5.2)

    if ((CENTRAL_HYDROGEN < 1e-5) and (CENTRAL_HELIUM < 1e-5)
                                  and (CENTRAL_CARBON < 0.02)
                                  and (CO_CORE_MASS <= 1.35)
                                  and (HE_CORE_MASS <= 1.35)
                                  and (MASS < 10.0)):
        # AGB star - WHITE DWARF
        plausible_outcomes.append(4.5)

    if ((CENTRAL_HYDROGEN < 1e-5) and (CENTRAL_HELIUM > 0.9)
                                  and (LUMINOSITY_HE > 1)
                                  and (MASS < 1.4)):
        # Helium White Dwarf
        plausible_outcomes.append(6.1)

    if ((CENTRAL_HYDROGEN < 0.001) and (CENTRAL_HELIUM > 0.9)
                                   and (CENTRAL_TEMP < 7.0)):
        # White Dwarf
        plausible_outcomes.append(6.2)

    if ((CENTRAL_HYDROGEN < 0.001) and (CENTRAL_HELIUM > 0.9)
                                   and (MASS < 1.5)):
        # White Dwarf
        plausible_outcomes.append(6.3)

    if ((CENTRAL_HYDROGEN < 0.001) and (CENTRAL_HELIUM < 0.9)
                                   and (MASS < 1.5)):
        # White Dwarf
        plausible_outcomes.append(6.4)

    if ((CENTRAL_HYDROGEN < 0.001) and (CENTRAL_HELIUM < 1e-5)
                                   and (HE_CORE_MASS < 1.30)
                                   and (CO_CORE_MASS < 1.30)):
        # WHITE DWARF
        plausible_outcomes.append(6.5)

    if ((MASS < 2.5) and (CENTRAL_HYDROGEN < 1e-5)
                     and (CENTRAL_HELIUM < 0.1)
                     and (MAXIMUM_TEMP > 8.8)
                     and (CO_CORE_MASS < 1.3)):
        # WHITE DWARF
        plausible_outcomes.append(6.6)

    if ((MASS < 2.1) and (MASS > 1.9)
                     and (HE_CORE_MASS < 1.3)
                     and (CO_CORE_MASS > 0.2)):
        plausible_outcomes.append(6.7)

    if ((CENTRAL_HYDROGEN < 1e-5) and (CENTRAL_HELIUM < 1e-5)
                                  and (CENTRAL_CARBON < 0.2)
                                  and (CO_CORE_MASS > 1.35)
                                  and (CENTRAL_TEMP > 8.6)
                                  and (HE_CORE_MASS > 1.35)):
        #SNe - will explode but needs further pushing
        plausible_outcomes.append(1.4)

    if ((MASS < 3.0) and (MASS > 1.4)
                     and (HE_CORE_MASS > 1.3)
                     and (CO_CORE_MASS > 1.3)
                     and (CENTRAL_HYDROGEN < 1e-5)
                     and (CENTRAL_HELIUM < 1e-5)
                     and (MAXIMUM_TEMP > 8.8)):
        # SNE probably
        plausible_outcomes.append(1.3)

    if ((MASS < 2.0) and (MASS > 1.4)
                     and (CENTRAL_HYDROGEN < 1e-5)
                     and (CENTRAL_HELIUM < 0.1)
                     and (MAXIMUM_TEMP > 8.8)
                     and (CO_CORE_MASS > 1.3)):
        # ODD TRANSIENTS - Ia's?
        plausible_outcomes.append(2.0)

    if ((CENTRAL_HYDROGEN < 1e-5) and (CENTRAL_HELIUM < 1e-5)
                                  and (CENTRAL_CARBON < 0.1)
                                  and (CO_CORE_MASS > 1.35)
                                  and (CENTRAL_TEMP > 8.8)
                                  and (HE_CORE_MASS > 1.35)):
        #SNe - not quite close enough to end but okay
        outcome=1.2

    if ((MASS > 2.0) and (HE_CORE_MASS > 1.4)
                     and (CO_CORE_MASS > 1.3)
                     and (ONE_CORE_MASS > 0.1)):
        #SNe - all big enough cores so definitely okay
        plausible_outcomes.append(1.1)

    if (out.nan_num() > 8):
        plausible_outcomes.append(-2.0)

    if (RADIUS > sep and AGE < 1e6):
        plausible_outcomes.append(-1.0)

    if (MODEL_NUM < 1000 and BINARY_PERIOD < 1.0):
        plausible_outcomes.append(-1.0)

    if (plausible_outcomes[-1] < 0.1 and AGE > 14e9):
        #Universe ain't old enough.
        plausible_outcomes.append(7.0)

    if (0.1 < HE_CORE_MASS < 0.55) and (CENTRAL_TEMP > 7.9):
        # Helium flash
        plausible_outcomes.append(8.0)

        NMOD_of_first_peak_Helium = (LHE > 2).argmin()
        LHe_post_HF = LHE[NMOD_of_first_peak_Helium+1:]

        if (LHe_post_HF.max()/LHE[LHE > 2].min() > 2):
            # Ratio of peak LHe to LHe of Helium Flash is > 2
            # so a thermal pulse is probably occuring...
            plausible_outcomes.append(9.0)

    outcome = plausible_outcomes[-1]

    if explain:
        explain_result(explainer, plausible_outcomes, outcome)

    if as_string: outcome = to_str(outcome)
    return outcome

def to_str(code):
    return strings()[code]

def strings():
    return {
       -2.0: 'Too Many NaNs',
       -1.0: 'Too young',
        0.0: 'Unclassifiable',
        1.1: 'SNe',
        1.2: 'SNe',
        1.3: 'SNe',
        1.4: 'SNe (requires a little more push)',
        2.0: 'Type Ia?',
        4.1: "White Dwarf",
        4.2: 'White Dwarf',
        4.3: 'CO White Dwarf',
        4.4: 'White Dwarf',
        4.5: 'AGB Star (White Dwarf)',
        5.1: 'sAGB Star (ONe White Dwarf)',
        5.2: 'ONe White Dwarf',
        6.1: "He White Dwarf",
        6.2: 'White Dwarf',
        6.3: 'White Dwarf',
        6.4: 'White Dwarf',
        6.5: 'White Dwarf',
        6.6: 'White Dwarf',
        6.7: 'White Dwarf',
        7.0: 'Too old',
        8.0: "Helium Flash",
        9.0: "Thermally Pulsing"
    }

def explain_result(explainer, plausible_outcomes, outcome):
    max_line_length = 78

    header = ["" for _ in range(4)]
    header[0] = " Model metadata "
    header[1] = " Inferred from out file "
    header[2] = " Inferred from plot file "
    header[3] = " Results "

    exclude_from_table = ['RADIUS', 'outfile_loc', 'plotfile_loc']

    table = tabulate.tabulate([[k, v] for k, v in explainer.items() if k not in exclude_from_table], ['Parameter', 'Value'], tablefmt='fancy_grid')

    for i in range(4):
        header[i] = "=" * int((max_line_length - len(header[i])) / 2) + header[i] + "=" * int((max_line_length - len(header[i])) / 2)

    outstr = f"""
==============================================================================
|  _  __     _ _   _       _    _    _____ _               _  __ _           |
| | |/ /    (_) | (_)     | |  (_)  / ____| |             (_)/ _(_)          |
| | ' / __ _ _| |_ _  __ _| | ___  | |    | | __ _ ___ ___ _| |_ _  ___ _ __ |
| |  < / _` | | __| |/ _` | |/ / | | |    | |/ _` / __/ __| |  _| |/ _ \ '__||
| | . \ (_| | | |_| | (_| |   <| | | |____| | (_| \__ \__ \ | | | |  __/ |   |
| |_|\_\__,_|_|\__|_|\__,_|_|\_\_|  \_____|_|\__,_|___/___/_|_| |_|\___|_|   |
==============================================================================

{header[0]}
out file: {explainer['outfile_loc']}
plot file: {explainer['plotfile_loc']}

{header[1]}
{table}

{header[2]}
RADIUS = {explainer['RADIUS']} Solar Radii

{header[3]}
I classify as {outcome} ({to_str(outcome)})"""

    multiple_outstr = ""

    if len(plausible_outcomes) > 2:
        outstr += "\n"
        multiple_outstr = "Degenerate result. Alternatives:\n"
        subtr = 0
        for i, oc in enumerate(plausible_outcomes):
            if oc == 0.0:
                subtr = 1
                continue
            multiple_outstr += f"    {i+1-subtr}. {oc} ({to_str(oc)})\n"

    kaitiaki.debug('info', outstr)
    if multiple_outstr != "":
        kaitiaki.debug('warning', multiple_outstr.rstrip())