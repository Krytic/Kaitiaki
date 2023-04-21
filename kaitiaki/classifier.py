import os

import numpy as np
import matplotlib.pyplot as plt
import tabulate

import kaitiaki

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

def diagnosis_plot(files_loc='.', time_axis='age', fname=None, fext='', temp=None, lum=None):
    """Generates a diagnosis plot for a model.

    This is a 3x3 plot, with the following axes (left to right, top to bottom):
        - Mass
        - Kippenhahn Diagram (Star 1)
        - log(L)

        - Radius
        - HR Diagram
        - log(T)

        - Helium Luminosity
        - Kippenhahn Diagram (Star 2) OR Chemical Composition Plot (Star 1)
        - Mass Loss Rate

    The secondary, if it exists, is shown in dashed lines. If the secondary does not exist, the Kippenhahn for Star 2 will instead be the Chemical Composition Plot for Star 1.

    Args:
        files_loc (str): the location of the files to load (default: `'.'`)
        time_axis (str): whether the x-axis should be age or timestep (default: `'age'`)
        fname (str): The filename to save the figure as. If None, no figure is saved. (default: `None`)
        fext (str): The file extension of the files to be plotted, with no leading period (e.g. 'bak' for 'plot.bak'). If an empty string is passed, it is assumed the files have no extension. (default: '')
    """

    # time_axis should accept modelnum not timestep -- or the kippenhahn
    # plotter should change.
    assert time_axis in ['age', 'timestep'], ("time_axis must be "
                                              "age/timestep, not "
                                             f"{time_axis}")

    if fext != '':
        fext = f'.{fext}'

    # Utility function to share X axis between plots after creation.
    def conjoin_axes(*ax):
        ax = ax[0]
        ax[0].get_shared_x_axes().join(ax[0], *ax)
        for axi in ax[:-1]: axi.set_xticklabels([])

    fig, axes = plt.subplots(3,3,figsize=(15,15))

    period_axis = axes[1,0].twinx()

    conjoin_axes(axes[:,0])
    conjoin_axes(axes[:,2])

    ls = {'': '-', '2': '--'}

    if os.path.exists(f"{files_loc}/plot2{fext}"):
        stars = ['', '2']
        secondary_exists = True
    else:
        stars = ['']
        secondary_exists = False

    for itr, star in enumerate(stars):
        pf = kaitiaki.file.plot(f'{files_loc}/plot{star}{fext}')

        # Column 1
        pf.plot(time_axis, 'M', ax=axes[0,0], ls=ls[star])

        pf.plot(time_axis, 'log(R)',
                           transform=kaitiaki.utils.transforms.unlog_y,
                           ax=axes[1,0],
                           ls=ls[star])

        pf.plot(time_axis, 'L_(He)', ax=axes[2,0], ls=ls[star])

        pf.plot(time_axis, 'P_bin', ax=period_axis, ls='dotted', c='magenta')

        pf.plot(time_axis, 'a', ax=axes[1,0], ls='--', c='k')

        # Column 2
        axis = 'modelnum' if time_axis == 'timestep' else time_axis
        if time_axis == 'age': axis = 'collapsetime'

        ax = {'': axes[0,1], '2': axes[2,1]}
        pf.kippenhahn_diagram(ax=ax[star], x_axis=axis, legend=False)
        pf.hr_diagram(ax=axes[1,1], ls=ls[star])

        if temp != None:
            axes[1,1].errorbar(temp[itr].n, lum[itr].n, xerr=temp[itr].s, yerr=lum[itr].s)

        if not secondary_exists:
            cs = 'rgbkm'

            for i, species in enumerate('XYCNO'):
                pf.plot(time_axis, species,
                                   ax=axes[2,1],
                                   c=cs[i],
                                   label=species)

            axes[2,1].legend()

        # Column 3
        pf.plot(time_axis, 'log(L)', ax=axes[0,2], ls=ls[star])
        pf.plot(time_axis, 'log(T)', ax=axes[1,2], ls=ls[star])
        pf.plot(time_axis, 'DM1W', ax=axes[2,2], ls=ls[star])

        axes[0,0].set_title("mass")
        axes[1,0].set_title("radius")
        axes[2,0].set_title("Helium luminosity")

        axes[0,2].set_title('log(L)')
        axes[1,2].set_title('log(T)')
        axes[2,2].set_title('mass loss rate')

    if fname is not None:
        if fname.split(".")[-1] != 'png':
            fname = fname + '.png'
        plt.savefig(fname, dpi=600, facecolor='white', bbox_inches='tight')

    explainer1, plausible_outcomes1, outcome1 = go(f'{files_loc}/out',
                                                   f'{files_loc}/plot',
                                                   detailed_return=True)

    explainer2, plausible_outcomes2, outcome2 = go(f'{files_loc}/out2',
                                                   f'{files_loc}/plot2',
                                                   detailed_return=True)

    dual_explain((explainer1, plausible_outcomes1, outcome1),
                 (explainer2, plausible_outcomes2, outcome2))

def go(outfile_loc: str='out', plotfile_loc: str='plot', as_string: bool=False, explain: bool=False, detailed_return: bool=False, has_he_flash: bool=False):
    """Classifies a model according to the modelcheck criteria used by Jan.

    Performs a rudimentary classification based on manually inspected values. It is my wish to replace this one day with a CNN.

    Args:
        outfile_loc (str): The outfile to load. (default: `'out'`)
        plotfile_loc (str): The plotfile to load (default: `'plot'`)
        as_string (bool): Whether to return the result as a string like 'Too Old' or the numeric code like 8.0 (default: `False`)
        explain (bool): Whether to use kaitiaki's debugger to print out an explanation as to how the result was obtained. (default: `False`)
        detailed_return (bool): If True, instead of returning anything else, this will return all the arguments that can be passed to :code:`explain`. Note that this is a more complicated way of going explain=True, and is only really useful to the :code:`dual_explain()` function. (default: `False`)
        has_he_flash (bool): Whether or not the model has a helium flash you want it to consider. (default: `False`)

    Returns:
        mixed: If :code:`detailed_return` is True, a 3-tuple containing the arguments to :code:`explain()`. Else, if :code:`as_string` is True, the string representing the classification (e.g. "SNe" or "White Dwarf"). Else, the numeric code representing the result -- the whole part being the type (e.g. 1 is SNe) and the decimal part being the confidence, where a LOWER confidence is better. e.g. 1.1 means SNe with a better confidence than 1.3.
    """

    if has_he_flash and not outfile_loc.endswith('.posthf'):
        outfile_loc = outfile_loc + ".posthf"

    out = kaitiaki.file.out(outfile_loc)

    out = out.last()

    plot = kaitiaki.file.plot(plotfile_loc)

    # Fetch the necessary variables.

    MAXIMUM_TEMP = out.get('temp', 'Tmax')
    MASS = out.get('mass')
    AGE = out.get('age')
    RADIUS = 10e0**plot.access()['log(R)'].to_numpy()[-1]
    LHE = plot.access()['L_(He)'].to_numpy()

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

    exclude = ['as_string', 'explain', 'out', 'plot', 'exclude', 'LHE', 'detailed_return']
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

    if detailed_return:
        return explainer, plausible_outcomes, outcome

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

def explain_result(explainer, plausible_outcomes, outcome, speak=True):
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

    if speak:
        kaitiaki.debug('info', outstr)
        if multiple_outstr != "":
            kaitiaki.debug('warning', multiple_outstr.rstrip())
    else:
        return outstr, multiple_outstr.rstrip()

def dual_explain(obj1, obj2):
    outstr1, multiple_outstr1 = explain_result(*obj1, speak=False)
    outstr2, multiple_outstr2 = explain_result(*obj2, speak=False)

    dual_outstr = []

    outstr2 = outstr2.split("\n")
    outstr1 = outstr1.split("\n")

    while len(outstr2) < len(outstr1):
        outstr2.append("")

    while len(outstr1) < len(outstr2):
        outstr1.append("")

    for i, line in enumerate(outstr1):
        dual_outstr.append(f"{line.rstrip().ljust(78)} | {outstr2[i].rstrip().ljust(78)}")

    dual_outstr = "\n".join(dual_outstr)
    kaitiaki.debug('info', dual_outstr)

    if multiple_outstr1 != "":
        kaitiaki.debug('warning', "Star 1:\n" + multiple_outstr1)

    if multiple_outstr2 != "":
        kaitiaki.debug('warning', "Star 2:\n" + multiple_outstr2)
