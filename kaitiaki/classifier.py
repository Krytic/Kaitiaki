import numpy as np

import kaitiaki
import takahe

def go(outfile_loc, plotfile_loc):
    """
    Classifies a model according to the modelcheck criteria used by Jan.

    Interpretation:
        result[0] == category
        result[1] == confidence level bounded between 0 and 1.
                     LOWER is better.
    """

    outcome = 0.0

    out = kaitiaki.out.outfile(outfile_loc).last()
    plot = kaitiaki.plotfile.Plotfile(plotfile_loc).last()

    """
    TBD_* means I need to implement this (probably from plotfile)
    """

    TBD_ONE_CORE_MASS = 0.0

    if out.get('dens', 'cntr') > 8.0:
        # White Dwarf - Central density is high
        outcome = 4.1

    if (out.get('mass') < 2.0) and (out.get('temp', 'Tmax') > 8.1):
        # White Dwarf - total mass is less than 2Msun and max temp is highish
        outcome = 6.6

    if ((out.get('H1', 'cntr') < 1e-5) and (out.get('He4', 'cntr') < 1e-5)
                                       and (out.get('mass'))):
        # WHITE DWARF - mass less than Chandrasekhar mass and no central hydrogen and helium
        outcome = 4.2

    if ((out.get('H1', 'cntr') < 1e-5) and (out.get('He4', 'cntr') < 1e-5)
                                       and (out.get("MHe") > 0.3) # actually MCO.
                                       and (out.get('MHe') < out.get('mass')+0.1)
                                       and (out.get('MH')-out.get('MHe') < 2.0)
                                       and (out.get('MHe') < 1.3)):
        # CO WHITE DWARF similar to above but when envelope hasn't been removed, also requires that the hydrogen and helium burning shells are close together
        outcome = 4.3

    if ((out.get('H1', 'cntr') < 1e-5) and (out.get('He4', 'cntr') < 1e-5)
                                       and (out.get("MHe") > 0.3) # actually MCO.
                                       and (out.get('MHe') < out.get('mass')+0.1)):
        # WHITE DWARF - similar to above but without requiring low enevelope mass and low CO core mass
        outcome = 4.4

    if ((out.get('H1', 'cntr') < 1e-5) and (out.get('He4', 'cntr') < 1e-5)
                                       and (out.get("MHe") > 0.5) # actually MCO.
                                       and (out.get("MHe") < 1.3) # actually MCO.
                                       and (out.get('mass') < 1.45)):
        # WHITE DWARF
        outcome = 4.5

    if ((out.get('H1', 'cntr') < 1e-5) and (out.get('He4', 'cntr') < 1e-5)
                                       and (TBD_ONE_CORE_MASS > 1.3)
                                       and (out.get("MHe") > 1.3)
                                       and (out.get("MH")) > 1.3):
        # SAGB STAR - ONe WD
        outcome = 5.1

    if ((out.get('H1', 'cntr') < 1e-5) and (out.get('He4', 'cntr') < 1e-5)
                                       and (out.get("MH") > 1.3)
                                       and (out.get("MHe") > 1.3)
                                       and (out.get("MH") < 1.4)
                                       and (out.get("MHe") < 1.4)):
        # ONe WHITE DWARF
        outcome = 5.2

    if ((out.get('H1', 'cntr') < 1e-5) and (out.get('He4', 'cntr') < 1e-5)
                                       and (out.get("C12", 'cntr') < 0.02)
                                       and (out.get("MHe") <= 1.35)
                                       and (out.get("MH") <= 1.35)
                                       and (out.get('mass') < 10.0)):
        # AGB star - WHITE DWARF
        outcome = 4.5

    if(out.get('H1', 'cntr') < 0.001 and out.get('He4', 'cntr') > 0.9
                                     and out.get('temp', 'cntr') < 7.0):
        # White Dwarf
        outcome = 6.2

    if(out.get('H1', 'cntr') < 0.001 and out.get('He4', 'cntr') > 0.9
                                     and out.get('mass') < 1.5):
        # White Dwarf
        outcome = 6.3

    if(out.get('H1', 'cntr') < 0.001 and out.get('He4', 'cntr') < 0.9
                                     and out.get('mass') < 1.5):
        # White Dwarf
        outcome = 6.4

    if(out.get('H1', 'cntr') < 0.001 and out.get('He4', 'cntr') < 1e-5
                                     and out.get('MH') < 1.30
                                     and out.get('MHe') < 1.30):
        # WHITE DWARF
        outcome = 6.5

    if(out.get('mass') < 2.5 and out.get('H1', 'cntr') < 1e-5
                          and out.get('He4', 'cntr') < 0.1
                          and out.get('temp', 'Tmax') > 8.8
                          and out.get('MHe') < 1.3):
        # WHITE DWARF
        outcome = 6.6

    if(out.get('mass') < 2.1 and out.get('mass') > 1.9
                          and out.get('MH') < 1.3
                          and out.get('MHe') > 0.2):
        outcome = 6.7

    if(out.get('H1', 'cntr') < 1e-5 and out.get('He4', 'cntr') < 1e-5
                                    and out.get('C12', 'cntr') < 0.2
                                    and out.get('MHe') > 1.35
                                    and out.get('temp', 'cntr') > 8.6
                                    and out.get('MH') > 1.35):
        #SNe - will explode but needs further pushing
        outcome = 1.4

    if(out.get('mass') < 3.0 and out.get('mass') > 1.4
                          and out.get('MH') > 1.3
                          and out.get('MHe') > 1.3
                          and out.get('H1', 'cntr') < 1e-5
                          and out.get('He4', 'cntr') < 1e-5
                          and out.get('temp', 'Tmax') > 8.8):
        # SNE probably
        outcome = 1.3

    if(out.get('mass') < 2.0 and out.get('mass') > 1.4
                          and out.get('H1', 'cntr') < 1e-5
                          and out.get('He4', 'cntr') < 0.1
                          and out.get('temp', 'Tmax') > 8.8
                          and out.get('MHe') > 1.3):
        # ODD TRANSIENTS Ia's?
        outcome = 2.0

    if(out.get('H1', 'cntr') < 1e-5 and out.get('He4', 'cntr') < 1e-5
                                    and out.get('C12', 'cntr') < 0.1
                                    and out.get('MHe') > 1.35
                                    and out.get('temp', 'cntr') > 8.8
                                    and out.get('MH') > 1.35):
        #SNe - not quite close enough to end but okay
        outcome=1.2

    if(out.get('mass') > 2.0 and out.get('MH') > 1.4
                          and out.get('MHe') > 1.3
                          and TBD_ONE_CORE_MASS > 0.1):
        #SNe - all big enough cores so definitely okay
        outcome = 1.1

    if(out.nan_num() > 8):
        outcome = 0.0

    sep = takahe.helpers.compute_separation(out.get('P'),
                                            out.get('mass'),
                                            out.get('Mb')-out.get('mass'))

    if(10e0**plot['logR'].to_numpy()[0] > sep and out.get('age') < 1e6):
        outcome = -1.0;

    if(out.get('modelnum') < 1000 and out.get('P') < 1.0):
        outcome = -1.0

    if(outcome < 0.1 and out.get('age') > 14e9):
        #Universe ain't old enough.
        outcome = 7.0


    if (0.1 < out.get('MHe') < 0.55) and (out.get('temp', 'cntr') > 7.9):
        # Helium flash
        outcome = 8.0

    return outcome

def to_str(code):
    return strings()[code]

def strings():
    return {
        4.1: "White Dwarf",
        6.6: "White Dwarf",
        4.2: 'White Dwarf',
        4.3: 'CO White Dwarf',
        4.4: 'White Dwarf',
        4.5: 'AGB Star (White Dwarf)',
        5.1: 'sAGB Star (ONe White Dwarf)',
        5.2: 'ONe White Dwarf',
        6.2: 'White Dwarf',
        6.3: 'White Dwarf',
        6.4: 'White Dwarf',
        6.5: 'White Dwarf',
        6.6: 'White Dwarf',
        6.7: 'White Dwarf',
        1.4: 'SNe (requires a little more push)',
        1.3: 'SNe',
        2.0: 'Type Ia?',
        1.2: 'SNe',
        1.1: 'SNe',
        0.0: 'Unclassifiable',
       -1.0: 'Too young',
        7.0: 'Too old',
        8.0: "Helium Flash"
    }