import numpy as np

import kaitiaki
import takahe

## TODO: This is slow as hell!

def go(outfile_loc, plotfile_loc):
    """
    Classifies a model according to the modelcheck criteria used by Jan.

    Interpretation:
        result[0] == category
        result[1] == confidence level bounded between 0 and 1.
                     LOWER is better.
    """

    outcome = 0.0

    STARS = kaitiaki.STARS.STARSController()

    out = STARS.get_last_converged_model(outfile_loc, as_obj=True)
    plot = kaitiaki.plotfile.Plotfile(plotfile_loc).last()

    """
    Fetch the necessary variables.
    TBD_* means I need to implement this (probably from plotfile)
    """


    MAXIMUM_TEMP = out.get('temp', 'Tmax')
    MASS = out.get('mass')
    AGE = out.get('age')
    RADIUS = 10e0**plot['logR'].to_numpy()[0]

    MODEL_NUM = out.get('modelnum')

    CENTRAL_DENSITY = out.get('dens', 'cntr')
    CENTRAL_HYDROGEN = out.get('H1', 'cntr')
    CENTRAL_HELIUM = out.get('He4', 'cntr')
    CENTRAL_CARBON = out.get('C12', 'cntr')
    CENTRAL_TEMP = out.get('temp', 'cntr')

    LUMINOSITY_HE = out.get('LHe')

    CO_CORE_MASS = out.get('MHe')
    HE_CORE_MASS = out.get('MH')
    ONE_CORE_MASS = MASS - CO_CORE_MASS - HE_CORE_MASS
    # H_CORE_MASS = ...

    BINARY_MASS = out.get('Mb')
    BINARY_PERIOD = out.get('P')

    # $datsne[5] = MASS

    if (CENTRAL_DENSITY > 8.0):
        # White Dwarf - Central density is high
        outcome = 4.1

    if ((MASS < 2.0) and (MAXIMUM_TEMP > 8.1)):
        # White Dwarf - total mass is less than 2Msun and max temp is highish
        outcome = 6.6

    if ((CENTRAL_HYDROGEN < 1e-5) and (CENTRAL_HELIUM < 1e-5)
                                  and (MASS < 1.4)):
        # WHITE DWARF - mass less than Chandrasekhar mass and no central hydrogen and helium
        outcome = 4.2

    if ((CENTRAL_HYDROGEN < 1e-5) and (CENTRAL_HELIUM < 1e-5)
                                  and (CO_CORE_MASS > 0.3) # actually MCO.
                                  and (HE_CORE_MASS < CO_CORE_MASS+0.1)
                                  and (MASS - HE_CORE_MASS < 2.0)
                                  and (HE_CORE_MASS < 1.3)):
        # CO WHITE DWARF similar to above but when envelope hasn't been removed, also requires that the hydrogen and helium burning shells are close together
        outcome = 4.3

    if ((CENTRAL_HYDROGEN < 1e-5) and (CENTRAL_HELIUM < 1e-5)
                                  and (CO_CORE_MASS > 0.3)
                                  and (HE_CORE_MASS < CO_CORE_MASS+0.1)):
        # WHITE DWARF - similar to above but without requiring low enevelope mass and low CO core mass
        outcome = 4.4

    if ((CENTRAL_HYDROGEN < 1e-5) and (CENTRAL_HELIUM < 1e-5)
                                  and (CO_CORE_MASS > 0.5)
                                  and (CO_CORE_MASS < 1.3)
                                  and (HE_CORE_MASS < 1.45)):
        # WHITE DWARF
        outcome = 4.5

    if ((CENTRAL_HYDROGEN < 1e-5) and (CENTRAL_HELIUM < 1e-5)
                                  and (ONE_CORE_MASS > 1.3)
                                  and (CO_CORE_MASS > 1.3)
                                  and (HE_CORE_MASS) > 1.3):
        # SAGB STAR - ONe WD
        outcome = 5.1

    if ((CENTRAL_HYDROGEN < 1e-5) and (CENTRAL_HELIUM < 1e-5)
                                  and (HE_CORE_MASS > 1.3)
                                  and (CO_CORE_MASS > 1.3)
                                  and (HE_CORE_MASS < 1.4)
                                  and (CO_CORE_MASS < 1.4)):
        # ONe WHITE DWARF
        outcome = 5.2

    if ((CENTRAL_HYDROGEN < 1e-5) and (CENTRAL_HELIUM < 1e-5)
                                  and (CENTRAL_CARBON < 0.02)
                                  and (CO_CORE_MASS <= 1.35)
                                  and (HE_CORE_MASS <= 1.35)
                                  and (MASS < 10.0)):
        # AGB star - WHITE DWARF
        outcome = 4.5

    if ((CENTRAL_HYDROGEN < 1e-5) and (CENTRAL_HELIUM > 0.9)
                                  and (LUMINOSITY_HE > 1)
                                  and (MASS < 1.4)):
        # Helium White Dwarf
        outcome = 6.1

    if ((CENTRAL_HYDROGEN < 0.001) and (CENTRAL_HELIUM > 0.9)
                                   and (CENTRAL_TEMP < 7.0)):
        # White Dwarf
        outcome = 6.2

    if ((CENTRAL_HYDROGEN < 0.001) and (CENTRAL_HELIUM > 0.9)
                                   and (MASS < 1.5)):
        # White Dwarf
        outcome = 6.3

    if ((CENTRAL_HYDROGEN < 0.001) and (CENTRAL_HELIUM < 0.9)
                                   and (MASS < 1.5)):
        # White Dwarf
        outcome = 6.4

    if ((CENTRAL_HYDROGEN < 0.001) and (CENTRAL_HELIUM < 1e-5)
                                   and (HE_CORE_MASS < 1.30)
                                   and (CO_CORE_MASS < 1.30)):
        # WHITE DWARF
        outcome = 6.5

    if ((MASS < 2.5) and (CENTRAL_HYDROGEN < 1e-5)
                     and (CENTRAL_HELIUM < 0.1)
                     and (MAXIMUM_TEMP > 8.8)
                     and (CO_CORE_MASS < 1.3)):
        # WHITE DWARF
        outcome = 6.6

    if ((MASS < 2.1) and (MASS > 1.9)
                     and (HE_CORE_MASS < 1.3)
                     and (CO_CORE_MASS > 0.2)):
        outcome = 6.7

    if ((CENTRAL_HYDROGEN < 1e-5) and (CENTRAL_HELIUM < 1e-5)
                                  and (CENTRAL_CARBON < 0.2)
                                  and (CO_CORE_MASS > 1.35)
                                  and (CENTRAL_TEMP > 8.6)
                                  and (HE_CORE_MASS > 1.35)):
        #SNe - will explode but needs further pushing
        outcome = 1.4

    if ((MASS < 3.0) and (MASS > 1.4)
                     and (HE_CORE_MASS > 1.3)
                     and (CO_CORE_MASS > 1.3)
                     and (CENTRAL_HYDROGEN < 1e-5)
                     and (CENTRAL_HELIUM < 1e-5)
                     and (MAXIMUM_TEMP > 8.8)):
        # SNE probably
        outcome = 1.3

    if ((MASS < 2.0) and (MASS > 1.4)
                     and (CENTRAL_HYDROGEN < 1e-5)
                     and (CENTRAL_HELIUM < 0.1)
                     and (MAXIMUM_TEMP > 8.8)
                     and (CO_CORE_MASS > 1.3)):
        # ODD TRANSIENTS - Ia's?
        outcome = 2.0

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
        outcome = 1.1

    if (out.nan_num() > 8):
        outcome = -2.0

    sep = takahe.helpers.compute_separation(BINARY_PERIOD,
                                            MASS,
                                            BINARY_MASS - MASS)

    if (RADIUS > sep and AGE < 1e6):
        outcome = -1.0

    if (MODEL_NUM < 1000 and BINARY_PERIOD < 1.0):
        outcome = -1.0

    if (outcome < 0.1 and AGE > 14e9):
        #Universe ain't old enough.
        outcome = 7.0


    if (0.1 < HE_CORE_MASS < 0.55) and (CENTRAL_TEMP > 7.9):
        # Helium flash
        outcome = 8.0

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
        8.0: "Helium Flash"
    }