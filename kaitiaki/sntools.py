import kaitiaki
import numpy as np
from scipy.integrate import quad
from scipy.optimize import fminbound, fmin
import warnings


from uncertainties import ufloat

# from takahe.constants
G = 6.67430e-11
c = 299792458

# Solar Units
SOLAR_RADIUS = 696340000
SOLAR_MASS = 1.989e30


# TODO: busted :(
def compute_remnant_mass(path_to_file: str):
    plot = kaitiaki.file.plot(path_to_file)
    # sneplot = kaitiaki.file.sneplot(path_to_file.replace('plot', 'sneplot'))

    mass = plot.get('M')
    radius = 10**plot.get('log(R)')

    def r(M):
        return radius[np.where(mass < M)[0][-1]] * SOLAR_RADIUS  # m

    M_star = mass[-1] * SOLAR_MASS  # kg
    U = G*M_star/radius[-1]**2  # Joules

    energy_budget = 1e44  # Joules

    # print(f"BE: {U}")

    def f(M_rem):
        M_rem *= SOLAR_MASS
        return quad(lambda M: (U-G*M/r(M)), M_rem, M_star)[0]-energy_budget

    # print(f"f(1.4): {f(1.4)}")

    best_value, _, _, _ = fminbound(f, 0.16, 1000,  # MSun
                                    maxfun=10000,
                                    full_output=1)

    # best_value = fmin(func=f, x0=1.4*SOLAR_MASS) / SOLAR_MASS

    return best_value  # solar masses


def post_supernova_parameters(path_to_file, kick_prescription, fidelity=1000):
    plot = kaitiaki.file.plot(path_to_file)
    plot2 = kaitiaki.file.plot(path_to_file.replace('plot', 'plot2'))

    M2_initial = plot2.get('M')[0]
    v_kick = kaitiaki.kicks.sample(kick_prescription, fidelity)

    v_constant = np.sqrt(G * SOLAR_MASS / SOLAR_RADIUS)

    M1 = plot.get('M')[-1]
    M2 = plot2.get('M')[-1]

    Mrem_SN = compute_remnant_mass(path_to_file)

    a = plot.get('a')[-1]

    v = v_constant*1e-3*np.sqrt((M1+M2)/a)

    # Sample the direction of the kick uniformly across the 4pi steradians
    # of the orbit.
    theta = np.arccos(2*np.random.uniform()-1)  # Between 0 and pi
    phi = 2*np.pi*np.random.uniform()           # Between 0 and 2pi

    mtilde = (M1+M2) / (Mrem_SN+M2)

    P = 1-2*mtilde+(v_kick**2)/(v**2)+2*(v_kick/v)*np.cos(theta)
    Q = 1+P/mtilde - ((v_kick*np.sin(theta)*np.cos(phi))**2)/(mtilde * v**2)
    eta = 2+P/mtilde
    newa1 = (a/(2-eta))

    # K3L:
    tosqrt = (((newa1)/215.0)**3)/(Mrem_SN + M2_initial)  # P^2 (in years)

    with warnings.catch_warnings(action="ignore"):
        newP = np.log10(365.25*np.sqrt(tosqrt))  # in days

    ecc = np.sqrt(1+(eta-2)*(Q+1))

    newP = np.ma.masked_invalid(newP)
    ecc = np.ma.masked_invalid(ecc)

    newP = ufloat(newP.mean(), newP.std())
    ecc = ufloat(ecc.mean(), ecc.std())

    return newP, ecc, Mrem_SN