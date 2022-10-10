from copy import deepcopy
from decimal import Decimal
from os import path
# from imageio_ffmpeg import get_ffmpeg_exe

import pandas as pd
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.cm as cm
from matplotlib.patches import Wedge, Patch
import subprocess
from tqdm import tqdm
import kaitiaki

# Use ffmpeg executable from imageio_ffmpeg to avoid installation issues
# plt.rcParams['animation.ffmpeg_path'] = get_ffmpeg_exe()

from matplotlib.animation import FuncAnimation

"""


WWWWWWWW                           WWWWWWWWIIIIIIIIIIPPPPPPPPPPPPPPPPP
W::::::W                           W::::::WI::::::::IP::::::::::::::::P
W::::::W                           W::::::WI::::::::IP::::::PPPPPP:::::P
W::::::W                           W::::::WII::::::IIPP:::::P     P:::::P
 W:::::W           WWWWW           W:::::W   I::::I    P::::P     P:::::P
  W:::::W         W:::::W         W:::::W    I::::I    P::::P     P:::::P
   W:::::W       W:::::::W       W:::::W     I::::I    P::::PPPPPP:::::P
    W:::::W     W:::::::::W     W:::::W      I::::I    P:::::::::::::PP
     W:::::W   W:::::W:::::W   W:::::W       I::::I    P::::PPPPPPPPP
      W:::::W W:::::W W:::::W W:::::W        I::::I    P::::P
       W:::::W:::::W   W:::::W:::::W         I::::I    P::::P
        W:::::::::W     W:::::::::W          I::::I    P::::P
         W:::::::W       W:::::::W         II::::::IIPP::::::PP
          W:::::W         W:::::W          I::::::::IP::::::::P
           W:::W           W:::W           I::::::::IP::::::::P
            WWW             WWW            IIIIIIIIIIPPPPPPPPPP

This package is a work in progress. Comments, structure, and methods
may change without warning whilst this warning is in place.
"""

class ModelFile:
    def __init__(self, file_location, datafile_location='data'):
        self.__models = []
        with open(file_location, 'r') as f:
            contents = f.readlines()

        first_line = contents[0].replace('-',' -').replace('E -','E-').split()
        nmesh = int(first_line[6])

        number_of_lines_per_model = 2 * nmesh + 1

        self.__length = 0

        for i in range(0, len(contents), number_of_lines_per_model):
            data = contents[i:i+number_of_lines_per_model]
            header = data[0].strip().split()
            model = data[1:nmesh+1]
            deltas = data[nmesh+1:]

            this_model = Model(*header, model, deltas)

            self.__models.append(this_model)

            self.__length += 1

    def get_bounds(self, parameter):
        global_minimum = np.infty
        global_maximum = -np.infty

        for model in self.__models:
            data = model.get(parameter)

            amin = np.amin(data)
            amax = np.amax(data)

            global_maximum = max(amax, global_maximum)
            global_minimum = min(amin, global_minimum)

        return global_minimum, global_maximum


    def request(self, model_num):
        return self.__models[model_num]

    def animate(self, parameter, fps=15, output_fname='animation', get_every=1):
        fig = plt.figure()

        amin, amax = self.get_bounds(parameter)

        def _init():
            pass
        def _animate(i):
            self.request(i).plot(parameter, fig=fig, amin=amin, amax=amax)

        ani = FuncAnimation(fig, _animate, tqdm(np.arange(0, self.__length, get_every)), init_func=_init, interval=1000 / fps, blit=False, repeat=False)

        ani.save(output_fname + ".mp4", writer="ffmpeg", extra_args=['-vcodec', 'libx264'])

class Model:
    def __init__(self, mass, dt,
                 age, period, total_mass,
                 energy_generation, nmesh, nmod,
                 nmod_start, nstar, h_pressure, he_pressure,
                 model_lines,
                 delta_lines):

        self.__mass = float(mass)
        self.__dt = float(dt)
        self.__age = float(age)
        self.__period = float(period)
        self.__total_mass = float(total_mass)
        self.__energy_generation = float(energy_generation)
        self.__nmesh = int(nmesh)
        self.__nmod = int(nmod)
        self.__nmod_start = int(nmod_start)
        self.__nstar = int(nstar)
        self.__h_pressure = float(h_pressure)
        self.__he_pressure = float(he_pressure)

        cleaned_models = []

        for i in range(self.__nmesh):
            this_model = model_lines[i].replace('-',' -').replace('E -','E-').split()
            for j in range(len(this_model)):
                this_model[j] = float(this_model[j])

            cleaned_models.append(this_model)


            # delta_lines[i] = [delta_lines[i][j:j+15].strip() for j in range(0, len(delta_lines[i]), 15)]
            # if delta_lines[i][-1] in ["\n", '']:
            #     delta_lines[i] = delta_lines[i][:-1]
            # for j in range(len(delta_lines[i])):
            #     delta_lines[i][j] = float(delta_lines[i][j])

        df = pd.DataFrame(cleaned_models, columns=['logf', 'logT', 'X_O',
                                                'logm', 'X_H', 'dQ/dK',
                                                'logr', 'L', 'X_He',
                                                'X_C', "X_Ne", 'X_N',
                                                'H_orb', 'H_spin', 'X_3He'])

        # df_deltas = pd.DataFrame(delta_lines, columns=['logf', 'logT', 'X_O',
        #                                                'logm', 'X_H', 'dQ/dK',
        #                                                'logr', 'L', 'X_He',
        #                                                'X_C', "X_Ne", 'X_N',
        #                                                'H_orb', 'H_spin', 'X_3He'])

        # These four parameters are log_e whereas other log parameters
        # are log_10. Change their base:
        df['logf'] /= np.log(10)
        df['logT'] /= np.log(10)
        df['logm'] /= np.log(10)
        df['logr'] /= np.log(10)

        # df_deltas['logf'] /= np.log(10)
        # df_deltas['logT'] /= np.log(10)
        # df_deltas['logm'] /= np.log(10)
        # df_deltas['logr'] /= np.log(10)

        self.__data = df
        # self.__deltas = df_deltas

    def get(self, column):
        return self.__data[column].to_numpy()

    def plot(self, parameter, fig=None, amin=0, amax=20):
        data = self.get(parameter)
        mass = 10**self.get('logm')

        ax_min = np.amin(mass)
        ax_max = np.amax(mass)

        abs_max = max(np.abs(ax_max), np.abs(ax_min))

        if fig == None:
            fig = plt.figure()
        else:
            plt.clf()
        ax = plt.gca()

        norm = matplotlib.colors.Normalize(vmin=amin,
                                           vmax=amax)
        cmap = cm.plasma

        m = cm.ScalarMappable(norm=norm, cmap=cmap)

        for i in range(0, self.__nmesh):
            r1 = i * mass[i] / self.__nmesh
            r2 = (i+1) * mass[i] / self.__nmesh
            w = (r2-r1)

            cval = data[i]

            wedge = Wedge((0, 0), r2, theta1=0,
                                      theta2=360,
                                      width=w,
                                      color=m.to_rgba(cval))

            ax.add_artist(wedge)

        plt.xlim(-(abs_max+abs_max / 10), abs_max+abs_max / 10)
        plt.ylim(-(abs_max+abs_max / 10), abs_max+abs_max / 10)

        plt.xlabel(r"Mass [$M_\odot$]")
        plt.ylabel(r"Mass [$M_\odot$]")

        cbar = fig.colorbar(m)

        cbar.set_label(parameter)

        plt.title(f"{parameter} of star, t={self.__age/1e6} Myr")

        for spine in ['top', 'left', 'bottom', 'right']:
            ax.spines[spine].set_visible(True)

        plt.grid(False)