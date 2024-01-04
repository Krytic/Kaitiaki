import kaitiaki

import numpy as np
import matplotlib
import matplotlib.pyplot as plt

from mpl_toolkits.axes_grid1 import make_axes_locatable


def radial_plot(values,
                NM2=199,
                labelname='',
                ax=None,
                cax=None,
                normal_value=None):

    rlist = np.arange(0, NM2, 1)

    vmax = normal_value
    vmin = None

    if normal_value is not None:
        vmin = 0

    if values.ndim == 1:
        # 1d array
        thetalist = np.radians(np.arange(0, 361, 1))
        rmesh, thetamesh = np.meshgrid(rlist, thetalist)
        results = np.column_stack((values for _ in range(len(thetalist)))).T
    else:
        # not a 1d array
        thetalist = np.radians(np.linspace(0, 360, len(values)))
        rmesh, thetamesh = np.meshgrid(rlist, thetalist)
        results = values

    if ax is None:
        _, ax = plt.subplots(dpi=120,
                             subplot_kw=dict(projection='polar'))
    else:
        if not isinstance(ax, matplotlib.projections.polar.PolarAxes):
            raise TypeError(("You need to pass radial_plot a polar axis "
                             "for this to work."))

    ax.set_theta_direction(-1)
    ax.set_theta_offset(np.pi / 2.0)

    im = ax.contourf(thetamesh, rmesh, results, NM2+1,
                     cmap='inferno',
                     vmin=vmin,
                     vmax=vmax)

    ax.tick_params(
        axis='x',           # changes apply to the x-axis
        which='both',       # both major and minor ticks are affected
        bottom=False,       # ticks along the bottom edge are off
        top=False,          # ticks along the top edge are off
        labelbottom=False)  # labels along the bottom edge are off

    ax.xaxis.grid(False)

    fig = ax.get_figure()

    if cax is None:
        cbar = fig.colorbar(im, orientation='horizontal', pad=0.2)
    else:
        cbar = fig.colorbar(im, cax=cax, pad=0.2)

    if vmax is not None:
        cbar.ax.set_ylim(vmin, vmax)

    cbar.set_label(labelname)


def radial_and_rectilinear(radial_values, NM2=199, labelname=''):
    fig = plt.figure()
    radial_axis = fig.add_subplot(121, projection='polar')
    rectilinear_axis = fig.add_subplot(122)

    cax = fig.add_axes([0.05, 0.05, 0.03, 0.8])

    radial_plot(radial_values, NM2, labelname, radial_axis, cax=cax)

    cax.yaxis.tick_left()
    cax.yaxis.set_label_position('left')

    return radial_axis, rectilinear_axis


def biradial(values_left, values_right,
             NM2=199,
             labelname='',
             normalize=False):

    fig = plt.figure()

    radial_left = fig.add_subplot(121, projection='polar')
    radial_right = fig.add_subplot(122, projection='polar')

    cax_left = fig.add_axes([0.05, 0.05, 0.03, 0.8])
    cax_right = fig.add_axes([0.92, 0.05, 0.03, 0.8])

    if normalize:
        normal_value = max(values_left.max(), values_right.max())
    else:
        normal_value = None

    radial_plot(values_left, NM2, labelname,
                radial_left, cax_left, normal_value)

    radial_plot(values_right, NM2, labelname,
                radial_right, cax_right, normal_value)

    cax_left.yaxis.set_ticks_position('left')
    cax_left.yaxis.set_label_position('left')

    return fig
