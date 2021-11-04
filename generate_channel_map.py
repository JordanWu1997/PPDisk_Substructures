#!/usr/bin/env python
# -*- coding:utf-8 -*-

# ############################################################################
# :Author: Sheng-Jun Lin, Kuan-Hsien Wu
# :Email: jordankhwu@gapp.nthu.edu.tw
# :Date: 2021-11-04
# :Description: This code is to generate channel map from ALMA data
# ############################################################################

import aplpy
import astropy.io.fits as pyfits
import matplotlib.pyplot as plt
import numpy as np
import pyregion
from astropy.wcs import WCS
from matplotlib import rc, rcParams, style

from ALMA_data_manipulation import ObsData

rc('text', usetex=True)
rcParams.update({'mathtext.default': 'regular'})
style.use('classic')


def plot_single_channel(fitsfile,
                        fig,
                        sub_tuple,
                        title,
                        colorbarlabel,
                        center_ra_deg=225.9325,
                        center_dec_deg=-37.93194,
                        dist_pc=148.,
                        width_asec=5.,
                        height_asec=5.,
                        cmap='inferno',
                        m=None,
                        M=None,
                        alpha=1,
                        dim=[0, 1],
                        north=False,
                        bmaj_deg=0,
                        bmin_deg=0,
                        bpa_deg=0,
                        sbar=False,
                        hide_Dec=False,
                        hide_RA=False,
                        **kwargs):

    sub = aplpy.FITSFigure(fitsfile,
                           figure=fig,
                           subplot=sub_tuple,
                           dimensions=dim,
                           north=north,
                           **kwargs)
    sub.show_colorscale(cmap=cmap, vmin=m, vmax=M)

    try:
        sub.recenter(center_ra_deg,
                     center_dec_deg,
                     width=width_asec / 3600,
                     height=height_asec / 3600)
    except Exception as msg:
        print(msg)

    sub.set_nan_color('w')
    sub.tick_labels.set_xformat('dd:mm:ss')
    sub.tick_labels.set_yformat('dd:mm:ss')
    sub.axis_labels.hide_x()
    sub.axis_labels.hide_y()
    if hide_Dec:
        sub.tick_labels.hide_y()
    if hide_RA:
        sub.tick_labels.hide_x()
    sub.ticks.set_color('white')

    if sbar:
        sub.add_scalebar(10. / dist_pc / 3600.)  #deg
        sub.scalebar.set_label('1000 au')
        sub.scalebar.set_linewidth(4)  #pt
        sub.scalebar.set_color('k')
    sub.image.set_alpha(alpha)

    if bmaj_deg != 0:
        pass
        # sub._header.update({'BMAJ': bmaj_deg})
        # sub._header.update({'BMIN': bmin_deg})
        # sub._header.update({'BPA': bpa_deg})
        # sub.add_beam()

    return sub


def plot_4x4_channel(PPDisk, chan_start=0):

    fig = plt.figure(figsize=(20, 20))

    vel_ms_axis = PPDisk.freqax2velax()
    for i in range(chan_start, chan_start + 16):
        print('channel = {:d}'.format(i))

        hide_Ra_flag, hide_Dec_flag = True, True
        if ((i - chan_start) % 4 == 0):
            hide_Dec_flag = False
        if (i - chan_start) > 11:
            hide_Ra_flag = False

        sub0 = plot_single_channel(PPDisk.fitsfile,
                                   fig, (4, 4, i - chan_start + 1),
                                   'vel km/s',
                                   '',
                                   hide_RA=hide_Ra_flag,
                                   hide_Dec=hide_Dec_flag,
                                   slices=[i + 1, 0],
                                   dim=[0, 1],
                                   M=0.013,
                                   m=-0.006,
                                   bmaj_deg=PPDisk.bmaj_deg,
                                   bmin_deg=PPDisk.bmin_deg,
                                   bpa_deg=PPDisk.bpa_deg)

        sub0.set_title('{:.2f} km/s'.format(vel_ms_axis[i] / 1000.),
                       size='x-large')

    plt.subplots_adjust(wspace=0.05, hspace=0.1)
    plt.suptitle("IMLup 12CO J=2-1 Channel Map", fontsize=20)
    fig.savefig('{}_chan_map_{:2d}_{:2d}.pdf'.format(PPDisk.name, chan_start,
                                                     chan_start + 16),
                bbox_inches='tight')


def main():
    IM_Lup = ObsData(
        '/mazu/users/jordan/PPDisk_Project/DSHARP_DR/IMLup/IMLup_CO.fits',
        158.,
        name='IM_Lup')
    IM_Lup.stellar_property(1.12, 4250)
    IM_Lup.disk_property(47.5,
                         144.5,
                         np.array([117.]),
                         np.array([117. * 0.13]),
                         offset_x=-1.5,
                         offset_y=1.0)
    IM_Lup.planet_property(np.array([69.]), np.array([0.77]))
    plot_4x4_channel(IM_Lup, chan_start=22)


if __name__ == '__main__':
    main()
