#!/bin/env python
# -*- coding:utf-8 -*-
# #########################################################
# Author : Sheng-Jun Lin, Kuan-Hsien Wu
# Email : sj.lin@gapp.nthu.edu.tw
#         jordankhwu@gapp.nthu.edu.tw
# Description : This is code to generate channel map
# Date : 2021-03-28
# #########################################################

from matplotlib import rcParams
from astropy.wcs import WCS
from matplotlib import rc, rcParams, style

import astropy.constants as const
import matplotlib as mpl
import matplotlib.pyplot as plt
import astropy.io.fits as pyfits
import numpy as np
import aplpy
import pyregion
import sys
import copy

# Matplotlib setup
rc('text', usetex=True)
rcParams.update({'mathtext.default': 'regular'})
style.use('classic')

# Data input and output
cube_file = '/mazu/users/jordan/DSHARP_Project/DSHARP_DR/IMLup/IMLup_CO.fits'
save_name = 'IM_Lup_Channel_Maps.png'

# Physical property
d_pc = 158
w_asec = 2.5
h_asec = 2.5

# Functions
def plot_a_map(fitsfile, fig, sub_tuple, title, colorbarlabel, \
               inp_cmap='Blues_r', factor=1, \
               m=None, M=None, alpha=1, dim=[0,1], **kwargs):
    '''
    This code is modified from Seb's work
    This is to plot channel map of ALMA cube data
    '''
    sub = aplpy.FITSFigure(fitsfile, figure=fig, subplot=sub_tuple, \
                           dimensions=dim, **kwargs)
    if factor != 1:
        sub._data *= factor
    sub.show_colorscale(cmap=inp_cmap, vmin=m, vmax=M)
    sub.set_nan_color('w')
    sub.tick_labels.set_xformat('hh:mm:ss')
    sub.tick_labels.set_yformat('dd:mm:ss')
    sub.ticks.set_color('black')
    return sub

def main():
    '''
    Main code is to plot channel map
    '''
    # Plot channel map
    fig = plt.figure(figsize=(45*0.7,60*0.7))
    for i in range(15, 45):
        print(i)
        sub0 = plot_a_map(cube_file,
                          fig, (6,5, i-14), 'vel km/s', '',
                          m=-0.0005, M=0.03,
                          slices=[i, 0], dim=[0,1])
        frq_rpix = sub0._header['CRPIX3']
        frq_delt = sub0._header['CDELT3']
        frq_rval = sub0._header['CRVAL3']
        restfrq  = sub0._header['RESTFRQ']
        frq = frq_rval + (i + 1 - frq_rpix) * frq_delt
        vel_kms = (restfrq - frq)/restfrq * const.c.value*1e-3
    plt.savefig(save_name)

if __name__ == '__main__':
    main()
