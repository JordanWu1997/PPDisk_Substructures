#!/usr/bin/env python
# -*- coding:utf-8 -*-

# ############################################################################
# :Author: Kuan-Hsien Wu
# :Email: jordankhwu@gapp.nthu.edu.tw
# :Date: 2021-11-04
# :Description: This code is to generate moment map from ALMA data by using
#               bettermoment CLI [https://bettermoments.readthedocs.io/en/\
#               latest/user/command_line.html]
# ############################################################################

import aplpy
import matplotlib.pyplot as plt
import numpy as np
from astropy.io import fits
from os import system, path

from ALMA_data_manipulation import ObsData
from azimuthal_PV_diagram import AzimuthalPVDiagram


def generate_disk_mask(PPDisk, radius_pix, save_mask=True):
    """Generate disk mask with given radius and disk PA, INC

    Args:
      PPDisk(ObsData class): disk for generating mask
      radius_pix(int): radius to mask (unit: pixel)
      save_mask: save output mask (Default value = True)

    Returns:
      4D float array: mask in 4D fr ALMA cube data
      str: output mask fitsfile filename

    """
    azimuthal_PV_diagram = AzimuthalPVDiagram()
    radii_sky_pix, _ = azimuthal_PV_diagram.set_disk_polar_to_2D_map(
        PPDisk.data[0], PPDisk.disk_pos, PPDisk.disk_inc, PPDisk.offset_x_pix,
        PPDisk.offset_y_pix)

    n_vel_chan = len(PPDisk.freqax2velax())
    mask_array = np.ones(np.shape(radii_sky_pix))
    mask_array[radii_sky_pix > radius_pix] = 0

    mask_3D_list = []
    for i in range(len(PPDisk.freqax2velax())):
        mask_3D_list.append(mask_array)
    mask_3D_array = np.array(mask_3D_list)
    mask_4D_array = np.array([mask_3D_array])

    mask_name = '{}_disk_mask.fits'.format(PPDisk.name)
    fits.writeto(mask_name,
                 mask_4D_array,
                 header=PPDisk.header,
                 overwrite=True)

    return mask_4D_array, mask_name


def generate_mom0_map(fitsfile, chan_start, chan_end, mask=None):
    """

    Args:
      fitsfile: cube data fitsfile
      chan_start: channel that starts to integrate
      chan_end: channel that ends
      mask: user-defined cube mask (Default value = None)

    Returns:

    """
    if mask == None:
        system(
            "bettermoments {} -firstchannel {:d} -lastchannel {:d} -method zeroth"
            .format(fitsfile, chan_start, chan_end))
    else:
        system(
            "bettermoments {} -firstchannel {:d} -lastchannel {:d} -method zeroth -mask {}"
            .format(fitsfile, chan_start, chan_end, mask))


def generate_mom1_map(fitsfile, chan_start, chan_end, mask=None):
    """

    Args:
      fitsfile: cube data fitsfile
      chan_start: channel that starts to integrate
      chan_end: channel that ends
      mask: user-defined cube mask (Default value = None)

    Returns:

    """
    if mask == None:
        system(
            "bettermoments {} -firstchannel {:d} -lastchannel {:d} -method first"
            .format(fitsfile, chan_start, chan_end))
    else:
        system(
            "bettermoments {} -firstchannel {:d} -lastchannel {:d} -method first -mask {}"
            .format(fitsfile, chan_start, chan_end, mask))


def plot_mom0_map(fitsfile):
    """Plot moment 0 map by Aplpy

    Args:
      fitsfile(str): moment 0 map file

    Returns:

    """
    fig = plt.figure(figsize=(8, 8))
    sub = aplpy.FITSFigure(fitsfile, figure=fig)
    sub.show_colorscale(cmap="inferno")
    sub.add_colorbar()
    plt.show()


def plot_mom1_map(fitsfile):
    """

    Args:
      fitsfile(str): moment 1 map file

    Returns:

    """
    fig = plt.figure(figsize=(8, 8))
    sub = aplpy.FITSFigure(fitsfile, figure=fig)
    sub.show_colorscale(cmap="coolwarm")
    sub.add_colorbar()
    sub.colorbar.set_axis_label_text('km/s')
    sub.colorbar.set_axis_label_rotation(270)
    plt.show()


def main():
    """Main function in generate_moment_map.py"""

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

    disk_mask, mask_name = generate_disk_mask(IM_Lup, 400)

    if not path.isfile(
            '/mazu/users/jordan/PPDisk_Project/DSHARP_DR/IMLup/IMLup_CO_M0.fits'
    ):
        generate_mom0_map(IM_Lup.fitsfile, 12, 45)
    if not path.isfile(
            '/mazu/users/jordan/PPDisk_Project/DSHARP_DR/IMLup/IMLup_CO_M1.fits'
    ):
        generate_mom1_map(IM_Lup.fitsfile,
                          12,
                          45,
                          mask='./{}'.format(mask_name))

    plot_mom0_map(
        '/mazu/users/jordan/PPDisk_Project/DSHARP_DR/IMLup/IMLup_CO_M0.fits')
    plot_mom1_map(
        '/mazu/users/jordan/PPDisk_Project/DSHARP_DR/IMLup/IMLup_CO_M1.fits')


if __name__ == '__main__':
    main()
