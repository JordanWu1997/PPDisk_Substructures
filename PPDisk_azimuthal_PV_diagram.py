#!/usr/bin/env python
# -*- coding:utf-8 -*-

# ############################################################################
# :Author: Kuan-Hsien Wu
# :Email: jordankhwu@gapp.nthu.edu.tw
# :Date: 2021-10-31
# :Description: This code is for ALMA protoplanetary disk data PV cut diagram
# ############################################################################

import matplotlib
import matplotlib.pyplot as plt
import numpy as np
from astropy.io import fits

from ALMA_data_manipulation import ObsData
from azimuthal_PV_diagram import generate_azimuthal_PV_diagram
from generate_Keplerian_disk_PV import plot_LOS_Kepleran_PV_cut
from azimuthal_PV_diagram import AzimuthalPVDiagram


def plot_ring_PV_diagram(
        PPDisk,
        r_center,
        r_width,
        chan_start=0,
        chan_end=-1,
        plot_planet_loc=True,
        color_list=['lime', 'cyan', 'red', 'magenta', 'yellow', 'blue'],
        aspect_ratio_list=[0.0, 0.2, 0.4],
        new_figure=True,
        show_sub_title=True,
        show_sub_x_axis=False,
        show_sub_y_axis=False,
        show_sub_legend=False,
        show_plot=True):
    """

    Args:
      PPDisk(ObsData class): Protoplanetary disk
      r_center(int): center of radius (unit: pixel)
      r_width(int): width of ring (unit: pixel)
      chan_start(int, optional): index of chan to start plotting (Default value = 1)
      chan_end(int, optional): index of chan to end plotting (Default value = 1)
      plot_planet_loc(bool, optional): whethere plot planet (Default value = True)
      color_list(1D str list, optional): color list for Keplerian motion plot
      aspect_ratio_list(1D float list, optional): aspect ratio for disk model
      new_figure(bool, optional): whether use new figure (Default value = True)
      show_sub_title(bool, optional): whether show title (Default value = True)
      show_sub_x_axis(bool, optional): whether show x axis (Default value = False)
      show_sub_y_axis(bool, optional): whether show y axis (Default value = False)
      show_sub_legend(bool, optional): whether show legend (Default value = False)
      show_plot(bool, optional): whether show plot (Default value = True)

    Returns:

    """

    if new_figure:
        matplotlib.rcParams.update({'font.size': 10})
        plt.figure(figsize=(18, 8))

    generate_azimuthal_PV_diagram(PPDisk.data,
                                  PPDisk.freqax2velax(),
                                  r_center,
                                  r_width,
                                  PPDisk.disk_pos,
                                  PPDisk.disk_inc,
                                  PPDisk.offset_x_pix,
                                  PPDisk.offset_y_pix,
                                  new_figure=new_figure,
                                  show_PV_cut_data_point=False,
                                  show_title=show_sub_title,
                                  show_x_axis=show_sub_x_axis,
                                  show_y_axis=show_sub_y_axis,
                                  show_legend=show_sub_legend,
                                  show_plot=show_plot)

    if plot_planet_loc:
        for planet_pos, planet_sep in zip(PPDisk.planet_pos,
                                          PPDisk.planet_sep):
            planet_pos_on_disk_axis = 360. - (PPDisk.disk_pos +
                                              PPDisk.planet_pos)

            if abs(r_center -
                   PPDisk.au2pix(PPDisk.asec2au(planet_sep, PPDisk.dist_pc))
                   ) <= r_width / 2:
                plt.axvline(
                    planet_pos_on_disk_axis,
                    c='r',
                    ls='-',
                    label=
                    'planet candidate PA\n@ planet location\n(r={:.1f} au ({:d} pix))'
                    .format(
                        PPDisk.asec2au(PPDisk.planet_sep[0], PPDisk.dist_pc),
                        int(PPDisk.asec2pix(PPDisk.planet_sep[0]))))
            else:
                plt.axvline(
                    planet_pos_on_disk_axis,
                    c='k',
                    ls='-',
                    label='planet candidate PA\n(r={:.1f} au ({:d} pix))'.
                    format(
                        PPDisk.asec2au(PPDisk.planet_sep[0], PPDisk.dist_pc),
                        int(PPDisk.asec2pix(PPDisk.planet_sep[0]))))

    for aspect_ratio, color in zip(aspect_ratio_list, color_list):
        plot_LOS_Kepleran_PV_cut(PPDisk,
                                 r_center,
                                 aspect_ratio,
                                 theta_n=361,
                                 plot_SI=True,
                                 plot_chan=False,
                                 plot_thin_disk_r=False,
                                 plot_thick_disk_r=True,
                                 plot_thick_disk_r_z=False,
                                 color_thin_disk_r=color,
                                 color_thick_disk_r=color,
                                 color_thick_disk_r_z=color,
                                 new_figure=False,
                                 show_title=False,
                                 show_x_axis=False,
                                 show_y_axis=False,
                                 show_legend=False,
                                 show_plot=False)

    if show_sub_title:
        plt.title(
            'Azimuthal Position-Velocity Diagram\n'\
            '(r_center={:.1f} au ({:d} pix), r_width={:.1f} au ({:d} pix),'\
            ' star={:.1f} solar_mass)'
            .format(PPDisk.pix2au(r_center), r_center, PPDisk.pix2au(r_width),
                    r_width, PPDisk.stellar_mass))

    axes = plt.gca()
    axes.set_xlim([0., 360.])
    axes.set_ylim(
        [PPDisk.freqax2velax()[chan_start],
         PPDisk.freqax2velax()[chan_end]])

    if show_sub_legend:
        leg = plt.legend()
        for lh in leg.legendHandles:
            lh.set_alpha(1)

    if show_plot:
        plt.show()


def plot_cut_ring(PPDisk,
                  r_center_list,
                  r_width_list,
                  image_size=400,
                  ref_2D_data=None,
                  ref_cmap='coolwarm',
                  show_plot=True):
    """

    Args:
      PPDisk(ObsData class): Protoplanetary disk
      r_center_list(1D int list): list of the center of the ring (unit: pixel)
      r_width_list(1D int list): list of the width of the ring (unit: pixel)
      image_size(int, optional): image size to plot (unit: pixel) (Default value = 400)
      ref_2D_data(2D float array, optional): data to plot reference map (Default value = None)
      ref_cmap(str. optional): colormap for plotting reference map (Default value = 'coolwarm')
      show_plot(bool, optional): whether show plot (Default value = True)

    Returns:

    """
    # Plot cut location
    matplotlib.rcParams.update({'font.size': 12})
    if ref_2D_data is None:
        cut_ring_array = np.zeros(PPDisk.data[0].shape)
        azimuthal_PV_diagram = AzimuthalPVDiagram()
        for i, (r_center, r_width) in enumerate(zip(r_center_list, r_width_list)):
            print('{:d}, r_center={:d} pix, r_width={:d} pix'.format(
                i, r_center, r_width))
            radii, _ = azimuthal_PV_diagram.set_disk_polar_to_2D_map(
                PPDisk.data[0], PPDisk.disk_pos, PPDisk.disk_inc,
                PPDisk.offset_x_pix, PPDisk.offset_y_pix)
            cut_ring_array[abs(radii - r_center) <= r_width / 2.] = r_center

        plt.figure(figsize=(8, 8))
        plt.imshow(cut_ring_array, origin='lower', cmap="gray")
        plt.xlabel('pixel')
        plt.ylabel('pixel')
        ax = plt.gca()
        ax.set_xlim(
            np.shape(PPDisk.data[0])[0] / 2 - image_size / 2,
            np.shape(PPDisk.data[0])[0] / 2 + image_size / 2)
        ax.set_ylim(
            np.shape(PPDisk.data[0])[1] / 2 - image_size / 2,
            np.shape(PPDisk.data[0])[1] / 2 + image_size / 2)
        cbar = plt.colorbar()
        cbar.set_label('ring radius (pixel)')
    else:
        cut_ring_array = np.zeros(ref_2D_data.shape)
        azimuthal_PV_diagram = AzimuthalPVDiagram()
        for i, (r_center, r_width) in enumerate(zip(r_center_list, r_width_list)):
            print('{:d}, r_center={:d} pix, r_width={:d} pix'.format(
                i, r_center, r_width))
            radii, _ = azimuthal_PV_diagram.set_disk_polar_to_2D_map(
                ref_2D_data, PPDisk.disk_pos, PPDisk.disk_inc,
                PPDisk.offset_x_pix, PPDisk.offset_y_pix)
            cut_ring_array[abs(radii - r_center) <= r_width / 2.] = r_center
        plt.subplots(1, 2, figsize=(16, 8))
        plt.subplot(1, 2, 1)
        plt.imshow(cut_ring_array, origin='lower', cmap="gray")
        plt.xlabel('pixel')
        plt.ylabel('pixel')
        ax = plt.gca()
        ax.set_xlim(
            np.shape(ref_2D_data)[0] / 2 - image_size / 2,
            np.shape(ref_2D_data)[0] / 2 + image_size / 2)
        ax.set_ylim(
            np.shape(ref_2D_data)[1] / 2 - image_size / 2,
            np.shape(ref_2D_data)[1] / 2 + image_size / 2)
        cbar = plt.colorbar()
        cbar.set_label('ring radius (pixel)')
        plt.subplot(1, 2, 2)
        plt.imshow(ref_2D_data, cmap=ref_cmap)
        ax = plt.gca()
        ax.set_xlim(
            np.shape(ref_2D_data)[0] / 2 - image_size / 2,
            np.shape(ref_2D_data)[0] / 2 + image_size / 2)
        ax.set_ylim(
            np.shape(ref_2D_data)[1] / 2 - image_size / 2,
            np.shape(ref_2D_data)[1] / 2 + image_size / 2)
        plt.colorbar()

    if show_plot:
        plt.show()


def plot_sequential_ring_PV_diagram(PPDisk,
                                    r_center_list,
                                    r_width_list,
                                    chan_start=0,
                                    chan_end=-1,
                                    save=True,
                                    show_plot=True):
    """

    Args:
      PPDisk(ObsData class): Protoplanetary disk
      r_center_list(1D int list): list of the center of the ring (unit: pixel)
      r_width_list(1D int list): list of the width of the ring (unit: pixel)
      chan_start(int, optional): start channel for velocity axis (Default value = 0)
      chan_end(int, optional): end channel for velocity axis (Default value = -1)
      save(bool, optional): save as pdf file (Default value = True)
      show_plot(bool, optional): whether show plot (Default value = True)

    Returns:

    """
    matplotlib.rcParams.update({'font.size': 12})
    plt.subplots(len(r_center_list), 1, figsize=(18, 2.1 * len(r_center_list)))
    for i, (r_center, r_width) in enumerate(zip(r_center_list, r_width_list)):
        print('{:d}, r_center={:d} pix, r_width={:d} pix'.format(
            i, r_center, r_width))

        if i != len(r_center_list) - 1:
            show_sub_axis_flag = False
        else:
            show_sub_axis_flag = True

        plt.subplot(len(r_center_list), 1, i + 1)
        plot_ring_PV_diagram(PPDisk,
                             r_center,
                             r_width,
                             chan_start=chan_start,
                             chan_end=chan_end,
                             new_figure=False,
                             show_sub_x_axis=show_sub_axis_flag,
                             show_sub_y_axis=True,
                             show_sub_title=False,
                             show_plot=False)

        for gap_location, gap_width in zip(PPDisk.gap_location,
                                           PPDisk.gap_width):
            if abs(r_center - PPDisk.au2pix(gap_location)) <= r_width / 2:
                title_obj = plt.title(
                    'r_center={:.1f} au ({:d} pix), r_width={:.1f} au ({:d} pix)'\
                    '@ gap (center={:.1f} au ({:d} pix), width={:.1f} au ({:d} pix))'
                    .format(PPDisk.pix2au(r_center), r_center,
                            PPDisk.pix2au(r_width), r_width,
                            gap_location, ObsData.round_to_pix(PPDisk.au2pix(gap_location)),
                            gap_width, ObsData.round_to_pix(PPDisk.au2pix(gap_width))
                            ))
                plt.setp(title_obj, color='r')
            else:
                plt.title(
                    'r_center={:.1f} au ({:d} pix), r_width={:.1f} au ({:d} pix)'
                    .format(PPDisk.pix2au(r_center), r_center,
                            PPDisk.pix2au(r_width), r_width))

    plt.suptitle(
        '{} Azimuthal Position-Velocity Diagram\n'\
        '(beam={:.1f} au x {:.1f} au, channel_width={:.1f} m/s)\n'\
        '(systematic_velocity={:.1f} m/s)'\
        .format(
            PPDisk.name,
            PPDisk.asec2au(PPDisk.deg2asec(PPDisk.bmaj_deg), PPDisk.dist_pc),
            PPDisk.asec2au(PPDisk.deg2asec(PPDisk.bmin_deg), PPDisk.dist_pc),
            PPDisk.get_channel_width(), PPDisk.sys_vel))

    figure, axe = plt.gcf(), plt.gca()
    lines, labels = axe.get_legend_handles_labels()
    plt.figlegend(lines, labels, ncol=1, bbox_to_anchor=(0.98, 0.98))

    if save:
        plt.savefig('{}_azimuthal_PV_diagram.pdf'.format(PPDisk.name))

    if show_plot:
        plt.show()


def main():
    """Main function in PPDisk_azimuthal_PV_diagram.py"""

    # IM Lup (Load ALMA_observation data)
    PPDisk = ObsData(
        '/mazu/users/jordan/PPDisk_Project/DSHARP_DR/IMLup/IMLup_CO.fits',
        158.,
        name='IM_Lup')
    PPDisk.stellar_property(1.12, 4250)
    PPDisk.disk_property(47.5,
                         144.5,
                         np.array([117.]),
                         np.array([117. * 0.13]),
                         offset_x=-1.5,
                         offset_y=1.0)
    PPDisk.planet_property(np.array([69.]), np.array([0.77]))

    # Generate r_center, r_width of ring for azimuthal PV diagram
    r_center_list = [
        ObsData.round_to_pix(PPDisk.au2pix(PPDisk.gap_location))[0] +
        ObsData.round_to_pix(PPDisk.au2pix(PPDisk.gap_width))[0] * i
        for i in range(-5, 10)
    ]
    r_width_list = [ObsData.round_to_pix(PPDisk.au2pix(PPDisk.gap_width))[0]
                    ] * len(r_center_list)

    # PLot cut ring for visualization
    ref_2D_data = fits.getdata(
        '/mazu/users/jordan/PPDisk_Project/DSHARP_DR/IMLup/IMLup_CO_M1.fits')
    plot_cut_ring(PPDisk, r_center_list, r_width_list, ref_2D_data=ref_2D_data)

    # Plot sequential ring PV diagram
    plot_sequential_ring_PV_diagram(PPDisk,
                                    r_center_list,
                                    r_width_list,
                                    chan_start=16,
                                    chan_end=-16,
                                    save=False)


if __name__ == '__main__':
    main()
