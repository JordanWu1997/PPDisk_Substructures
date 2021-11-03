#!/bin/env python
# -*- coding:utf-8 -*-

# ############################################################################
# :Author: Kuan-Hsien Wu
# :Email: jordankhwu@gapp.nthu.edu.tw
# :Date: 2021-10-31
# :Description: This code is for ALMA protoplanetary disk data PV cut diagram
# ############################################################################

import matplotlib
import matplotlib.pyplot as plt

from ALMA_data_manipulation import ObsData
from azimuthal_PV_diagram import generate_azimuthal_PV_diagram
from generate_Keplerian_disk_PV import *


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
        planet_pos_on_disk_axis = 360. - (PPDisk.disk_pos +
                                          PPDisk.planet_pos[0])
        plt.axvline(
            planet_pos_on_disk_axis,
            c='k',
            ls='-',
            label='planet candidate PA\n(r={:.1f} au ({:d} pix))'.format(
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
            'Azimuthal Position-Velocity Diagram\n(r_center={:.1f} au ({:d} pix), r_width={:.1f} au ({:d} pix), star={:.1f} solar_mass)'
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


def main():
    """Main function in PPDisk_azimuthal_PV_diagram.py"""

    # Load ALMA_observation data
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

    # plot_planet_ring_PV_diagram(PPDisk)
    r_center_list = [
        ObsData.round_to_pix(PPDisk.au2pix(PPDisk.gap_location))[0] +
        ObsData.round_to_pix(PPDisk.au2pix(PPDisk.gap_width))[0] * i
        for i in range(-5, 6)
    ]
    r_width_list = [ObsData.round_to_pix(PPDisk.au2pix(PPDisk.gap_width))[0]
                    ] * len(r_center_list)

    matplotlib.rcParams.update({'font.size': 14})
    plt.subplots(11, 1, figsize=(18, 28))
    for i, (r_center, r_width) in enumerate(zip(r_center_list, r_width_list)):
        print(i)
        if i != len(r_center_list) - 1:
            show_sub_axis_flag = False
        else:
            show_sub_axis_flag = True

        plt.subplot(len(r_center_list), 1, i + 1)
        # plt.subplots_adjust(left=0.05, right=0.82)
        plot_ring_PV_diagram(PPDisk,
                             r_center,
                             r_width,
                             chan_start=16,
                             chan_end=-16,
                             new_figure=False,
                             show_sub_x_axis=show_sub_axis_flag,
                             show_sub_y_axis=True,
                             show_sub_title=False,
                             show_plot=False)
        plt.title(
            'r_center={:.1f} au ({:d} pix), r_width={:.1f} au ({:d} pix)'.
            format(PPDisk.pix2au(r_center), r_center, PPDisk.pix2au(r_width),
                   r_width))

    plt.suptitle(
        "{} Azimuthal Position-Velocity Diagram\n(beam={:.1f} au x {:.1f} au, channel_width={:.1f} m/s)"
        .format(
            PPDisk.name,
            PPDisk.asec2au(PPDisk.deg2asec(PPDisk.bmaj_deg), PPDisk.dist_pc),
            PPDisk.asec2au(PPDisk.deg2asec(PPDisk.bmin_deg), PPDisk.dist_pc),
            PPDisk.get_channel_width()))

    figure, axe = plt.gcf(), plt.gca()
    lines, labels = axe.get_legend_handles_labels()
    plt.figlegend(lines, labels, ncol=1, bbox_to_anchor=(0.98, 0.98))

    plt.savefig('{}_azimuthal_PV_diagram.pdf'.format(PPDisk.name))


if __name__ == '__main__':
    main()
