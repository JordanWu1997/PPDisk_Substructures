#!/usr/bin/env python
# -*- coding:utf-8 -*-

# ############################################################################
# :Author: Kuan-Hsien Wu
# :Email: jordankhwu@gapp.nthu.edu.tw
# :Date: 2021-10-31
# :Description: This code is for HD 97048 protoplanetary disk ALMA data PV cut diagram
# ############################################################################

import matplotlib
import matplotlib.pyplot as plt
import numpy as np
from astropy.io import fits

from ALMA_data_manipulation import ObsData
from azimuthal_PV_diagram import AzimuthalPVDiagram
from PPDisk_azimuthal_PV_diagram import plot_ring_PV_diagram


def main():

    # HD 97048 (Load ALMA_observation data)
    PPDisk = ObsData(
        '/mazu/users/jordan/PPDisk_Project/Pinte_DR/HD_97048_13CO32_briggs_selfcal_nocontsub.image.fits',
        185.,
        name='HD_97048')
    PPDisk.stellar_property(2.4, 4750)
    PPDisk.disk_property(-38.,
                         185.,
                         np.array([130.]),
                         np.array([130. * 12 / 38]),
                         offset_x=0.0,
                         offset_y=0.0)
    PPDisk.planet_property(np.array([55.]), np.array([0.45]))

    # plot_planet_ring_PV_diagram(PPDisk)
    r_center_list = [
        ObsData.round_to_pix(PPDisk.au2pix(PPDisk.gap_location))[0] +
        ObsData.round_to_pix(PPDisk.au2pix(PPDisk.gap_width))[0] * i
        for i in range(-2, 10)
    ]
    r_width_list = [ObsData.round_to_pix(PPDisk.au2pix(PPDisk.gap_width))[0]
                    ] * len(r_center_list)

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
                             chan_start=55,
                             chan_end=-55,
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

    plt.savefig('{}_azimuthal_PV_diagram.pdf'.format(PPDisk.name))

    # Plot cut location
    cont_data = fits.getdata('/mazu/users/jordan/PPDisk_Project/Pinte_DR/'\
                             'HD_97048_b7_continuum_centered_selfcal_manmask_briggs.image.fits')[0,0]
    cut_ring_array = np.zeros(cont_data.shape)
    azimuthal_PV_diagram = AzimuthalPVDiagram()
    for i, (r_center, r_width) in enumerate(zip(r_center_list, r_width_list)):
        print('{:d}, r_center={:d} pix, r_width={:d} pix'.format(
            i, r_center, r_width))
        radii, _ = azimuthal_PV_diagram.set_disk_polar_to_2D_map(
            cont_data, PPDisk.disk_pos, PPDisk.disk_inc,
            PPDisk.offset_x_pix, PPDisk.offset_y_pix)
        cut_ring_array[abs(radii - r_center) <= r_width / 2.] = r_center

    matplotlib.rcParams.update({'font.size': 12})
    plt.subplots(1, 2, figsize=(16,8))
    plt.subplot(1,2,1)
    plt.imshow(cut_ring_array, origin='lower', cmap="gray")
    ax = plt.gca()
    ax.set_xlim(np.shape(cont_data)[0]/2 - 200, np.shape(cont_data)[0]/2 + 200)
    ax.set_ylim(np.shape(cont_data)[1]/2 - 200, np.shape(cont_data)[1]/2 + 200)
    plt.colorbar()
    plt.subplot(1,2,2)
    plt.imshow(cont_data)
    ax = plt.gca()
    ax.set_xlim(np.shape(cont_data)[0]/2 - 200, np.shape(cont_data)[0]/2 + 200)
    ax.set_ylim(np.shape(cont_data)[1]/2 - 200, np.shape(cont_data)[1]/2 + 200)
    plt.colorbar()
    plt.show()


if __name__ == "__main__":
    main()
