#!/usr/bin/env python
# -*- coding:utf-8 -*-

# ############################################################################
# :Author: Kuan-Hsien Wu
# :Email: jordankhwu@gapp.nthu.edu.tw
# :Date: 2021-10-31
# :Description: This code is for HD 97048 protoplanetary disk ALMA data PV cut diagram
# ############################################################################

from sys import path

path.append(
    '/mazu/users/jordan/PPDisk_Project/Python_Package/Source_Code/PPDisk_Substructures'
)

import matplotlib
import matplotlib.pyplot as plt
import numpy as np
from ALMA_data_manipulation import ObsData
from astropy.io import fits
from azimuthal_PV_diagram import AzimuthalPVDiagram, AzimuthalPVDiagramTest
from PPDisk_azimuthal_PV_diagram import (plot_cut_ring, plot_ring_PV_diagram,
                                         plot_sequential_ring_PV_diagram)


def main():
    """Main function in HD_97048_azimuthal_PV_diagram.py"""

    # HD 97048 (Load ALMA_observation data from Pinte)
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

    # PLot cut ring for visualization
    # ref_2D_data = fits.getdata('/mazu/users/jordan/PPDisk_Project/Pinte_DR/'\
                             # 'HD_97048_b7_continuum_centered_selfcal_manmask_briggs.image.fits')[0, 0]
    ref_2D_data = fits.getdata('/mazu/users/jordan/PPDisk_Project/Pinte_DR/'\
                             'HD_97048_13CO32_briggs_selfcal_nocontsub.image_M1.fits')

    azimuthal_PV_diagram_test = AzimuthalPVDiagramTest()
    azimuthal_PV_diagram_test.set_disk_polar_to_2D_map_test(
        test_2D_map=ref_2D_data,
        test_disk_pos=PPDisk.disk_pos,
        test_disk_inc=PPDisk.disk_inc,
        test_ring_radius_sky_pix=20,
        test_offset_x_pix=PPDisk.offset_x_pix,
        test_offset_y_pix=PPDisk.offset_y_pix)

    plot_cut_ring(PPDisk,
                  r_center_list,
                  r_width_list,
                  ref_2D_data=ref_2D_data,
                  ref_cmap='inferno')

    # Plot sequential ring PV diagram
    plot_sequential_ring_PV_diagram(PPDisk,
                                    r_center_list,
                                    r_width_list,
                                    chan_start=55,
                                    chan_end=-55,
                                    save=True)


if __name__ == "__main__":
    main()
