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
from azimuthal_PV_diagram import AzimuthalPVDiagram
from PPDisk_azimuthal_PV_diagram import (plot_cut_ring, plot_ring_PV_diagram,
                                         plot_sequential_ring_PV_diagram)


def main():
    """Main function in IM_Lup_azimuthal_PV_diagram.py"""

    # IM Lup (Load ALMA_observation data from DSHARP)
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
                                    save=True)


if __name__ == "__main__":
    main()
