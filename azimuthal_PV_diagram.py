#!/usr/bin/env python
# -*- coding:utf-8 -*-

# ############################################################################
# :Author: Kuan-Hsien Wu
# :Email: jordankhwu@gapp.nthu.edu.tw
# :Date: 2021-03-28
# :Description: This code is to generate azimuthal PV diagram from ALMA data
# ############################################################################

import matplotlib.pyplot as plt
import numpy as np
from scipy.interpolate import griddata

from ALMA_data_manipulation import ObsData
from coordinate_transformation import CoordinateTransformation


class AzimuthalPVDiagram():
    """Generate azimuthal position-velocity diagram"""
    def __init__(self):
        pass

    @staticmethod
    def round_to_pix(prop):
        """Round input float property to pixel number

        Args:
          prop(float): property to be rounded to pixel(integer)

        Returns:
          int: pixel number

        """
        return np.rint(prop).astype(int)

    def find_physical_in_pix(self, prop, r_physical_pix):
        """Find physical property corresponding pixel number

        Args:
          prop(float): physical property
          r_physical_pixel(float): ratio between physical property and pixel

        Returns:
          int: pixel number

        """
        return self.round_to_pix(prop / r_physical_pix)

    def set_disk_polar_to_vel_map(self, vel_map, disk_pos, disk_inc, offset_x,
                                  offset_y, r_physical_pix, return_cartesian=False):
        """Set disk corresponding radius and theta on sky plane for 2D map

        Args:
          vel_map(2D float array):
          disk_pos(float): position angle of disk
          disk_inc(float): inclination angle of disk
          offset_x(float): disk center x position offset
          offset_y(float): disk center y position offset
          r_physical_pix:
          r_physical_pixel(float): ratio of physical property and pixel
          return_cartesian: set True to return x,y position on sky plane

        Returns:
          2D float array: Disk corresponding radius on sky plane for 2D map
          2D float array: Disk corresponding theta on sky plane for 2D map
          2D float array: Disk corresponding x position on sky plane for 2D map
          2D float array: Disk corresponding y position on sky plane for 2D map

        """
        empty_map = np.empty_like(vel_map)
        range_y_pix, range_x_pix = np.shape(empty_map)
        offset_x_pix = self.find_physical_in_pix(offset_x, r_physical_pix)
        offset_y_pix = self.find_physical_in_pix(offset_y, r_physical_pix)
        disk_center_x_pix = int(range_x_pix / 2) + offset_x_pix
        disk_center_y_pix = int(range_y_pix / 2) + offset_y_pix

        x_pix = np.arange(0, range_x_pix, 1).astype(int) - disk_center_x_pix
        y_pix = np.arange(0, range_y_pix, 1).astype(int) - disk_center_y_pix
        xs_pix, ys_pix = np.meshgrid(x_pix, y_pix)
        xs_sky_pix, ys_sky_pix = CoordinateTransformation.transform_to_sky_coord(
            xs_pix, ys_pix, disk_pos, disk_inc)
        thetas_sky_deg = np.rad2deg(np.arctan2(ys_sky_pix, xs_sky_pix))
        radii_sky_pix = (xs_sky_pix**2 + ys_sky_pix**2)**0.5

        if return_cartesian:
            return radii_sky_pix, thetas_sky_deg, xs_sky_pix, ys_sky_pix
        else:
            return radii_sky_pix, thetas_sky_deg

    def set_disk_polar_to_vel_map_test(self):
        """Test set_disk_polar_to_vel_map"""

        test_vel_map = np.ones((100, 100))
        test_disk_pos = 120
        test_disk_inc = 47.5
        test_radius_sky_pix = 30
        test_offset_x = 3.0
        test_offset_y = 2.0
        test_r_physical_pix = 0.3

        test_radius_sky_pix, test_theta_sky, x_pixs, y_pixs = self.set_disk_polar_to_vel_map(
            test_vel_map, test_disk_pos, test_disk_inc, test_offset_x,
            test_offset_y, test_r_physical_pix, return_cartesian=True)

        plt.subplots(1, 4, figsize=(16, 4))
        plt.subplot(1, 4, 1)
        plt.imshow(test_radius_sky_pix, origin='lower')
        plt.title("radius (x^2+y^2(disk->sky))")
        plt.colorbar()
        plt.grid()
        plt.subplot(1, 4, 2)
        plt.imshow(test_theta_sky, origin='lower')
        plt.title("theta (arctan(y/x) sky)")
        plt.colorbar()
        plt.grid()
        plt.subplot(1, 4, 3)
        plt.imshow(x_pixs, origin='lower')
        plt.title("x (disk->sky)")
        plt.colorbar()
        plt.grid()
        plt.subplot(1, 4, 4)
        plt.imshow(y_pixs, origin='lower')
        plt.title("y (disk->sky)")
        plt.colorbar()
        plt.grid()
        plt.suptitle(
            'Transform from disk coord. to sky coord. (disk_inc={:.1f}°, disk_pos={:.1f}°)'
            .format(test_disk_inc, test_disk_pos))
        plt.show()


if __name__ == '__main__':
    azimuthal_PV_diagram = AzimuthalPVDiagram()
    azimuthal_PV_diagram.set_disk_polar_to_vel_map_test()
