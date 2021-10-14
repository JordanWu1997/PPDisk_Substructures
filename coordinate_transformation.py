#!/usr/bin/env python
# -*- coding:utf-8 -*-

# ############################################################################
# :Author: Kuan-Hsien Wu
# :Email: jordankhwu@gapp.nthu.edu.tw
# :Date: 2021-03-28
# :Description: This code is to transform between disk/sky coordinate
# ############################################################################

import matplotlib.pyplot as plt
import numpy as np


class CoordinateTransformation():
    """ """
    def __init__(self):
        pass

    @staticmethod
    def transform_to_disk_coord(x, y, disk_pos, disk_inc):
        """Transform sky coord to disk coord

        Args:
          x(float): x position on sky coord
          y(float): y position on sky coord
          disk_pos(float): position angle of disk
          disk_inc(float): inclination angle of disk

        Returns:
          float: x position on disk coord
          float: y position on disk coord

        """
        pass

    @staticmethod
    def transform_to_sky_coord(x, y, disk_pos, disk_inc):
        """Transform disk coord to sky coord

        Args:
          x(float): x position on disk coord
          y(float): y position on disk coord
          disk_pos(float): position angle of disk
          disk_inc(float): inclination angle of disk

        Returns:
          float: x position on sky coord
          float: y position on sky coord

        """
        disk_pos_rad, disk_inc_rad = np.deg2rad(disk_pos), np.deg2rad(disk_inc)
        sky_x = x * np.cos(disk_pos_rad) + y * np.sin(disk_pos_rad) * np.cos(
            disk_inc_rad)
        sky_y = x * np.sin(disk_pos_rad) - y * np.cos(disk_pos_rad) * np.cos(
            disk_inc_rad)
        return sky_x, sky_y

    def transform_to_disk_coord_test(self):
        """Test transform_to_disk_coord function"""
        pass

    def transform_to_sky_coord_test(self):
        """Test transform_to_sky_coord function"""

        # Test ring
        max_width = 50.
        test_disk_pos = 120.
        test_disk_inc = 47.5
        test_radius = 30.

        # Axis for both disk coord and sky coord
        disk_axis = np.array([[max_width, -max_width, 0., 0.],
                              [0., 0., max_width, -max_width]])
        sky_axis = np.array(
            self.transform_to_sky_coord(disk_axis[0], disk_axis[1],
                                        test_disk_pos, test_disk_inc))
        plt.plot(disk_axis[0][:2],
                 disk_axis[1][:2],
                 c='r',
                 label='disk coord. axes')
        plt.plot(disk_axis[0][2:], disk_axis[1][2:], c='r')
        plt.plot(sky_axis[0][:2],
                 sky_axis[1][:2],
                 c='g',
                 label='sky coord. axes')
        plt.plot(sky_axis[0][2:], sky_axis[1][2:], c='g')

        # Test circular disk for coordinate transformation
        theta = np.linspace(0, 359, 360, endpoint=True)
        test_ring_disk_x = test_radius * np.cos(theta)
        test_ring_disk_y = test_radius * np.sin(theta)
        test_ring_sky_x, test_ring_sky_y = self.transform_to_sky_coord(
            test_ring_disk_x, test_ring_disk_y, test_disk_pos, test_disk_inc)
        plt.scatter(test_ring_disk_x,
                    test_ring_disk_y,
                    label='disk coord. ring',
                    color='r',
                    alpha=0.1)
        plt.scatter(test_ring_sky_x,
                    test_ring_sky_y,
                    label='sky coord. ring',
                    color='g',
                    alpha=0.1)

        # Test origin of polar angle for coordinate transformation
        test_dot_disk_x, test_dot_disk_y = test_radius, 0
        test_dot_sky_x, test_dot_sky_y = self.transform_to_sky_coord(
            test_dot_disk_x, test_dot_disk_y, test_disk_pos, test_disk_inc)
        plt.scatter(test_dot_disk_x,
                    test_dot_disk_y,
                    label='disk coord. polar angle origin',
                    color='r',
                    alpha=1)
        plt.scatter(test_dot_sky_x,
                    test_dot_sky_y,
                    label='sky coord. polar angle origin',
                    color='g',
                    alpha=1)

        # Option for plot
        plt.title(
            'Transform from disk coord. to sky coord. (disk_inc={:.1f}°, disk_pos={:.1f}°)'
            .format(test_disk_inc, test_disk_pos))
        plt.axis('equal')
        plt.legend()
        plt.show()


if __name__ == '__main__':
    coordinate_transformation = CoordinateTransformation()
    coordinate_transformation.transform_to_disk_coord_test()
    coordinate_transformation.transform_to_sky_coord_test()
