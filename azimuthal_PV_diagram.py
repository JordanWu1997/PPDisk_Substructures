#!/usr/bin/env python
# -*- coding:utf-8 -*-

# ############################################################################
# :Author: Kuan-Hsien Wu
# :Email: jordankhwu@gapp.nthu.edu.tw
# :Date: 2021-03-28
# :Description: This code is to generate azimuthal PV diagram from ALMA data
# ############################################################################

import matplotlib
import matplotlib.pyplot as plt
import numpy as np

from ALMA_data_manipulation import ObsData


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
          r_physical_pix(float): ratio between physical property and pixel

        Returns:
          int: pixel number

        """
        return self.round_to_pix(prop / r_physical_pix)

    @staticmethod
    def rotate_disk_dot_to_sky_dot(x, y, disk_pos):
        """Rotate disk dot to sky dot (No shrink)

        Args:
          x(float): x position on disk coord
          y(float): y position on disk coord
          disk_pos(float): position angle of disk (unit: deg)

        Returns:
          float: x position on sky coord
          float: y position on sky coord

        """
        disk_pos_rad = np.deg2rad(disk_pos)

        # Rotate disk coordinate to sky coord
        sky_x = x * np.cos(-disk_pos_rad) + y * np.sin(-disk_pos_rad)
        sky_y = x * -np.sin(-disk_pos_rad) + y * np.cos(-disk_pos_rad)

        return sky_x, sky_y

    @staticmethod
    def rotate_disk_axis_to_sky_axis(x, y, disk_pos):
        """Rotate disk axis to sky axis (No shrink)

        Args:
          x(float): x position on disk coord
          y(float): y position on disk coord
          disk_pos(float): position angle of disk (unit: deg)

        Returns:
          float: x position on sky coord
          float: y position on sky coord

        """
        disk_pos_rad = np.deg2rad(disk_pos)

        # Rotate disk coordinate to sky coord
        sky_x = x * np.cos(-disk_pos_rad) + y * -np.sin(-disk_pos_rad)
        sky_y = x * np.sin(-disk_pos_rad) + y * np.cos(-disk_pos_rad)

        return sky_x, sky_y

    def set_disk_polar_to_2D_map(self,
                                 map_2D,
                                 disk_pos,
                                 disk_inc,
                                 offset_x_pix,
                                 offset_y_pix,
                                 return_cartesian=False):
        """Set disk corresponding radius and theta on sky plane for 2D map

        Args:
          map_2D(2D float array): 2D map
          disk_pos(float): disk position angle (unit: deg)
          disk_inc(float): disk inclination angle (unit: deg)
          offset_x_pix(int): x-offset of disk center (unit: pixel)
          offset_y_pix(int): y-offset of disk center (unit: pixel)
          return_cartesian(bool): return cartesian coordinates (Default value = False)

        Returns:
          2D float array: Disk corresponding radius on sky plane for 2D map 2D float array: Disk corresponding theta on sky plane for 2D map 2D float array: Disk corresponding x position on sky plane for 2D map 2D float array: Disk corresponding y position on sky plane for 2D map

        """
        empty_map = np.empty_like(map_2D)
        range_y_pix, range_x_pix = np.shape(empty_map)
        disk_center_x_pix = int(range_x_pix / 2) + offset_x_pix
        disk_center_y_pix = int(range_y_pix / 2) + offset_y_pix

        x_pix = np.arange(0, range_x_pix, 1).astype(int) - disk_center_x_pix
        y_pix = np.arange(0, range_y_pix, 1).astype(int) - disk_center_y_pix
        xs_pix, ys_pix = np.meshgrid(x_pix, y_pix)

        # Additional 90 deg is needed for 0 deg starts at x-axis (disk major axis)
        xs_sky_pix, ys_sky_pix = self.rotate_disk_axis_to_sky_axis(
            xs_pix, ys_pix, disk_pos + 90)

        thetas_sky_deg = np.rad2deg(np.arctan2(ys_sky_pix, xs_sky_pix))
        thetas_sky_deg[thetas_sky_deg < 0.] += 360.

        # The cos(disk_inc) factor is important to get real length on sky plane
        radii_sky_pix = ((ys_sky_pix / np.cos(np.deg2rad(disk_inc)))**2 +
                         (xs_sky_pix)**2)**0.5

        if return_cartesian:
            return radii_sky_pix, thetas_sky_deg, xs_sky_pix, ys_sky_pix
        else:
            return radii_sky_pix, thetas_sky_deg

    def azimuthal_PV_cut_from_cube(self, cube, vel_axis, r_center_pix,
                                   r_width_pix, disk_pos, disk_inc,
                                   offset_x_pix, offset_y_pix):
        """Generate azimuthal PV cut from 3D cube datat

        Args:
          cube(3D float array): 3D datacube map
          vel_axis(1D float array): velocity axis (unit: m/x)
          r_center_pix(int): central radius of ring to cut (unit: pixel)
          r_width_pix(int): width of ring to cut (unit: pixel)
          disk_pos(float): position angle of disk (unit: deg)
          disk_inc(float): inclination angle of disk (unit: deg)
          offset_x_pix(int): x-offset of disk center (unit: pixel)
          offset_y_pix(int): y-offset of disk center (unit: pixel)

        Returns:
          2D float array: theta array (shepe: #vel x #theta)
          2D float array: velocity array (shape: #vel x #theta)
          2D float array: intensity array (shape: #vel x #theta)

        """

        # Map polar coordinate to 2D map
        radii_sky_pix, thetas_sky_deg = self.set_disk_polar_to_2D_map(
            cube[0], disk_pos, disk_inc, offset_x_pix, offset_y_pix)
        x_rcw, y_rcw = np.where(
            abs(radii_sky_pix - r_center_pix) <= r_width_pix / 2.)

        # Cut azimuthally along velocity axis
        cut_theta_list, cut_velocity_list, cut_intensity_list = [], [], []
        for i, velocity in enumerate(vel_axis):
            velocity_axis_cut_theta_list = []
            velocity_axis_cut_velocity_list = []
            velocity_axis_cut_intensity_list = []
            for x, y in zip(x_rcw.T, y_rcw.T):
                velocity_axis_cut_theta_list.append(thetas_sky_deg[y, x])
                velocity_axis_cut_velocity_list.append(velocity)
                velocity_axis_cut_intensity_list.append(cube[i][y, x])
            cut_theta_list.append(velocity_axis_cut_theta_list)
            cut_velocity_list.append(velocity_axis_cut_velocity_list)
            cut_intensity_list.append(velocity_axis_cut_intensity_list)
        cut_theta_array = np.array(cut_theta_list)
        cut_velocity_array = np.array(cut_velocity_list)
        cut_intensity_array = np.array(cut_intensity_list)

        # Sort array according to theta from small to big
        sort_index = np.argsort(cut_theta_array[0])
        sort_theta_list, sort_velocity_list, sort_intensity_list = [], [], []
        for theta, velocity, intensity in zip(cut_theta_array,
                                              cut_velocity_array,
                                              cut_intensity_array):
            sort_theta, sort_velocity, sort_intensity = theta[
                sort_index], velocity[sort_index], intensity[sort_index]
            sort_theta_list.append(sort_theta)
            sort_velocity_list.append(sort_velocity)
            sort_intensity_list.append(sort_intensity)
        sort_theta_array = np.array(sort_theta_list)
        sort_velocity_array = np.array(sort_velocity_list)
        sort_intensity_array = np.array(sort_intensity_list)

        return sort_theta_array, sort_velocity_array, sort_intensity_array

    @staticmethod
    def get_peak_intensity_velocity_on_PV_diagram(velocity_array,
                                                  intensity_array):
        """Get peak intensity velocity along theta axis from PV cut array (shape: #vel,#theta)

        Args:
          velocity_array(2D float array): velocity array from PV cut result
          intensity_array(2D float array): intensity array form PV cut result

        Returns:
          peak_intensity_velocity_array(2D float array): peak intensity array

        """
        peak_intensity_velocity_list = []
        for i in range(len(intensity_array[0])):
            intensity_along_theta_array = intensity_array[:, i]
            peak_index = np.where(intensity_along_theta_array == np.max(
                intensity_along_theta_array))[0]
            peak_intensity_velocity_list.append(velocity_array[peak_index, i])
        peak_intensity_velocity_array = np.array(peak_intensity_velocity_list)
        return peak_intensity_velocity_array

    def plot_contour_PV_diagram(self,
                                theta_array,
                                velocity_array,
                                intensity_array,
                                show_peak_intensity_velocity=True,
                                show_PV_cut_data_point=True,
                                show_plot=True):
        """Plot contour azimuthal PV diagram from PV cut array (shape: #vel,#theta)

        Args:
          theta_array(2D float array): theta array from PV cut result
          velocity_array(2D float array): velocity array from PV cut result
          intensity_array(2D float array): intensity array form PV cut result
          show_peak_intensity_velocity(bool): whether show peak intensity velocity or not (Default value = True)
          show_PV_cut_data_point(bool): whether show PV cut data point or not (Default value = True)
          show_plot(bool): whether show plot or not (Default value = True)

        Returns:


        """
        X, Y, Z = theta_array, velocity_array, intensity_array

        matplotlib.rcParams.update({'font.size': 10})
        plt.figure(figsize=(18, 4))
        plt.contourf(X, Y, Z)

        peak_intensity_velocity_array = self.get_peak_intensity_velocity_on_PV_diagram(
            velocity_array, intensity_array)

        if show_PV_cut_data_point:
            plt.scatter(X,
                        Y,
                        c='black',
                        label='data_point_from_PV_cut',
                        alpha=0.1)

        if show_peak_intensity_velocity:
            plt.scatter(X[0],
                        peak_intensity_velocity_array,
                        c='red',
                        label='peak_intensity_velocity')

        # Label, Ticks and Title
        nxtick, nytick = 25, 10
        xtick_loc = np.linspace(np.amin(X[0]), np.amax(X[0]), nxtick)
        xtick_lab = np.around(
            np.linspace(np.amin(X[0]), np.amax(X[0]), nxtick, endpoint=True),
            1)
        ytick_loc = np.linspace(np.amin(Y[:, 0]),
                                np.amax(Y[:, 0]),
                                nytick,
                                endpoint=True)
        ytick_lab = np.around(
            np.linspace(np.amin(Y[:, 0]),
                        np.amax(Y[:, 0]),
                        nytick,
                        endpoint=True), 1)
        plt.xticks(xtick_loc, xtick_lab)
        plt.yticks(ytick_loc, ytick_lab)
        plt.xlabel('Theta (deg)')
        plt.ylabel('Velocity (m/s)')
        plt.title('Azimuthal Position-Velocity Cut Contour Diagram')
        plt.legend(loc=1)
        plt.grid(alpha=0.5)

        if show_plot:
            plt.show()


class AzimuthalPVDiagramTest(AzimuthalPVDiagram):
    def rotate_disk_dot_to_sky_dot_test(self,
                                        max_width=50.,
                                        test_disk_pos=120.,
                                        test_radius=30.):
        """Test rotate_disk_coord_to_sky_coord function"""

        # Axis for both disk coord and sky coord
        disk_axis = np.array([[max_width, -max_width, 0., 0.],
                              [0., 0., max_width, -max_width]])
        sky_axis = np.array(
            self.rotate_disk_dot_to_sky_dot(disk_axis[0], disk_axis[1],
                                            test_disk_pos))
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

        # Test origin of polar angle for coordinate transformation
        test_dot_disk_x, test_dot_disk_y = 0, test_radius
        test_dot_sky_x, test_dot_sky_y = self.rotate_disk_dot_to_sky_dot(
            test_dot_disk_x, test_dot_disk_y, test_disk_pos)

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
        plt.title('Rotate disk coord. to sky coord. ' +
                  '(disk_pos={:.1f}°)'.format(test_disk_pos))
        plt.axis('equal')
        plt.legend()
        plt.show()

    def set_disk_polar_to_2D_map_test(self,
                                      test_2D_map=np.ones((100, 100)),
                                      test_disk_pos=120.,
                                      test_disk_inc=47.5,
                                      test_ring_radius_sky_pix=30,
                                      test_offset_x_pix=0,
                                      test_offset_y_pix=0):
        """Test set_disk_polar_to_2D_map"""

        outputs = self.set_disk_polar_to_2D_map(test_2D_map,
                                                test_disk_pos,
                                                test_disk_inc,
                                                test_offset_x_pix,
                                                test_offset_y_pix,
                                                return_cartesian=True)
        test_radius_sky_pix = outputs[0]
        test_theta_sky = outputs[1]
        x_pixs = outputs[2]
        y_pixs = outputs[3]

        test_radius_sky_pix[
            test_radius_sky_pix <= test_ring_radius_sky_pix] = 0.

        plt.subplots(1, 4, figsize=(16, 5))
        plt.subplot(1, 4, 1)
        plt.imshow(test_radius_sky_pix, origin='lower')
        plt.title("radius\n((x/cos(inc))^2+(y)^2)^0.5\n(disk->sky)")
        plt.colorbar()
        plt.grid()
        plt.subplot(1, 4, 2)
        plt.imshow(test_theta_sky, origin='lower')
        plt.title("theta\narctan(y/x)\n(sky)")
        plt.colorbar()
        plt.grid()
        plt.subplot(1, 4, 3)
        plt.imshow(x_pixs, origin='lower')
        plt.title("x\n(disk->sky)")
        plt.colorbar()
        plt.grid()
        plt.subplot(1, 4, 4)
        plt.imshow(y_pixs, origin='lower')
        plt.title("y\n(disk->sky)")
        plt.colorbar()
        plt.grid()
        plt.suptitle(
            'Transform from disk coord. to sky coord.\n' +
            '(disk_inc={:.1f}°, disk_pos={:.1f}°, '.format(
                test_disk_inc, test_disk_pos) +
            'xoffset={:d} pixel, yoffset={:d} pixel, '.format(
                test_offset_x_pix, test_offset_y_pix) +
            'ring_radius={:d} pixel)'.format(test_ring_radius_sky_pix))
        plt.show()


def generate_azimuthal_PV_diagram(cube, vel_axis, r_center_pix, r_width_pix,
                                  disk_pos, disk_inc, offset_x_pix,
                                  offset_y_pix):
    """Generate azimuthal PV cut from 3D cube data

    Args:
      cube(3D float array): 3D datacube map
      vel_axis(1D float array): velocity axis (unit: m/x)
      r_center_pix(int): central radius of ring to cut (unit: pixel)
      r_width_pix(int): width of ring to cut (unit: pixel)
      disk_pos(float): position angle of disk (unit: deg)
      disk_inc(float): inclination angle of disk (unit: deg)
      offset_x_pix(int): x-offset of disk center (unit: pixel)
      offset_y_pix(int): y-offset of disk center (unit: pixel)

    Returns:


    """

    outputs = AzimuthalPVDiagram().azimuthal_PV_cut_from_cube(
        cube, vel_axis, r_center_pix, r_width_pix, disk_pos, disk_inc,
        offset_x_pix, offset_y_pix)

    AzimuthalPVDiagram().plot_contour_PV_diagram(
        outputs[0],
        outputs[1],
        outputs[2],
        show_peak_intensity_velocity=True,
        show_PV_cut_data_point=True,
        show_plot=True)


def main():
    """Main functions in AzimuthalPVDiagram.py"""

    # Load ALMA data
    IM_Lup = ObsData(
        '/mazu/users/jordan/PPDisk_Project/DSHARP_DR/IMLup/IMLup_CO.fits',
        158.)
    IM_Lup.stellar_property(1.12, 4250)
    IM_Lup.disk_property(47.5,
                         144.5,
                         np.array([117.]),
                         np.array([117. * 0.13]),
                         offset_x=-1.5,
                         offset_y=1.0)
    IM_Lup.planet_property(np.array([69.]), np.array([0.77]))

    # Test functions in azimuthal_PV_diagram module
    azimuthal_PV_diagram_test = AzimuthalPVDiagramTest()
    azimuthal_PV_diagram_test.rotate_disk_dot_to_sky_dot_test(
        test_disk_pos=IM_Lup.disk_pos)
    azimuthal_PV_diagram_test.set_disk_polar_to_2D_map_test(
        test_disk_pos=IM_Lup.disk_pos, test_disk_inc=IM_Lup.disk_inc)

    # Generate azimuthal PV diagram
    generate_azimuthal_PV_diagram(IM_Lup.data, IM_Lup.freqax2velax(), 20, 1,
                                  IM_Lup.disk_pos, IM_Lup.disk_inc,
                                  IM_Lup.offset_x_pix, IM_Lup.offset_y_pix)


if __name__ == '__main__':
    main()
