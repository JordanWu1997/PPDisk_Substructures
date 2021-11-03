#!/bin/env python
# -*- coding:utf-8 -*-

# ############################################################################
# :Author: Kuan-Hsien Wu
# :Email: jordankhwu@gapp.nthu.edu.tw
# :Date: 2021-10-30
# :Description: This code is for Keplerian disk PV cut diagram
# ############################################################################

import astropy.constants as const
import matplotlib.pyplot as plt
import numpy as np

from ALMA_data_manipulation import ObsData


def generate_LOS_Keplerian_velocity_thin_disk_r(radius_pix, theta,
                                                ratio_au_pix, inclination,
                                                mass_star, sys_vel,
                                                aspect_ratio):
    """Generate LOS Keplerian velocity for thin disk model.

    Args:
      radius_pix(int): radius to calculate Keplerian velocity (unit: pixel)
      theta(float): theta to calculate Keplerian velocity (unit: deg)
      ratio_au_pix(float): ratio of au and pixel
      inclination(float): inclination of disk (unit: deg)
      mass_star(float): central star mass (unit: solar_mass)
      sys_vel(float): systematic velocity (unit: float)
      aspect_ratio(float): ratio of height of disk and radius

    Returns:
      1D float array: LOS Keplerian velocity (unit: m/s)

    """

    radius_m = ObsData.au * ratio_au_pix * radius_pix
    kep_vel_r = (ObsData.G * mass_star * ObsData.M_s / radius_m)**0.5
    vel_r_LOS = sys_vel + kep_vel_r * np.sin(np.deg2rad(inclination)) * np.cos(
        np.deg2rad(theta))
    return vel_r_LOS


def generate_LOS_Keplerian_velocity_thick_disk_r(radius_pix, theta,
                                                 ratio_au_pix, inclination,
                                                 mass_star, sys_vel,
                                                 aspect_ratio):
    """Generate LOS Keplerian velocity for thick disk model

    Args:
      radius_pix(int): radius to calculate Keplerian velocity (unit: pixel)
      theta(float): theta to calculate Keplerian velocity (unit: deg)
      ratio_au_pix(float): ratio of au and pixel
      inclination(float): inclination of disk (unit: deg)
      mass_star(float): central star mass (unit: solar_mass)
      sys_vel(float): systematic velocity (unit: float)
      aspect_ratio(float): ratio of height of disk and radius

    Returns:
      1D float array: LOS Keplerian velocity (unit: m/s)

    """

    radius_m = ObsData.au * ratio_au_pix * radius_pix
    z_m = radius_m * aspect_ratio
    kep_vel_r = (ObsData.G * mass_star * ObsData.M_s /
                 (radius_m**2 + z_m**2)**0.5)**0.5 * (
                     radius_m / (radius_m**2 + z_m**2)**0.5)**0.5
    vel_r_LOS = sys_vel + kep_vel_r * np.sin(np.deg2rad(inclination)) * np.cos(
        np.deg2rad(theta))
    return vel_r_LOS


def generate_LOS_Keplerian_velocity_thick_disk_r_z(radius_pix, theta,
                                                   ratio_au_pix, inclination,
                                                   mass_star, sys_vel,
                                                   aspect_ratio):
    """Generate LOS Keplerian velocity for thick disk model

    Args:
      radius_pix(int): radius to calculate Keplerian velocity (unit: pixel)
      theta(float): theta to calculate Keplerian velocity (unit: deg)
      ratio_au_pix(float): ratio of au and pixel
      inclination(float): inclination of disk (unit: deg)
      mass_star(float): central star mass (unit: solar_mass)
      sys_vel(float): systematic velocity (unit: float)
      aspect_ratio(float): ratio of height of disk and radius

    Returns:
      1D float array: LOS Keplerian velocity (unit: m/s)

    """

    radius_m = ObsData.au * ratio_au_pix * radius_pix
    z_m = radius_m * aspect_ratio
    kep_vel_r = (ObsData.G * mass_star * ObsData.M_s /
                 (radius_m**2 + z_m**2)**0.5)**0.5
    vel_r_LOS = sys_vel + kep_vel_r * np.sin(np.deg2rad(inclination)) * np.cos(
        np.deg2rad(theta))
    return vel_r_LOS


def generate_LOS_Keplerian_PV_cut(disk_model,
                                  velax,
                                  radius_pix,
                                  thetas,
                                  ratio_au_pix,
                                  inclination,
                                  mass_star,
                                  sys_vel,
                                  aspect_ratio,
                                  output='SI'):
    """Generate LOS Keplerian PV cut of different disk models

    Args:
      disk_model(func): model to calculate Keplerian velocity
      velax(1D float array): velocty axis from cube data
      radius_pix(int): radius to calculate Keplerian velocity (unit: pixel)
      thetas(float): theta to calculate Keplerian velocity (unit: deg)
      ratio_au_pix(float): ratio of au and pixel
      inclination(float): inclination of disk (unit: deg)
      mass_star(float): central star mass (unit: solar_mass)
      sys_vel(float): systematic velocity (unit: float)
      aspect_ratio(float): ratio of height of disk and radius
      output(str, optional): result from velocity or channel (Default value = 'SI')

    Returns:
      1D float array: PV cut from Keplerian disk model

    """

    velocity_list, vel_channel_list = [], []
    for theta in thetas:
        vel_LOS = disk_model(radius_pix, theta, ratio_au_pix, inclination,
                             mass_star, sys_vel, aspect_ratio)
        vel_LOS_channel = ObsData.vel2chan(vel_LOS, velax)
        velocity_list.append(vel_LOS)
        vel_channel_list.append(vel_LOS_channel)

    if output == 'SI':
        velocity_array = np.array(velocity_list)
        return velocity_array
    if output == 'chan':
        vel_channel_velocity_array = np.array(
            [velax[vc] for vc in vel_channel_list])
        return vel_channel_velocity_array


def plot_LOS_Kepleran_PV_cut(PPDisk,
                             radius_pix,
                             aspect_ratio,
                             theta_n=361,
                             theta_start=0.,
                             theta_end=360.,
                             plot_SI=True,
                             plot_chan=True,
                             plot_thin_disk_r=True,
                             plot_thick_disk_r=True,
                             plot_thick_disk_r_z=True,
                             color_thin_disk_r='red',
                             color_thick_disk_r='lime',
                             color_thick_disk_r_z='blue',
                             new_figure=True,
                             show_title=True,
                             show_x_axis=True,
                             show_y_axis=True,
                             show_legend=True,
                             show_plot=True):
    """Plot LOS Keplerian PV cut diagram for different models

    Args:
      PPDisk(ObsData class): Protoplanetary disk for Keplerian disk
      radius_pix(int): radius to calculate Keplerian velocity (unit: pixel)
      aspect_ratio(float): aspect ratio of Keplerian disk
      theta_n(int, optional): number of slices in theta direction (Default value = 361)
      theta_start(float, optional): start of theta direction (Default value = 0.)
      theta_end(float, optional): end of theta direction (Default value = 360.)
      plot_SI(bool, optional): plot velocity in SI output (Default value = True)
      plot_chan(bool, optional): plot velocity in chan ouput (Default value = True)
      plot_thin_disk_r(bool, optional): plot thin_disk_r model (Default value = True)
      plot_thick_disk_r(bool, optional): plot thick_disk_r model (Default value = True)
      plot_thick_disk_r_z(bool, optional): plot thick disk_r_z model (Default value = True)
      color_thin_disk_r(str, optional): color for thin_disk_r model (Default value = 'red')
      color_thick_disk_r(str, optional): color for thick_disk_r model (Default value = 'lime')
      color_thick_disk_r_z(str, optional): color for thick disk_r_z model (Default value = 'blue')
      new_figure(bool, optional): whether use new figure or not (Default value = True)
      show_plot(bool, optional): whether show plot or not (Default value = True)
      show_title(bool, optional): whether show plot or not (Default value = True)
      show_x_axis(bool, optional): whether show x axis or not (Default value = True)
      show_y_axis(bool, optional): whether show y axis or not (Default value = True)
      show_legend(bool, optional): whether show legend or not (Default value = True)

    Returns:

    """

    # Theta list
    thetas = np.linspace(theta_start, theta_end, theta_n, endpoint=True)

    # Plot
    if new_figure:
        plt.figure(figsize=(12, 6))

    # Plot disk model PV diagram
    if plot_SI:
        if plot_thin_disk_r:
            LOS_vel_thin_vel = generate_LOS_Keplerian_PV_cut(
                generate_LOS_Keplerian_velocity_thin_disk_r,
                PPDisk.freqax2velax(),
                radius_pix,
                thetas,
                PPDisk.r_au_pix,
                PPDisk.disk_inc,
                PPDisk.stellar_mass,
                PPDisk.sys_vel,
                0.,
                output='SI')
            plt.plot(thetas,
                     LOS_vel_thin_vel,
                     ls='--',
                     c=color_thin_disk_r,
                     label=r'$V_{SI}=(\frac{GM}{r^2})^\frac{1}{2}$' +
                     ', z/r={:.2f}'.format(aspect_ratio))
        if plot_thick_disk_r:
            LOS_vel_thick_vel_r = generate_LOS_Keplerian_PV_cut(
                generate_LOS_Keplerian_velocity_thick_disk_r,
                PPDisk.freqax2velax(),
                radius_pix,
                thetas,
                PPDisk.r_au_pix,
                PPDisk.disk_inc,
                PPDisk.stellar_mass,
                PPDisk.sys_vel,
                aspect_ratio,
                output='SI')
            plt.plot(
                thetas,
                LOS_vel_thick_vel_r,
                ls='--',
                c=color_thick_disk_r,
                label=
                r'$V_{r\_SI}=(\frac{GM}{r^2+z^2})^\frac{1}{2}(\frac{r}{(r^2+z^2)^{\frac{1}{2}}})^\frac{1}{2}$'
                + ', z/r={:.2f}'.format(aspect_ratio))
        if plot_thick_disk_r_z:
            LOS_vel_thick_vel_r_z = generate_LOS_Keplerian_PV_cut(
                generate_LOS_Keplerian_velocity_thick_disk_r_z,
                PPDisk.freqax2velax(),
                radius_pix,
                thetas,
                PPDisk.r_au_pix,
                PPDisk.disk_inc,
                PPDisk.stellar_mass,
                PPDisk.sys_vel,
                aspect_ratio,
                output='SI')
            plt.plot(thetas,
                     LOS_vel_thick_vel_r_z,
                     ls='--',
                     c=color_thick_disk_r_z,
                     label=r'$V_{r\_z\_SI}=(\frac{GM}{r^2+z^2})^\frac{1}{2}$' +
                     ', z/r={:.2f}'.format(aspect_ratio))
    if plot_chan:
        if plot_thin_disk_r:
            LOS_vel_thin_chan = generate_LOS_Keplerian_PV_cut(
                generate_LOS_Keplerian_velocity_thin_disk_r,
                PPDisk.freqax2velax(),
                radius_pix,
                thetas,
                PPDisk.r_au_pix,
                PPDisk.disk_inc,
                PPDisk.stellar_mass,
                PPDisk.sys_vel,
                0.,
                output='chan')
            plt.plot(thetas,
                     LOS_vel_thin_chan,
                     ls='--',
                     c=color_thin_disk_r,
                     label=r'$V_{chan}=(\frac{GM}{r^2})^\frac{1}{2}$' +
                     ', z/r={:.2f}'.format(aspect_ratio))
        if plot_thick_disk_r:
            LOS_vel_thick_chan_r = generate_LOS_Keplerian_PV_cut(
                generate_LOS_Keplerian_velocity_thick_disk_r,
                PPDisk.freqax2velax(),
                radius_pix,
                thetas,
                PPDisk.r_au_pix,
                PPDisk.disk_inc,
                PPDisk.stellar_mass,
                PPDisk.sys_vel,
                aspect_ratio,
                output='chan')
            plt.plot(
                thetas,
                LOS_vel_thick_chan_r,
                ls='--',
                c=color_thick_disk_r,
                label=
                r'$V_{r\_chan}=(\frac{GM}{r^2+z^2})^\frac{1}{2}(\frac{r}{(r^2+z^2)^{\frac{1}{2}}})^\frac{1}{2}$'
                + ', z/r={:.2f}'.format(aspect_ratio))
        if plot_thick_disk_r_z:
            LOS_vel_thick_chan_r_z = generate_LOS_Keplerian_PV_cut(
                generate_LOS_Keplerian_velocity_thick_disk_r_z,
                PPDisk.freqax2velax(),
                radius_pix,
                thetas,
                PPDisk.r_au_pix,
                PPDisk.disk_inc,
                PPDisk.stellar_mass,
                PPDisk.sys_vel,
                aspect_ratio,
                output='chan')
            plt.plot(
                thetas,
                LOS_vel_thick_chan_r_z,
                ls='--',
                c=color_thick_disk_r_z,
                label=r'$V_{r\_z\_chan}=(\frac{GM}{r^2+z^2})^\frac{1}{2}$' +
                ', z/r={:.2f}'.format(aspect_ratio))

    if show_x_axis:
        plt.xlabel('Theta (deg)')

    if show_y_axis:
        plt.ylabel('Velocity (m/s)')

    if show_legend:
        plt.legend()

    if show_title:
        plt.title(
            'Keplerian Disk Azimuthal Position-Velocity Diagram\n(r={:.1f} au ({:d} pix), star={:.1f} solar_mass, inc={:.1f} deg)'
            .format(PPDisk.pix2au(radius_pix), radius_pix, PPDisk.stellar_mass,
                    PPDisk.disk_inc))

    if show_plot:
        plt.grid()
        plt.show()


def generate_Keplerian_disk_PV_and_plot_test():
    """Test plot_LOS_Kepleran_PV_cut for different disk models"""

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

    plot_LOS_Kepleran_PV_cut(IM_Lup, 100, 1.0)


def main():
    """Main functions in generate_Keplerian_disk_PV.py"""

    generate_Keplerian_disk_PV_and_plot_test()


if __name__ == '__main__':
    main()
