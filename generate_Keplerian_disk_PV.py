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
      vel_r_LOS(1D float array): LOS Keplerian velocity (unit: m/s)

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
      vel_r_LOS(1D float array): LOS Keplerian velocity (unit: m/s)

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
      vel_r_LOS(1D float array): LOS Keplerian velocity (unit: m/s)

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
                                  output='velocity'):
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
      output(str): result from velocity or channel (Default value = 'velocity')

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

    if output == 'velocity':
        velocity_array = np.array(velocity_list)
        return velocity_array
    if output == 'channel':
        vel_channel_velocity_array = np.array(
            [velax[vc] for vc in vel_channel_list])
        return vel_channel_velocity_array


def plot_LOS_Kepleran_PV_cut(PPDisk,
                             radius_pix,
                             theta_n,
                             aspect_ratio,
                             new_figure=True,
                             show_plot=True):
    """ Plot LOS Keplerian PV cut diagram for different models

    Args:
      PPDisk(ObsData class): Protoplanetary disk for Keplerian disk
      radius_pix(int): radius to calculate Keplerian velocity (unit: pixel)
      theta_n: number of slices in theta direction
      aspect_ratio: aspect ratio of Keplerian disk
      new_figure(bool): whether use new figure or not (Default value = True)
      show_plot(bool): whether show plot or not (Default value = True)

    Returns:

    """

    # Theta list
    thetas = np.linspace(0, 360, theta_n, endpoint=True)

    # Generate Keplerian PV of different model
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
        output='velocity')
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
        output='channel')
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
        output='velocity')
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
        output='channel')
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
        output='velocity')
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
        output='channel')


    # Plot
    if new_figure:
        plt.figure(figsize=(12, 6))

    plt.plot(thetas,
             LOS_vel_thin_vel,
             ls='-',
             c='r',
             label=r'$V_{SI}=(\frac{GM}{r^2})^\frac{1}{2}$')
    plt.plot(thetas,
             LOS_vel_thin_chan,
             ls='--',
             c='r',
             label=r'$V_{chan}=(\frac{GM}{r^2})^\frac{1}{2}$')
    plt.plot(
        thetas,
        LOS_vel_thick_vel_r,
        ls='-',
        c='g',
        label=
        r'$V_{r\_SI}=(\frac{GM}{r^2+z^2})^\frac{1}{2}(\frac{r}{(r^2+z^2)^{\frac{1}{2}}})^\frac{1}{2}$'
    )
    plt.plot(
        thetas,
        LOS_vel_thick_chan_r,
        ls='--',
        c='g',
        label=
        r'$V_{r\_chan}=(\frac{GM}{r^2+z^2})^\frac{1}{2}(\frac{r}{(r^2+z^2)^{\frac{1}{2}}})^\frac{1}{2}$'
    )
    plt.plot(thetas,
             LOS_vel_thick_vel_r_z,
             ls='-',
             c='b',
             label=r'$V_{r\_z\_SI}=(\frac{GM}{r^2+z^2})^\frac{1}{2}$')
    plt.plot(thetas,
             LOS_vel_thick_chan_r_z,
             ls='--',
             c='b',
             label=r'$V_{r\_z\_chan}=(\frac{GM}{r^2+z^2})^\frac{1}{2}$')
    plt.xlabel('Theta (deg)')
    plt.ylabel('Velocity (m/s)')
    plt.title(
        'Keplerian Disk Azimuthal Position-Velocity Diagram\n(r={:.1f} au ({:d} pix), z/r={:.2f}, star={:.1f} solar_mass)'
        .format(PPDisk.pix2au(radius_pix), radius_pix, aspect_ratio,
                PPDisk.stellar_mass))
    plt.legend()

    if show_plot:
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

    plot_LOS_Kepleran_PV_cut(IM_Lup, 100, 361, 1.0)


def main():
    """Main functions in generate_Keplerian_disk_PV.py"""

    generate_Keplerian_disk_PV_and_plot_test()


if __name__ == '__main__':
    main()
