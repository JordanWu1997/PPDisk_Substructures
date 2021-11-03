#!/usr/bin/env python
# -*- coding:utf-8 -*-

# ############################################################################
# :Author: Kuan-Hsien Wu
# :Email: jordankhwu@gapp.nthu.edu.tw
# :Date: 2021-03-28
# :Description: This code is for ALMA observation data manipulation
# ############################################################################

import astropy.constants as const
import numpy as np
from astropy.io import fits


class ObsData():
    """Observation data property loading, storage and manipulation"""
    # Constants (SI unit)
    G = const.G.value
    au = const.au.value
    c = const.c.value
    M_s = const.M_sun.value

    def __init__(self, fitsfile, dist_pc, name='PPDisk'):
        data = fits.getdata(fitsfile)[0, :, :, :]
        header = fits.getheader(fitsfile)
        bmaj_deg = float(abs(header['BMAJ']))
        bmin_deg = float(abs(header['BMIN']))
        r_deg_pix = float(abs(header['CDELT2']))
        r_asec_pix = self.deg2asec(r_deg_pix)
        r_au_pix = self.asec2au(r_asec_pix, dist_pc)

        self.fitsfile = fitsfile
        self.name = name
        self.data = data
        self.header = header
        self.dist_pc = dist_pc
        self.bmaj_deg = bmaj_deg
        self.bmin_deg = bmin_deg
        self.r_deg_pix = r_deg_pix
        self.r_asec_pix = r_asec_pix
        self.r_au_pix = r_au_pix

    def stellar_property(self, stellar_mass, sys_vel):
        """Input stellar property for observing object

        Args:
          stellar_mass(float): mass of stellar object (unit: solar mass)
          sys_vel(float): system velocity of stellar object (unit: m/s)

        Returns:

        """
        self.stellar_mass = stellar_mass
        self.sys_vel = sys_vel

    def disk_property(self,
                      disk_inc,
                      disk_pos,
                      gap_location,
                      gap_width,
                      offset_x=0.,
                      offset_y=0.):
        """Input disk property for observing object

        Args:
          disk_inc(float): inclination angle of disk (unit: deg)
          disk_pos(float): position angle of disk (unit: deg)
          gap_location(1D float array): location of gap (unit: au)
          gap_width(1D float array): width of gap (unit: au)
          offset_x(float, optional): x-dir offset of disk center (unit:masec) (Default value = 0.)
          offset_y(float, optional): y-dir offset of disk center (unit:masec) (Default value = 0.)

        Returns:

        """
        self.disk_inc = disk_inc
        self.disk_pos = disk_pos
        self.gap_location = gap_location
        self.gap_width = gap_width
        self.offset_x = offset_x
        self.offset_y = offset_y
        self.offset_x_pix = self.round_to_pix(offset_x * 1.e-3 /
                                              self.r_asec_pix)
        self.offset_y_pix = self.round_to_pix(offset_y * 1.e-3 /
                                              self.r_asec_pix)

    def planet_property(self, planet_pos, planet_sep):
        """Input planet property for observing object

        Args:
          planet_pos(1D float array): position angle of planet (unit: deg)
          planet_sep(1D float array): separation of planet (unit: asec)

        Returns:

        """
        self.planet_pos = planet_pos
        self.planet_sep = planet_sep

    @staticmethod
    def round_to_pix(prop):
        """Round input float property to pixel number

        Args:
          prop(float): property to be rounded to pixel(integer)

        Returns:
          int: pixel number

        """
        return np.rint(prop).astype(int)

    def freqax2velax(self, axis=3):
        """Transform from frequency axis to velocity axis

        Args:
          axis(int, optional): axis of frequency (Default value = 3)

        Returns:
          float: velocity axis (unit: m/s)

        """
        x0 = self.header['CRVAL{:d}'.format(axis)]
        dx = self.header['CDELT{:d}'.format(axis)]
        nu0 = self.header['RESTFRQ']
        # Get the velocity axis (for ALMA, it is the third axis).
        velax = np.arange(self.header['NAXIS{:d}'.format(axis)]) - self.header[
            'CRPIX{:d}'.format(axis)] + 1.0
        velax = x0 + velax * dx
        velax = self.c * (nu0 - velax) / nu0
        velax = velax.T

        return velax

    def get_channel_width(self, unit='velocity', axis=3):
        """

        Args:
          unit:  (Default value = 'velocity')
          axis:  (Default value = 3)

        Returns:
          float: channel width

        """
        if unit == 'velocity':
            return abs(self.freqax2velax()[1] - self.freqax2velax()[0])
        else:
            return abs(self.header['CDELT{:d}'.format(axis)])

    @staticmethod
    def vel2chan(vel, velax):
        """Transform from velocity to channel index

        Args:
          vel(float): velocity to transfom (unit: m/s)
          velax(1D float array): velocity axis (unit: m/s)

        Returns:
          int: velocity channel index

        """
        chan_index_list = np.where(velax >= vel)[0]
        chan_bd1 = chan_index_list[0]
        chan_bd2 = chan_bd1 - 1
        dist_bd1 = abs(vel - velax[chan_bd1])
        dist_bd2 = abs(vel - velax[chan_bd2])
        if dist_bd1 > dist_bd2:
            chan_index = chan_bd2
        else:
            chan_index = chan_bd1

        return chan_index

    def pix2asec(self, pix):
        """Transform from pixel number to arcsec

        Args:
          pix(int): pixel number (unit: count)

        Returns:
          float: arcsec

        """
        return pix * self.r_asec_pix

    def pix2au(self, pix):
        """Transform from pixel number to au

        Args:
          pix(int): pixel number (unit: count)

        Returns:
          float: au

        """
        return pix * self.r_au_pix

    def au2pix(self, au):
        """Transform from pixel number to au

        Args:
          au(float): au

        Returns:
          int: pixel number (unit: count)

        """
        return au / self.r_au_pix

    def asec2pix(self, asec):
        """Transform from arcsec to pixel number

        Args:
          asec(float): arcsec (unit: arcsec)

        Returns:
          int: pixel number

        """
        return self.round_to_pix(asec / self.r_asec_pix)

    @staticmethod
    def deg2asec(deg):
        """Transform from degree to arcsec

        Args:
          deg(float): degree (unit: deg)

        Returns:
          float: arcsec

        """
        return deg * 3600.

    @staticmethod
    def asec2au(asec, dist_pc):
        """Transform from arcsec to au

        Args:
          asec(float): arcsec (unit: arcsec)
          dist_pc(float): distance to observing object (unit: pc)

        Returns:
          float: au

        """
        return dist_pc * asec


def main():
    """Main functions in ALMA_data_manipulation.py"""

    # Load ALMA data
    IM_Lup = ObsData(
        '/mazu/users/jordan/PPDisk_Project/DSHARP_DR/IMLup/IMLup_CO.fits',
        158.,
        name='IM_Lup')
    IM_Lup.stellar_property(1.12, 4250)
    IM_Lup.disk_property(47.5,
                         144.5,
                         np.array([117.]),
                         np.array([117. * 0.13]),
                         offset_x=-1.5,
                         offset_y=1.0)
    IM_Lup.planet_property(np.array([69.]), np.array([0.77]))


if __name__ == '__main__':
    main()
