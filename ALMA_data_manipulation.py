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
    c = const.c.value
    M_s = const.M_sun.value

    def __init__(self, fitsfile, dist):
        data = fits.getdata(fitsfile)[0, :, :, :]
        header = fits.getheader(fitsfile)
        bmaj_deg = float(abs(header['BMAJ']))
        bmin_deg = float(abs(header['BMIN']))
        r_deg_pix = float(abs(header['CDELT2']))
        r_asec_pix = self.deg2asec(r_deg_pix)

        self.fitsfile = fitsfile
        self.data = data
        self.header = header
        self.dist = dist
        self.bmaj_deg = bmaj_deg
        self.bmin_deg = bmin_deg
        self.r_deg_pix = r_deg_pix
        self.r_asec_pix = r_asec_pix

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
        self.offset_x_pix = self.round_to_pix(offset_x * 1.e-3 / self.r_asec_pix)
        self.offset_y_pix = self.round_to_pix(offset_y * 1.e-3 / self.r_asec_pix)

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
          velocity axis (unit: m/s)

        """
        x0 = self.header['CRVAL{:d}'.format(axis)]
        dx = self.header['CDELT{:d}'.format(axis)]
        nu0 = self.header['RESTFRQ']
        # Get the velocity axis (for ALMA, it is the third axis).
        velax = np.arange(self.header['NAXIS{:d}'.format(axis)]) - self.header[
            'CRPIX{:d}'.format(axis)] + 1.0
        velax = x0 + velax * dx
        velax = self.c * (nu0 - velax) / nu0

        return velax

    @staticmethod
    def vel2chan(vel, velax):
        """Transform from velocity to channel index

        Args:
          vel(float): velocity to transfom (unit: m/s)
          velax(1D float array): velocity axis (unit: m/s)

        Returns:
          channel index

        """
        return np.where(min(velax[velax >= vel]))[0][0]

    def pix2asec(self, pix):
        """Transform from pixel number to arcsec

        Args:
          pix(int): pixel number (unit: count)

        Returns:
          arcsec

        """
        return pix * self.r_asec_pix

    def asec2pix(self, asec):
        """Transform from arcsec to pixel number

        Args:
          asec(float): arcsec (unit: arcsec)

        Returns:
          pixel number

        """
        return round(asec / self.r_asec_pix)

    @staticmethod
    def deg2asec(deg):
        """Transform from degree to arcsec

        Args:
          deg(float): degree (unit: deg)

        Returns:
          arcsec

        """
        return deg * 3600.

    @staticmethod
    def asec2au(asec, dist_pc):
        """Transform from arcsec to au

        Args:
          asec: arcsec (unit: arcsec)
          dist_pc(float): distance to observing object (unit: pc)

        Returns:
          au

        """
        return dist_pc * asec
