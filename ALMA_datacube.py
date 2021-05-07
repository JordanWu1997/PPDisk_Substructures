#!/bin/env python
# -*- coding:utf-8 -*-
# #########################################################
# Author : Kuan-Hsien Wu
# Email : jordankhwu@gapp.nthu.edu.tw
# Description : This is code to deal with alma data cube
# Date : 2021-03-28
# #########################################################

from astropy.io import fits

import astropy.constants as const
import numpy as np
import math
import copy

class ALMA_datacube():
    '''
    Load ALMA data and transform to 3D data and header
    '''
    # Constants (SI unit)
    G      = const.G.value
    c      = const.c.value
    M_s    = const.M_sun.value

    def __init__(self, fitsfile):
        '''
        fitsfile: fitsfile location
        '''
        data = fits.getdata(fitsfile)[0,:,:,:]
        header = fits.getheader(fitsfile)
        deg_pix = float(abs(header['CDELT2']))
        asec_pix = self.Deg_to_arcsec(deg_pix)
        Bmaj_deg = float(abs(header['BMAJ']))
        Bmin_deg = float(abs(header['BMIN']))

        self.fitsfile = fitsfile
        self.data = data
        self.header = header
        self.deg_pix = deg_pix
        self.asec_pix = asec_pix
        self.Bmaj_deg = Bmaj_deg
        self.Bmin_deg = Bmin_deg

    def Get_au_pix(self, d_pc):
        '''
        Get ratio of AU and pixel
        Input:
            d_pc: distance (unit: parsec)
        '''
        au_pix = d_pc * self.asec_pix
        return au_pix

    def Print_transformation(self, d_pc):
        '''
        This is to print unit transformation for observation data
        '''
        # Print out unit transformation table
        print('{:20} = {:.5f} au'.format('1 pix', self.Get_au_pix(d_pc)))
        print('{:20} = {:.5f} asec'.format('1 pix', self.asec_pix))
        # Print beam major axis information
        print('{:20} = {:.5f} au'.format('Bmaj', self.Deg_to_arcsec(self.Bmaj_deg)/self.asec_pix*self.Get_au_pix(d_pc)))
        print('{:20} = {:.5f} asec'.format('Bmaj', self.Deg_to_arcsec(self.Bmaj_deg)))
        # Print beam minor axis information
        print('{:20} = {:.5f} au'.format('Bmin', self.Deg_to_arcsec(self.Bmin_deg)/self.asec_pix*self.Get_au_pix(d_pc)))
        print('{:20} = {:.5f} asec'.format('Bmin', self.Deg_to_arcsec(self.Bmin_deg)))

    def Get_velax(self, axis=3):
        '''
        Convert frequency axis to get velocity axis
        '''
        # Load param from header
        x0 = self.header['CRVAL{:d}'.format(axis)]
        dx = self.header['CDELT{:d}'.format(axis)]
        nu = self.header['RESTFRQ']
        # Get the velocity axis (assuming it is the third axis).
        velax = np.arange(self.header['NAXIS{:d}'.format(axis)]) \
                - self.header['CRPIX{:d}'.format(axis)] + 1.0
        velax = x0 + velax * dx
        # Convert frequency to velocity.
        velax = self.c * (nu - velax) / nu
        return velax

    @staticmethod
    def Deg_to_arcsec(deg):
        '''
        From degree to arcsec
        '''
        return deg*3600.

    # @staticmethod
    # def Deg_to_rad(deg):
        # '''
        # From Degree to radian
        # '''
        # return deg/180.*np.pi

    # @staticmethod
    # def Rad_to_deg(rad):
        # '''
        # From Degree to radian
        # '''
        # return rad/np.pi*180.

    @staticmethod
    def Convert_velocity_to_channel(v, velax):
        '''
        Convert velocity to channel map
        Input:
            v     : velocity (in m/s)
            velax : velocity axis of channel (in m/s)
        Output:
            channel : channel number
        '''
        channel = velax[(velax >= v) & (velax < v)]
        return channel
