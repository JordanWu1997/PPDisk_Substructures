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

# Constants (SI unit)
G      = const.G.value
c      = const.c.value
M_s    = const.M_sun.value

def Get_ALMA_cube_data_head(data_path):
    '''
    Get ALMA data cube data and header
    Output:
        data: 3D data cube
        head: ALMA data header
    '''
    data = fits.getdata(data_path)[0,:,:,:]
    head = fits.getheader(data_path)
    return data, head

def Convert_freq_to_vel(data, head, axis=3):
    '''
    Convert frequency axis to velocity axis
    Input:
        data: 3D data cube (WI frequency axis)
        head: header of input data cube
        unit: velocity (m/s)
    '''
    # Load param from header
    x0 = head['CRVAL{:d}'.format(axis)]
    dx = head['CDELT{:d}'.format(axis)]
    nu = head['RESTFRQ']

    # Get the velocity axis (assuming it is the third axis).
    velax = np.arange(head['NAXIS{:d}'.format(axis)]) - head['CRPIX{:d}'.format(axis)] + 1.0
    velax = x0 + velax * dx

    # Convert frequency to velocity.
    velax = c * (nu - velax) / nu
    return velax

def Get_au_pix(head, d_pc):
    '''
    Get ratio of AU and pixel
    Input:
        head: header of input data cube (rval unit:deg)
        d_pc: distance (parsec)
    '''
    deg      = head['CDELT2']
    asec_deg = deg * 3600
    au_pix   = d_pc * asec_deg
    au_pix   = abs(au_pix)
    return au_pix

def Convert_pix_to_au(data_1D_pix, head, d_pc):
    '''
    Convert pixel to AU
    '''
    au_pix     = Get_au_pix(head, d_pc)
    data_1D_au = np.array([data_pix*au_pix for data_pix in data_1D_pix])
    return data_1D_au

def Convert_au_to_pixel(data_1D_au, head, d_pc):
    '''
    Convert AU to pixel
    '''
    au_pix      = Get_au_pix(head, d_pc)
    data_1D_pix = np.rint([data_au/au_pix for data_au in data_1D_au]).astype(int)
    return data_1D_pix

def Rad_to_deg(rad):
    '''
    Convert radian to degree
    '''
    deg = rad/np.pi*180
    return deg

def Deg_to_rad(deg):
    '''
    Convert degree to radian
    '''
    rad = deg/180*np.pi
    return rad

def Convert_velocity_to_channel(v, velax):
    '''
    Convert velocity to channel map
    Input:
        v     : velocity (in m/s)
        velax : velocity axis of channel (in m/s)
    Output:
        channel : channel number
    '''
    channel = np.nan
    for i in range(len(velax)-1):
        if (v < velax[i+1]) and (v > velax[i]):
            channel = i
    return channel

def Print_unit_transformation(distance_pc, au_pix, beamsize_au):
    '''
    This is to print unit transformation for observation data
    '''
    # Print out unit transformation table
    print('{:20} = {:.5f} au'.format('1 pix', au_pix))
    print('{:20} = {:.5f} asec'.format('1 pix', au_pix/distance_pc))
    # Print beamsize information
    print('{:20} = {:.5f} au'.format('beamsize', beamsize_au))
    print('{:20} = {:.5f} asec'.format('beamsize', beamsize_au/distance_pc))
    # Return ratio
    asec_pix = au_pix/distance_pc
    return au_pix, asec_pix
