#!/bin/env python

from load_alma_data import *
from generate_keplerian_pv import *
import astropy.constants as const
import numpy as np

G    = const.G.value
au   = const.au.value
c    = const.c.value
M_s  = const.M_sun.value

def Generate_LOS_Keplerian_vt_thin_disk(r_pix, t, au_pix, inc, m_star, v_sys, t_origin_shift=0.):
    '''
    Generate velocity of Keplerian motion in line of sight (SI unit)
    Input:
        r_pix  : Radius (in pixel)
        t      : Theta (in deg)
        au_pix : Ratio of A.U. and pixel
        inc    : Inclination angle (in deg)
        m_star : Mass of central star (in kg)
        v_sys  : Systematic velocity (in m/s)
        t_origin_shift: shift of origin of theta direction
    Output:
        vt_obs : unit (in m/s)
    '''
    r      = au * au_pix * r_pix
    v_r    = (G * M_s / r) ** 0.5
    vt_obs = v_sys + v_r * np.sin(Deg_to_rad(inc)) * np.cos(Deg_to_rad(t - t_origin_shift))
    return vt_obs

def Generate_LOS_Keplerian_vt_thick_disk(r_pix, t, au_pix, inc, m_star, v_sys, aspect_ratio=0.1, t_origin_shift=0.):
    '''
    Generate velocity of Keplerian motion in line of sight (SI unit)
    This model consider thick gas disk (z-component is needed)
    Centrifugal force is from radial component of graviety from central star
    Z-component is derived from aspect ratio (where aspect ratio = z/r)

    Input:
        r_pix  : Radius (in pixel)
        t      : Theta (in deg)
        au_pix : Ratio of A.U. and pixel
        inc    : Inclination angle (in deg)
        m_star : Mass of central star (in kg)
        v_sys  : Systematic velocity (in m/s)
        t_origin_shift: shift of origin of theta direction
        aspect_ratio: z/r
    Output:
        vt_obs : unit (in m/s)
    '''
    r      = au * au_pix * r_pix
    z      = r * aspect_ratio
    v_r    = (G * M_s * r**2 / (r**2+z**2)**(3/2)) ** 0.5
    vt_obs = v_sys + v_r * np.sin(Deg_to_rad(inc)) * np.cos(Deg_to_rad(t - t_origin_shift))
    return vt_obs

def Generate_Keplerian_PV_cut(disk_model, velax, r_pix, ts, au_pix, pa, inc, m_star, v_sys, output='velocity',  **kwargs):
    '''
    Generate PV cut based on Keplerian motion (SI unit)
    Input:
        velax  : velocity axis of channel (in m/s)
        r_pix  : Radius (in pixel)
        ts     : Polar angle array (in deg)
        au_pix : Ratio of A.U. and pixel
        pa     : Position angle (in deg)
        inc    : Inclination angle (in deg)
        m_star : Mass of central star (in kg)
        v_sys  : Systematic velocity (in m/s)
        t_origin_shift: shift of origin of theta direction
    Output:
        t_list : velcocity channel index list correspoding to polar angle (in #channel)
    '''
    vv_list, vc_list = [], []
    for t in ts:
        vv = disk_model(r_pix, t, au_pix, inc, m_star, v_sys, **kwargs)
        vv_list.append(vv)
        vc = Convert_velocity_to_channel(vv, velax)
        vc_list.append(vc)
    if output == 'velocity':
        return vv_list
    elif output == 'channel':
        vcv_list = [velax[vc] for vc in vc_list]
        return vcv_list
