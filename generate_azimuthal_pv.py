#!/bin/env python



from load_alma_data import *
from generate_keplerian_pv import *
from scipy.interpolate import griddata

import numpy as np

def Project_to_sky_plane(x, y, pa, inc):
    '''
    Project x,y (disk coordinate) to xp,yp (sky-plane coordinate)
    Input:
        pa      : Position angle (in deg)
        inc     : Inclination angle (in deg)
        x, y    : input disk coordinate (in pixel)
    Output:
        xp, yp  : Ouput sky-plane coordinate　(in pixel)
    '''
    pa_rad, inc_rad = Deg_to_rad(pa), Deg_to_rad(inc)
    xp = x * np.cos(pa_rad) + y * np.sin(pa_rad) * np.cos(inc_rad)
    yp = x * np.sin(pa_rad) - y * np.cos(inc_rad) * np.cos(pa_rad)
    return xp, yp

def Transform_radius_to_XY(data_2D, r, pa, inc, slices=360, **kwarg):
    '''
    Transform radius to x,y coordinate (sky-plane)
    Input:
        data_2D  : 2D array
        r        : Radius (in pixel)
        pa       : Position angle (in deg)
        inc      : Inclination angle (in deg)
    Output:
        xpr, ypr : x,y  coordinate int array(in pixel)
    '''
    # Theta list
    th_list = np.linspace(0, 360, slices)

    # Center of input data (in pixel)
    empty = np.empty_like(data_2D)
    if ('cx' in kwarg) and ('cy' in kwarg):
        cx, cy = kwarg['cx'], kwarg['cy']
    else:
        cx, cy = int(len(empty[0])/2), int(len(empty)/2)

    # xpr, ypr
    xpr_list, ypr_list = [], []
    for th in th_list:
        xr, yr = r * np.cos(Deg_to_rad(th)), r * np.sin(Deg_to_rad(th))
        xpr, ypr = Project_to_sky_plane(xr, yr, pa, inc)
        xpr, ypr = xpr + cx, ypr + cy
        xpr_list.append(xpr)
        ypr_list.append(ypr)
    xpr_array = np.array(np.rint(xpr_list), dtype=int)
    ypr_array = np.array(np.rint(ypr_list), dtype=int)
    return xpr_array, ypr_array

def Generate_inner_outer_boundary(data_2D, rc, rw, pa, inc):
    '''
    Generate inner/outer boundary xy coordinate according to radius
    Input:
        data_2D : 2D array (for shape reference)
        rc      : center of radius of ring (in pixel)
        rw      : width of ring (in pixel)
        pa      : position angle (in deg)
        inc     : inclination angle (in deg)
    Output:
        bound_inner/bound_outer
    '''
    r_inner = rc - rw/2.
    r_outer = rc + rw/2.
    bound_inner = Transform_radius_to_XY(data_2D, r_inner, pa, inc)
    bound_outer = Transform_radius_to_XY(data_2D, r_outer, pa, inc)
    return bound_inner, bound_outer

def Generate_radii_and_thetas_of_2D_data(data_2D, pa, inc, **kwarg):
    '''
    Generate radius and polar angle profile of inoyt 2D data array
    Input:
        data_2D : 2D array
        pa      : Position angle (in deg)
        inc     : Inclination angle (in deg)
        xc, yc  : Center of data (in pixel)
    Output:
        rs_pix: 2D array with corresponding radius (in pixel)
        theta : 2D array with corresponding theta angle (in deg)
    '''
    # Initialization
    empty = np.empty_like(data_2D)
    # Center of input data (in pixel)
    if ('cx' in kwarg) and ('cy' in kwarg):
        cx, cy = kwarg['cx'], kwarg['cy']
    else:
        cx, cy = int(len(empty[0])/2), int(len(empty)/2)

    # Generate radius, polar angle
    x_pix  = np.array([np.arange(0, len(empty[0]), 1).astype(int) - cx])
    y_pix  = np.array([np.arange(0, len(empty), 1).astype(int) - cy])

    # Radius after project from circular to ellptical ring
    x_pixs, y_pixs   = np.meshgrid(x_pix, y_pix)
    xp_pixs, yp_pixs = Project_to_sky_plane(x_pixs, y_pixs, pa, inc)
    rp_pixs = (xp_pixs ** 2 + yp_pixs ** 2) ** 0.5
    thetas  = Rad_to_deg(np.arctan2(xp_pixs, yp_pixs))
    return rp_pixs, thetas

def Generate_non_average_PV_cut(data_3D, velax, rc, rw, pa, inc, interpolation=False, **kwarg):
    '''
    Generate non average PV_cut (Exact data pixel by pixel + Interpolation)
    Input:
        data_3D : 3D PPV data cube
        velax   : Velocity on velocity axis (frequency axis)
        rc      : center of radius of ring (in pixel)
        rw      : width of ring (in pixel)
        pa      : position angle (in deg)
        inc     : inclination angle (in deg)
    Output:
        X : 2D array of theta
        Y : 2D array of velocity
        Z : 2D array of intensity
    '''
    # Generate polar angle and velocity position coordinate in specific radius and width
    rp_pixs, thetas = Generate_radii_and_thetas_of_2D_data(data_3D[0], pa, inc, **kwarg)
    x_rcw, y_rcw = np.where(abs(rp_pixs - rc) <= rw/2.)

    # Get points from input data
    theta_list = []
    intensity_list = []
    velocity_list  = []
    for i, velocity in enumerate(velax):
        for x, y in zip(x_rcw.T, y_rcw.T):
            theta_list.append(thetas[y, x])
            intensity_list.append(data_3D[i][y, x])
            velocity_list.append(velocity)
    x = np.array(theta_list)
    y = np.array(velocity_list)
    z = np.array(intensity_list)

    if interpolation:
        # Start inteprolation from 3D sparse to 3D continuum
        ny, nx = len(velax), 720
        xmin, xmax = min(theta_list), max(theta_list)
        ymin, ymax = min(velocity_list), max(velocity_list)
        xi   = np.linspace(xmin, xmax, nx)
        yi   = np.linspace(ymin, ymax, ny)
        X, Y = np.meshgrid(xi, yi)
        Z    = griddata((x, y), z, (X, Y), method='nearest')
        return X, Y, Z
    else:
        # Reshape and Sort
        xr = x.reshape(len(velax), int(len(x)/len(velax)))
        yr = y.reshape(len(velax), int(len(y)/len(velax)))
        zr = z.reshape(len(velax), int(len(z)/len(velax)))
        xrs, yrs, zrs = [], [], []
        sort_ind = np.argsort(xr[0])
        for x, y, z in zip(xr, yr, zr):
            xs, ys, zs = x[sort_ind], y[sort_ind], z[sort_ind]
            xrs.append(xs)
            yrs.append(ys)
            zrs.append(zs)
        xrs, yrs, zrs = np.array(xrs), np.array(yrs), np.array(zrs)
        return xrs, yrs, zrs

def Get_peak_intensity_velocity(X, Y, Z):
    '''
    Input:
        X : 2D array of theta
        Y : 2D array of velocity
        Z : 2D array of intensity
    Output:
        peak_list (velocity along polar angle (theta))
    '''
    # Plot peak
    peak_list = []
    for i in range(len(Z[0])):
        peak = Y[list(Z[:,i]).index(max(Z[:,i])), i]
        peak_list.append(peak)
    return peak_list
