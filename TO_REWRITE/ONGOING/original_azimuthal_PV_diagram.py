def Generate_radii_and_thetas_of_2D_data(data_2D, pa, inc, **kwarg):
    """Generate radius and polar angle profile of input 2D data array
    Input:
        data_2D : 2D array
        pa      : Position angle (in deg)
        inc     : Inclination angle (in deg)
        xc, yc  : Center of data (in pixel)
    Output:
        rs_pix: 2D array with corresponding radius (in pixel)
        theta : 2D array with corresponding theta angle (in deg)

    Args:
      data_2D: param pa:
      inc: param **kwarg:
      pa:
      **kwarg:

    Returns:

    """
    # Initialization
    empty = np.empty_like(data_2D)
    # Center of input data (in pixel)
    if ('cx' in kwarg) and ('cy' in kwarg):
        cx, cy = kwarg['cx'], kwarg['cy']
    else:
        cx, cy = int(len(empty[0])/2), int(len(empty)/2)

    # Generate radius, polar angle
    x_pix = np.array([np.arange(0, len(empty[0]), 1).astype(int) - cx])
    y_pix = np.array([np.arange(0, len(empty), 1).astype(int) - cy])

    # Radius after project from circular to ellptical ring
    x_pixs, y_pixs   = np.meshgrid(x_pix, y_pix)
    xp_pixs, yp_pixs = Project_to_sky_plane(x_pixs, y_pixs, pa, inc)
    rp_pixs = (xp_pixs ** 2 + yp_pixs ** 2) ** 0.5

    # Note: For normal arctan, it should be y/x, but here I want
    thetas  = np.rad2deg(np.arctan2(xp_pixs, yp_pixs))
    return rp_pixs, thetas

def Generate_azimuthal_cut(data_2D, rc_cut_pix, rw_cut_pix, pa, inc, **kwarg):
    """Generate azimuthal cut of 2D data array
    Input:
        data_2D : 2D array
        pa      : Position angle (in deg)
        inc     : Inclination angle (in deg)
        xc, yc  : Center of data (in pixel)
    Output:
        theta_list : 1D list with corresponding theta angle (in deg)
        cut_list   : 1D list with corresponding cut value

    Args:
      data_2D: param rc_cut_pix:
      rw_cut_pix: param pa:
      inc: param **kwarg:
      rc_cut_pix:
      pa:
      **kwarg:

    Returns:

    """
    # Initialization
    empty = np.empty_like(data_2D)
    # Center of input data (in pixel)
    if ('cx' in kwarg) and ('cy' in kwarg):
        cx, cy = kwarg['cx'], kwarg['cy']
    else:
        cx, cy = int(len(empty[0])/2), int(len(empty)/2)
    # Generate rs, thetas
    rp_pixs, thetas = Generate_radii_and_thetas_of_2D_data(data_2D, pa, inc,
                                                           cx=cx, cy=cy)
    # Azimuthal cut
    theta_list, cut_list = [], []
    x_rcw, y_rcw = np.where(abs(rp_pixs - rc_cut_pix) <= rw_cut_pix/2.)
    for x, y in zip(x_rcw.T, y_rcw.T):
        theta_list.append(thetas[y, x])
        cut_list.append(data_2D[y, x])
    return theta_list, cut_list

def Generate_non_average_PV_cut(data_3D, velax, rc, rw, pa, inc, interpolation=False, **kwarg):
    """Generate non average PV_cut (Exact data pixel by pixel + Interpolation)
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

    Args:
      data_3D: param velax:
      rc: param rw:
      pa: param inc:
      interpolation: Default value = False)
      velax:
      rw:
      inc:
      **kwarg:

    Returns:

    """
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
    """Input:
        X : 2D array of theta
        Y : 2D array of velocity
        Z : 2D array of intensity
    Output:
        peak_list (velocity along polar angle (theta))

    Args:
      X: param Y:
      Z:
      Y:

    Returns:

    """
    # Plot peak
    peak_list = []
    for i in range(len(Z[0])):
        peak = Y[list(Z[:, i]).index(max(Z[:, i])), i]
        peak_list.append(peak)
    return peak_list
