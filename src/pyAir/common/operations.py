# -*- coding: utf-8 -*-
"""
Created on Tue Feb 17 10:56:26 2026

@author: U69432
"""

import numpy as np
from math import cos
from math import sin


def rotate_2Dstress(sigmax: float, sigmay: float, tau: float, angle_rad: float) -> np.array:
    '''
    Given a rotation angle in radians, it transforms the stress vector.

    Parameters
    ----------
    sigmax : float
        Value of the axial stress in the X direction.
    sigmay : float
        Value of the axial stress in the Y direction.
    tau : float
        Value of the shear stress.
    angle_rad : float
        Angle of roation in radians.

    Returns
    -------
    out_vec : array
        Array with the rotated (sigmax, sigmay, tau)

    Validation
    ---------
    Example 1 : Rotation of -90 deg of (0.0,1.0,0.0). It must return the 
    vector (1.0, 0.0 , 0.0)

    Example 2 : Rotation of 45 deg of (0.0,0.0,1.0). It is known that a pure
    shear state rotated 45 deg should generate a tensor with one axial load in
    tension and another in compression, in other words, it should return
    (1.0, -1.0, 0.0)

    '''
    ang = angle_rad
    c = cos(ang)
    s = sin(ang)
    matrix = np.zeros((3, 3))
    matrix[0][0] = c**2
    matrix[0][1] = (s**2)
    matrix[0][2] = 2 * c * s
    matrix[1][0] = s**2
    matrix[1][1] = c**2
    matrix[1][2] = -2 * c * s
    matrix[2][0] = -1 * s * c
    matrix[2][1] = s * c
    matrix[2][2] = c**2 - s**2

    inp_vec = np.array((sigmax, sigmay, tau))

    out_vec = np.dot(matrix, inp_vec)

    return out_vec
