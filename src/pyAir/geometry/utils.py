from geometry.constants import CSYS_TYPE, ANGLE_UNITS
from geometry.vector import Vector
import math


def angle_between_vectors(vector1: Vector, vector2: Vector, opt=ANGLE_UNITS.RAD) -> float:
    '''
    Returns the angle between two vectors.

    opt = ANGLE_UNITS.RAD. Returns the angle in radiands. Default value
    opt = ANGLE_UNITS.DEG. Returns the angle in degrees.
    '''
    cos = (vector1 @ vector2) / (vector1.norm * vector2.norm)
    theta_rad = math.acos(cos)
    if opt == ANGLE_UNITS.DEG:
        theta_deg = math.degrees(theta_rad)
        return theta_deg
    else:
        return theta_rad
