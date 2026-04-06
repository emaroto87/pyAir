from geometry.constants import CSYS_TYPE, ANGLE_UNITS
from geometry.vector import Vector
from geometry.point import Point
import numpy as np
import math


class Csys:
    def __init__(self,
                 origin: Point,
                 dir1: Vector,
                 dir2: Vector,
                 dir3: Vector,
                 csys_type: CSYS_TYPE = CSYS_TYPE.REC):

        self.origin = origin
        self.dir1 = dir1
        self.dir2 = dir2
        self.dir3 = dir3
        self.type = csys_type

    @classmethod
    def create_REC_by_3points(cls,
                              origin: Point,
                              point_on_axis: Point,
                              point_on_plane: Point) -> "Csys":
        '''
        Creates a rectangular coordinate system using 3 points.

        Parameters
        ----------
        origin : Point
            Point of the origin of the Csys.
        point_on_axis : Point
            Point on the X-axis of the Csys
        point_on_plane : Point
            Point on the XY plane.

        '''
        v1 = Vector(origin, point_on_axis).to_array()
        v2 = Vector(origin, point_on_plane).to_array()
        n = np.cross(v1, v2)
        v3 = np.cross(n, v1)
        dir1 = v1
        dir2 = Vector(Point(0, 0, 0), Point(*v3))
        dir3 = Vector(Point(0, 0, 0), Point(*n))

        return cls(origin, dir1, dir2, dir3, csys_type=CSYS_TYPE.REC)

    def rotate_by_vector(self, vector: Vector) -> "Csys":
        '''
        Given a vector (v), projects the vector into the XY plane (vp) and 
        rotates the system along the Z axis so the the X-axis matches the 
        projected vector vp.

        '''
        o = self.origin
        v = vector
        n = self.dir3
        vp = v - (v@n)*n
        v1 = vp.to_array()
        v2 = np.cross(n, v1)
        dir1 = Vector(Point(0, 0, 0), Point(*v1))
        dir2 = Vector(Point(0, 0, 0), Point(*v2))
        dir3 = n

        return Csys(o, dir1, dir2, dir3, csys_type=CSYS_TYPE.REC)

    def rotate_by_csys(self, csys: "Csys") -> "Csys":
        v = csys.dir1
        return self.rotate_by_vector(v)

    def rotate_by_angle(self, angle_deg: float) -> "Csys":
        angle = math.radians(angle_deg)
        c = math.cos(angle)
        s = math.sin(angle)
        rot_matrix = np.array([
            [c, -s, 0.0],
            [s, c, 0.0],
            [0.0, 0.0, 1.0]
        ])
        origin = self.origin
        dir1_rot = rot_matrix @ self.dir1
        dir2_rot = rot_matrix @ self.dir2
        dir3_rot = rot_matrix @ self.dir3
        dir1_rot_v = Vector(Point(0, 0, 0), Point(*dir1_rot))
        dir2_rot_v = Vector(Point(0, 0, 0), Point(*dir2_rot))
        dir3_rot_v = Vector(Point(0, 0, 0), Point(*dir3_rot))
        csys_type = self.type

        return Csys(origin, dir1_rot_v, dir2_rot_v, dir3_rot_v, csys_type)
