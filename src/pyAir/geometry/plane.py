from geometry.point import Point
from geometry.vector import Vector
import numpy as np


class Plane:
    def __init__(self):
        pass


class Plane_by3points(Plane):
    def __init__(self,
                 point1: Point,
                 point2: Point,
                 point3: Point,
                 plane_id=None,
                 inside_tol: None | float = None):
        self.point1 = point1
        self.point2 = point2
        self.point3 = point3
        self.vecdir1 = Vector(point1, point2)
        self.vecdir2 = Vector(point1, point3)
        self.plane_id = plane_id
        self.inside_tol = inside_tol

    @ property
    def normal(self):
        arr1 = self.vecdir1.to_array()
        arr2 = self.vecdir2.to_array()
        n = np.cross(arr1, arr2)
        v = Vector(Point(0, 0, 0), Point(n[0], n[1], n[2]))
        v_norm = v.unitary
        return v_norm

    @ property
    def implict_eq_coeff(self):
        n = self.normal
        a = n[0]
        b = n[1]
        c = n[2]
        x1 = self.point1[0]
        x2 = self.point1[1]
        x3 = self.point1[2]
        d = - a * x1 - b * x2 - c * x3
        return (a, b, c, d)

    def point_inplane(self, point: Point):
        v1 = Vector(self.point1, point)
        v2 = Vector(self.point1, self.point2)
        v3 = Vector(self.point2, self.point3)
        # print(v1)
        # print(v2)
        # print(v3)

        matrix = [
            v1.to_list(),
            v2.to_list(),
            v3.to_list()
        ]

        det = np.linalg.det(np.array(matrix))
        print(det)
        if det == 0.0:
            return True
        else:
            return False

    def proyect_point(self, point, direction: Vector = None) -> Point:
        if direction is None:
            n = self.normal
            d = self.implict_eq_coeff[3]
            nu = ((n @ point + d) / (n @ n)) * n
            p_prime = point - nu
            return p_prime
        else:
            n = self.normal
            v = direction
            nom = (point - self.point1) @ n
            den = n @ v
            if den == 0:
                # print('Unable to compute the proyection. The direction is '
                #       'perpendicular to the normal of the plane.')
                return None
            else:
                nu = (nom / den) * v
                p_prime = point - nu
                return p_prime

    def intriangle(self, point, verbose=False):
        """
        Check if the point complies with the following requirements.
            1. Belongs to the plane
            2. It is contained within the triangle defined by the points
               P1, P2 and P3
        """
        tol = self.inside_tol
        # edge1 = Vector(self.point1, self.point2)
        # edge2 = Vector(self.point2, self.point3)
        # edge3 = Vector(self.point3, self.point1)
        # l_edge1 = edge1.norm
        # l_edge2 = edge2.norm
        # l_edge3 = edge3.norm
        # L = max([l_edge1, l_edge2, l_edge3])
        if verbose:
            print(f'Tol : {tol}')

        v0 = Vector(self.point1, self.point2)
        v1 = Vector(self.point1, self.point3)
        v2 = Vector(self.point1, point)
        d00 = v0 @ v0
        d01 = v0 @ v1
        d11 = v1 @ v1
        d20 = v2 @ v0
        d21 = v2 @ v1
        den = d00 * d11 - d01**2

        if den == 0:
            return False  # degenerated triangle
        else:
            u = (d11 * d20 - d01 * d21) / den
            v = (d00 * d21 - d01 * d20) / den
            w = 1 - u - v
            if verbose:
                print(u)
            if verbose:
                print(v)
            if verbose:
                print(w)
            inside = (u >= -tol) and (v >= -tol) and (w >= -tol)
            on_edge = inside and (abs(u) <= tol or abs(v)
                                  <= tol or abs(w) <= tol)
            if inside:
                if verbose:
                    print('Inside')
                return True
            elif on_edge:
                if verbose:
                    print('On edge')
                return True
            else:
                if verbose:
                    print('Outside')
                return False
