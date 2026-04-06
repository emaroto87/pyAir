# -*- coding: utf-8 -*-
"""
Created on Thu Feb 19 09:42:54 2026

@author: U69432
"""
from geometry.constants import CSYS_TYPE, ANGLE_UNITS
import numpy as np
import math

__version__ = 1.0


class Point:
    def __init__(self,
                 coord1: float,
                 coord2: float,
                 coord3: float,
                 csys_type: CSYS_TYPE = CSYS_TYPE.REC,
                 point_id=None):

        self.__coords = (coord1, coord2, coord3)
        if csys_type is CSYS_TYPE.REC:
            self.x = self.__coords[0]
            self.y = self.__coords[1]
            self.z = self.__coords[2]

        elif csys_type is CSYS_TYPE.CYL:
            self.rho = self.__coords[0]
            self.phi = self.__coords[1]
            self.z = self.__coords[2]

        elif csys_type is CSYS_TYPE.SPH:
            self.rho = self.__coords[0]
            self.phi = self.__coords[1]
            self.theta = self.__coords[2]

        else:
            raise ValueError(
                'Incorrect value of csys_type argument.'
            )
        self.csys_type = csys_type
        self.point_id = point_id

    def __repr__(self):
        if self.csys_type == CSYS_TYPE.REC:
            msg = f'[x: {self.x}, y: {self.y}, z: {self.z}]'
        elif self.csys_type == CSYS_TYPE.CYL:
            msg = f'[rho: {self.rho}, phi : {self.phi}, z: {self.z}]'
        elif self.csys_type == CSYS_TYPE.SPH:
            msg = f'[rho: {self.rho}, phi : {self.phi}, theta: {self.theta}]'
        return msg

    def to_list(self):
        return list(self.__coords)

    def to_array(self):
        return np.array(self.__coords)

    def __add__(self, other):
        if isinstance(other, (list, tuple, np.ndarray, Point)):
            if len(other) == 3:
                out = Point(
                    self.__coords[0] + other[0],
                    self.__coords[1] + other[1],
                    self.__coords[2] + other[2],
                    self.csys_type
                )
                return out
            else:
                raise ValueError(
                    f'Unable to add {other} becaues the length is not equal to 3'
                )

    def __sub__(self, other):
        if isinstance(other, (list, tuple, np.ndarray, Point)):
            if len(other) == 3:
                out = Point(
                    self.__coords[0] - other[0],
                    self.__coords[1] - other[1],
                    self.__coords[2] - other[2],
                    self.csys_type
                )
                return out
            else:
                raise ValueError(
                    f'Unable to add {other} becaues the length is not equal to 3'
                )
        else:
            raise ValueError(
                f'It nos possible to substract the point with a type {type(other)}'
            )

    def __len__(self):
        return 3

    def __eq__(self, point):
        cond = [True if self.__coords[i] ==
                point.__coords[i] else False for i in range(3)]
        return all(cond)

    def __getitem__(self, index):
        return self.__coords[index]

    # def __truediv__(self, number):
    #     if n

    def __matmul__(self, other):
        if isinstance(other, (list, tuple, np.ndarray, Point)):
            if len(other) == 3:
                out = self.__coords[0] * other[0] + self.__coords[1] * \
                    other[1] + self.__coords[2] * other[2]
                return out
            else:
                raise ValueError(
                    f'Unable to add {other} becaues the length is not equal to 3'
                )
        else:
            raise ValueError(
                f'It nos possible to multiply the point with a type {type(other)}'
            )

    def __mul__(self, scalar):
        if isinstance(scalar, (int, float)):
            return Point(
                self.__coords[0] * scalar,
                self.__coords[1] * scalar,
                self.__coords[2] * scalar,
            )
        return NotImplemented

    def __rmul__(self, scalar):
        return self.__mul__(scalar)

    def __truediv__(self, scalar):
        if isinstance(scalar, (int, float)):
            return self.__mul__(1/scalar)
        return NotImplemented

    @property
    def norm(self):
        return self.__matmul__(self.__coords)**0.5

    @property
    def unitary(self):
        mod = self.norm
        return Point(
            self.__coords[0]/mod,
            self.__coords[1]/mod,
            self.__coords[2]/mod,
        )
