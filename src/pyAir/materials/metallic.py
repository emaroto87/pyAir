# -*- coding: utf-8 -*-

__version__ = '0.0.1'

'''
Cosas a pulir : 
    1. crear una función __repr__ y __str__
'''

import math
from dataclasses import dataclass


@dataclass
class Metallic:
    name: str = None
    norm: str = None
    shape: str = None
    temp_treat: str = None
    area_range: tuple = None
    thick_range: tuple = None
    basis: str = None
    Et: float = None
    Ec: float = None
    Ftu: float = None
    Fty: float = None
    Fcy: float = None
    Fsu: float = None
    Fbru_ed_1p5: float = None
    Fbry_ed_1p5: float = None
    Fbru_ed_2p0: float = None
    Fbry_ed_2p0: float = None
    e: float = None

    def ramberg_osgood(self):
        return

    # def __repr__(self):
    #     return None

    # def __str__(self):
    #     return None

    @property
    def ramberg_osgood_Nt(self) -> float:
        '''
        Ramberg-Osgood shape factor N in tension

        Returns
        -------
        n : float

        Misecellanius
        -------------
        - Ref : C.E. Part III eq. 1.21
        - Checked : Yes
        '''
        num = math.log10(self.e/0.002)
        den = math.log10(self.Ftu/self.Fty)
        n = num / den
        return n

    @property
    def ramberg_osgood_Nc(self) -> float:
        '''
        Ramberg-Osgood shape factor N in compresion

        Returns
        -------
        n : float

        Misecellanius
        -------------
        - Ref : C.E. Part III eq. 1.21
        - Checked : Yes
        '''
        num = math.log10(self.e/0.002)
        den = math.log10(self.Ftu/self.Fcy)
        n = num / den
        return n

    @property
    def ramberg_osgood_Ns(self) -> float:
        '''
        Ramberg-Osgood shape factor N in shear

        Returns
        -------
        n : float

        Misecellanius
        -------------
        - Ref : C.E. Part III eq. 1.29
        - Checked : Yes
        '''

        nt = self.ramberg_osgood_Nt
        nc = self.ramberg_osgood_Nc
        ns = (nt + nc)/2

        return ns

    def strain(self, stress: float) -> float:
        '''
        Strain taking into account the inelastic behaviour according to
        Ramberg-Osgood formulation.

        Returns
        -------
        n : float

        Misecellanius
        -------------
        - Ref : C.E. Part III eq. 1.18 - 1.22
        - Checked : Yes
        '''
        s = stress
        Et = self.Et
        Ec = self.Ec
        fty = self.Fty
        fcy = self.Fcy
        nt = self.ramberg_osgood_Nt
        nc = self.ramberg_osgood_Nc

        if s > 0:
            eps = (s/Et)+0.002*(s/fty)**nt
        elif s < 0:
            eps = (s/Ec)+0.002*(s/fcy)**nc
        else:
            eps = 0.0

        return eps

    def E_tan(self, stress: float) -> float:
        '''
        Tangential Young's Moduli taking into account the inelastic behaviour 
        according to Ramberg-Osgood formulation.

        Returns
        -------
        n : float

        Misecellanius
        -------------
        - Ref : C.E. Part III eq. 1.20 - 1.24
        - Checked : Yes
        '''

        s = stress
        Et = self.Et
        Ec = self.Ec
        fty = self.Fty
        fcy = self.Fcy
        nt = self.ramberg_osgood_Nt
        nc = self.ramberg_osgood_Nc

        if s > 0:
            a = Et
            b = 1 + ((0.002*nt*Et)/fty) * ((s/fty)**(nt-1))
            E_tan = a / b
        elif s < 0:
            a = Ec
            b = 1 + ((0.002*nc*Ec)/fcy) * ((s/fcy)**(nc-1))
            E_tan = a / b
        else:
            E_tan = 0.0

        return E_tan

    def E_sec(self, stress: float) -> float:
        '''
        Secant Young's Moduli taking into account the inelastic behaviour 
        according to Ramberg-Osgood formulation.

        Returns
        -------
        n : float

        Misecellanius
        -------------
        - Ref : C.E. Part III eq. 1.19 - 1.23
        - Checked : Yes
        '''

        s = stress
        Et = self.Et
        Ec = self.Ec
        fty = self.Fty
        fcy = self.Fcy
        nt = self.ramberg_osgood_Nt
        nc = self.ramberg_osgood_Nc

        if s > 0:
            a = Et
            b = 1 + ((0.002*Et)/fty) * ((s/fty)**(nt-1))
            E_sec = a / b
        elif s < 0:
            a = Ec
            b = 1 + ((0.002*Ec)/fcy) * ((s/fcy)**(nc-1))
            E_sec = a / b
        else:
            E_sec = 0.0

        return E_sec

    def plastic_region(self, stress: float) -> bool:
        s = stress
        strain = self.strain(s)

        if strain > 0.002:
            return True
        else:
            return False


def __ramberg_osgood_test():
    # Prueba Ramberg-Osgood

    Al_7075 = Metallic(name='Al_7075',
                       norm=None,
                       shape='Sheet',
                       temp_treat='T3',
                       area_range=None,
                       thick_range='0.53-1.57',
                       basis='A',
                       Et=72395,
                       Ec=73744,
                       Ftu=413.68,
                       Fty=303.37,
                       Fcy=248.21,
                       Fsu=255.10,
                       Fbru_ed_1p5=None,
                       Fbry_ed_1p5=None,
                       Fbru_ed_2p0=None,
                       Fbry_ed_2p0=None,
                       e=0.12
                       )

    return Al_7075
