#!/usr/bin/python
# -*- coding: utf-8 -*-

'''
Numerosos ensayos realizados con perfiles de pared delgada de poca longitud
sometidos a carga de compresion axial, han demostrado que el fallo de los
mismos no se produce al pandear localmente sus tramos, sino que el perfil
tiene capacidad para soportar posteriores incrementos de carga.El fallo de
un perfil de poca longitud sometido a cargas de compresion axial se denomina
 "crippling".

Los ensayos realizados muestran el siguiente comportamiento de la sección
transversal:

* Inicialmente, la carga de compresion se distribuye uniformemente sobre toda
la sección transversal.

* Aumentando la carga de compresion, se llega a un nivel de esfuerzos en el que
se produce el pandeo local de algun tramo de la sección del perfil, con lo cual
se produce una concentración de carga en las zonas adyacentes más estables de
la sección transversal, tales como intersecciones o esquinas.

* Posteriores aumentos de carga de compresioón hacen que se repita el proceso
indicado en el párrafo anterior, hata que se llegue al fallo local de la seccion
transversal cuando la int4ensidad del esfuerzo alcance un valor lo suficiente-
mente alto como para provocar la total deformación de dicha sección.

La capacidad de post-pandeo en compresión del perfil, es decir, el incremento de
carga que es capaz de soportar con posterioridad al pandeo de sus tramos, depende
básicamente de las relaciones b/t de los mismos. En efecto, si las relaciones b/t
son pequeñas, los pandeos locales se producirán a unos niveles de esfuerzos grandes,
tales que el crippling del perfil se porduce prácticamente de forma inmediata tras
el pandeo. Por el contrario, en perfiles con relaciones b/t de los tramos grandes,
el crippling aparece muy posteriormente al pandeo.

No hay una solución teorica para el calculo del esfuerzo de cripping para cualquier
tipo de perfil. Los metodods desarrollados tienen un caracter semiempirico.
'''

__version__ = '0.1'


import math
from collections import namedtuple
from collections import defaultdict
from dataclasses import dataclass, asdict

cripp_segment = namedtuple('cripp_segment', ['width', 'thick','r1','r2','folded','rect'])


@dataclass
class Metallic:
    name : str
    norm : str
    shape : str
    temp_treat : str
    area_range : tuple
    thick_range : tuple
    basis : str
    Et : float
    Ec : float
    Ftu : float
    Fty : float
    Fcy : float
    Fsu : float
    Fbru_ed_1p5 : float
    Fbry_ed_1p5 : float
    Fbru_ed_2p0 : float
    Fbry_ed_2p0 : float
    e : float
    clad : bool = False
    
    def ramberg_osgood(self):
        return
    
    @property
    def ramberg_osgood_Nt(self)->float:
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
    def ramberg_osgood_Nc(self)->float:
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
    def ramberg_osgood_Nc(self)->float:
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
        ns = ( nt + nc )/2
        
        return ns        
    
    def strain(self, stress : float)->float:
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
        else :
            eps = 0.0
        
        return eps
    
    def E_tan(self, stress : float)->float:
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
    
    def E_sec(self, stress : float)->float:
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
    
    def plastic_region(self, stress : float)->bool:
        s = stress
        strain = self.strain(s)
        
        if strain > 0.002 : 
            return True
        else: 
            return False
        



# from Materials import Metallic

def mat_correction_factor_eq(material : Metallic)->float:
    
    Fty = material.Fty
    E = material.Ec
    x = Fty/E
    coeff0 =  10.9595770137
    coeff1 = -46.3417236368
    coeff2 =  137.629157242
    coeff3 = -240.701823736
    coeff4 =  250.573444286
    coeff5 = -152.366927844
    coeff6 =  49.8640284887
    coeff7 = -6.77533482838

    K = coeff0 + (coeff1 * x) + (coeff2 * x **2) + (coeff3 * x **3) \
    + (coeff4 * x **4) + (coeff5 * x **5) + (coeff6 * x **6) + (coeff7 * x **7)
    
    return K


def mat_corrector_factor(material : Metallic ,material_base : Metallic)->float:
    Kx = mat_correction_factor_eq(material)
    Ko = mat_correction_factor_eq(material_base)

def Fcci_seg(
        Fty : float,
        E : float,
        y0 : float,
        m : float,
        b : float,
        t : float,
        K : float,
        xs : float,
        alpha : float,
        n : float
        ) -> float:
    
    beta = b/t * math.sqrt(Fty/(K*E))

    if beta <= xs:
        Fcci_seg = Fty *(y0 - m * beta)
    else:
        Fcci_seg = Fty * alpha / (beta**n)
    
    print('Beta :', beta)
    print('Fcci_seg: ',Fcci_seg)

    
    return Fcci_seg

def w_eff(segment : cripp_segment):
    '''
    Cálculo de los anchos efectivos según el manual de cálculo
    
    Parameters
    ----------
    segment : cripp_segment
        DESCRIPTION.
    type_code : 'str'
        DESCRIPTION.

    Returns
    -------
    None.

    '''
    
    w = segment.width
    t = segment.thick
    r1 = segment.r1
    r2 = segment.r2
    folded_flag = segment.folded
    rect_flag = segment.rect
    
    if r1!=0 and r2==0: 
        r = r1
    elif r1==0 and r2!=0:
        r = r2
    elif r1!=0 and r2!=0:
        r = None
    else:
        r = 0
   
    # Tramos circulares de 180 grados 
    if rect_flag == False:
        b8 = w - 0.5 * (t + r/2)
        b9 = w - (t + r/2)
        b10 = w - 0.5 * (t + r/2)
    
    # Tramos rectos
    else:
        if folded_flag:
            if r is None:

                print('Folded Supported Edge')
                'Tramo intermedio con dos radios internos r1 y r2'
                b = w - 0.5 * (t + r1/2) - 0.5 * (t + r2/2)
                K = 3.6
            else:
                print('Folded Free Edge')
                'Tramo con borde un borde libre'
                b = w - 0.5 * (t + r/2)
                K = 0.41
        else:
    
            b2 = w - 0.5 * (1 - 0.2*(r/t)**2)*t
            b3 = w - 0.5 * (1 - 0.2*(r/t)**2)*t
            b4 = w - 0.5 * (1 - 0.2*((r1 + r2)/(2*t))**2)*t
        
            b6 = w - 0.5 * (1 - 0.2*(r1/t)**2)*t - 0.5 * (1 - 0.2*(r2/t)**2)*t
            b7 = w - 0.5 * (t + r1/2) - 0.5 * (1 - 0.2*(r2/t)**2)*t
            b8 = w - 0.5 * (t + r/2)
            b9 = w - (t + r/2)
            b10 = w - 0.5 * (t + r/2)

    return b,K


def crippling_strength(material : Metallic, cripp_segments : list[cripp_segment]) -> float:
    Fty = material.Fty
    E = material.Ec
    clad = material.clad
    name = material.name
    shape = material.shape
    temp_treat = material.temp_treat
    material_base = ''

    if ('2024' in name) and (temp_treat.upper() in ['T3','T42','T3511','T42']):
        y0 = 1.4
        m = 0.628
        alpha = 0.780
        n = 1
        xs = 1.114
    elif ('7075' in name) :
        print('7075 option')
        y0 = 1.3
        m = 0.572
        alpha = 0.724
        n = 0.8
        xs = 1.007

    Pcc = 0
    Acc = 0    
    for cripp_seg in cripp_segments:
        b, K  = w_eff(cripp_seg)
        t = cripp_seg.thick
    
        fcci = Fcci_seg(Fty, E, y0, m, b, t, K, xs, alpha, n)
        Pcc += fcci * b * t
        Acc += b * t
    
    res = dict(
        Pcc = Pcc,
        Acc = Acc,
        sigma_cc = Pcc/Acc
        )

    return res


Al_7075_T6 = Metallic(
    name = 'AL_7075',
    norm = None,
    shape = 'Sheet',
    temp_treat = 'T6', 
    area_range = (None,None),
    thick_range = (None,None),
    basis = 'A',
    Et = 71020,
    Ec = 72390,
    Ftu = 489.5,
    Fty = 427.5,
    Fcy = 420.6,
    Fsu = None,
    Fbru_ed_1p5 = None,
    Fbry_ed_1p5 = None,
    Fbru_ed_2p0 = None,
    Fbry_ed_2p0 = None,
    e = 8.0
    )

seg1 = [cripp_segment(30.0,1.0,2.5,2.5,True,True)]
print(crippling_strength(Al_7075_T6, seg1))