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
import math
from collections import namedtuple
from collections import defaultdict
cripp_segment = namedtuple('cripp_segment', ['width', 'thick','K_supp'])


from Materials import Metallic

def mat_correction_factor_eq(material : Metallic)->float:
    Fty = material.Fty
    E = material.Ec
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

def Fcci(material : Metallic, width : float, thick : float, K:float) -> float:
    Fty = material.Fty
    E = material.Ec
    clad = material.clad
    name = material.name
    shape = material.shape
    temp_treat = material.temp_treat
    material_base = ''

    if ('2024' in name) and (temp_treat.upper() in ['T3','T42','T3511','T42']:
        y0 = 1.4
        m = 0.628
        alpha = 0.780
        n = 1
        xs = 1.114
    elif ('7075' in name) :
        y0 = 1.3
        m = 0.572
        alpha = 0.724
        n = 0.8
        xs = 1.007

    b = width
    t = thick
    beta = b/t * math.sqrt(Fty/(K*E))
    if beta <= xs:
        Fcci = Fty *(y0 -m * beta)
    else:
        Fcci = Fty * alpha / (beta**n)

