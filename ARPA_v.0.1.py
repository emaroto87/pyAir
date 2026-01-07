# -*- coding: utf-8 -*-
"""
Created on Wed Dec 17 11:24:51 2025

@author: U69432
"""
# ARPA

import numpy as np
import scipy as sci
from math import sin
from math import cos
from math import pi

'''
Notas:
    
    
    Y
    
    ^
    |
    |
    ++++++++++++++++++++++++++++++++++
    +                                +
    +                                +
    +                                + B
    +                                +
    +                                +
    ++++++++++++++++++++++++++++++++++  ---> X
                    L
    
    a) Dimensiones: Longitud L y anchura B
    
    b) Apoyado en sus cuatro bordes, con rigideces GJx en los bordes longitu-
    dinales y GJy en los transversales.
    
    c) Con cambios transversales (dir. y) de espesor y/o material, es decir, 
    compuestos por M tramos de anchura bm de rigidez constante (las anchuras
    y las rigideces pueden cambiar de un tramo a otro). La matriz geenralizada
    de rigidez de cada tramo vale [D*]m. El tramo "m" está situado entre
    las distancias ym1 e ym2 con respecto al borde izquierdo del panel
    
    d) Cada tramo longitudinal puede estar compuesto por capas de materiales
    diferentes, cada una de las cuales tiene un espesor tij (capa j del tramo i)
                                                
    e) El material de cada una de las capas anteriores ha de ser ortotropico 
    y elastico
    
    f) Reforda por medio de "Nr" elementos longitudinales de rigidez axial
    constante (que puede variar de un rigidizador a otro), los cuales están
    situados a distancias yrn del borde longitudilan del panel. Estos rigidi-
    zadores tienen unos modulos elasticos Ern,Grn. Sus propiedades geometricas
    son Arn, Irn, y constantes de torsion y warping Jrn y Cwrn.

    g) Sometido a cargas en el plano 
         - Longitudinal : Px
         - Carga transversal de cada tramo: Pym
         - Flujo de cortadura en cada tramo: Pxym
'''

L = 100
B = 100

i = 1  # tramo longitudinal 1
j = 1  # capa 1
m = 1  # tramo tranversal 1

d1y_list = []


def alphai(i: int, L: float) -> float:
    return i * pi / L


def alphak(k: int, L: float) -> float:
    return k * pi / L


def betaj(j: int, B: float) -> float:
    return j * pi / B


def betal(l: int, B: float) -> float:
    return l * pi / B


def integral_ss(A: float, B: float, x: float) -> float:
    if A == B:
        out = x/2 - sin(2*A*x) / (4*A)
    else:
        num1 = sin((A-B)*x)
        den1 = 2*(A-B)
        num2 = sin((A+B)*x)
        den2 = 2*(A+B)
        out = num1/den1 - num2/den2
    return out


def integral_cc(A: float, B: float, x: float) -> float:
    if A == B:
        out = x/2 + sin(2*A*x) / (4*A)
    else:
        num1 = sin((A-B)*x)
        den1 = 2*(A-B)
        num2 = sin((A+B)*x)
        den2 = 2*(A+B)
        out = num1/den1 + num2/den2
    return out


def integral_sc(A: float, B: float, x: float) -> float:
    if A == B:
        out = - ((cos(A*x))**2)/(2*A)
    else:
        num1 = cos((A-B)*x)
        den1 = 2*(A-B)
        num2 = cos((A+B)*x)
        den2 = 2*(A+B)
        out = -num1/den1 - num2/den2
    return out


def fss(A: float, B: float, x_ini: float, x_end: float) -> float:
    return integral_ss(A, B, x_end) - integral_ss(A, B, x_ini)


def fcc(A: float, B: float, x_ini: float, x_end: float) -> float:
    return integral_cc(A, B, x_end) - integral_cc(A, B, x_ini)


def fsc(A: float, B: float, x_ini: float, x_end: float) -> float:
    return integral_sc(A, B, x_end) - integral_sc(A, B, x_ini)


def fss_ik(i: int, k: int) -> float:
    a = i * pi / L
    b = k * pi / L
    return fss(a, b, 0, L)


def fcc_ik(i: int, k: int) -> float:
    a = i * pi / L
    b = k * pi / L
    return fcc(a, b, 0, L)


def fsc_ik(i: int, k: int) -> float:
    a = i * pi / L
    b = k * pi / L
    return fsc(a, b, 0, L)


def fss_jl(j: int, l: int) -> float:
    a = j * pi / L
    b = l * pi / L
    return fss(a, b, 0, B)


def fcc_jl(j: int, l: int) -> float:
    a = j * pi / B
    b = l * pi / B
    return fcc(a, b, 0, B)


def fsc_jl(j: int, l: int) -> float:
    a = j * pi / L
    b = l * pi / L
    return fsc(a, b, 0, B)


def f1jl_m(D11m_star: float, ym1: float, ym2: float, j: float, l: float) -> float:
    return D11m_star * fss(a, b, ym1, ym2)


def f1jl(d1y_list: list[tuple], j: float, l: float) -> float:
    f1jl = 0
    for D11m_star, ym1, ym2 in d1y_list:
        f1jl += f1jl_m(D11m_star, ym1, ym2, j*pi/L, l*pi/L)
    return f1jl


def f2jl_m(D22m_star: float, ym1: float, ym2: float, j: float, l: float) -> float:
    a = j * pi / B
    b = l * pi / B
    return D22m_star * fss(a, b, ym1, ym2)


def f2jl(d2y_list: list[tuple], j: float, l: float) -> float:
    f2jl = 0
    for D22m_star, ym1, ym2 in d2y_list:
        f2jl += f2jl_m(D22m_star, ym1, ym2, j, l)
    return f2jl


def f3jl_m(D33m_star: float, ym1: float, ym2: float, j: float, l: float) -> float:
    a = j * pi / B
    b = l * pi / B
    return D33m_star * fcc(a, b, ym1, ym2)


def f3jl(d3y_list: list[tuple], j: float, l: float) -> float:
    f3jl = 0
    for D33m_star, ym1, ym2 in d3y_list:
        f3jl += f3jl_m(D33m_star, ym1, ym2, j, l)
    return f3jl


def f4jl_m(D12m_star: float, ym1: float, ym2: float, j: float, l: float) -> float:
    a = j * pi / B
    b = l * pi / B
    return D12m_star * fss(a, b, ym1, ym2)


def f4jl(d4y_list: list[tuple], j: float, l: float) -> float:
    f4jl = 0
    for D12m_star, ym1, ym2 in d4y_list:
        f4jl += f4jl_m(D12m_star, ym1, ym2, j, l)
    return f4jl


def f5jl_m(D13m_star: float, ym1: float, ym2: float, j: float, l: float) -> float:
    a = j * pi / B
    b = l * pi / B
    return D13m_star * fsc(a, b, ym1, ym2)


def f5jl(d5y_list: list[tuple], j: float, l: float) -> float:
    f5jl = 0
    for D13m_star, ym1, ym2 in d5y_list:
        f5jl += f4jl_m(D13m_star, ym1, ym2, j, l)
    return f5jl


def f6jl_m(D23m_star: float, ym1: float, ym2: float, j: float, l: float) -> float:
    a = j * pi / B
    b = l * pi / B
    return D23m_star * fsc(a, b, ym1, ym2)


def f6jl(d6y_list: list[tuple], j: float, l: float) -> float:
    f6jl = 0
    for D23m_star, ym1, ym2 in d6y_list:
        f6jl += f6jl_m(D23m_star, ym1, ym2, j, l)
    return f6jl


def du1dwij_kl(i: int, j: int, k: int, l: int, nx: int, ny: int):
    alpha_i = alphai(i, L)
    alpha_k = alphak(k, L)
    alpha2_i = alpha_i**2
    alpha2_k = alpha_k**2
    beta_j = betaj(j, B)
    beta_l = betal(l, B)
    beta2_j = betaj**2
    beta2_l = betal**2
    fss_ik = fss(alpha_i, alpha_k, 0, L)
    f1_jl = f1jl(d1y_list, j, l)
    return alpha2_i * alpha2_k * fss_ik * f1_jl


def du2dwij_kl(i: int, j: int, k: int, l: int, nx: int, ny: int):
    alpha_i = alphai(i, L)
    alpha_k = alphak(k, L)
    alpha2_i = alpha_i**2
    alpha2_k = alpha_k**2
    beta_j = betaj(j, B)
    beta_l = betal(l, B)
    beta2_j = betaj**2
    beta2_l = betal**2
    fss_ik = fss(alpha_i, alpha_k, 0, L)
    f2_jl = f2jl(d2y_list, j, l)
    return beta2_j * beta2_l * beta_j * beta_l * fss_ik * f1_jl
