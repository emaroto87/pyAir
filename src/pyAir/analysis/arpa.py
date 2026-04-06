# -*- coding: utf-8 -*-
"""
Created on Wed Dec 17 11:24:51 2025

@author: U69432
"""
# ARPA

import numpy as np
from scipy.linalg import eig
from math import sin
from math import cos
from math import pi

import matplotlib.pyplot as plt
from matplotlib import cm

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

'''
CONTROL DE VERSIONES
=====================
v.0.6:
-----
    - Añadimos parte de visualización en 3D con matplotlib
    - Añadidos ejemplos Ex1 al Ex4. El Ex4 que es cortadura pura no soy capaz
    de correlar los resultados.
v.0.7:
-----
    - He generado las funciones para poder ver la deformada de pandeo y estoy 
    viendo que las ondas de pandeo no están coincidiendo ni por asomo.
v.0.8:
------
    - Sigo sin poder hacer que correlen las ondas de pandeo. Además he corrobo
    rado con un FEM que lo que saca PANPA y NASTRAN es lo mismo. Por tanto el 
    problema debe estar en algo que está mal en código. Problablemente sea en
    la interpretación de la posición de los indices "i" e "j" a la hora de 
    ensamblar la matrix Q y T y por tanto la posición de los miembros del 
    vector Wij.
    
    - También he podido comprobar que el script obtiene bien los 5 primeros 
    autovalores pero a partir de ahí los resultados que proporciona empieza
    diferenciarese con respecto a los obtenidos en Nastran.
    
    - Arreglado el problema que obligaba a introducir el signo negativo
'''

__version__ = '0.8'
__release_date__ = '06/02/2026'

# Inputs
# nx = 3
# ny = 3

L = 300.0
B = 100.0

# i = 1  # tramo longitudinal 1
# j = 1  # capa 1
# m = 1  # tramo tranversal 1

# d1y_list = [(143154.44, 0, B)]
# d2y_list = [(125020.68, 0, B)]
# d3y_list = [(59793.6, 0, B)]
# d4y_list = [(56899.24, 0, B)]
# d5y_list = [(6346.81, 0, B)]
# d6y_list = [(6346.81, 0, B)]

# fxm_list = [(1000000.0, 0, B)]
# fym_list = [(0, 0, B)]
# fxym_list = [(0.0, 0, B)]


def alpha(n: int) -> float:
    return n * pi / L


def alpha2(n: int) -> float:
    x = alpha(n)
    return x**2


def beta(n: int,) -> float:
    return n * pi / B


def beta2(n: int,) -> float:
    x = beta(n)
    return x**2


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
    a = j * pi / B
    b = l * pi / B
    return fss(a, b, 0, B)


def fcc_jl(j: int, l: int) -> float:
    a = j * pi / B
    b = l * pi / B
    return fcc(a, b, 0, B)


def fsc_jl(j: int, l: int) -> float:
    a = j * pi / B
    b = l * pi / B
    return fsc(a, b, 0, B)


# PARTE I : ENERGIA INTERNA U
# =============================================================================


def f1jl_m(D11m_star: float, ym1: float, ym2: float, j: float, l: float) -> float:
    a = j * pi / B
    b = l * pi / B
    return D11m_star * fss(a, b, ym1, ym2)


def f1jl(d1y_list: list[tuple], j: float, l: float) -> float:
    f1jl = 0
    for D11m_star, ym1, ym2 in d1y_list:
        f1jl += f1jl_m(D11m_star, ym1, ym2, j, l)
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
        f5jl += f5jl_m(D13m_star, ym1, ym2, j, l)
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


def du1dwij_kl(i: int, j: int, k: int, l: int) -> float:
    return alpha2(i) * alpha2(k) * fss_ik(i, k) * f1jl(d1y_list, j, l)


def du2dwij_kl(i: int, j: int, k: int, l: int) -> float:
    return beta2(j) * beta2(l) * fss_ik(i, k) * f2jl(d2y_list, j, l)


def du3dwij_kl(i: int, j: int, k: int, l: int) -> float:
    return 4 * alpha(i) * alpha(k) * beta(j) * beta(l) * fcc_ik(i, k) * f3jl(d3y_list, j, l)


def du4dwij_kl(i: int, j: int, k: int, l: int) -> float:
    return (alpha2(i) * beta2(l) + alpha2(k) * beta2(j)) * fss_ik(i, k) * f4jl(d4y_list, j, l)


def du5dwij_kl(i: int, j: int, k: int, l: int) -> float:
    first = alpha2(i) * alpha(k) * beta(l) * \
        fsc_ik(i, k) * f5jl(d5y_list, j, l)
    # Nota: en la doc, aparece fsc_ik(k,i), puede que sea una errata
    secnd = alpha(i) * alpha2(k) * beta(j) * \
        fsc_ik(k, i) * f5jl(d5y_list, l, j)
    return -2 * (first + secnd)


def du6dwij_kl(i: int, j: int, k: int, l: int) -> float:
    first = alpha(k) * beta2(j) * beta(l) * fsc_ik(i, k) * f6jl(d6y_list, j, l)
    secnd = alpha(i) * beta(j) * beta2(l) * fsc_ik(k, i) * f6jl(d6y_list, l, j)
    return -2 * (first + secnd)


def du_dwij_kl(i: int, j: int, k: int, l: int) -> float:
    return du1dwij_kl(i, j, k, l) + du2dwij_kl(i, j, k, l) + \
        du3dwij_kl(i, j, k, l) + du4dwij_kl(i, j, k, l) + \
        du5dwij_kl(i, j, k, l) + du6dwij_kl(i, j, k, l)

# PARTE II : ENERGIA POTENCIAL V
# =============================================================================


def g1jl_m(fxm: float, ym1: float, ym2: float, j: float, l: float) -> float:
    a = j * pi / B
    b = l * pi / B
    return fxm * fss(a, b, ym1, ym2)


def g1jl(fxm_list: list[tuple], j: float, l: float) -> float:
    g1jl = 0
    for fxm, ym1, ym2 in fxm_list:
        g1jl += g1jl_m(fxm, ym1, ym2, j, l)
    return g1jl


def g2jl_m(fym: float, ym1: float, ym2: float, j: float, l: float) -> float:
    a = j * pi / B
    b = l * pi / B
    return fym * fcc(a, b, ym1, ym2)


def g2jl(fym_list: list[tuple], j: float, l: float) -> float:
    g2jl = 0
    for fym, ym1, ym2 in fym_list:
        g2jl += g2jl_m(fym, ym1, ym2, j, l)
    return g2jl


def g3jl_m(fxym: float, ym1: float, ym2: float, j: float, l: float) -> float:
    a = j * pi / B
    b = l * pi / B
    return fxym * fsc(a, b, ym1, ym2)


def g3jl(fxym_list: list[tuple], j: float, l: float) -> float:
    g3jl = 0
    for fxym, ym1, ym2 in fxym_list:
        g3jl += g3jl_m(fxym, ym1, ym2, j, l)
    return g3jl


def dv1dwij_kl(i: int, j: int, k: int, l: int) -> float:
    return alpha(i) * alpha(k) * fss_ik(i, k) * g1jl(fxm_list, j, l)


def dv2dwij_kl(i: int, j: int, k: int, l: int) -> float:
    return beta(j) * beta(l) * fss_ik(i, k) * g2jl(fym_list, j, l)


def dv3dwij_kl(i: int, j: int, k: int, l: int) -> float:
    first = alpha(i) * beta(l) * fsc_ik(k, i) * g3jl(fxym_list, j, l)
    secnd = alpha(k) * beta(j) * fsc_ik(i, k) * g3jl(fxym_list, l, j)
    return (first + secnd)


def dv_dwij_kl(i: int, j: int, k: int, l: int) -> float:
    return dv1dwij_kl(i, j, k, l) + dv2dwij_kl(i, j, k, l) + \
        dv3dwij_kl(i, j, k, l)


# PARTE III : RESOLUCION DEL PROBLEMA DE AUTOVALORES
# =============================================================================
def idx(i, j, ny):
    return (i-1)*ny + (j-1)


def QT_matrices(tol=1e-6):
    q = np.zeros((nx*ny, nx*ny))
    t = np.zeros((nx*ny, nx*ny))
    for i in range(1, nx+1):
        for j in range(1, ny+1):
            row = idx(i, j, ny)
            for k in range(1, nx+1):
                for l in range(1, ny+1):
                    col = idx(k, l, ny)
                    qij_kl = du_dwij_kl(i, j, k, l)
                    tij_kl = dv_dwij_kl(i, j, k, l)
                    if qij_kl < tol:
                        qij_kl = 0.0
                    if tij_kl < tol:
                        tik_kl = 0.0
                    q[row, col] = qij_kl
                    t[row, col] = tij_kl
                    # print(
                    #     f" {i} {j} {k} {l} Q{row}{col}={qij_kl}")
    return q, t


def solve_eig():
    q, t = QT_matrices()
    eigvalues, eigvectors = eig(q, -t)
    return eigvalues, eigvectors


def solve_eig_tikhonov(lambda_reg=1e-8):
    """
    Solve a eigen problem using the Tikhonov stabilization for ill-conditioned
    problems.

    Recommended values of lambda_reg:
        - 1e-6 -> High stability but introduces greater smoothinf
        - 1e-8 -> Balanced stability (recommended)
        - 1e-10 -> Low stability improvement but less numeric alteration.

    Returns (eigenvalues, eigenvectors)
    """
    # Obtenemos las matrices Q y T
    q, t = QT_matrices()
    # Creamos las matrices Q y T regularizadas
    n = q.shape[0]
    q_reg = q + lambda_reg * np.eye(n)
    t_reg = t + lambda_reg * np.eye(n)
    eigvalues, eigvectors = eig(q_reg, -t_reg)
    return eigvalues, eigvectors


# PARTE IV : VISUALIZACION DE RESULTADOS
# =============================================================================


def eigmode(nx: int, ny: int, wij_coeff_array: np.array, tol=1e-3):
    wij_coeff_array_real = wij_coeff_array.real

    w_max = wij_coeff_array_real[wij_coeff_array_real.argmax()]
    w_min = wij_coeff_array_real[wij_coeff_array_real.argmin()]

    if abs(w_max) >= abs(w_min):
        w_maxabs = w_max
    else:
        w_maxabs = w_min

    # print(w_maxabs)

    for i in range(0, nx+1):
        for j in range(0, ny+1):
            index = idx(i, j, ny)
            wij_coeff = wij_coeff_array_real[index]
            if abs(wij_coeff) < tol:
                wij_coeff = 0.0
            # print(index, wij_coeff)
            if wij_coeff == w_maxabs:
                mode_i = i+1
                mode_j = j+1
                main_w = wij_coeff
            else:
                pass
            index += 1
    return main_w, mode_i, mode_j


def XYZ(i: int, j: int, wij: float, L: float, B: float, points_wave=9):

    # Está comprobado, junto con la función draw. Que dibuja lo que tiene
    # dibujar.

    step_x = L/(i*(points_wave-1))
    step_y = B/(j*(points_wave-1))

    x = np.arange(0, L + step_x, step_x)
    y = np.arange(0, B + step_y, step_y)
    x, y = np.meshgrid(x, y)

    z = wij * np.sin(alpha(i)*x) * np.sin(beta(j)*y)
    return x, y, z


def draw(x: np.array, y: np. array, z: np.array, scale_factor=1.0):

    # Plot the surface
    x = scale_factor * x
    y = scale_factor * y
    z = scale_factor * z
    fig, ax = plt.subplots(subplot_kw={"projection": "3d"})
    ax.plot_surface(x, y, z, vmin=z.min() * 2, cmap=cm.Blues)
    # Alter fig size
    # fig.set_figwidth(3)
    # fig.set_figheight(1)

    # ax.set(xticklabels=[],
    #        yticklabels=[],
    #        zticklabels=[])

    plt.show()


def eigendata(nx: int, ny: int, eigenvalues: np.array, eigenvectors: np.array):

    eigenvalues_list = eigenvalues.tolist()

    for k in range(len(eigenvalues)):
        eigenvalue = eigenvalues[k]
        eigenvector = eigenvectors[k]
        main_w, mode_i, mode_j = eigmode(nx, ny, eigenvector)


# PARTE V : EJEMPLOS Y VALIDACION
# =============================================================================


def check_dus_dwij():
    for i in range(1, nx+1):
        for j in range(1, ny+1):
            wij_du1 = 0
            wij_du2 = 0
            wij_du3 = 0
            wij_du4 = 0
            wij_du5 = 0
            wij_du6 = 0
            for k in range(1, nx+1):
                for l in range(1, ny+1):
                    wij_du1 += du1dwij_kl(i, j, k, l)
                    wij_du2 += du2dwij_kl(i, j, k, l)
                    wij_du3 += du3dwij_kl(i, j, k, l)
                    wij_du4 += du4dwij_kl(i, j, k, l)
                    wij_du5 += du5dwij_kl(i, j, k, l)
                    wij_du6 += du6dwij_kl(i, j, k, l)
            wij_coef = wij_du1 + wij_du2 + wij_du3 + wij_du4 + wij_du5 + wij_du6
            print(f'Coef for term W_{(i-1)+(j-1)},{(k-1)+(l-1)} = {wij_coef}')
            print(f'\t dU1 term : {wij_du1}')
            print(f'\t dU2 term : {wij_du2}')
            print(f'\t dU3 term : {wij_du3}')
            print(f'\t dU4 term : {wij_du4}')
            print(f'\t dU5 term : {wij_du5}')
            print(f'\t dU6 term : {wij_du6}')
            _ = input('Press Enter to contine...')


def ex1_check():
    global L, B, nx, ny
    global d1y_list, d2y_list, d3y_list, d4y_list, d5y_list, d6y_list
    global fxm_list, fym_list, fxym_list

    L = 3.14
    B = 3.14
    nx = 7
    ny = 7

    d1y_list = [(143154.44, 0, B)]
    d2y_list = [(125020.68, 0, B)]
    d3y_list = [(59793.6, 0, B)]
    d4y_list = [(56899.24, 0, B)]
    d5y_list = [(6346.81, 0, B)]
    d6y_list = [(6346.81, 0, B)]

    fxm_list = [(1000000.0, 0, B)]
    fym_list = [(0, 0, B)]
    fxym_list = [(0.0, 0, B)]

    print('EXAMPLE #1:')
    print('================================================================')
    print('')
    print('Material Data :')
    print('---------------')
    print('Type : Orthotropic')
    print('Mechanical Properties:')
    print('    t = 0.125')
    print('    E11 = 160000')
    print('    E12 = 5760')
    print('    G12 = 3020')
    print('    nu12 = 0.33')
    print('')
    print('Laminate data :')
    print('---------------')
    print('Type : Monolotic')
    print('Lay-up : 45/-45/0/90/45/-45/0/90/45/-45/0/90/90/0/-45/45/90/0/-45/45/90/0/-45/45')
    print('Matrix_D = |  143154.44  56899.24  6346.81  |')
    print('           |   56899.24 125020.68  6346.81  |')
    print('           |    6346.81   6346.81  59793.6  |')
    print('')
    print('Plate Data :')
    print('---------------')
    print('thickness (t) = 3.0')
    print(f'Length (L) = {L}')
    print(f'Width (B) = {B}')
    print('')
    print('Loads:')
    print('---------------')
    print('Nxx = 1000000.0')
    print('Nyy = 0.0')
    print('Nxy = 0.0')
    print('')
    print('Solution terms:')
    print('---------------')
    print(f'nx = {nx}')
    print(f'ny = {ny}')
    print('')

    eigenvalues, eigenvectors = solve_eig()
    eigenvalues_list = [e.real for e in list(eigenvalues)]
    eigenvalues_list.sort(reverse=False)

    # Transforming from array to list
    print('')
    print('First Five minimum Eigencalues:')
    print('-------------------------------')
    print(eigenvalues_list[0:5])
    return eigenvalues_list


def ex2_check(mode_index):
    global L, B, nx, ny
    global d1y_list, d2y_list, d3y_list, d4y_list, d5y_list, d6y_list
    global fxm_list, fym_list, fxym_list

    L = 300.0
    B = 100.0
    nx = 9
    ny = 7

    d1y_list = [(143154.44, 0, B)]
    d2y_list = [(125020.68, 0, B)]
    d3y_list = [(59793.6, 0, B)]
    d4y_list = [(56899.24, 0, B)]
    d5y_list = [(6346.81, 0, B)]
    d6y_list = [(6346.81, 0, B)]

    fxm_list = [(-100.0, 0, B)]
    fym_list = [(0, 0, B)]
    fxym_list = [(0.0, 0, B)]

    print('EXAMPLE #2:')
    print('================================================================')
    print('')
    print('Material Data :')
    print('---------------')
    print('Type : Orthotropic')
    print('Mechanical Properties:')
    print('    t = 0.125')
    print('    E11 = 160000')
    print('    E12 = 5760')
    print('    G12 = 3020')
    print('    nu12 = 0.33')
    print('')
    print('Laminate data :')
    print('---------------')
    print('Type : Monolotic')
    print('Lay-up : 45/-45/0/90/45/-45/0/90/45/-45/0/90/90/0/-45/45/90/0/-45/45/90/0/-45/45')
    print('Matrix_D = |  143154.44  56899.24  6346.81  |')
    print('           |   56899.24 125020.68  6346.81  |')
    print('           |    6346.81   6346.81  59793.6  |')
    print('')
    print('Plate Data :')
    print('---------------')
    print('thickness (t) = 3.0')
    print(f'Length (L) = {L}')
    print(f'Width (B) = {B}')
    print('')
    print('Loads:')
    print('---------------')
    print('Nxx = 100.0')
    print('Nyy = 0.0')
    print('Nxy = 0.0')
    print('')
    print('Solution terms:')
    print('---------------')
    print(f'nx = {nx}')
    print(f'ny = {ny}')
    print('')

    eigenvalues, eigenvectors = solve_eig_tikhonov(1e-8)
    idx = np.argsort(eigenvalues.real)
    eigenvalues = eigenvalues[idx]
    eigenvectors = eigenvectors[:, idx]
    vec = eigenvectors / np.max(np.abs(eigenvectors))
    # eigenvalues_list = [e.real for e in list(eigenvalues)]
    # eigenvalues_list.sort(reverse=False)

    # Transforming from array to list
    print('')
    print('First Five minimum Eigencalues:')
    print('-------------------------------')
    # print(eigenvalues_list[0:5])

    # Valores Con Tikhonov lambda = 1e-8
    # Eigen1 pos = 56
    # Eigen2 pos = 48
    # Eigen3 pos = 49

    wij_array_mode = eigenvectors[:mode_index]
    print(nx, ny)
    main_w_mode, mode_i, mode_j = eigmode(nx, ny, wij_array_mode)
    print(mode_i, mode_j)
    x, y, z = XYZ(mode_i, mode_j, main_w_mode, L, B)
    draw(x, y, z)

    return eigenvalues, eigenvectors


def ex3_check():
    global L, B, nx, ny
    global d1y_list, d2y_list, d3y_list, d4y_list, d5y_list, d6y_list
    global fxm_list, fym_list, fxym_list

    L = 300.0
    B = 100.0
    nx = 9
    ny = 7

    d1y_list = [(143154.44, 0, B)]
    d2y_list = [(125020.68, 0, B)]
    d3y_list = [(59793.6, 0, B)]
    d4y_list = [(56899.24, 0, B)]
    d5y_list = [(6346.81, 0, B)]
    d6y_list = [(6346.81, 0, B)]

    fxm_list = [(0.0, 0, B)]
    fym_list = [(100.0, 0, B)]
    fxym_list = [(0.0, 0, B)]

    print('EXAMPLE #3:')
    print('================================================================')
    print('')
    print('Material Data :')
    print('---------------')
    print('Type : Orthotropic')
    print('Mechanical Properties:')
    print('    t = 0.125')
    print('    E11 = 160000')
    print('    E12 = 5760')
    print('    G12 = 3020')
    print('    nu12 = 0.33')
    print('')
    print('Laminate data :')
    print('---------------')
    print('Type : Monolotic')
    print('Lay-up : 45/-45/0/90/45/-45/0/90/45/-45/0/90/90/0/-45/45/90/0/-45/45/90/0/-45/45')
    print('Matrix_D = |  143154.44  56899.24  6346.81  |')
    print('           |   56899.24 125020.68  6346.81  |')
    print('           |    6346.81   6346.81  59793.6  |')
    print('')
    print('Plate Data :')
    print('---------------')
    print('thickness (t) = 3.0')
    print(f'Length (L) = {L}')
    print(f'Width (B) = {B}')
    print('')
    print('Loads:')
    print('---------------')
    print('Nxx = 0.0')
    print('Nyy = -100.0')
    print('Nxy = 0.0')
    print('')
    print('Solution terms:')
    print('---------------')
    print(f'nx = {nx}')
    print(f'ny = {ny}')
    print('')

    eigenvalues, eigenvectors = solve_eig()
    eigenvalues_list = [e.real for e in list(eigenvalues)]
    eigenvalues_list.sort(reverse=False)

    # Transforming from array to list
    print('')
    print('First Five minimum Eigenvalues:')
    print('-------------------------------')
    print(eigenvalues_list[0:5])
    return eigenvalues, eigenvectors, eigenvalues_list


# def ex4_check():
#     global L, B, nx, ny
#     global d1y_list, d2y_list, d3y_list, d4y_list, d5y_list, d6y_list
#     global fxm_list, fym_list, fxym_list

#     L = 300.0
#     B = 100.0
#     nx = 9
#     ny = 7

#     d1y_list = [(143154.44, 0, B)]
#     d2y_list = [(125020.68, 0, B)]
#     d3y_list = [(59793.6, 0, B)]
#     d4y_list = [(56899.24, 0, B)]
#     d5y_list = [(6346.81, 0, B)]
#     d6y_list = [(6346.81, 0, B)]

#     fxm_list = [(0.0, 0, B)]
#     fym_list = [(0.0, 0, B)]
#     fxym_list = [(100.0, 0, B)]

#     print('EXAMPLE #4:')
#     print('================================================================')
#     print('')
#     print('Material Data :')
#     print('---------------')
#     print('Type : Orthotropic')
#     print('Mechanical Properties:')
#     print('    t = 0.125')
#     print('    E11 = 160000')
#     print('    E12 = 5760')
#     print('    G12 = 3020')
#     print('    nu12 = 0.33')
#     print('')
#     print('Laminate data :')
#     print('---------------')
#     print('Type : Monolotic')
#     print('Lay-up : 45/-45/0/90/45/-45/0/90/45/-45/0/90/90/0/-45/45/90/0/-45/45/90/0/-45/45')
#     print('Matrix_D = |  143154.44  56899.24  6346.81  |')
#     print('           |   56899.24 125020.68  6346.81  |')
#     print('           |    6346.81   6346.81  59793.6  |')
#     print('')
#     print('Plate Data :')
#     print('---------------')
#     print('thickness (t) = 3.0')
#     print(f'Length (L) = {L}')
#     print(f'Width (B) = {B}')
#     print('')
#     print('Loads:')
#     print('---------------')
#     print('Nxx = 0.0')
#     print('Nyy = 0.0')
#     print('Nxy = 100.0')
#     print('')
#     print('Solution terms:')
#     print('---------------')
#     print(f'nx = {nx}')
#     print(f'ny = {ny}')
#     print('')

#     eigenvalues, eigenvectors = solve_eig()
#     eigenvalues_list = [e.real for e in list(eigenvalues)]
#     eigenvalues_list.sort(reverse=False)

#     # Transforming from array to list
#     print('')
#     print('First Five minimum Eigenvalues:')
#     print('-------------------------------')
#     print(eigenvalues_list[0:5])
#     return eigenvalues_list
