#!/usr/bin/python
# -*- coding: utf-8 -*-
import numpy as np
import math

from matplotlib import pyplot as plt
from matplotlib.patches import Polygon

__version__ = "1.0.0"
__author__ = "e.maroto"
__email__ = "e.maroto.py@gmail.com"

__notes__ = "El codigo funciona usando el algoritmo <<shoelace>>. Para que algoritmo fucnione \
los segmentos deben ser definidos en el sentido antihorario. Y todos los puntos deben estar \
situados en el primer cuadrante de los ejes X-Y"

# ---------------------------------------------------------
# Funciones auxiliares para arcos
# ---------------------------------------------------------


def arc_center(p1, p2, R):
    sign = math.copysign(1, R)
    R = abs(R)
    x1, y1 = p1
    x2, y2 = p2

    # Computing length of the chord and check
    dx, dy = x2-x1, y2-y1
    c = math.hypot(dx, dy)
    if c > 2*abs(R):
        raise ValueError(
            f'There is no arc for radius {(R)} defined for section between {p1} and {p2}. '
            f'The chord is greater than the twice radius : {c} > {2*R}'
        )
    # Computing the center of the chord (M) and arc origin (O)
    mx, my = (x1+x2)/2, (y1+y2)/2
    h = math.sqrt(R**2 - (c/2)**2)
    # Computing normal direction to the chord
    nx, ny = -dy/c, dx/c
    # Computing arc center
    cx, cy = mx + sign*h*nx, my + sign*h*ny
    return cx, cy


def point_angle(p, center):
    x, y = p
    cx, cy = center
    angle = math.atan2(y-cy, x-cx)
    return angle


def normaliza_delta_t(t1, t2, R):
    sign = math.copysign(1, R)
    dt = t2 - t1
    if sign == 1 and dt < 0:
        dt += 2*math.pi
    if sign == -1 and dt > 0:
        dt -= 2*math.pi
    t1_norm = t1
    t2_norm = t1 + dt
    return t1_norm, t2_norm

# ---------------------------------------------------------
# Arc Contour Integration
# ---------------------------------------------------------


def arc_section(p1, p2, R, npoints=100):

    cx, cy = arc_center(p1, p2, R)
    t1 = point_angle(p1, (cx, cy))
    t2 = point_angle(p2, (cx, cy))
    t1_norm, t2_norm = normaliza_delta_t(t1, t2, R)

    ts = np.linspace(t1_norm, t2_norm, npoints)
    # print(ts)
    xs = cx+abs(R)*np.cos(ts)
    ys = cy+abs(R)*np.sin(ts)

    A = Qx = Qy = Ix = Iy = Ixy = 0
    points = []

    for i in range(npoints - 1):
        p1 = float(xs[i]), float(ys[i])
        p2 = float(xs[i + 1]), float(ys[i + 1])
        A_, Qx_, Qy_, Ix_, Iy_, Ixy_, line_points = line_section(p1, p2)
        A += A_
        Qx += Qx_
        Qy += Qy_
        Ix += Ix_
        Iy += Iy_
        Ixy += Ixy_
        points += line_points

    return A, Qx, Qy, Ix, Iy, Ixy, points

# ---------------------------------------------------------
# Line Contour Integration
# ---------------------------------------------------------


def line_section(p1, p2):
    x1, y1 = p1
    x2, y2 = p2
    cross = x1*y2 - x2*y1
    A = 0.5*cross
    Qx = (1/6)*((x1+x2)*cross)
    Qy = (1/6)*((y1+y2)*cross)
    Ix = (1/12)*((y1**2+y1*y2+y2**2)*cross)
    Iy = (1/12)*((x1**2+x1*x2+x2**2)*cross)
    Ixy = (1/24)*((x1*y2+2*x1*y1+2*x2*y2+x2*y1)*cross)
    points = [(p1, p2)]
    return A, Qx, Qy, Ix, Iy, Ixy, points

# ---------------------------------------------------------
# Line Contour Integration
# ---------------------------------------------------------


def plot_polygon(sections: list[tuple]):
    # Note:
    # The items of the list should be a pair of points
    # stored into a tuple.
    # Example:
    # sections = [
    #       ( (1.0,2.0) ,  (2.0,2.0) ),
    #       ( (2.0,2.0) ,  (2.0,3.0) ),
    # ]

    # Vertex extraction in counter-clockwise order
    # --------------------------------------------

    # Inital vertex of the serie
    vertices = [sections[0][0]]
    # Rest of vertex
    for s in sections:
        vertices += [s[1]]
    # Crear figura
    fig, ax = plt.subplots()

    # Dibujar polígono
    polygon = Polygon(vertices, closed=True, edgecolor="black",
                      facecolor="skyblue", alpha=0.6)
    ax.add_patch(polygon)

    # Ajustar ejes
    ax.set_aspect('equal', adjustable='box')
    ax.set_xlim(min(x for x, y in vertices) - 1,
                max(x for x, y in vertices) + 1)
    ax.set_ylim(min(y for x, y in vertices) - 1,
                max(y for x, y in vertices) + 1)

    # Mostrar
    plt.show()

# ---------------------------------------------------------
# Main function
# ---------------------------------------------------------


def section_properties(segments, npoints=100):

    A = Qx = Qy = Ix = Iy = Ixy = 0
    section_segments = []
    for s in segments:
        if len(s) == 2:
            a, qx, qy, ix, iy, ixy, points = line_section(p1=s[0], p2=s[1])
        elif len(s) == 3:
            a, qx, qy, ix, iy, ixy, points = arc_section(
                p1=s[0], p2=s[1], R=s[2], npoints=npoints)
        else:
            raise ValueError(f"Unable to identify type of section for : {s}")
        A += a
        Qx += qx
        Qy += qy
        Ix += ix
        Iy += iy
        Ixy += ixy
        section_segments += points

    Cx = Qx/A
    Cy = Qy/A
    Ix_cg = Ix - A*Cy**2
    Iy_cg = Iy - A*Cx**2
    Ixy_cg = Ixy - A*Cx*Cy
    J = Ix_cg + Iy_cg

    res = {
        "Area": A,
        "Cx": Cx,
        "Cy": Cy,
        "Ix_cg": Ix_cg,
        "Iy_cg": Iy_cg,
        "Ixy_cg": Ixy_cg,
        "J_aprox": J
    }
    for k, v in res.items():
        if abs(v) < 1e-6:
            res[k] = 0.0
        else:
            res[k] = v

    return res, section_segments


# ---------------------------------------------------------
# Ejemplos de uso
# ---------------------------------------------------------
if __name__ == "__main__":

    # # Ejemplo 1: Semicírculo de radio 1 con base
    # # --------------------------------------
    # # Ref: BDSM 6010 - Pag 22
    # # xcg = R = 1.0
    # # ycg = 4*R / (3*pi) = 0.424413
    # # Ix =  (pi/8-8/(9*pi))*R⁴ = 0.10975
    # # Iy = pi*R⁴/8 = 0.39269
    # # J = (pi/4-8/(9*pi))*R⁴= 0.5024
    # R=1.0
    # p1=(2,0)
    # p2=(0,0)
    # tramos=[
    #     (p1,p2,R),
    #     (p2,p1)
    # ]
    # print("=== Semicírculo + base ===")
    # res,seg=section_properties(tramos)
    # for k,v in res.items():
    #     print(f"{k}: {v}")
    #
    # # Ejemplo 2: Círculo completo de radio 1
    # # --------------------------------------
    # # Ref: BDSM 6010 - Pag 21
    # # Ix = Iy = pi/32 = 0.098125 (ejes de inercia)
    # # Lo definimos con dos arcos de 180°
    # p1=(2,0)
    # p2=(2,2)
    # tramos=[
    #     (p1,p2,1.0),
    #     (p2,p1,1.0)
    # ]
    # print("\n=== Círculo completo ===")
    # res,seg=section_properties(tramos)
    # for k,v in res.items():
    #     print(f"{k}: {v}")
    #
    # # Ejemplo 6: Triangulo base 2, altura2
    # # Ref: BDSM 6010 - Pag 25
    # # --------------------------------------------
    # # cx = B/2 = 1
    # # cy = H/3 = 2/3 = 0.6666
    # # Ix = BH³/36 = 16/36 = 0.4444
    # # Iy = HB³/48 = 16/48 = 0.3333
    # # J = (4BH³ + 3B³H)/124 = 112/144 = 0.777
    # # --------------------------------------------
    # # Verificado y correcto
    # p1 = (0.0,0.0)
    # p2 = (2.0,0.0)
    # p3 = (1.0,2.0)
    # tramos=[
    #     (p1,p2),
    #     (p2,p3),
    #     (p3,p1),
    # ]
    # print("\n=== Triangulo base 2, altura2 ===")
    # res,seg=section_properties(tramos)
    # for k,v in res.items():
    #     print(f"{k}: {v}")

    # Ejemplo 4 : T-section with radius 6
    # # Ref: AIDA example
    # # --------------------------------------------
    # # cx = 46.0
    # # cy = 44.5
    # # Ix = 222803.7
    # # Iy = 650045.1
    # # Ixy = 0
    # # --------------------------------------------
    # # Verificado y correcto

    p1 = (0, 56)
    p2 = (92, 56)
    p3 = (92, 46)
    p4 = (55, 46)
    p5 = (49, 40)
    p6 = (49, 0)
    p7 = (43, 0)
    p8 = (43, 40)
    p9 = (37, 46)
    p10 = (0, 46)

    tramos = [
        (p10, p9),
        (p9, p8, -6),
        (p8, p7),
        (p7, p6),
        (p6, p5),
        (p5, p4, -6),
        (p4, p3),
        (p3, p2),
        (p2, p1),
        (p1, p10)
    ]
    print("\n=== Aida Example 3===")
    res, seg = section_properties(tramos)
    for k, v in res.items():
        print(f"{k}: {v}")
    plot_polygon(seg)
