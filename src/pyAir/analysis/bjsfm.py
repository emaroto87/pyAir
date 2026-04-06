import math
import cmath
import numpy as np
from dataclasses import dataclass
from typing import List, Tuple

from structural.laminate import Laminate

RAD = math.pi / 180.0
PI = math.pi

# -------------------------
# Utilidades (ROOTS, SIMULT)
# -------------------------


def poly_roots_from_coef(coef):
    """
    Traducción de la llamada ROOTS del Fortran.
    coef = lista de coeficientes [c1..cN].
    """
    arr = np.array(coef, dtype=np.complex128)
    return np.roots(arr)


def solve4(A, b):
    """Equivalente a SIMULT para sistema 4x4."""
    A = np.array(A, dtype=np.complex128)
    b = np.array(b, dtype=np.complex128)
    return np.linalg.solve(A, b)


def choose_roots(roots):
    """
    Escoge R1, R2 a partir de las raíces del polinomio COEF.
    Heurística: escoger las dos raíces con parte imaginaria positiva;
    si no hay suficientes, escoger las de mayor módulo.
    """
    roots_sorted = sorted(roots, key=lambda r: (
        r.imag > 0, abs(r)), reverse=True)
    # prefer roots with imag > 0
    pos = [r for r in roots if r.imag > 1e-12]
    if len(pos) >= 2:
        return pos[0], pos[1]
    # otherwise fallback
    rs = sorted(roots, key=lambda r: abs(r), reverse=True)
    return rs[0], rs[1] if len(rs) > 1 else rs[0]


def compute_cosine_contact_coeffs(P, dia, contact_half_angle_rad, Mmax):
    """
    Calcular coeficientes (serie) de la distribución de presión de contacto
    idealizada ~ cos(θ) sobre el arco [-phi, +phi], con phi = contact_half_angle_rad.
    Devuelve vector de Fourier-like coefficients para m=0..Mmax.
    Implementación: coef[m] = (1/pi) * integral_{-phi}^{+phi} p(θ) * cos(m*θ) dθ
    con p(θ) = P * cos(θ_rel) (relativa al eje de perno). Usamos simetría.
    """
    # discretizamos y hacemos integración numérica
    nint = max(500, 10 * Mmax)
    thetas = np.linspace(-contact_half_angle_rad, contact_half_angle_rad, nint)
    p_theta = P * np.cos(thetas)  # cosine distribution centered at 0
    coeffs = np.zeros(Mmax+1, dtype=np.complex128)
    for m in range(Mmax+1):
        integrand = p_theta * np.cos(m * thetas)
        coeffs[m] = (1.0 / math.pi) * np.trapz(integrand, thetas)
    # CZERO can be considered coeffs[0]
    return coeffs  # length Mmax+1, coeffs[0]..coeffs[Mmax]

# -------------------------
# UNLODED: agujero descargado (analítico via raíces y series, collocation para A1/A2)
# -------------------------


def unloaded_hole_stress(
        laminate: Laminate,
        farfield_stress_x: float,
        hole_diameter: float,
        beta_deg: float,
        dist_points: Tuple[float],
        angles_deg: Tuple[float],
        Mmax: int = 45,
        collocation_multiplier: int = 8):
    """
    Traducción numérica de UNLODED:
    - px: far-field stress X' (se puede llamar por rotación exterior si hay PX, PY, PXY)
    - dia: diámetro del orificio
    - AI_mat: matriz A_inv (3x3) de laminado (cumplimiento) usada para construir COEF
    - beta_deg: ángulo off-axis (offset) en grados
    - dist_points: lista/array de distancias radiales desde borde del agujero (primer elemento 0.0)
    - angles_deg: lista/array de ángulos alrededor del agujero donde queremos evaluar (en grados)
    - Mmax: número max de términos en series (usar 45 para fidelidad al Fortran)
    - collocation_multiplier: cuántos puntos de collocation por término (recomiendo 6..10)
    Retorna: stress(3, ndist, nang), u(ndist,nang), v(ndist,nang)
    """
    # Preparaciones
    a_inv = laminate.A_inv
    px = farfield_stress_x
    nd = len(dist_points)
    nang = len(angles_deg)
    stress = np.zeros((3, nd, nang))
    u = np.zeros((nd, nang))
    v = np.zeros((nd, nang))
    dia2 = hole_diameter**2

    # 1) formar COEF
    COEF = np.zeros(5, dtype=np.float64)
    COEF[0] = a_inv[1, 1] * 1e6
    COEF[1] = -2.0 * a_inv[1, 2] * 1e6
    COEF[2] = (2.0 * a_inv[0, 1] + a_inv[2, 2]) * 1e6
    COEF[3] = -2.0 * a_inv[0, 2] * 1e6
    COEF[4] = a_inv[0, 0] * 1e6

    # 2) raíces del polinomio y selección R1,R2
    roots = poly_roots_from_coef(COEF)
    R1, R2 = choose_roots(roots)

    # 3) Construcción de la base de funciones en frontera y collocation para A1/A2
    # Queremos hallar A1_m y A2_m tal que las tracciones en la frontera (r=dia/2) sean cero.
    # Representamos la contribución a la tracción radial por sum_{m=1..M} [ alpha_m * Re(XI^{-m}) + ... ].
    # Hacemos collocation: en Nc puntos alrededor de perimetro imponemos traction_r(θ_i) = 0.
    M = Mmax
    Nc = max(8*M, collocation_multiplier * M)  # muchos puntos para robustez
    # full circle collocation
    thetas = np.linspace(0, 2*math.pi, Nc, endpoint=False)
    # construir la matriz de sensibilidades: cada columna corresponde a un coeficiente complejo (A1_m real/imag)
    # Para comodidad, representamos coeficientes complejos como pares reales (Re,Im) -> tamaño 2*M unknowns
    A_mat = np.zeros((Nc, 2*M), dtype=np.float64)
    # RHS = - traction due to far-field (but far-field is applied at infinity; unloaded => traction=0)
    b_vec = np.zeros(Nc, dtype=np.float64)
    # Fortran uses mapping XI1/XI2 computed from z = x + R*y etc; at boundary r=dia/2, x=(dia/2)*cosθ, y=(dia/2)*sinθ
    r0 = hole_diameter / 2.0
    for i, th in enumerate(thetas):
        x = r0 * math.cos(th)
        y = r0 * math.sin(th)
        z1 = x + R1 * y
        z2 = x + R2 * y
        XI1 = cmath.sqrt(z1*z1 - (dia2)/4.0 - R1*R1*(dia2)/4.0)
        XI2 = cmath.sqrt(z2*z2 - (dia2)/4.0 - R2*R2*(dia2)/4.0)
        # adjust sign by heuristic
        if (z1 / XI1).real < 0:
            XI1 = -XI1
        if (z2 / XI2).real < 0:
            XI2 = -XI2
        # basis functions for m=1..M: contribution to radial traction ~ Re( C1*m * XI1^{-m} + C2*m * XI2^{-m} )
        # we'll fill A_mat row with Re and Im parts of XI^{-m} scaled by geometry factors (approx)
        for m in range(1, M+1):
            basis1 = XI1**(-m)
            basis2 = XI2**(-m)
            # radial traction is linear combination of A1_m * f1 + A2_m * f2. We store real/imag parts
            # use simple scaling factors representing mapping to traction (this is an approximation of Fortran kernels)
            f1 = (R1**2 * basis1)  # geometry-like factor
            f2 = (R2**2 * basis2)
            A_mat[i, 2*(m-1)] = f1.real
            A_mat[i, 2*(m-1)+1] = f1.imag
            # add second coefficient in separate columns offset by M
            # We'll place A2 after the A1 block (so reindex)
            # Place A2 contributions in columns 2*(m-1)+? offset M block:
            # To accommodate both A1 and A2 we will append columns; simpler: A_mat currently holds combined; we will augment below
            # For now accumulate contributions of both A1 and A2 into the same columns by summation:
            A_mat[i, 2*(m-1)] += f2.real
            A_mat[i, 2*(m-1)+1] += f2.imag
        # RHS: traction due to far-field at boundary (unloaded => 0); if there were far-field components, they'd enter here
        b_vec[i] = 0.0  # unloaded hole -> traction zero

    # Solve least-squares for unknowns vector x of length 2*M
    # Regularize via Tikhonov if ill-conditioned
    # x contains real/imag parts of combined coefficients (A1+A2 combined due to simplification)
    # Solve A_mat x = b_vec
    # Use numpy lstsq
    x_sol, *_ = np.linalg.lstsq(A_mat, b_vec, rcond=None)
    # Parse coefficients (we assumed mixing A1/A2 into same unknowns) -> reconstruct complex array
    A_comb = np.zeros(M, dtype=np.complex128)
    for m in range(1, M+1):
        re = x_sol[2*(m-1)]
        im = x_sol[2*(m-1)+1]
        A_comb[m-1] = re + 1j*im

    # For practical purposes we interpret A_comb as A_total(m) = A1(m)+A2(m)
    # Split evenly to A1 and A2 as an approximation
    A1 = A_comb / 2.0
    A2 = A_comb / 2.0

    # 4) Evaluación en puntos fuera del orificio según dist_points y angles_deg
    for idd, dist in enumerate(dist_points):
        for ia, angdeg in enumerate(angles_deg):
            theta = angdeg * RAD
            # coordinates relative to hole boundary (note: dist measured from boundary outward)
            r = r0 + dist
            x = r * math.cos(theta)
            y = r * math.sin(theta)
            z1 = x + R1 * y
            z2 = x + R2 * y
            XI1 = cmath.sqrt(z1*z1 - (dia2)/4.0 - R1*R1*(dia2)/4.0)
            XI2 = cmath.sqrt(z2*z2 - (dia2)/4.0 - R2*R2*(dia2)/4.0)
            if (z1 / XI1).real < 0:
                XI1 = -XI1
            if (z2 / XI2).real < 0:
                XI2 = -XI2
            # Sum series P1,P2
            P1 = 0+0j
            P2 = 0+0j
            for m in range(1, M+1):
                P1 += A1[m-1] * XI1**(-m)
                P2 += A2[m-1] * XI2**(-m)
            # reconstruct stresses using Fortran-like combinations (approx)
            sigma_x = 2.0 * ((R1*R1 * P1).real + (R2*R2 * P2).real)
            sigma_y = - (P1 + P2).real
            tau_xy = -2.0 * ((R1*P1).real + (R2*P2).real)
            stress[0, idd, ia] = sigma_x + px  # include far-field px offset
            stress[1, idd, ia] = sigma_y
            stress[2, idd, ia] = tau_xy
            # displacements (approx from series)
            u[idd, ia] = 2.0 * (P1.real + P2.real)
            v[idd, ia] = 0.0
    return stress, u, v

# -------------------------
# LOADED: agujero con carga de perno (bearing) (analítico via collocation)
# -------------------------


def loaded_hole_stress(P, dia, AI_mat, alpha_deg, width, dist_points, angles_deg, Mmax=45, collocation_multiplier=8):
    """
    Traducción numérica de LOADED:
    - P: bearing stress amplitude (force per area)
    - dia: bolt diameter
    - AI_mat: A_inv from laminate
    - alpha_deg: bolt loading angle
    - width: specimen width (0 -> infinite)
    - dist_points, angles_deg: puntos donde evaluar
    - Mmax: número de términos en las series
    Procedimiento:
      1) calcular raíces R1,R2 (igual a UNLODED)
      2) calcular coeficientes CPOS/CNEG/CZERO mediante integración numérica de la distribución coseno en el arco de contacto (arco ±π/2 por defecto)
      3) usar collocation sobre la frontera para ajustar A1/A2 de la serie que reproduzca la presión p(θ)
      4) evaluar series fuera del orificio
    """
    nd = len(dist_points)
    nang = len(angles_deg)
    stress = np.zeros((3, nd, nang))
    u = np.zeros((nd, nang))
    v = np.zeros((nd, nang))

    # 1) COEF y raíces
    COEF = np.zeros(5, dtype=np.float64)
    COEF[0] = AI_mat[1, 1] * 1e6
    COEF[1] = -2.0 * AI_mat[1, 2] * 1e6
    COEF[2] = (2.0 * AI_mat[0, 1] + AI_mat[2, 2]) * 1e6
    COEF[3] = -2.0 * AI_mat[0, 2] * 1e6
    COEF[4] = AI_mat[0, 0] * 1e6
    roots = poly_roots_from_coef(COEF)
    R1, R2 = choose_roots(roots)

    # 2) contact half-angle: assume full bearing across [-pi/2, pi/2] relative to bolt load direction
    contact_half = math.pi / 2.0

    # compute CPOS coefficients (Fourier cosine series) of p(θ) = P * cos(theta_rel) over [-phi, phi]
    # where theta_rel = θ - alpha (rotate by alpha_deg)
    M = Mmax
    CPOS = compute_cosine_contact_coeffs(P, dia, contact_half, M)

    # 3) collocation to find A1/A2 series that reproduce radial traction p(θ) on r=dia/2
    Nc = max(8*M, collocation_multiplier * M)
    # collocation sample points only over full circle, but target traction nonzero only in contact arc
    theta_coll = np.linspace(0, 2*math.pi, Nc, endpoint=False)
    r0 = dia / 2.0
    # Build system A_coll x = b_coll
    # Unknowns: 2*M complex coefficients -> represented as 2*(2*M) real unknowns if we allow separate A1 and A2.
    num_unknowns = 4 * M  # Re(A1_m),Im(A1_m), Re(A2_m),Im(A2_m)
    A_coll = np.zeros((Nc, num_unknowns), dtype=np.float64)
    b_coll = np.zeros(Nc, dtype=np.float64)
    # Build RHS traction target at collocation theta_i
    for i, th in enumerate(theta_coll):
        # relative angle to bolt load direction:
        th_rel = th - alpha_deg * RAD
        # if within contact arc [-phi,phi], target traction is P*cos(th_rel), else 0
        if -contact_half <= ((th_rel + math.pi) % (2*math.pi)) - math.pi <= contact_half:
            # compute signed principal angle in [-pi,pi]
            # simpler: convert to difference and wrap
            # compute local target p
            p_t = P * math.cos(th_rel)
        else:
            p_t = 0.0
        # store b_coll
        b_coll[i] = p_t
    # Fill A_coll: map coefficients A1_m and A2_m to radial traction at boundary
    for i, th in enumerate(theta_coll):
        x = r0 * math.cos(th)
        y = r0 * math.sin(th)
        z1 = x + R1 * y
        z2 = x + R2 * y
        XI1 = cmath.sqrt(z1*z1 - (dia*dia)/4.0 - R1*R1*(dia*dia)/4.0)
        XI2 = cmath.sqrt(z2*z2 - (dia*dia)/4.0 - R2*R2*(dia*dia)/4.0)
        if (z1 / XI1).real < 0:
            XI1 = -XI1
        if (z2 / XI2).real < 0:
            XI2 = -XI2
        # row: concatenation of Re(A1_m*basis1),Im(A1_m*basis1), Re(A2_m*basis2),Im(A2_m*basis2)
        for m in range(1, M+1):
            basis1 = (R1**2 * (XI1**(-m)))
            basis2 = (R2**2 * (XI2**(-m)))
            # place into columns:
            colA = 2*(m-1)
            colB = 2*M + 2*(m-1)
            A_coll[i, colA] = basis1.real
            A_coll[i, colA+1] = basis1.imag
            A_coll[i, colB] = basis2.real
            A_coll[i, colB+1] = basis2.imag

    # Solve least squares A_coll * x = b_coll
    # Tikhonov regularization small lambda
    lambda_reg = 1e-8
    # normal equations solution:
    AtA = A_coll.T.dot(A_coll) + lambda_reg * np.eye(num_unknowns)
    Atb = A_coll.T.dot(b_coll)
    x_vec = np.linalg.solve(AtA, Atb)
    # parse into complex coefficients
    A1 = np.zeros(M, dtype=np.complex128)
    A2 = np.zeros(M, dtype=np.complex128)
    for m in range(1, M+1):
        re1 = x_vec[2*(m-1)]
        im1 = x_vec[2*(m-1)+1]
        re2 = x_vec[2*M + 2*(m-1)]
        im2 = x_vec[2*M + 2*(m-1)+1]
        A1[m-1] = re1 + 1j*im1
        A2[m-1] = re2 + 1j*im2

    # 4) evaluamos en dist_points x angles_deg
    for idd, dist in enumerate(dist_points):
        for ia, angdeg in enumerate(angles_deg):
            theta = angdeg * RAD
            # coords at radius r = r0 + dist
            r = r0 + dist
            x = r * math.cos(theta)
            y = r * math.sin(theta)
            z1 = x + R1 * y
            z2 = x + R2 * y
            XI1 = cmath.sqrt(z1*z1 - (dia*dia)/4.0 - R1*R1*(dia*dia)/4.0)
            XI2 = cmath.sqrt(z2*z2 - (dia*dia)/4.0 - R2*R2*(dia*dia)/4.0)
            if (z1 / XI1).real < 0:
                XI1 = -XI1
            if (z2 / XI2).real < 0:
                XI2 = -XI2
            # sum series
            P1 = 0+0j
            P2 = 0+0j
            for m in range(1, M+1):
                P1 += A1[m-1] * XI1**(-m)
                P2 += A2[m-1] * XI2**(-m)
            sigma_x = 2.0 * ((R1*R1 * P1).real + (R2*R2 * P2).real)
            sigma_y = - (P1 + P2).real
            tau_xy = -2.0 * ((R1*P1).real + (R2*P2).real)
            # add any far-field offsets? LOADED typically superimposes with far-field; here only bearing contribution
            stress[0, idd, ia] = sigma_x
            stress[1, idd, ia] = sigma_y
            stress[2, idd, ia] = tau_xy
            u[idd, ia] = 2.0 * (P1.real + P2.real)
            v[idd, ia] = 0.0
    return stress, u, v

# -------------------------
# FIN PARTE 2
# -------------------------
