from typing import Optional
import numpy as np
from structural.laminate import Laminate
from structural.ply import Ply
from utils.utils import transform_strain_global_to_local
from analysis.rfs import Kind_of_load, Type_parameter, RF


def solve_midplane(
        laminate: Laminate,
        N: np.ndarray,
        M: Optional[np.ndarray] = None,
        load_case: Optional[str] = None) -> np.ndarray:
    rhs = np.concatenate([N, M])
    sol = np.linalg.solve(laminate.ABD_matrices, rhs)
    eps = sol[:3]
    kappa = sol[3:]
    return eps, kappa


def laminate_plain_strength(
        laminate: Laminate,
        N: Optional[np.ndarray] = None,
        M: Optional[np.ndarray] = None,
        safety_factor=1.0,
        load_case=None,
        testing_mode=False,
):

    if N is None and M is None:
        print('[Error]: No input load.')
        return

    if N is None:
        N = np.zeros(3,)
    if M is None:
        M = np.zeros(3,)

    eps0, kappa0 = solve_midplane(laminate=laminate, N=N, M=M)
    stacking = laminate.stacking
    z = laminate.z_interfaces()

    RF_critical = float('inf')
    RF_table = []
    for k, ply in enumerate(stacking, start=1):
        z_bot = z[k-1]
        z_top = z[k]
        z_mid = (z_top + z_bot) / 2

        eps_global = eps0 + z_mid * kappa0

        if isinstance(ply, Ply):
            print(
                f'Analysing Ply {k} : {ply.material.name} / {ply.theta_deg} / zmid {z_mid}')
            eps_local = transform_strain_global_to_local(
                eps_global, ply.theta_deg)
            u_e1 = float(eps_local[0]) * 1.0e+6
            u_e2 = float(eps_local[1]) * 1.0e+6
            u_g12 = float(eps_local[2]) * 1.0e+6

            u_et_all = ply.material.u_et_all
            u_ec_all = ply.material.u_ec_all

            # Nota importante:
            # ----------------
            # El programa CLA, mete un valor de poison de 0.3 por defecto. Por
            # tanto a la hora de hacer comprobaciones con resultados de este
            # programa hay que activar la opcion "testing_mode = True"
            nu12 = 0.3 if testing_mode else ply.material.nu12
            nu12 = 0.3

            # Computing RF in direction 1
            if u_e1 > 0:
                u_all_1 = u_et_all
                mode1 = Kind_of_load.TENSION
                RF1 = abs((u_all_1 / (u_e1 * safety_factor)))
            elif u_e1 < 0:
                u_all_1 = u_ec_all
                mode1 = Kind_of_load.COMPRESSION
                RF1 = abs((u_all_1 / (u_e1 * safety_factor)))
            else:
                RF1 = float('inf')
                mode1 = 'none'
                u_all_1 = 'none'
            # Computing RF in direction 2
            if u_e2 > 0:
                u_all_2 = u_et_all
                mode2 = Kind_of_load.TENSION
                RF2 = abs((u_all_2 / (u_e2 * safety_factor)))
            elif u_e2 < 0:
                u_all_2 = u_ec_all
                mode2 = Kind_of_load.COMPRESSION
                RF2 = abs((u_all_2 / (u_e2 * safety_factor)))
            else:
                RF2 = float('inf')
                mode2 = 'none'
                u_all_1 = 'none'

            # Computing RF in direction 12
            u_all_12 = (1 + nu12) * max(u_et_all, u_ec_all)
            u_e12 = abs(u_e1 - u_e2)
            mode12 = Kind_of_load.SHEAR
            if u_e12 == 0.0:
                RF12 = float('inf')
            else:
                RF12 = abs((u_all_12 / (u_e12 * safety_factor)))

            RFs = [RF1, RF2, RF12]
            u_e = [u_e1, u_e2, u_e12]
            u_all = [u_all_1, u_all_2, u_all_12]
            modes = [mode1, mode2, mode12]

            RFmin = min(RFs)
            u_e_min = u_e[RFs.index(RFmin)]
            u_all_min = u_all[RFs.index(RFmin)]
            mode_min = modes[RFs.index(RFmin)]

            ply_data = RF(
                location=str(laminate.name) + "-Ply_" + str(k),
                material=ply.material.name,
                load_case=load_case,
                type_parameter=Type_parameter.STRAIN,
                units_parameter='micro-strains',
                kind_of_load=modes,
                ultimate_value=u_e,
                allowable_value=u_all,
                RF=RFs,
                remarks=''
            )

            if RFmin < RF_critical:
                RF_critical = RFmin
                crit_ply_data = ply_data

            RF_table += [ply_data]

    return crit_ply_data, RF_table
