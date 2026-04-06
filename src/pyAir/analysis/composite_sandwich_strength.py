# -*- coding: utf-8 -*-
import numpy as np
import math
from dataclasses import dataclass
from structural.laminate import Laminate
from structural.panel_core import Core


@dataclass
class Sandwich:

    top_facesheet: Laminate
    core: Core
    bot_facesheet: Laminate
    name: str | None = None

    @property
    def stacking(self):
        stack_b = self.bot_facesheet.stacking
        stack_t = self.top_facesheet.stacking
        stacking = stack_b + [self.core] + stack_t
        return stacking

    @property
    def thickness(self):
        top_facesheet_t = self.top_facesheet.thickness
        bot_facesheet_t = self.bot_facesheet.thickness
        core_t = self.core.material.thickness
        return bot_facesheet_t + core_t + top_facesheet_t

    def z_interfaces(self, z0: float = 0.0) -> list:
        if z0 == 0.0:
            z = [-1*self.thickness/2]
        # Computing the midplane height z
        else:
            z = [z0]

        for ply in self.stacking:
            z.append(z[-1] + ply.material.thickness)

        return z

    def stiffness_matrices(self):

        core_t = self.core.material.thickness
        fs_b = self.bot_facesheet
        fs_t = self.top_facesheet
        A_top, B_top, D_top = fs_t.stiffness_matrices(z0=-self.thickness/2)
        A_bot, B_bot, D_bot = fs_b.stiffness_matrices(z0=-self.thickness/2)

        A = A_top + A_bot
        B = B_top + B_bot
        D = D_top + D_bot

        K_sx = self.core.material.G13 * core_t
        K_sy = self.core.material.G23 * core_t

        return A, B, D, K_sx, K_sy

    @property
    def ABD_matrices(self):
        A, B, D, K_sx, K_sy = self.stiffness_matrices()
        return np.block([[A, B], [B, D]])

    def core_shear_strength(self,
                            Qxz: float,
                            Qyz: float,
                            N: float,
                            safety_factor=1.0,
                            thick_approach='conservative'
                            ):
        """
        Compute the Core Shear Strentgh according to CMH-17.Vol 16 eqs 4.6.2(b)

        Notes regarding N
        -----------------
        N values according to the different types of Honeycombs
        See CMH-17. Vol6, table 4.6.2

        Core Material   Cell Geometry   Cell Size   Density   N
        -------------   -------------   ---------   -------   ----
        Fiberglass      Hexagonal       3/8         3.5       1.34
        Fiberglass      Hexagonal       3/8         4.5       1.32
        Nomex           Flex            F50         5.5       1.45
        Nomex           Over.expand     3/16        3.0       1.50
        Kevlar          Hexagonal       1/8         3.0       1.23
        Korex           Hexagonal       1/8         6.0       1.20

        """
        if thick_approach.upper() == 'CONSERVATIVE':
            t = self.core.material.thickness
        else:
            t = self.core.material.thickness + 0.5 * \
                (self.top_facesheet.thickness + self.bot_facesheet.thickness)

        # Stresses
        tau_xz = Qxz / t
        tau_yz = Qyz / t
        tau_xz_all = self.core.material.F13
        tau_yz_all = self.core.material.F23

        # Computation of Reserve Factors
        RF_tau_xz = tau_xz_all / (safety_factor * tau_xz)
        RF_tau_yz = tau_yz_all / (safety_factor * tau_yz)

        # Core Shear Interaction R.F
        RF_shear = 1 / ((1/RF_tau_xz)**N + (1/RF_tau_yz)**N)**(1/N)

        return RF_tau_xz, RF_tau_yz, RF_shear

    def core_crushing_strength(self,
                               Mx: float = 0.0,
                               My: float = 0.0,
                               Pressure: float = 0.0,
                               safety_factor=1.0):

        # Formula according to CMH-17, vol6 eq. 4.6.4
        # -------------------------------------------
        core_t = self.core.material.thickness
        top_f_t = self.top_facesheet.thickness
        bot_f_t = self.bot_facesheet.thickness
        d = core_t + (top_f_t + bot_f_t)/2.0

        A, B, D, K_sx, K_sy = self.stiffness_matrices()

        D = np.linalg.inv(D)
        Dx = 1/D[0, 0]
        Dy = 1/D[1, 1]

        sigma_bn_x = (Mx**2)/(d * Dx)
        sigma_bn_y = (My**2)/(d * Dy)
        sigma_press = Pressure

        sigma = sigma_bn_x + sigma_bn_y + sigma_press
        sigma_all = self.core.material.F3c

        RF_core_flex_comp = sigma_all / (safety_factor * sigma)

        print('stress_app : ', -1*sigma)
        print('stress_all : ', sigma_all)

        return RF_core_flex_comp

    def core_shear_crimping_strength(self):

        Kb = self.core.material.Kbasis
        t_u = self.top_facesheet.thickness
        t_l = self.top_facesheet.thickness
        t_c = self.core.material.thickness
        Gxz = self.core.material.G13
        Gyz = self.core.material.G23
        d = t_c + (t_u + t_l) / 2
        denom = (t_u + t_l) * t_c
        Fc_x = Kb * d**2 * Gxz / denom
        Fc_y = Kb * d**2 * Gyz / denom
        Fs_xy = Kb * d**2 * math.sqrt(Gxz * Gyz) / denom

        return Fc_x, Fc_y, Fs_xy

    def wrinkling_strength(self,
                           wrink_coeff1: float = 0.247,
                           wrink_coeff2: float = 0.078,
                           wrink_coeff3: float = 0.33,
                           wrink_coeff4: float = 0.0,
                           tc_override: bool = False
                           ):

        # No me cuadra con lo que me sale de SANDRES

        # Computation of the facesheet_wrinkling (Fw), CMH-17. Vol6
        # Eq. 4.6.6.3 (a), (b) and (c)
        # ----------------------------------------------------------------

        C1 = wrink_coeff1
        C2 = wrink_coeff2
        C3 = wrink_coeff3
        C4 = wrink_coeff4

        facesheets = [self.top_facesheet, self.bot_facesheet]
        tc = self.core.material.thickness
        Gxz = self.core.material.G13
        Gyz = self.core.material.G23
        Ec = self.core.material.E3

        wrink_all = []

        for fs in facesheets:

            tf = fs.thickness
            A, B, D = fs.A_B_D

            if fs.is_symmetric:
                D = D
            else:
                D = fs.D_star

            # Usando D_inv
            D_inv = np.linalg.inv(D)
            D11_inv = D_inv[0, 0]
            D22_inv = D_inv[1, 1]
            Ex = 12 / (D11_inv * tf**3)
            Ey = 12 / (D22_inv * tf**3)
            # Usando D
            # D11 = D[0,0]
            # D22 = D[1,1]
            # #
            # Ex = (12 / (tf**3) ) * D11
            # Ey = (12 / (tf**3) ) * D22
            #
            # eqv_moduli = fs.laminate_apparent_moduli()
            # Ex = eqv_moduli['Ex']
            # Ey = eqv_moduli['Ey']
            print('Ex', Ex)
            tc_x_lim = 1.82 * tf * (Ex*Ec/(Gxz**2))**(1/3)
            tc_y_lim = 1.82 * tf * (Ey*Ec/(Gyz**2))**(1/3)
            tc_lim = min(tc_x_lim, tc_y_lim)

            for E, G in [(Ex, Gxz), (Ey, Gyz)]:

                if tc_override:
                    tc_lim = 1.82 * tf * (E*Ec/(G**2))**(1/3)
                print(tc_lim)
                if tc >= tc_lim:  # Thick core
                    Fw = C1 * (Ec * E * G)**(1/3) + C2 * G * (tc / tf)
                    print('Thick')
                else:  # Thin core
                    Fw = C3 * (Ec * E * (tf / tc))**(1/2) + C4 * G * (tc / tf)
                wrink_all += [Fw]

        return wrink_all

    def wrinkling_RFs(self,
                      N: np.array = None,
                      M: np.array = None,
                      wrink_coeffs: tuple = (0.247, 0.078, 0.33, 0.0),
                      safety_factor=1.0
                      ):

        # Nota:
        # No me sale igual que en SANDRES porque igual hay que computarlo en
        # ejes principales taly como indica el CMH-17 Vol.16

        C1, C2, C3, C4 = wrink_coeffs

        wrink_all = self.wrinkling_strength(C1, C2, C3, C4)

        Fwr_x_top = wrink_all[0]
        Fwr_y_top = wrink_all[1]

        Fwr_x_bot = wrink_all[2]
        Fwr_y_bot = wrink_all[3]

        loads = self.facesheet_stresses(N, M)
        sigma_x_top = loads['sigma_x_top']
        sigma_y_top = loads['sigma_y_top']
        sigma_xy_top = loads['sigma_xy_top']

        sigma_x_bot = loads['sigma_x_bot']
        sigma_y_bot = loads['sigma_y_bot']
        sigma_xy_bot = loads['sigma_xy_bot']

        RF_wr_x_top = Fwr_x_top / \
            (safety_factor * abs(sigma_x_top)) if sigma_x_top < 0.0 else float('inf')
        RF_wr_y_top = Fwr_y_top / \
            (safety_factor * abs(sigma_y_top)) if sigma_y_top < 0.0 else float('inf')

        RF_wr_x_bot = Fwr_x_bot / \
            (safety_factor * abs(sigma_x_bot)) if sigma_x_bot < 0.0 else float('inf')
        RF_wr_y_bot = Fwr_y_bot / \
            (safety_factor * abs(sigma_y_bot)) if sigma_y_bot < 0.0 else float('inf')

        res = dict(
            RF_wr_x_top=RF_wr_x_top,
            RF_wr_x_bot=RF_wr_x_bot,
            RF_wr_y_top=RF_wr_y_top,
            RF_wr_y_bot=RF_wr_y_bot
        )

        return res

    def facesheet_stresses(self, N: np.array = None, M: np.array = None):
        '''
        Computes the stresses of the top and bottom facesheets of a composite
        sandwich panel with a honeycomb core.

        Parameters
        ----------
        N : np.array, optional
            Mid-plane  forces per length unit {Nx, Ny and Nxy}
        M : np.array, optional
            Mid-plane flexural moment per length unit {Mx ,My and Mxy}

        Returns
        -------
        res : TYPE
            DESCRIPTION.

        '''

        if N is None:
            N = np.zeros((3,))
        if M is None:
            M = np.zeros((3,))

        core_t = self.core.material.thickness
        top_f_t = self.top_facesheet.thickness
        bot_f_t = self.bot_facesheet.thickness
        d = core_t + (top_f_t + bot_f_t) / 2

        Nx, Ny, Nxy = N
        Mx, My, Mxy = M

        sigma_x_Nx = Nx / (top_f_t + bot_f_t)
        sigma_y_Ny = Ny / (top_f_t + bot_f_t)
        sigma_xy_Nxy = Nxy / (top_f_t + bot_f_t)

        sigma_x_Mx = Mx / d
        sigma_y_My = My / d
        sigma_xy_Mxy = Mxy / d

        sigma_x_top = sigma_x_Nx + (sigma_x_Mx / top_f_t)
        sigma_x_bot = sigma_x_Nx - (sigma_x_Mx / bot_f_t)

        sigma_y_top = sigma_y_Ny + (sigma_y_My / top_f_t)
        sigma_y_bot = sigma_y_Ny - (sigma_y_My / bot_f_t)

        sigma_xy_top = sigma_xy_Nxy + (sigma_xy_Mxy / top_f_t)
        sigma_xy_bot = sigma_xy_Nxy - (sigma_xy_Mxy / bot_f_t)

        res = dict(
            sigma_x_top=sigma_x_top,
            sigma_y_top=sigma_y_top,
            sigma_xy_top=sigma_xy_top,
            sigma_x_bot=sigma_x_bot,
            sigma_y_bot=sigma_y_bot,
            sigma_xy_bot=sigma_xy_bot,
        )

        return res

    def dimpling_RFs(self, N: np.array = None, M: np.array = None, Kb: float = 1.0):

        # Computation according to CMH-17. Vol6
        # Eqs 4.6.5.4
        # -------------------------------------
        cs = self.core.material.cell_size
        Fc_dimp_top, Fs_dimp_top = self.top_facesheet.dimpling_strength(
            cs, Kb=Kb)
        Fc_dimp_bot, Fs_dimp_bot = self.bot_facesheet.dimpling_strength(
            cs, Kb=Kb)

        loads = self.facesheet_stresses(N, M)
        sigma_x_top = loads['sigma_x_top']
        sigma_y_top = loads['sigma_y_top']
        sigma_xy_top = loads['sigma_xy_top']
        sigma_x_bot = loads['sigma_x_bot']
        sigma_y_bot = loads['sigma_y_bot']
        sigma_xy_bot = loads['sigma_xy_bot']

        # top facesheet
        fx_top = sigma_x_top if sigma_x_top < 0.0 else 0.0
        fy_top = sigma_y_top if sigma_y_top < 0.0 else 0.0
        Rc_dimp_top = abs((fx_top + fy_top) / Fc_dimp_top)
        Rs_dimp_top = abs(sigma_xy_top / Fs_dimp_top)
        R_comb_den = (2 * Rs_dimp_top)
        R_comb_num = (-1 * Rc_dimp_top +
                      math.sqrt(Rc_dimp_top**2 + 4 * Rs_dimp_top**2))

        if (Rc_dimp_top != 0.0 and Rs_dimp_top != 0.0):
            R_comb_top = (R_comb_num / R_comb_den)

        elif Rc_dimp_top == 0.0:
            R_comb_top = Rs_dimp_top

        else:
            R_comb_top = Rc_dimp_top

        # bot facesheet
        fx_bot = sigma_x_bot if sigma_x_bot < 0.0 else 0.0
        fy_bot = sigma_y_bot if sigma_y_bot < 0.0 else 0.0
        Rc_dimp_bot = abs((fx_bot + fy_bot) / Fc_dimp_bot)
        Rs_dimp_bot = abs(sigma_xy_bot / Fs_dimp_bot)
        R_comb_den = (2 * Rs_dimp_bot)
        R_comb_num = (-1 * Rc_dimp_bot +
                      math.sqrt(Rc_dimp_bot**2 + 4 * Rs_dimp_bot**2))

        if (Rc_dimp_bot != 0.0 and Rs_dimp_bot != 0.0):
            R_comb_bot = (R_comb_num / R_comb_den)

        elif Rc_dimp_bot == 0.0:
            R_comb_bot = 1/Rs_dimp_bot

        else:
            R_comb_bot = 1/Rc_dimp_bot

        return R_comb_top, R_comb_bot
