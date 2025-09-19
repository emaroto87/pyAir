# -*- coding: utf-8 -*-
"""
Created on Fri Sep  5 08:11:24 2025

@author: U69432
"""

import math

import numpy as np
from typing import Optional
from dataclasses import dataclass, asdict
import abc

from enum import Enum

__version__= '0.0.1'
__release_date__ = '10/09/2025'

# ---->>> BEGIN SECTION
# ==========================================================================
# UTILITIES
# =========================================================================


def transform_stress_global_to_local(stress_global: np.ndarray, theta_deg: float) -> np.ndarray:
    '''
    Performs a rotation transformation of the components of a plane-stress vector
    {sigma_x, sigma_y, sigma_xy} from a global cordinate system into a local.

    Parameters
    ----------
    stress_global : np.ndarray
    
    theta_deg : float
        Angle of rotation about the 3-axis

    Returns
    -------
    stress_local : np.ndarray
        DESCRIPTION.

    '''
    
    theta = math.radians(theta_deg)
    c = math.cos(theta)
    s = math.sin(theta)
    c2 = c**2
    s2 = s**2
    cs = c*s
    sx , sy , txy = stress_global
    s1 = c2 * sx + s2 * sy + 2 * cs * txy
    s2 = s2 * sx + c2 * sy - 2 * cs * txy
    t12 = -cs * sx + cs * sy + (c2 - s2) * txy
    stress_local = np.array([s1,s2,t12])
    return stress_local

def transform_strain_global_to_local(strain_global: np.ndarray, theta_deg: float) -> np.ndarray:
    theta = math.radians(theta_deg)
    c = math.cos(theta)
    s = math.sin(theta)
    c2 = c**2
    s2 = s**2
    cs = c*s
    ex , ey , gxy = strain_global
    e1 = c2 * ex + s2 * ey + cs * gxy
    e2 = s2 * ex + c2 * ey - cs * gxy
    g12 = -2.0 * cs * ex + 2.0 * cs * ey + (c2 - s2) * gxy
    strain_local = np.array([e1,e2,g12])
    return strain_local


### <<<---- END SECTION

# ------------------------------------------------->[RF_ACD class and auxiliar]
@dataclass
class Type_parameter:
    FORCE = 'F'
    MOMENT = 'M'
    STRESS = 'St'
    STRAIN = 'Str'
    RUNNING_LOAD = 'FI'
    DEFLECTION = 'D'
    PRESSURE = 'Pr'
    LOAD_FACTOR = 'g'
    
    def __repr__(self):
        return str(asdict(self))
    def __str__(self):
        return str(asdict(self))
        
    
@dataclass
class Kind_of_load:
    SHEAR = 'S'
    TENSION = 'T'
    COMPRESSION = 'C'
    BENDING = 'Bn'
    BEARING = 'B'
    COMPLEX = 'Clpx'
    

@dataclass
class RF_ACD:
    location : None | str = None
    material : None | str = None
    load_case : None | str = None
    type_parameter : None | Type_parameter = None
    units_parameter : None | str = None
    kind_of_load : None | Kind_of_load = None
    ultimate_value : None | float = None
    allowable_value : None | float = None
    RF : None | float = None
    remarks : None | str = None
    
    # [Posibles Mejoras:]
    # - Una función que te transforme en un DataFrame   
# <------------------------------------------------[ RF_ACD class and auxiliar]


# ------------------------------------------------->[RF_ACD class and auxiliar]
@dataclass
class Forpla:
    load_case : str
    Nx : float = 0.0
    Ny : float = 0.0 
    Nxy : float = 0.0 
    Mx : float = 0.0
    My : float = 0.0
    Mxy : float = 0.0
    Qx : float = 0.0
    Qy : float = 0.0
    Pr : float = 0.0

    # def __str__(self):
    #     return vars(self).items()
    def to_array(self):
        vec = [ v for k,v in vars(self).items() if not k == 'load_case' ]
        return np.array(vec)
    
    

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

@dataclass
class Orthotropic:
    name : str
    thickness: float
    E1 : float
    E2 : float
    G12 : float
    nu12 : float
    norm : None | str = None # Norm or spec
    cond : None | str = None # Condition Dry/Wet
    temp : None | str = None # Normally, -55C , RT, +55C, +70C, +90C and +120C
    Kbm : None | float = 1.0 # Mean to B-basis Knock-down factor
    E1b : None | float = None 
    E2b : None | float = None
    G12b : None | float = None
    CTE1 : None | float = None # Coeff. Thermal Exp. at direction 1
    CTE2 : None | float = None # Coeff. Therma Exp.  at direction 2
    # Strength allowables
    u_et_all : None | float = None # rupture tension micro-strains allowable
    u_ec_all : float | float = None# rupture compresion micro-strains allowable
    F1t : float | None = None #
    F1c : float | None = None #
    # Damage torelance allowables
    uTAI : float | None = None # damage torlerance micro-strains tension after impact
    uCAI : float | None = None # damage tolerance micro-strains compresion after impact
    uBAI : float | None = None # damage tolerance micro-strains bending after impact
    iCDT: float | None = None
    
    
    
    @property
    def nu21(self) -> float:
        nuyx = self.nu12* self.E2 / self.E1
        return nuyx
    
@dataclass
class Honeycomb:
    name : str
    thickness : float 
    cell_size : float
    G13 : float
    G23 : float
    E3 : float
    F13 : float # Tau_13 stress allowable 
    F23 : float # Tau_23 stress allowable
    F3t : float # Sigma_33 traction stress allowable
    F3c : float # Sigma_33 compression stress allowable
    Kbasis : float # Mean to K-basis factor
    ekdf : float # Environmental Knock-Down Factor


@dataclass    
class Ply:
    material : Orthotropic
    theta_deg : float
    
    def q_plane_stress(self) ->tuple:
        d = 1 - self.material.nu12 * self.material.nu21
        Q11 = self.material.E1/d
        Q12 = self.material.nu12*self.material.E2/d
        Q22 = self.material.E2/d
        Q66 =  self.material.G12
        
        return Q11, Q12, Q22, Q66
    
    def Qbar(self):
        Q11, Q12, Q22, Q66 = self.q_plane_stress()
        theta = math.radians(self.theta_deg)
        # Sins and Cosins
        c = math.cos(theta)
        s = math.sin(theta)
        c2 = c*c
        s2 = s*s
        c4 = c2*c2
        s4 = s2*s2
        s3c = s2*s*c
        c3s = c2*c*s
        
        # Coefficients
        Q11b = Q11*c4 + 2*(Q12 + 2*Q66)*s2*c2 + Q22*s4
        Q22b = Q11*s4 + 2*(Q12 + 2*Q66)*s2*c2 + Q22*c4
        Q12b = (Q11 + Q22 - 4*Q66)*s2*c2 + Q12*(s4+c4)
        Q16b = (Q11 - Q12 - 2*Q66)*c3s-(Q22-Q12-2*Q66)*s3c
        Q26b = (Q11 - Q12 - 2*Q66)*s3c-(Q22-Q12-2*Q66)*c3s
        Q66b = (Q11 + Q22 - 2*Q12 - 2*Q66)*s2*c2 + Q66*(s4+c4)
        
        return [
            [Q11b,Q12b,Q16b],
            [Q12b,Q22b,Q26b],
            [Q16b,Q26b,Q66b]
            ]

@dataclass
class Core:
    material : Honeycomb
    theta_deg : float
    
    def core_shear_strength(self,
                            Q13 : float,
                            Q23 : float,
                            N : float, # Core shear interaction factor
                            safety_factor=1.0):
       
        # Notes regarding N
        # -----------------
        # N values according to the different types of Honeycombs
        # See CMH-17. Vol6, table 4.6.2
        # 
        # Core Material   Cell Geometry   Cell Size   Density   N
        # Fiberglass      Hexagonal       3/8         3.5       1.34
        # Fiberglass      Hexagonal       3/8         4.5       1.32
        # Nomex           Flex            F50         5.5       1.45
        # Nomex           Over.expand     3/16        3.0       1.50
        # Kevlar          Hexagonal       1/8         3.0       1.23
        # Korex           Hexagonal       1/8         6.0       1.20
        
        thickness = self.material.thickness
        
        # According to CMH-17.Vol 6 eqs 4.6.2(b)
        # --------------------------------------
        # Stresses
        tau_xz = Q13 / thickness
        tau_yz = Q23 / thickness
        tau_xz_all = self.material.F13
        tau_yz_all = self.material.F23
        
        # Computation of Reserve Factors
        RF_tau_xz = tau_xz_all / (safety_factor * tau_xz)
        RF_tau_yz = tau_yz_all / (safety_factor * tau_yz)
        
        # Core Shear Interaction R.F
        RF_shear = 1 / ((1/RF_tau_xz)**N + (1/RF_tau_yz)**N)**(1/N)
        
        # Results
        
        return RF_tau_xz, RF_tau_yz , RF_shear
        
def is_symmetric(stacking):
    n = len(stacking)
    j = int(n/2) if n % 2 == 0.0 else int((n-1)/2)
    
    _ = []
    for i in range(j):
        cond1 = stacking[i].material == stacking[-(i+1)].material
        cond2 = stacking[i].theta_deg == stacking[-(i+1)].theta_deg
        cond = cond1 and cond2
        _ += [cond]
    is_symm = True if all(_) else False

    return is_symm     
        
        
        
    
    
          
@dataclass
class Laminate:                                                                #### -> LAMINATE CLASS
    stacking : list[Ply]
    name : str | None = None
    '''
    [Idea] Formas de dar el stacking:
        - Lista : [Tupla(Ply,angulo)]
        - Dict  : {Ply, secuencia)
                   
    '''
        
    @property
    def thickness(self):
        t = sum(ply.material.thickness for ply in self.stacking)
        return t
    
    def z_interfaces(self, z0 : float = 0.0)->list:
        if z0==0.0:
            z = [-1*self.thickness/2]
        # Computing the midplane height z
        else :
            z = [z0]

        for ply in self.stacking:
            z.append(z[-1] + ply.material.thickness)
    
        return z
    
    def stiffness_matrices(self , z0 : float = 0.0):
        
        # Computing ply interfaces
        z = self.z_interfaces(z0 = z0)
        
        # Computing A, B and D matrices 
        A = np.zeros((3,3)) # Pure Membrane stiffness Matrix
        B = np.zeros((3,3)) # Pure Bending stiffness Matrix
        D = np.zeros((3,3)) # Membrane-Bending Coupling stiffness Matrix
            
        for k, ply in enumerate(self.stacking ,start=1):
            Qbar = np.array(ply.Qbar())
            z_bot = z[k-1]
            z_top = z[k]
            A += Qbar * (z_top - z_bot)
            B += (1/2) * Qbar * (z_top**2 - z_bot**2)
            D += (1/3) * Qbar * (z_top**3 - z_bot**3)
        
        return np.round(A,2), np.round(B,2), np.round(D,2)
   
    @property
    def ABD_matrices(self):
        A,B,D = self.stiffness_matrices()
        return np.block([[A,B],[B,D]])


    def laminate_apparent_moduli(self):

        A,B,D = self.stiffness_matrices()
        Abar = A / self.thickness
        S = np.linalg.inv(Abar)
        S11 = S[0,0]
        S22 = S[1,1]
        S66 = S[2,2]
        Ex = 1.0 / S11
        Ey = 1.0 / S22
        Gxy = 1.0 / S66
        nuxy = - S[0,1] / S[0,0]

        return dict(Ex = Ex, Ey = Ey, Gxy =Gxy, nuxy = nuxy)

    @property
    def D_star(self):
        A, B, D = self.stiffness_matrices()
        A_inv = np.linalg.inv(A)
        D_star = D - B @ A_inv @ B                                             # @ Means product of matrices
        return D_star

    def dimpling_strength(self, core_cell_size: float , Kb : float = 1.0):
        
        tf = self.thickness
        cs = core_cell_size
        A,B,D = self.stiffness_matrices()
        pi = math.pi
        
        # Computation of compression strength in dimpling (Fc), CMH-17.Vol6
        # Eq 4.6.5.1 (c):
        # -----------------------------------------------------------------
        if is_symmetric(self.stacking):
            D_prime = D
            print (True)
        else:
            D_prime = self.D_star
            
        D_prime = D

        D11 = D_prime[0,0]
        D12 = D_prime[0,1]
        D22 = D_prime[1,1]
        D66 = D_prime[2,2]
        
        Fc = Kb * (1 / tf) * ((pi/cs)**2) * ( D11 + 2 * (D12 + 2 * D66) + D22 )
        
        # Computation of the shear strength in dimpling (Fs), CMH-17. Vol6
        # Eq. 4.6.5.3
        # ----------------------------------------------------------------
        eqv_moduli = self.laminate_apparent_moduli()
        Ex = eqv_moduli['Ex']
        Ey = eqv_moduli['Ey']
        Fs = Kb * 0.6 * min (Ex,Ey) * (tf/cs)**(1.5)
        
        return Fc, Fs
    

# --> CLASS SANDWICH    
@dataclass
class Sandwich:
    
    top_facesheet : Laminate
    core : Core
    bot_facesheet : Laminate
    name : str | None = None
    
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
    
    def z_interfaces(self, z0 : float = 0.0)->list:
        if z0==0.0:
            z = [-1*self.thickness/2]
        # Computing the midplane height z
        else :
            z = [z0]

        for ply in self.stacking:
            z.append(z[-1] + ply.material.thickness)
    
        return z    

    def stiffness_matrices(self):
        
        core_t = self.core.material.thickness
        fs_b = self.bot_facesheet
        fs_t = self.top_facesheet
        A_top, B_top, D_top = fs_t.stiffness_matrices(z0 = -self.thickness/2)
        A_bot, B_bot, D_bot = fs_b.stiffness_matrices(z0 = -self.thickness/2)

        A = A_top + A_bot
        B = B_top + B_bot
        D = D_top + D_bot

        K_sx = self.core.material.G13 * core_t
        K_sy = self.core.material.G23 * core_t
        
        return A, B , D , K_sx, K_sy
    
    @property
    def ABD_matrices(self):
        A, B, D, K_sx, K_sy= self.stiffness_matrices()
        return np.block([[A,B],[B,D]])
     
    def core_shear_strength(self,
                      Qx : float,
                      Qy : float,
                      N : float,
                      safety_factor = 1.0,
                      thick_approach = 'conservative'
                      ):
        
        if thick_approach.upper() == 'CONSERVATIVE':
            t = self.core.material.thickness
        else:
            t = self.core.material.thickness + 0.5*(self.top_facesheet.thickness + self.bot_facesheet.thickness)
            
        RF_tau_xz, RF_tau_yz, RF_shear = self.core.core_shear_strength(
                                        Q13 = Qx,
                                        Q23 = Qy,
                                        N = N,
                                        safety_factor = safety_factor)
        
        return RF_tau_xz, RF_tau_yz , RF_shear
                                        
    def core_crushing_strength(self,
                                  Mx : float = 0.0,
                                  My : float = 0.0,
                                  Pressure : float = 0.0,
                                  safety_factor = 1.0):
    
        # Formula according to CMH-17, vol6 eq. 4.6.4
        # -------------------------------------------
        core_t = self.core.material.thickness
        top_f_t = self.top_facesheet.thickness
        bot_f_t = self.bot_facesheet.thickness
        d = core_t + (top_f_t + bot_f_t)/2.0
        
        A, B , D , K_sx, K_sy = self.stiffness_matrices()

        D = np.linalg.inv(D)
        Dx = 1/D[0,0]
        Dy = 1/D[1,1]
        
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
        Fs_xy =Kb * d**2 * math.sqrt(Gxz * Gyz) / denom
        
        return Fc_x, Fc_y ,Fs_xy
        
    def wrinkling_strength(self,
                           wrink_coeff1 : float = 0.247,
                           wrink_coeff2 : float = 0.078,
                           wrink_coeff3 : float = 0.33,
                           wrink_coeff4 : float = 0.0,
                           tc_override : bool = False
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
            A , B , D = fs.stiffness_matrices()
            
            if is_symmetric(fs.stacking):
                D = D
            else:
                D = fs.D_star
            
            # Usando D_inv
            D_inv = np.linalg.inv(D)
            D11_inv = D_inv[0,0]
            D22_inv = D_inv[1,1]
            Ex = 12 / ( D11_inv * tf**3)
            Ey = 12 / ( D22_inv * tf**3)
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
            tc_lim = min(tc_x_lim,tc_y_lim)
            
            for E,G in [(Ex,Gxz),(Ey,Gyz)]:
                
                if tc_override: tc_lim = 1.82 * tf * (E*Ec/(G**2))**(1/3)
                print(tc_lim)
                if tc >= tc_lim : # Thick core
                    Fw = C1 * (Ec * E * G)**(1/3) + C2 * G * (tc / tf)
                    print('Thick')
                else : # Thin core
                    Fw = C3 * (Ec * E * (tf / tc))**(1/2) + C4 * G * (tc / tf)
                wrink_all += [Fw]
        
        return wrink_all
    
    def wrinkling_RFs(self, 
                      N : np.array = None, 
                      M : np.array = None,
                      wrink_coeffs : tuple = (0.247,0.078,0.33,0.0),
                      safety_factor = 1.0
                      ):
        
        ### Nota: 
        # No me sale igual que en SANDRES porque igual hay que computarlo en 
        # ejes principales taly como indica el CMH-17 Vol.16
        
        C1,C2,C3,C4 = wrink_coeffs
        
        wrink_all = self.wrinkling_strength(C1,C2,C3,C4)
        
        Fwr_x_top = wrink_all[0]
        Fwr_y_top = wrink_all[1]
        
        Fwr_x_bot = wrink_all[2]
        Fwr_y_bot = wrink_all[3]
        
        loads = self.facesheet_stresses(N,M)
        sigma_x_top = loads['sigma_x_top']
        sigma_y_top = loads['sigma_y_top']
        sigma_xy_top = loads['sigma_xy_top']
        
        sigma_x_bot = loads['sigma_x_bot']
        sigma_y_bot = loads['sigma_y_bot']
        sigma_xy_bot = loads['sigma_xy_bot']
            
        RF_wr_x_top = Fwr_x_top / (safety_factor * abs(sigma_x_top)) if sigma_x_top < 0.0 else float('inf')
        RF_wr_y_top = Fwr_y_top / (safety_factor * abs(sigma_y_top)) if sigma_y_top < 0.0 else float('inf')
        
        RF_wr_x_bot = Fwr_x_bot / (safety_factor * abs(sigma_x_bot)) if sigma_x_bot < 0.0 else float('inf')
        RF_wr_y_bot = Fwr_y_bot / (safety_factor * abs(sigma_y_bot)) if sigma_y_bot < 0.0 else float('inf')
        
        res = dict(
            RF_wr_x_top = RF_wr_x_top,
            RF_wr_x_bot = RF_wr_x_bot,
            RF_wr_y_top = RF_wr_y_top,
            RF_wr_y_bot = RF_wr_y_bot
            )
        
        return res
    
    def facesheet_stresses(self, N:np.array = None, M : np.array = None):
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
        
        
        if N is None : N = np.zeros((3,))
        if M is None : M = np.zeros((3,))
        
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
        
        sigma_x_top = sigma_x_Nx + ( sigma_x_Mx / top_f_t)
        sigma_x_bot = sigma_x_Nx - ( sigma_x_Mx / bot_f_t)
        
        sigma_y_top = sigma_y_Ny + ( sigma_y_My / top_f_t )
        sigma_y_bot = sigma_y_Ny - ( sigma_y_My / bot_f_t )
        
        sigma_xy_top = sigma_xy_Nxy + ( sigma_xy_Mxy / top_f_t )
        sigma_xy_bot = sigma_xy_Nxy - ( sigma_xy_Mxy / bot_f_t )
        
        res = dict(
            sigma_x_top = sigma_x_top,
            sigma_y_top = sigma_y_top,
            sigma_xy_top = sigma_xy_top,
            sigma_x_bot = sigma_x_bot,
            sigma_y_bot = sigma_y_bot,
            sigma_xy_bot = sigma_xy_bot,            
            )
        
        return res
        
    def dimpling_RFs(self, N: np.array = None, M: np.array = None, Kb : float = 1.0):
        
        # Computation according to CMH-17. Vol6 
        # Eqs 4.6.5.4
        # -------------------------------------
        cs = self.core.material.cell_size
        Fc_dimp_top, Fs_dimp_top = self.top_facesheet.dimpling_strength(cs,Kb=Kb)
        Fc_dimp_bot, Fs_dimp_bot = self.bot_facesheet.dimpling_strength(cs,Kb=Kb)
        
        
        loads = self.facesheet_stresses(N,M)
        sigma_x_top = loads['sigma_x_top']
        sigma_y_top = loads['sigma_y_top']
        sigma_xy_top = loads['sigma_xy_top']
        sigma_x_bot = loads['sigma_x_bot']
        sigma_y_bot = loads['sigma_y_bot']
        sigma_xy_bot = loads['sigma_xy_bot']

        # top facesheet
        fx_top = sigma_x_top if sigma_x_top < 0.0 else 0.0
        fy_top = sigma_y_top if sigma_y_top < 0.0 else 0.0
        Rc_dimp_top = abs(( fx_top + fy_top) / Fc_dimp_top)
        Rs_dimp_top = abs(sigma_xy_top / Fs_dimp_top)
        R_comb_den = (2 * Rs_dimp_top)
        R_comb_num = (-1* Rc_dimp_top + math.sqrt(Rc_dimp_top**2 + 4 * Rs_dimp_top**2))
        
        if ( Rc_dimp_top != 0.0 and Rs_dimp_top != 0.0):
            R_comb_top = (R_comb_num / R_comb_den) 
            
        elif Rc_dimp_top == 0.0:
            R_comb_top = Rs_dimp_top
        
        else:
            R_comb_top = Rc_dimp_top
        
        # bot facesheet
        fx_bot = sigma_x_bot if sigma_x_bot < 0.0 else 0.0
        fy_bot = sigma_y_bot if sigma_y_bot < 0.0 else 0.0
        Rc_dimp_bot = abs(( fx_bot + fy_bot) / Fc_dimp_bot)
        Rs_dimp_bot = abs(sigma_xy_bot / Fs_dimp_bot)
        R_comb_den = (2 * Rs_dimp_bot)
        R_comb_num = (-1* Rc_dimp_bot + math.sqrt(Rc_dimp_bot**2 + 4 * Rs_dimp_bot**2))
        
        if ( Rc_dimp_bot != 0.0 and Rs_dimp_bot != 0.0):
            R_comb_bot = (R_comb_num / R_comb_den) 
        
        elif Rc_dimp_bot == 0.0:
            R_comb_bot = 1/Rs_dimp_bot
       
        else:
            R_comb_bot = 1/Rc_dimp_bot
        
        return R_comb_top, R_comb_bot
        
    
# <-- CLASS SANDWICH         

       

def solve_midplane( laminate : Laminate , N: np.ndarray, M: Optional[np.ndarray] = None, load_case : Optional[str] = None) -> np.ndarray:
    rhs = np.concatenate([N,M])
    sol = np.linalg.solve(laminate.ABD_matrices,rhs)
    eps = sol[:3]
    kappa = sol[3:]
    return eps, kappa

def plain_strength( 
        laminate : Laminate | Sandwich ,
        N: Optional[np.ndarray] = None ,
        M: Optional[np.ndarray] = None ,
        safety_factor = 1.0 ,
        load_case = None,
        testing_mode = False,
        ):
    
    if N is None and M is None:
        print('[Error]: No input load.')
        return
    
    if N is None : N = np.zeros(3,)
    if M is None : M = np.zeros(3,)
    
    eps0 , kappa0 = solve_midplane(laminate = laminate, N = N, M = M)
    stacking = laminate.stacking
    z = laminate.z_interfaces()
    
    RF_critical = float('inf')
    RF_table = []
    for k,ply in enumerate(stacking,start=1):
        z_bot = z[k-1]
        z_top = z[k]
        z_mid = ( z_top + z_bot ) / 2

        eps_global = eps0 + z_mid * kappa0
                
        if isinstance(ply, Ply):
            print(f'Analysing Ply {k} : {ply.material.name} / {ply.theta_deg} / zmid {z_mid}')
            eps_local = transform_strain_global_to_local(eps_global, ply.theta_deg)
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
                RF1 = abs((u_all_1 / ( u_e1 * safety_factor)))
            elif u_e1 < 0:
                u_all_1 = u_ec_all
                mode1 = Kind_of_load.COMPRESSION
                RF1 = abs((u_all_1/ ( u_e1 * safety_factor)))
            else:
                RF1 = float('inf')
                mode1 = 'none'
                u_all_1 = 'none'
            # Computing RF in direction 2 
            if u_e2 > 0: 
                u_all_2 = u_et_all
                mode2 = Kind_of_load.TENSION
                RF2 = abs((u_all_2 / ( u_e2 * safety_factor)))
            elif u_e2 < 0:
                u_all_2 = u_ec_all
                mode2 = Kind_of_load.COMPRESSION
                RF2 = abs((u_all_2 / ( u_e2 * safety_factor)))
            else:
                RF2 = float('inf')
                mode2 = 'none'
                u_all_1 = 'none'
                
            # Computing RF in direction 12 
            u_all_12 = (1 + nu12) * max(u_et_all,u_ec_all)
            u_e12 = abs(u_e1 - u_e2)
            mode12 = Kind_of_load.SHEAR
            if u_e12 == 0.0:
                RF12 = float('inf')
            else:
                RF12 = abs((u_all_12/ ( u_e12 * safety_factor)))
            
            RFs =  [RF1,RF2,RF12]
            u_e = [u_e1,u_e2,u_e12]
            u_all = [u_all_1,u_all_2,u_all_12]
            modes = [mode1,mode2,mode12]
            
            RFmin = min(RFs)
            u_e_min = u_e[RFs.index(RFmin)]
            u_all_min = u_all[RFs.index(RFmin)]
            mode_min = modes[RFs.index(RFmin)]
            
            ply_data = RF_ACD(
                location = str(laminate.name) + "-Ply_" + str(k),
                material = ply.material.name,
                load_case = load_case,
                type_parameter = Type_parameter.STRAIN,
                units_parameter = 'micro-strains',
                kind_of_load = modes,
                ultimate_value = u_e,
                allowable_value = u_all,
                RF = RFs,
                remarks = ''
                )
            
            if RFmin < RF_critical:
                RF_critical = RFmin
                crit_ply_data = ply_data
            
            RF_table += [ply_data]
            
    return crit_ply_data, RF_table
                
                
    
    
    
    

'''
Comparativa de resultados con los obtenidos a través de SANDRES
Ver archivo SANDRES_Panel_Mono_worksapce.out

'''
IMA21E_HW_70 = Orthotropic(
    name = 'IMA-M21E-UD-70',
    thickness = 0.184,
    E1 = 154000,
    E2 = 8500,
    G12 = 4200,
    nu12 = 0.350,
    u_et_all = 15500,
    u_ec_all = 7130
    )

T300PW_HW_70 = Orthotropic(
    name = 'T300PW_HW_70',
    thickness = 0.237,
    E1 = 50000,
    E2 = 50000,
    G12 = 1000,
    nu12 = 0.05,
    u_et_all = 9100,
    u_ec_all = 8160
    )

Honeycomb_N636_70 = Honeycomb(
    name = "Honeycomb_N636_70",
    thickness = 6.0,
    cell_size = 4.8,
    G13 = 26.2,
    G23 = 66.6,
    E3 = 25.8,
    F13 = 0.495,
    F23 = 0.664,
    F3t = 2.565,
    F3c = 0.865,
    Kbasis = 0.740, # Mean to K-basis factor
    ekdf = 1.0
    )

Core_N636_0 = Core(Honeycomb_N636_70 ,0)

stacking_0 = [
    Ply(T300PW_HW_70,0),
    Ply(T300PW_HW_70,45),
    Ply(T300PW_HW_70,0)
    ]

#Test_1
stacking_1 = [
    Ply(T300PW_HW_70,45),
    Ply(T300PW_HW_70,0),
    Ply(T300PW_HW_70,45),
    Ply(T300PW_HW_70,0),
    Ply(T300PW_HW_70,45),
    # SYM
    Ply(T300PW_HW_70,0),
    # SYM
    Ply(T300PW_HW_70,45),
    Ply(T300PW_HW_70,0),
    Ply(T300PW_HW_70,45),
    Ply(T300PW_HW_70,0),
    Ply(T300PW_HW_70,45)
    ]

#Test_2
stacking_2 = [
    Ply(IMA21E_HW_70,0),
    Ply(IMA21E_HW_70,45),
    Ply(IMA21E_HW_70,-45),
    Ply(IMA21E_HW_70,0),
    Ply(IMA21E_HW_70,45),
    Ply(IMA21E_HW_70,-45),
    Ply(IMA21E_HW_70,0),
    Ply(IMA21E_HW_70,-45),
    Ply(IMA21E_HW_70,45),
    Ply(IMA21E_HW_70,0),
    Ply(IMA21E_HW_70,-45),
    Ply(IMA21E_HW_70,45),
    Ply(IMA21E_HW_70,0)
    ]


# Test_3
stacking_3 = [
    Ply(IMA21E_HW_70,0),
    Ply(IMA21E_HW_70,90),
    Ply(IMA21E_HW_70,45),
    Ply(IMA21E_HW_70,-45),
    Ply(IMA21E_HW_70,0),
    Ply(IMA21E_HW_70,90),
    Ply(IMA21E_HW_70,45),
    Ply(IMA21E_HW_70,-45),
    #---SYM
    Ply(IMA21E_HW_70,0),
    #---SYM
    Ply(IMA21E_HW_70,-45),
    Ply(IMA21E_HW_70,45),
    Ply(IMA21E_HW_70,90),
    Ply(IMA21E_HW_70,0),
    Ply(IMA21E_HW_70,-45),
    Ply(IMA21E_HW_70,45),
    Ply(IMA21E_HW_70,90),
    Ply(IMA21E_HW_70,0),
    ]

# Test_4
stacking_4 = [
    Ply(IMA21E_HW_70,0),
    Ply(IMA21E_HW_70,0),
    Ply(IMA21E_HW_70,0),
    Ply(IMA21E_HW_70,0),
    Ply(IMA21E_HW_70,0),
    ]

l0 = Laminate(stacking_0)
l1 = Laminate(stacking_1)
l2 = Laminate(stacking_2)
l3 = Laminate(stacking_3)
l4 = Laminate(stacking_4)
l5 = Laminate(
    [
     Ply(T300PW_HW_70,45),
     Ply(T300PW_HW_70,45),
     Ply(T300PW_HW_70,45)
     ]
    )
l6 = Laminate(
    [
     Ply(T300PW_HW_70,0),
     Ply(T300PW_HW_70,0),
     Ply(T300PW_HW_70,0)
     ]
    )

a0,b0,d0 = l0.stiffness_matrices()
a1,b1,d1 = l0.stiffness_matrices()
a3,b3,d3 = l3.stiffness_matrices()
a4,b4,d4 = l4.stiffness_matrices()


sandwich1 = Sandwich(
    top_facesheet = l0,
    core = Core_N636_0,
    bot_facesheet= l0,
)
sandwich_unsym = Sandwich(
    top_facesheet = l5,
    core = Core_N636_0,
    bot_facesheet= l6,
)
        
# IMA21E_HW_70_UP
IMA21E_HW_70_UP = Laminate([
    Ply(IMA21E_HW_70,0),
    Ply(IMA21E_HW_70,45),
    Ply(IMA21E_HW_70,90),
    Ply(IMA21E_HW_70,0),
    Ply(IMA21E_HW_70,45)
])
IMA21E_HW_70_LO = Laminate([
    Ply(IMA21E_HW_70,90),
    Ply(IMA21E_HW_70,90),
    Ply(IMA21E_HW_70,45),
    Ply(IMA21E_HW_70,45),
    Ply(IMA21E_HW_70,0)
])

SANDRES_Panel_ASym_Thin_Dimpling = Sandwich(
    top_facesheet = IMA21E_HW_70_UP,
    core = Core_N636_0,
    bot_facesheet= IMA21E_HW_70_LO
)
M = np.array([1000.0,0.0,0.0])
N = np.array([0.0,0.0,0.0])
SANDRES_Panel_ASym_Thin_Dimpling.dimpling_RFs(
    N = N,
    M = M,
    Kb = 0.74
)
SANDRES_Panel_ASym_Thin_Dimpling.top_facesheet.dimpling_strength(
    core_cell_size = SANDRES_Panel_ASym_Thin_Dimpling.core.material.cell_size,
    Kb = 0.74
)