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

# ---->>> BEGIN SECTION
# ==========================================================================
# UTILITIES
# =========================================================================


def transform_stress_global_to_local(stress_global: np.ndarray, theta_deg: float) -> np.ndarray:
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
    stress_local = np.ndarray([s1,s2,t12])
    return stress_local

### <<<---- END SECTION
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
    
    def __post_init__(self):
        if not isinstance(self.type_parameter,(type(None),type(Type_parameter))):
            raise ValueError(
                f"Value of type_parameter must be one value of {list(Type_parameter)}"
                )
        if not isinstance(self.kind_of_load,(type(None),type(Kind_of_load))):
            raise ValueError(
                f"Value of type_parameter must be one value of {list(Type_parameter)}"
                )
    
    
    def __str__(self):
        return asdict(self)
    def __repr__(self):
        return asdict(self)
    def __call__(self):
        pass
    
    
    


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
    Kbm : None | float = None # Mean to B-basis Knock-down factor
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
    F3c : float # Sigma_33 traction stress allowable
    Kbasis : float # Mean to K-basis factor
    ekdf : float # Environmental Knock-Down Factor
    
    def core_shear_strength(self,
                            thickness : float,
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
        
        
        # According to CMH-17.Vol 6 eqs 4.6.2(b)
        # --------------------------------------
        # Stresses
        tau_xz = Q13 / thickness
        tau_yz = Q23 / thickness
        tau_xz_all = self.F13
        tau_yz_all = self.F23
        
        # Computation of Reserve Factors
        RF_tau_xz = tau_xz_all / (safety_factor * tau_xz)
        RF_tau_yz = tau_yz_all / (safety_factor * tau_yz)
        
        # Core Shear Interaction R.F
        RF_shear = 1 / math.sqrt((1/RF_tau_xz)**N + (1/RF_tau_yz)**N)
        
        # Results
        
        return RF_tau_xz, RF_tau_yz , RF_shear
    
        
@dataclass
class Laminate:
    stacking : list[Ply]
    
    '''
    [Idea] Formas de dar el stacking:
        - Lista : [Tupla(Ply,angulo)]
        - Dict  : {Ply, secuencia)
                   
                  
    
    '''
    @property
    def thickness(self):
        t = sum(ply.material.thickness for ply in self.stacking)
        return t
    
    def stiffness_matrices(self , z0 : float = 0.0):
        if z0==0.0:
            z = [-1*self.thickness/2]
        # Computing the midplane height z
        else :
            z = [z0]

        for ply in self.stacking:
            z.append(z[-1] + ply.material.thickness)

        
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

@dataclass
class Sandwich:
    top_facesheet : Laminate
    core : Honeycomb
    bot_facesheet : Laminate

    @property
    def thickness(self):
        top_facesheet_t = self.top_facesheet.thickness
        bot_facesheet_t = self.bot_facesheet.thickness
        core_t = self.core.thickness
        return bot_facesheet_t + core_t + top_facesheet_t
        

    def stiffness_matrices(self):
        core_t = self.core.thickness
        A_top, B_top, D_top= self.top_facesheet.stiffness_matrices(z0 = -self.thickness/2)
        A_bot, B_bot, D_bot = self.bot_facesheet.stiffness_matrices(z0 = -self.thickness/2)

        A = A_top + A_bot
        B = B_top + B_bot
        D = D_top + D_bot

        K_sx = self.core.G13 * core_t
        K_sy = self.core.G23 * core_t
        return A, B , D , K_sx, K_sy
    
    def core_strength(self,
                      Qx : float,
                      Qy : float,
                      N : float,
                      safety_factor = 1.0,
                      thick_approach = 'conservative'
                      ):
        
        if thick_approach.upper == 'CONSERVATIVE':
            t = self.core.thickness
        else:
            t = self.core.thickness + 0.5*(self.top_facesheet.thickness + self.bot_facesheet.thickness)
            
        RF_tau_xz, RF_tau_yz, RF_shear = self.core.core_shear_strength(
                                        thickness = t,
                                        Q13 = Qx,
                                        Q23 = Qy,
                                        N = N,
                                        safety_factor = safety_factor)
        
        return RF_tau_xz, RF_tau_yz , RF_shear
                                        




def solve_midplane( laminate : Laminate ,N: np.ndarray, M: Optional[np.ndarray] = None) -> np.ndarray:
    rhs = np.concatenate([N,M])
    sol = np.linalg.solve(laminate.ABD_matrices,rhs)
    eps = sol[:3]
    kappa = sol[3:]
    return eps, kappa


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

Core_N636 = Honeycomb(
    name = "Honeycomb_N636_70",
    thickness = 6.0,
    cell_size = 0.0,
    G13 = 26.2,
    G23 = 66.6,
    E3 = 25.8,
    F13 = 0.495,
    F23 = 0.664,
    F3t = 2.565,
    F3c = 0.865,
    Kbasis =0.740, # Mean to K-basis factor
    ekdf =1.0
    )

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

a0,b0,d0 = l0.stiffness_matrices()
a1,b1,d1 = l0.stiffness_matrices()
a3,b3,d3 = l3.stiffness_matrices()
a4,b4,d4 = l4.stiffness_matrices()


sandwich1 = Sandwich(
    top_facesheet = l0,
    core = Core_N636,
    bot_facesheet= l0,
)
        

