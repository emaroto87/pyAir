from materials.composite import Orthotropic, Honeycomb
from structural.ply import Ply
from structural.laminate import Laminate
from structural.panel_core import Core
from analysis.composite_laminate_plain_strength import Sandwich
import numpy as np

'''
Comparativa de resultados con los obtenidos a través de SANDRES
Ver archivo SANDRES_Panel_Mono_worksapce.out

'''
IMA21E_HW_70 = Orthotropic(
    name='IMA-M21E-UD-70',
    thickness=0.184,
    E1=154000,
    E2=8500,
    G12=4200,
    nu12=0.350,
    u_et_all=15500,
    u_ec_all=7130
)

T300PW_HW_70 = Orthotropic(
    name='T300PW_HW_70',
    thickness=0.237,
    E1=50000,
    E2=50000,
    G12=1000,
    nu12=0.05,
    u_et_all=9100,
    u_ec_all=8160
)

Honeycomb_N636_70 = Honeycomb(
    name="Honeycomb_N636_70",
    thickness=6.0,
    cell_size=4.8,
    G13=26.2,
    G23=66.6,
    E3=25.8,
    F13=0.495,
    F23=0.664,
    F3t=2.565,
    F3c=0.865,
    Kbasis=0.740,  # Mean to K-basis factor
    ekdf=1.0
)

Core_N636_0 = Core(Honeycomb_N636_70, 0)

stacking_0 = [
    Ply(T300PW_HW_70, 0),
    Ply(T300PW_HW_70, 45),
    Ply(T300PW_HW_70, 0)
]

# Test_1
stacking_1 = [
    Ply(T300PW_HW_70, 45),
    Ply(T300PW_HW_70, 0),
    Ply(T300PW_HW_70, 45),
    Ply(T300PW_HW_70, 0),
    Ply(T300PW_HW_70, 45),
    # SYM
    Ply(T300PW_HW_70, 0),
    # SYM
    Ply(T300PW_HW_70, 45),
    Ply(T300PW_HW_70, 0),
    Ply(T300PW_HW_70, 45),
    Ply(T300PW_HW_70, 0),
    Ply(T300PW_HW_70, 45)
]

# Test_2
stacking_2 = [
    Ply(IMA21E_HW_70, 0),
    Ply(IMA21E_HW_70, 45),
    Ply(IMA21E_HW_70, -45),
    Ply(IMA21E_HW_70, 0),
    Ply(IMA21E_HW_70, 45),
    Ply(IMA21E_HW_70, -45),
    Ply(IMA21E_HW_70, 0),
    Ply(IMA21E_HW_70, -45),
    Ply(IMA21E_HW_70, 45),
    Ply(IMA21E_HW_70, 0),
    Ply(IMA21E_HW_70, -45),
    Ply(IMA21E_HW_70, 45),
    Ply(IMA21E_HW_70, 0)
]


# Test_3
stacking_3 = [
    Ply(IMA21E_HW_70, 0),
    Ply(IMA21E_HW_70, 90),
    Ply(IMA21E_HW_70, 45),
    Ply(IMA21E_HW_70, -45),
    Ply(IMA21E_HW_70, 0),
    Ply(IMA21E_HW_70, 90),
    Ply(IMA21E_HW_70, 45),
    Ply(IMA21E_HW_70, -45),
    # ---SYM
    Ply(IMA21E_HW_70, 0),
    # ---SYM
    Ply(IMA21E_HW_70, -45),
    Ply(IMA21E_HW_70, 45),
    Ply(IMA21E_HW_70, 90),
    Ply(IMA21E_HW_70, 0),
    Ply(IMA21E_HW_70, -45),
    Ply(IMA21E_HW_70, 45),
    Ply(IMA21E_HW_70, 90),
    Ply(IMA21E_HW_70, 0),
]

# Test_4
stacking_4 = [
    Ply(IMA21E_HW_70, 0),
    Ply(IMA21E_HW_70, 0),
    Ply(IMA21E_HW_70, 0),
    Ply(IMA21E_HW_70, 0),
    Ply(IMA21E_HW_70, 0),
]

l0 = Laminate(stacking_0)
l1 = Laminate(stacking_1)
l2 = Laminate(stacking_2)
l3 = Laminate(stacking_3)
l4 = Laminate(stacking_4)
l5 = Laminate(
    [
        Ply(T300PW_HW_70, 45),
        Ply(T300PW_HW_70, 45),
        Ply(T300PW_HW_70, 45)
    ]
)
l6 = Laminate(
    [
        Ply(T300PW_HW_70, 0),
        Ply(T300PW_HW_70, 0),
        Ply(T300PW_HW_70, 0)
    ]
)


sandwich1 = Sandwich(
    top_facesheet=l0,
    core=Core_N636_0,
    bot_facesheet=l0,
)
sandwich_unsym = Sandwich(
    top_facesheet=l5,
    core=Core_N636_0,
    bot_facesheet=l6,
)

# IMA21E_HW_70_UP
IMA21E_HW_70_UP = Laminate([
    Ply(IMA21E_HW_70, 0),
    Ply(IMA21E_HW_70, 45),
    Ply(IMA21E_HW_70, 90),
    Ply(IMA21E_HW_70, 0),
    Ply(IMA21E_HW_70, 45)
])
IMA21E_HW_70_LO = Laminate([
    Ply(IMA21E_HW_70, 90),
    Ply(IMA21E_HW_70, 90),
    Ply(IMA21E_HW_70, 45),
    Ply(IMA21E_HW_70, 45),
    Ply(IMA21E_HW_70, 0)
])

SANDRES_Panel_ASym_Thin_Dimpling = Sandwich(
    top_facesheet=IMA21E_HW_70_UP,
    core=Core_N636_0,
    bot_facesheet=IMA21E_HW_70_LO
)
M = np.array([1000.0, 0.0, 0.0])
N = np.array([0.0, 0.0, 0.0])
SANDRES_Panel_ASym_Thin_Dimpling.dimpling_RFs(
    N=N,
    M=M,
    Kb=0.74
)
SANDRES_Panel_ASym_Thin_Dimpling.top_facesheet.dimpling_strength(
    core_cell_size=SANDRES_Panel_ASym_Thin_Dimpling.core.material.cell_size,
    Kb=0.74
)
