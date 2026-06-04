# Nos aseguramos que cargue los modulos
import pytest
import numpy as np
from analysis.composite_laminate_plain_strength import laminate_plain_strength, solve_midplane
from analysis.composite_sandwich_strength import Sandwich
from structural.panel_core import Core
from structural.laminate import Laminate
from structural.ply import Ply
from materials.composite import Orthotropic, Honeycomb
import sys
from pathlib import Path

PROJECT_ROOT = Path(__file__).resolve().parents[3]
SRC = PROJECT_ROOT / "src"
sys.path.insert(0, str(PROJECT_ROOT))


"""
v.0.2.1
-------
    - Comentarios:
        * Se han añadido test para verificar que los cálculos son correctos. 
        Se ha comprobado que efectivamente está todo bien
"""


@pytest.fixture
def IMA21E_HW_70():
    return Orthotropic(
        name='IMA-M21E-UD-70',
        thickness=0.184,
        E1=154000,
        E2=8500,
        G12=4200,
        nu12=0.350,
        u_et_all=15500,
        u_ec_all=7130
    )


@pytest.fixture
def MATERIAL_01():
    return Orthotropic(
        name='MATERIAL-01',
        thickness=0.180,
        E1=138000,
        E2=8900,
        G12=4700,
        nu12=0.3,
    )


@pytest.fixture
def laminate_num2(IMA21E_HW_70):
    """ Laminate IMA21E_HW_70 with stacking
       [0/45/-45/0/45/-45/0/-45/45/0/-45/45/0]
    """
    return Laminate([
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
    ])


@pytest.fixture
def ancolab_example1(MATERIAL_01):
    """ Laminate MATERIAL-01 with stacking 
       [45/-45/0/90/0/90/-45/45]s
    """
    return [
        Ply(MATERIAL_01, 45),
        Ply(MATERIAL_01, -45),
        Ply(MATERIAL_01, 0),
        Ply(MATERIAL_01, 90),
        Ply(MATERIAL_01, 0),
        Ply(MATERIAL_01, 90),
        Ply(MATERIAL_01, -45),
        Ply(MATERIAL_01, 45),
        # ---
        Ply(MATERIAL_01, 45),
        Ply(MATERIAL_01, -45),
        Ply(MATERIAL_01, 90),
        Ply(MATERIAL_01, 0),
        Ply(MATERIAL_01, 90),
        Ply(MATERIAL_01, 0),
        Ply(MATERIAL_01, -45),
        Ply(MATERIAL_01, 45),
    ]


def test_solve_midplane_onlyN(ancolab_example1):
    N = np.array([100.0, 0.0, 0.0])
    M = None
    eps0, kappa0 = solve_midplane(ancolab_example1, N, M)
    u_eps0 = eps0*1e+6
    u_kappa0 = kappa0*1e+6
    assert u_eps0[0] == pytest.approx(657.5)
    assert u_eps0[1] == pytest.approx(-20.4)
    assert u_eps0[2] == pytest.approx(0.0)
    assert u_kappa0[0] == pytest.approx(0.0)
    assert u_kappa0[1] == pytest.approx(0.0)
    assert u_kappa0[2] == pytest.approx(0.0)


def test_solve_midplane_ancolab_ex1(ancolab_example1):
    N = np.array([100.0, 0.0, 0.0])
    M = np.array([100.0, 0.0, 0.0])
    eps0, kappa0 = solve_midplane(ancolab_example1, N, M)
    u_eps0 = eps0*1e+6
    u_kappa0 = kappa0*1e+6
    assert u_eps0[0] == pytest.approx(657.5)
    assert u_eps0[1] == pytest.approx(-20.4)
    assert u_eps0[2] == pytest.approx(0.0)
    assert u_kappa0[0] == pytest.approx(95.5)
    assert u_kappa0[1] == pytest.approx(40.3)
    assert u_kappa0[2] == pytest.approx(0.0)


def test_solve_midplane_ancolab_ex3(ancolab_example1):
    N = np.array([1318.389, 1063.07, 179.841])
    M = np.array([0.0, 0.0, 0.0])
    eps0, kappa0 = solve_midplane(ancolab_example1, N, M)
    u_eps0 = eps0*1e+6
    u_kappa0 = kappa0*1e+6
    assert u_eps0[0] == pytest.approx(650.0)
    assert u_eps0[1] == pytest.approx(430.0)
    assert u_eps0[2] == pytest.approx(310.)
    assert u_kappa0[0] == pytest.approx(0.0)
    assert u_kappa0[1] == pytest.approx(0.0)
    assert u_kappa0[2] == pytest.approx(0.0)


MATERIAL_01 = Orthotropic(
    name='MATERIAL-01',
    thickness=0.180,
    E1=138000,
    E2=8900,
    G12=4700,
    nu12=0.3,
)

stacking_ex_ancolab_1 = [
    Ply(MATERIAL_01, 45),
    Ply(MATERIAL_01, -45),
    Ply(MATERIAL_01, 0),
    Ply(MATERIAL_01, 90),
    Ply(MATERIAL_01, 0),
    Ply(MATERIAL_01, 90),
    Ply(MATERIAL_01, -45),
    Ply(MATERIAL_01, 45),
    # ---
    Ply(MATERIAL_01, 45),
    Ply(MATERIAL_01, -45),
    Ply(MATERIAL_01, 90),
    Ply(MATERIAL_01, 0),
    Ply(MATERIAL_01, 90),
    Ply(MATERIAL_01, 0),
    Ply(MATERIAL_01, -45),
    Ply(MATERIAL_01, 45),
]
N = np.array([1318.389, 1063.07, 179.841])
M = np.array([0.0, 0.0, 0.0])
laminate_ex_ancolab_1 = Laminate(stacking_ex_ancolab_1)
eps0, kappa0 = solve_midplane(laminate_ex_ancolab_1, N, M)
