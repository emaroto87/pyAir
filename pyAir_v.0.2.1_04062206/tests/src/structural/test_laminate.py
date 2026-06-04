from src.structural.laminate import Laminate
import pytest
import numpy as np
from dataclasses import dataclass


@dataclass
class MockMaterial:
    name: str
    thickness: float


@dataclass
class MockPly:
    material: MockMaterial
    theta_deg: float
    Qbar: np.ndarray


# ============================== FIXTURES =====================================


@pytest.fixture
def simple_laminate():
    """Symmetric Laminate with 2 plies"""
    E = 150e9
    nu = 0.25
    E_star = E/(1-nu**2)
    Q = np.array(
        [E_star, nu*E_star, 0.0],
        [nu*E_star, E_star, 0.0],
        [0.0, 0.0, E/(2*(1+nu))]
    )
    ply1 = MockPly(MockMaterial('Mat1', 0.125), 0.0, Q)
    ply2 = MockPly(MockMaterial('Mat1', 0.125), 0.0, Q)
    return Laminate([ply1, ply2], name="TestLamSym")


@pytest.fixture
def nonsymmetric_laminate():
    """ Non-symmetric Laminate in order to see if B!=0"""
    E = 150e9
    nu = 0.25
    E_star = E/(1-nu**2)
    Q = np.array(
        [E_star, nu*E_star, 0.0],
        [nu*E_star, E_star, 0.0],
        [0.0, 0.0, E/(2*(1+nu))]
    )
    ply1 = MockPly(MockMaterial('Mat1', 0.125), 0.0, Q)
    ply2 = MockPly(MockMaterial('Mat1', 0.125), 45.0, Q)
    return Laminate([ply1, ply2], name="TestLamAsym")


# ================================ TESTS =====================================


def test_initialization(simple_laminate):
    lam = simple_laminate
    assert lam.thickness == pytest.approx(0.25)
    assert lam.name == "TestLamSym"
    assert lam.is_symmetric is True


def test_z_interfaces(simple_laminate):
    lam = simple_laminate
    z = lam.z_interfaces()
    # It should be 3 interfaces, (-t/2), 0, (t/2)
    assert len(z) == 3
    assert z[0] == pytest.approx(-0.125)
    assert z[-1] == pytest.approx(0.125)


def test_A_B_D_shapes(simple_laminate):
    A, B, D = simple_laminate.A_B_D
    assert A.shape == (3, 3)
    assert B.shape == (3, 3)
    assert D.shape == (3, 3)


def test_symmetric_laminate_has_zero_B(simple_laminate):
    _, B, _ = simple_laminate.A_B_D
    assert np.allclose(B, 0, atol=1e-9)


def test_nonsymmetric_laminate_has_nonzero_B(nonsymmetric_laminate):
    _, B, _ = nonsymmetric_laminate.A_B_D
    assert not np.allclose(B, 0, atol=1e-9)


def test_compliance_matrix(simple_laminate):
    lam = simple_laminate
    S = lam.compliance_matrix
    assert S.shape == (3, 3)
    # S should be the inverse of Abar
    assert np.allclose(S @ lam.Abar, np.eye(3), atol=1e-8)


def test_laminate_moduli(simple_laminate):
    mod = simple_laminate.laminate_apparent_moduli()
    assert "Ex" in mod
    assert "Ey" in mod
    assert "Gxy" in mod
    assert "nuxy" in mod
    assert mod["Ex"] > 0
    assert mod["Ey"] > 0
    assert -1.0 < mod["nuxy"] < 1.0


def test_dimpling_strength(simple_laminate):
    s = simple_laminate.summary()
    assert isinstance(s, str)
    assert "Laminate" in s
    assert "thickness" in s
