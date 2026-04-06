from __future__ import annotations
from dataclasses import dataclass
from typing import Optional, Final
import numpy as np

# Nota: Typing Final hace que la variable se convierta en "constante" para que
# que en un checkeo "estatico" (antes de correrlo) los analísis de código den


@dataclass(frozen=True, slots=True)
class Orthotropic:
    """
    Model of Orthotropic material for composite laminates under the hypothesis
    of plane stress.

    Notes:
    ------
    - nu21 is computed by reciprocity
        nu21 = nu12 * E2 / E1
    - It is validated the stability of the material:
        1 - nu12 * nu21 >0
    """

    name: str
    thickness: float
    E1: float
    E2: float
    G12: float
    nu12: float
    norm: Optional[str] = None  # Norm or spec
    cond: Optional[str] = None  # Condition Dry/Wet
    # Normally, -55C , RT, +55C, +70C, +90C and +120C
    temp: Optional[str] = None
    Kbm: Optional[float] = 1.0  # Mean to B-basis Knock-down factor
    E1b: Optional[float] = None
    E2b: Optional[float] = None
    G12b: Optional[float] = None
    CTE1: Optional[float] = None  # Coeff. Thermal Exp. at direction 1
    CTE2: Optional[float] = None  # Coeff. Therma Exp.  at direction 2

    # Strength allowables
    u_et_all: Optional[float] = None  # rupture tension micro-strains allowable
    # rupture compresion micro-strains allowable
    u_ec_all: Optional[float] = None
    F1t: Optional[float] = None
    F1c: Optional[float] = None

    # Damage torelance allowables
    # damage torlerance micro-strains tension after impact
    uTAI: Optional[float] = None
    # damage tolerance micro-strains compresion after impact
    uCAI: Optional[float] = None
    # damage tolerance micro-strains bending after impact
    uBAI: Optional[float] = None
    iCDT: Optional[float] = None

# -------------------------- Initialization -------------------------------
    def __post_init__(self):
        # Cheking Modulus
        if self.E1 <= 0 or self.E2 <= 0 or self.G12 <= 0:
            raise ValueError(
                f"Invalid modulus : E1={self.E1}, E2={self.E2}, G12={self.G12}"
                f"(all should be greater than 0.0)."
            )
        if self.thickness <= 0:
            raise ValueError(
                "Thickness should be greater then 0.0 (t = {self.thickness}")

        d = 1-self.nu12 * self.nu21
        if d <= 0:
            raise ValueError(
                f"Stability condition not meet: 1 - nu12*nu21 = {d} <=0"
            )

# ------------------------- Basic Properties ------------------------------
    @property
    def nu21(self) -> float:
        return self.nu12 * self.E2 / self.E1

    @property
    def Q(self) -> np.ndarray:
        """
        Returns the four independent terms of the reduced Q plane-stress.
        No rotation applied.
            (Q11, Q12, Q22, Q66)
        """
        d: Final[float] = 1 - self.material.nu12 * self.material.nu21
        Q11 = self.material.E1/d
        Q12 = self.material.nu12*self.material.E2/d
        Q22 = self.material.E2/d
        Q66 = self.material.G12

        return np.array([
            [Q11, Q12, 0.0],
            [Q12, Q22, 0.0],
            [0.0, 0.0, Q66]
        ], dtype=float)


@dataclass(frozen=True, slots=True)
class Honeycomb:
    name: str
    thickness: float
    cell_size: float
    G13: float
    G23: float
    E3: float
    F13: float  # Tau_13 stress allowable
    F23: float  # Tau_23 stress allowable
    F3t: float  # Sigma_33 traction stress allowable
    F3c: float  # Sigma_33 compression stress allowable
    Kbasis: float  # Mean to K-basis factor
    ekdf: float  # Environmental Knock-Down Factor
