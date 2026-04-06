from __future__ import annotations

from dataclasses import dataclass, field
from functools import cached_property
from typing import Tuple, Iterable, Dict, List, Union, Sequence, Mapping, TypeAlias
from structural.ply import Ply
import pandas as pd
import numpy as np
import math

StackingLike: TypeAlias = (
    Sequence[Ply]
)


@dataclass
class Laminate:  # -> LAMINATE CLASS
    """
    Composite Laminate using the Classical Laminate Theory (CLT) defined by 
    stacking with the following valid format:
        - stacking_tuple = ((Ply1,angle1),(Ply2,angle2),etc)
        - stacking_list  = [((Ply1,angle1),(Ply2,angle2),etc)]
        - stacking_dict = {Ply1: angle1, Ply2: angle2, etc}

    The convenction of axis is as follows:
        - The miplane is located at z = 0.0
        - By default, the computation of the properties are done at midplane

    """

    stacking: StackingLike
    name: str | None = None
    numeric_tol: float = 1e-12

    _thickness: float = field(init=False, repr=False)
    _is_symmetric: bool = field(init=False, repr=False)
    _A_B_D: tuple = field(init=False, repr=False)
    _D_star: np.array = field(init=False, repr=False)

    # -------------------------- Initialization -------------------------------
    def __post_init__(self) -> None:
        self._validate_inputs()
        self._thickness = sum(ply.material.thickness for ply in self.stacking)
        self._is_symmetric = self._compute_stacking_symmetry()

    def __str__(self):
        return self.summary()

    def __repr__(self):
        return self.summary()

    # ------------------------- Basic Properties ------------------------------

    @property
    def thickness(self) -> float:
        return self._thickness

    @property
    def is_symmetric(self):
        return self._is_symmetric

    @property
    def stacking_table(self):
        data = [
            (ply.material.name, ply.material.thickness, ply.theta_deg)
            for ply in self.stacking
        ]
        table = pd.DataFrame(
            data=data,
            columns=['Material', 'Thickness', 'Angle(deg)']
        )
        table.index.name = 'Ply'
        return table

    @cached_property
    def A_B_D(self):
        return self._compute_A_B_D()

    @cached_property
    def ABD_matrices(self):
        """Matriz ABD in 6x6 block"""
        A, B, D = self.A_B_D
        return np.block([[A, B], [B, D]])

    @cached_property
    def A_inv(self):
        A, B, D = self.A_B_D
        return np.linalg.solve(A, np.eye(3))

    @cached_property
    def Abar(self):
        A, _, _ = self.A_B_D
        return A / self.thickness

    @cached_property
    def compliance_matrix(self):
        return np.linalg.solve(self.Abar, np.eye(3))

    @cached_property
    def D_star(self):
        return self._compute_D_star()

    # ------------------------- Public Functions  -----------------------------

    def z_interfaces(self, z0: float = 0.0) -> list:
        z = [z0 - self.thickness/2] if z0 == 0 else [z0]
        for ply in self.stacking:
            z.append(z[-1] + ply.material.thickness)
        return z

    def laminate_apparent_moduli(self, angle_deg=0):
        S = self.compliance_matrix
        S11 = S[0, 0]
        S12 = S[0, 1]
        S22 = S[1, 1]
        S66 = S[2, 2]

        # Tranformation in other coordinates
        theta = math.radians(angle_deg)
        c = math.cos(theta)
        s = math.sin(theta)
        c4 = c**4
        c2 = c**2
        s4 = s**4
        s2 = s**2

        S11p = S11*c4 + S22*s4 + (2*S12 + S66)*c2*s2
        S12p = (S11 + S22 - S66)*c2*s2 + S12*(c4 + s4)
        S22p = S11*s4 + S22*c4 + (2*S12 + S66)*c2*s2
        S66p = (S11 + S22 - 2*S12 - S66)*c2*s2 + S66*(c4 + s4)

        Ex = 1.0 / S11p
        Ey = 1.0 / S22p
        Gxy = 1.0 / S66p
        nuxy = - S12p / S11p

        return dict(Ex=Ex, Ey=Ey, Gxy=Gxy, nuxy=nuxy)

    def dimpling_strength(self, core_cell_size: float, Kb: float = 1.0):

        tf = self.thickness
        cs = core_cell_size
        A, B, D = self.A_B_D
        pi = math.pi

        # Computation of compression strength in dimpling (Fc), CMH-17.Vol6
        # Eq 4.6.5.1 (c):
        # -----------------------------------------------------------------
        if self.is_symmetric:
            D_prime = D
        else:
            D_prime = self.D_star

        D11 = D_prime[0, 0]
        D12 = D_prime[0, 1]
        D22 = D_prime[1, 1]
        D66 = D_prime[2, 2]

        Fc = Kb * (1 / tf) * ((pi/cs)**2) * (D11 + 2 * (D12 + 2 * D66) + D22)

        # Computation of the shear strength in dimpling (Fs), CMH-17. Vol6
        # Eq. 4.6.5.3
        # ----------------------------------------------------------------
        eqv_moduli = self.laminate_apparent_moduli()
        Ex = eqv_moduli['Ex']
        Ey = eqv_moduli['Ey']
        Fs = Kb * 0.6 * min(Ex, Ey) * (tf/cs)**(1.5)

        return Fc, Fs

    def summary(self):

        eqv = self.laminate_apparent_moduli()
        lines = [
            f"Laminate: {self.name or 'No name'}",
            f"    - Total thickness (t) : {self.thickness:.2f}",
            f"    - Symmetric : {self.is_symmetric}",
            f"    - Aparent modulus (theta=0.0):",
            f"        E1 = {eqv['Ex'] : .0f} ",
            f"        E2 = {eqv['Ey'] : .0f} ",
            f"        E12 = {eqv['Gxy'] : .0f} ",
            f"        nu12 = {eqv['nuxy'] : .2f} ",
        ]
        return "\n".join(lines)

    # ------------------------- Private Functions  ----------------------------
    def _compute_D_star(self):
        """
        [D*] = [D] - [B][A^{-1}][B] for buckling/dimpling of non-symmetric
        laminates (B!=0)
        """
        A, B, D = self.A_B_D
        A_inv = self.A_inv
        return D - B @ A_inv @ B

    def _compute_stacking_symmetry(self):
        stacking = self.stacking
        tol = self.numeric_tol
        n = len(stacking)
        mid = n // 2
        return all(
            stacking[i].material == stacking[-1-i].material
            and abs(stacking[i].theta_deg - stacking[-1-i].theta_deg) < tol
            for i in range(mid)
        )

    def _compute_A_B_D(self,
                       z0: float = 0.0,
                       decimals: float | None = None) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
        """
        Compute (and cache) A,B,D matrices using the midplane as reference
        (z=0). For a symmetric laminate, B->0 within numeric tolerances.
        """
        # Computing ply interfaces
        z = self.z_interfaces(z0=z0)

        # Computing A, B and D matrices
        A = np.zeros((3, 3), dtype=float)  # Pure Membrane stiffness Matrix
        B = np.zeros((3, 3), dtype=float)  # Pure Bending stiffness Matrix
        # Membrane-Bending Coupling stiffness Matrix
        D = np.zeros((3, 3), dtype=float)

        for k, ply in enumerate(self.stacking, start=1):
            Qbar = np.asarray(ply.Qbar, dtype=float)
            z_bot = z[k-1]
            z_top = z[k]
            A += Qbar * (z_top - z_bot)
            B += (1/2) * Qbar * (z_top**2 - z_bot**2)
            D += (1/3) * Qbar * (z_top**3 - z_bot**3)
        if decimals:
            return np.round(A, decimals), np.round(B, decimals), np.round(D, decimals)
        else:
            return A, B, D

    def _validate_inputs(self) -> None:

        if isinstance(self.stacking, (tuple, list)):
            if isinstance(self.stacking, list):
                self.stacking = tuple(self.stacking)
        else:
            raise TypeError(
                "Stacking is not a sequence (tuple or list) of Ply objects."
            )

        if len(self.stacking) == 0:
            raise ValueError(
                "Laminate stacking must contain at least one ply"
            )
        # if not all(
        #         isinstance(ply, Ply)
        #         and isinstance(angle, float)
        #         for ply, angle in self.stacking
        # ):
        #     raise ValueError(
        #         "Laminate stacking member is not (Ply,float) type entry"
        #     )
