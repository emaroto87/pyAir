from dataclasses import dataclass
from functools import cached_property
from typing import Final
import numpy as np
import math
from materials.composite import Orthotropic


@dataclass(frozen=True, slots=True)
class Ply:
    material: Orthotropic
    theta_deg: float
    # -------------------------- Initialization -------------------------------

    def __post_init__(self):
        self._validate_inputs()

    def __str__(self):
        return self._summary()

    def __repr__(self):
        return self._summary()

    # ------------------------- Basic Properties ------------------------------
    def Qbar(self):
        Q = self.material.Q
        Q11, Q12, _, Q22, _, Q66 = Q[[0, 0, 0, 1, 1, 2], [0, 1, 2, 1, 2, 2]]
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

        return np.asarray([
            [Q11b, Q12b, Q16b],
            [Q12b, Q22b, Q26b],
            [Q16b, Q26b, Q66b]
        ], dtype=float)

    def _validate_inputs(self) -> None:

        if not isinstance(self.material, Orthotropic):
            raise TypeError(
                "Material is not an Orthotropic object."
            )
        if not isinstance(self.theta_deg, (int, float)):
            raise TypeError(
                "theta_deg is neither a float or integer type.")

        E1 = float(self.material.E1)
        E2 = float(self.material.E2)
        G12 = float(self.material.G12)
        nu12 = float(self.material.nu12)
        nu21 = float(self.material.nu21)

        if E1 <= 0 or E2 <= 0 or G12 <= 0:
            raise ValueError(
                f"Invalid modulus : E1={E1}, E2={E2}, G12={G12}"
                f"(all should be greater than 0.0)."
            )
        d = 1 - nu12 * nu21
        if d <= 0:
            raise ValueError(
                f"Stability condition not meet: 1 - nu12*nu21 = {d} <=0"
            )

    def _summary(self):
        lines = [
            f"Ply{{Material: {self.material.name or 'No name'}, "
            f"thick : {self.material.thickness:.2f}, "
            f"angle : {self.theta_deg:.0f}}}"
        ]
        return "\n".join(lines)
