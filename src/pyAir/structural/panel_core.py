from dataclasses import dataclass
from materials.composite import Honeycomb


@dataclass(frozen=True, slots=True)
class Core:
    material: Honeycomb
    theta_deg: float
