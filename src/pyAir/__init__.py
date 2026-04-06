"""
pyAir - Advanced FEM Toolkit for Structural Analysis (c) 2026 E.Maroto
"""
from importlib.metadata import version

try:
    __version__ = version("pyAir")
except:
    __version__ = "0.2.0-dev"


from .analysis import aida
from .analysis import arpa
from .analysis import composite_laminate_plain_strength
from .analysis import composite_sandwich_strength
from .analysis import metallic_crippling
from .fem.post import op2_reader
from .fem.pre import create_fasteners

__all__ = [
    "aida",
    "arpa",
    "composite_laminate_plain_strength",
    "composite_sandwich_strength",
    "metallic_crippling",
    "op2_reader",
    "create_fasteners",
    "__version__",
]
