# -*- coding: utf-8 -*-
from dataclasses import dataclass, asdict

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
class RF:
    location: None | str = None
    material: None | str = None
    load_case: None | str = None
    type_parameter: None | Type_parameter = None
    units_parameter: None | str = None
    kind_of_load: None | Kind_of_load = None
    ultimate_value: None | float = None
    allowable_value: None | float = None
    RF: None | float = None
    remarks: None | str = None

    # [Posibles Mejoras:]
    # - Una función que te transforme en un DataFrame
# <------------------------------------------------[ RF_ACD class and auxiliar]
