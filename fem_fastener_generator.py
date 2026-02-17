# -*- coding: utf-8 -*-
"""
Created on Tue Feb 17 15:14:59 2026

@author: U69432
"""

'''
Script para generar remaches de manera rápida en NASTRAN

Inputs:
    - Datos de remaches (localización, tipo, vector normal, radio, direción)
    - Malla 2D con propiedades aplicadas
    -

Pasos:
    1. Coger el punto de refencia (P) y su vector de dirección (Vd) y obtener
    los puntos proyectados en las superficies.
    2. Generar CBUSH
    3. Calcular Huth
        3.1 Obtener el apilado del elemento
        3.2 Calcular E1,E2 del apilado en la dirección que se desee
        3.3 Calular Kshear en esa dirección.
    4. Proceso recursivo

'''

__HEAD_TYPE_PTR = 'PTR'
__HEAD_TYPE_CSK = 'CSK'


class Material:
    pass


class Fastener:
    def __init__(self, material: Material, diameter: float, head_type: str):
        self.material = material
        self.diameter = diameter
        self.head_type = head_type


class MetallicPlate:
    def __init__(self, material: Material, thichness: float):
        self.material = material
        self.thickness = thickness


class CompositePlate:
    def __


def huth(plate1: Plate, plate2: Plate, fastener: Fastener)
