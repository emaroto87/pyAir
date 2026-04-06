from materials.metallic import Metallic
from structural.laminate import Laminate


class Plate:
    pass


class MetallicPlate:
    def __init__(self,
                 material: Metallic,
                 thickness: float,
                 plate_id: int = None,
                 linked_FEM_eid: int = None,
                 linked_FEM_pid: int = None):

        self.material = material
        self.thickness = thickness
        self.linked_FEM_eid = linked_FEM_eid
        self.linked_FEM_pid = linked_FEM_pid
        self.plate_id = plate_id

    def __repr__(self):
        msg = f'{self.material.name} Metallic Plate ID {self.plate_id}' \
              f'with thick = {self.thickness}. Linked ElemID :' \
              f'{self.linked_FEM_eid}'
        return msg


class CompositePlate:
    def __init__(self,
                 laminate: Laminate = None,
                 plate_id: int = None,
                 linked_FEM_eid: int = None,
                 linked_FEM_pid: int = None):

        self.laminate = laminate
        self.linked_FEM_eid = linked_FEM_eid
        self.linked_FEM_pid = linked_FEM_pid
        self.plate_id = plate_id
        self.thickness = laminate.thickness

    def get_equivalent_moduli(self, angle_deg: float = 0.0):
        '''
        Returns the equivalent values of Ex, Ey, Gxy and nuxy into another
        coordinate system which is defined as a rotation from the stacking
        coordinate system
        '''
        data = self.laminate.laminate_apparent_moduli(angle_deg)
        Ex_local, Ey_local, Gxy_local, nuxy_local = data.values()

        return Ex_local, Ey_local, Gxy_local, nuxy_local
