from geometry.vector import Vector
from geometry.point import Point
from enum import Enum
__HEAD_TYPE_PTR = 'PTR'
__HEAD_TYPE_CSK = 'CSK'


class HeadType(Enum):
    PTR = 'protruding'
    CSK = 'countersunk'


class Fastener:
    def __init__(self,
                 Et: float,
                 diameter: float,
                 rivet=True,
                 head_type: HeadType | None = None,
                 ref_point: Point | None = None,
                 shank_vector: Vector | None = None,
                 orientation_vector: Vector | None = None,
                 fastener_id: int | None = None,
                 fastener_label: str | None = None,
                 ):
        self.Et = Et
        self.diameter = diameter
        self.head_type = head_type
        self.ref_point = ref_point
        self.shank_vector = shank_vector
        self.orientation_vector = orientation_vector
        self.fastener_id = fastener_id
        self.fastener_label = fastener_label

        if rivet:
            self.is_bolt = False
            self.is_rivet = True
        else:
            self.is_bolt = True
            self.is_rivet = False

    def __str__(self):
        msg = f'Fastener['
        if not self.fastener_id:
            msg += f' ID:{self.fastener_id}'
        if not self.fastener_label:
            msg += f' Label:{self.fastener_label}'
        msg += f' D: {self.diameter} ]'
        return msg
