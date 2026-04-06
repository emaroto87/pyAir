
from geometry.point import Point


class Vector(Point):
    def __init__(self, A: Point, B: Point):
        # checking
        # if A.csys_type == B.csys_type:
        # self.__coords = B - A
        # self.csys_type = A.csys_type
        super().__init__(
            (B - A)[0],
            (B - A)[1],
            (B - A)[2],
            A.csys_type
        )
        # else:
        #     raise ValueError(
        #         f'Csys of point A does not math with point B.'
        #     )
