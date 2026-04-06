@dataclass
class Forpla:
    load_case: str
    Nx: float = 0.0
    Ny: float = 0.0
    Nxy: float = 0.0
    Mx: float = 0.0
    My: float = 0.0
    Mxy: float = 0.0
    Qx: float = 0.0
    Qy: float = 0.0
    Pr: float = 0.0

    # def __str__(self):
    #     return vars(self).items()
    def to_array(self):
        vec = [v for k, v in vars(self).items() if not k == 'load_case']
        return np.array(vec)
