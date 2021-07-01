from numba.experimental import jitclass
from numba import njit, types


members = [
    ('radius', types.double)
]


@jitclass(members)
class Circle:
    def __init__(self, radius: int):
        self.radius = radius

    def compute_boundary(self, theta):
        return self.radius

    def compute_dtheta_boundary(self, theta):
        return 0

    def compute_ddtheta_boundary(self, theta):
        return 0

    def compute_alpha_normal(self, theta):
        return math.pi / 2  # (?)