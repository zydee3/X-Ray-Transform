from numba.experimental import jitclass
from numba import njit, types
from numpy import sin, cos


members = [
    ('radius', types.double),
    ('amplitude', types.double),
    ('cycles', types.double)
]


@jitclass(members)
class Cosine:
    def __init__(self, radius, amplitude, cycles):
        self.radius = radius
        self.amplitude = amplitude
        self.cycles = cycles

    def compute_boundary(self, theta):
        return self.radius + self.amplitude * cos(self.cycles * theta)

    def compute_dtheta_boundary(self, theta):
        return -1 * self.cycles * self.amplitude * sin(self.cycles * theta)

    def compute_ddtheta_boundary(self, theta):
        return -1 * (self.cycles ** 2) * self.amplitude * cos(self.cycles * theta)
