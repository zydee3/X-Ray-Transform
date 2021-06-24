from .base_domain import BaseDomain
from math import floor, sin, cos


class Cosine(BaseDomain):
    def __init__(self, radius: int, amplitude: float, cycles: int):
        super().__init__(radius + amplitude)
        self.radius = radius
        self.amplitude = amplitude
        self.cycles = floor(cycles)

    def compute_boundary(self, theta):
        return self.radius + self.amplitude * cos(self.cycles * theta)

    def compute_dtheta_boundary(self, theta):
        return -1 * self.cycles * self.amplitude * sin(self.cycles * theta)

    def compute_ddtheta_boundary(self, theta):
        return -1 * (self.cycles ** 2) * self.amplitude * cos(self.cycles * theta)

    def compute_alpha_normal(self, theta):
        return super().compute_alpha_normal(theta)
