from .base_domain import BaseDomain
import math


class Circle(BaseDomain):
    def __init__(self, radius: int):
        super().__init__(radius)
        self.radius = radius

    def compute_boundary(self, theta):
        return radius

    def compute_dtheta_boundary(self, theta):
        return 0

    def compute_ddtheta_boundary(self, theta):
        return 0

    def compute_alpha_normal(self, theta):
        return math.pi / 2  # (?)
