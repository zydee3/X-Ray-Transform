from .base_domain import BaseDomain
from math import sqrt, sin, cos


class Ellipse(BaseDomain):
    def __init__(self, minor_radius, major_radius, theta_offset):
        super().__init__(major_radius if major_radius > minor_radius else minor_radius)
        self.minor_radius = minor_radius
        self.major_radius = major_radius
        self.theta_offset = theta_offset

    def compute_boundary(self, theta):
        delta_theta = theta - self.theta_offset
        distance = sqrt((self.major_radius * cos(delta_theta)) ** 2 + (self.minor_radius * sin(delta_theta)) ** 2)
        return self.minor_radius * self.major_radius / distance

    def compute_dtheta_boundary(self, theta):
        boundary = self.compute_boundary(theta)
        return (boundary ** 3) * sin(2 * (theta - self.theta_offset)) * (self.major_radius ** 2 - self.minor_radius ** 2) / 2 / ((self.minor_radius * self.major_radius) ** 2)

    def compute_ddtheta_boundary(self, theta):
        boundary = self.compute_boundary(theta)
        dtheta_boundary = self.compute_dtheta_boundary(theta)
        return ((3 * dtheta_boundary) ** 2) / boundary + (boundary ** 3) * cos(2 * (theta - self.theta_offset)) * (self.major_radius - self.minor_radius) / ((self.minor_radius * self.major_radius) ** 2)

    def compute_alpha_normal(self, theta):
        return 0
