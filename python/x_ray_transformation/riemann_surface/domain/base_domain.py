from abc import ABC, abstractmethod


class BaseDomain(ABC):
    def __init__(self, max_radius: int):
        self.max_radius = max_radius

    @abstractmethod
    def compute_boundary(self, theta):
        return 1

    @abstractmethod
    def compute_dtheta_boundary(self, theta):
        return 0

    @abstractmethod
    def compute_ddtheta_boundary(self, theta):
        return 0

    @abstractmethod
    def compute_alpha_normal(self, theta):
        return 0
