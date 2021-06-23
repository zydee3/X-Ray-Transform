from abc import ABC, abstractmethod


class BaseDomain(ABC):
    def __init__(self, max_radius: int):
        self.boundaries = []
        self.dx_boundaries = []
        self.ddx_boundaries = []
        self.max_radius = max_radius
