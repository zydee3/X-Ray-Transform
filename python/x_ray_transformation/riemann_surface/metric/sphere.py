from .base_metric import BaseMetric
from .metric_type import MetricType
import math


class Sphere(BaseMetric):
    def __init__(self, radius):
        super().__init__()
        self.radius = radius
        self.metric_type = MetricType.sphere

    def compute_values(self, x_values: list, y_values: list, precomputed_values: list = None):
        super().compute_values(x_values, y_values)

    def compute_log_g(self, x: int, y: int, precomputed_values: list = None) -> float:
        return math.log(self.radius ** 4 * 4 / (precomputed_values[1] ** 2))

    def compute_dx_log_g(self, x: int, y: int, precomputed_values: list = None) -> float:
        return -4 * x / precomputed_values[1]

    def compute_dy_log_g(self, x: int, y: int, precomputed_values: list = None) -> float:
        return -4 * y / precomputed_values[1]

    def compute_curvature(self, x: int, y: int, precomputed_values: list = None) -> float:
        return 1 / precomputed_values[0]

    def display(self, x_values: list, y_values: list):
        super().display()
