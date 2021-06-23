from .base_metric import BaseMetric
import math


class Sphere(BaseMetric):
    def __init__(self, radius):
        super().__init__()
        self.radius = radius

    def compute_values(self, x_values: list, y_values: list):
        super().compute_values(x_values, y_values)

    def compute_log_g(self, x: int, y: int) -> float:
        return math.log(self.radius ** 4 * 4) - 2 * math.log(x * x + y * y )

    def compute_dx_log_g(self, x: int, y: int) -> float:
        return -4 * x / (sum_of_squares(x, y, self.radius))

    def compute_dy_log_g(self, x: int, y: int) -> float:
        return -4 * y / (sum_of_squares(x, y, self.radius))

    def compute_curvature(self, x: int, y: int) -> float:
        return 1 / (self.radius ** 2)

    def display(self, x_values: list, y_values: list):
        super().display()


def sum_of_squares(a: int, b: int, c: int) -> int:
    return (a * a) + (b * b) + (c * c)
