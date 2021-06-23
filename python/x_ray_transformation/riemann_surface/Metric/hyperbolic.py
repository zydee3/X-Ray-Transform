from .base_metric import Metric
import math


class Hyperbolic(Metric):
    def __init__(self, radius):
        super().__init__()
        self.radius = radius

    def compute_values(self, x_values, y_values):
        super().compute_values(x_values, y_values)

    def compute_log_g(self, x, y) -> float:
        difference = difference_of_squares(self.radius, x, y)
        return math.log(self.radius ** 4 * 4) - math.log(difference * difference)

    def compute_dx_log_g(self, x, y) -> float:
        return 4 * x / (difference_of_squares(self.radius, x, y))

    def compute_dy_log_g(self, x, y) -> float:
        return 4 * y / (difference_of_squares(self.radius, x, y))

    def compute_curvature(self, x, y) -> float:
        return -1 / (self.radius ** 2)

    def display(self, x_values, y_values):
        super().display()


def difference_of_squares(a, b, c):
    return (a * a) - (b * b) - (c * c)
