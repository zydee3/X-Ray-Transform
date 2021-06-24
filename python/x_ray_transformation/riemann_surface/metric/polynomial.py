from .base_metric import BaseMetric
from .metric_type import MetricType
import math


class Polynomial(BaseMetric):
    def __init__(self, coefficients):
        super().__init__()
        if len(coefficients) != 6:
            raise Exception(f'Invalid number of coefficients (Expect: 6, Actual: {len(coefficients)}).')
        self.coefficients = coefficients
        self.metric_type = MetricType.polynomial

    def compute_values(self, x_values: list, y_values: list):
        super().compute_values(x_values, y_values)

    def compute_log_g(self, x: int, y: int, precomputed_values: list = None) -> float:
        cfs = self.coefficients
        return ((cfs[0] * x) ** 2) + (cfs[1] * x * y) + ((cfs[2] * y) ** 2) + (cfs[3] * x) + (cfs[4] * y) + cfs[5]

    def compute_dx_log_g(self, x: int, y: int, precomputed_values: list = None) -> float:
        cfs = self.coefficients
        return (2 * cfs[0] * x) + (cfs[1] * y) + cfs[3]

    def compute_dy_log_g(self, x: int, y: int, precomputed_values: list = None) -> float:
        cfs = self.coefficients
        return (cfs[1] * x) + (2 * cfs[2] * y) + cfs[4]

    def compute_curvature(self, x: int, y: int, precomputed_values: list = None) -> float:
        return -1 * (math.e ** self.compute_log_g(x, y) * (self.coefficients[0] + self.coefficients[1]))

    def display(self, x_values: list, y_values: list):
        super().display()
