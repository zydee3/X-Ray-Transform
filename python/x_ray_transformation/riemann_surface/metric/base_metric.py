from abc import ABC, abstractmethod
from x_ray_transformation.meta.derivative import ddx, dddx, ddy, dddy
from x_ray_transformation.riemann_surface.metric.metric_type import MetricType
import math


class BaseMetric(ABC):

    def __init__(self):
        self.log_g_of_t = None
        self.dx_log_g_of_t = None
        self.dy_log_g_of_t = None
        self.curvature = None

    def __str__(self):
        return f'Metric: {self.__class__.__name__}\n' \
               f'\tlog: {self.log_g_of_t}\n' \
               f'\tddx: {self.dx_log_g_of_t}\n' \
               f'\tddy: {self.dy_log_g_of_t}\n' \
               f'\tcur: {self.curvature}\n'

    @abstractmethod
    def compute_values(self, x_values: list, y_values: list):
        if len(x_values) != len(y_values):
            raise Exception(f'Number of x_values ({len(x_values)}) and y_values ({len(y_values)}) must be the same.')
        if hasattr(self, "widths") and len(x_values) != len(self.widths):
            raise Exception(f'Number of x,y ({len(x_values)}) values and coefficients ({len(self.widths)}) must be the same.')

        num_elements = len(x_values)
        self.log_g_of_t = [None] * num_elements
        self.dx_log_g_of_t = [None] * num_elements
        self.dy_log_g_of_t = [None] * num_elements
        self.curvature = [None] * num_elements
        precomputed_values = [None] * 4

        for index in range(len(x_values)):
            x = x_values[index]
            y = y_values[index]

            if self.metric_type == MetricType.sphere or self.metric_type == MetricType.hyperbolic:
                precomputed_values[0] = self.radius * self.radius
                precomputed_values[1] = x * x + y * y + precomputed_values[0]
            elif self.metric_type == MetricType.gaussian:
                # pre-compute these values so they are only computed once per x,y
                precomputed_values[0] = index
                precomputed_values[1] = self.compute_x_bar(x, index)
                precomputed_values[2] = self.compute_y_bar(y, index)
                precomputed_values[3] = self.compute_fj(index, self.x_bar, self.y_bar)

            self.log_g_of_t[index] = self.compute_log_g(x, y, precomputed_values)
            self.dx_log_g_of_t[index] = self.compute_dx_log_g(x, y, precomputed_values)
            self.dy_log_g_of_t[index] = self.compute_dy_log_g(x, y, precomputed_values)
            self.curvature[index] = self.compute_curvature(x, y, precomputed_values)

    @abstractmethod
    def compute_log_g(self, x: int, y: int, precomputed_values: list = None) -> float:
        raise Exception("Undefined function for log_g")

    @abstractmethod
    def compute_dx_log_g(self, x, y, precomputed_values: list = None) -> float:
        return ddx(self.log_g, x, y)

    @abstractmethod
    def compute_dy_log_g(self, x: int, y: int, precomputed_values: list = None) -> float:
        return ddy(self.log_g, x, y)

    @abstractmethod
    def compute_curvature(self, x: int, y: int, precomputed_values: list = None) -> float:
        return -(dddx(self.log_g, x, y) + dddy(self.log_g, x, y)) / (2 * math.exp(self.log_g(x, y)))

    @abstractmethod
    def display(self, x_values: list, y_values: list):
        if len(x_values) != len(y_values):
            raise Exception("Set of x and y points are not of the same length.")
        else:
            pass
