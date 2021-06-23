from abc import ABC, abstractmethod
from x_ray_transformation.meta.derivative import ddx, dddx, ddy, dddy
import math


class BaseMetric(ABC):

    def __init__(self):
        self.log_g_of_t = []
        self.dx_log_g_of_t = []
        self.dy_log_g_of_t = []
        self.curvature = []

        self.index = None
        self.x_bar = None
        self.y_bar = None
        self.fj = None

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

        for index in range(len(x_values)):
            x = x_values[index]
            y = y_values[index]

            if hasattr(self, "widths"):
                # pre-compute these values so they are only computed once per x,y
                self.index = index
                self.x_bar = self.compute_x_bar(x, index)
                self.y_bar = self.compute_y_bar(y, index)
                self.fj = self.compute_fj(index, self.x_bar, self.y_bar)

            self.log_g_of_t.append(self.compute_log_g(x, y))
            self.dx_log_g_of_t.append(self.compute_dx_log_g(x, y))
            self.dy_log_g_of_t.append(self.compute_dy_log_g(x, y))
            self.curvature.append(self.compute_curvature(x, y))

        # no longer need these, clear them for garbage collector
        self.index = None
        self.x_bar = None
        self.y_bar = None
        self.fj = None

    @abstractmethod
    def compute_log_g(self, x: int, y: int) -> float:
        raise Exception("Undefined function for log_g")

    @abstractmethod
    def compute_dx_log_g(self, x, y) -> float:
        return ddx(self.log_g, x, y)

    @abstractmethod
    def compute_dy_log_g(self, x: int, y: int) -> float:
        return ddy(self.log_g, x, y)

    @abstractmethod
    def compute_curvature(self, x: int, y: int) -> float:
        return -(dddx(self.log_g, x, y) + dddy(self.log_g, x, y)) / (2 * math.exp(self.log_g(x, y)))

    @abstractmethod
    def display(self, x_values: list, y_values: list):
        if len(x_values) != len(y_values):
            raise Exception("Set of x and y points are not of the same length.")
        else:
            pass
