from ..metric import Metric
import math


class Gaussian(Metric):
    def __init__(self, widths: list, weights: list, center_x: list, center_y: list):
        """ sig = width, k = weights """
        super().__init__()

        len_widths = len(widths)
        len_weights = len(weights)
        len_center_x = len(center_x)
        len_center_y = len(center_y)

        if len_widths != len_weights or len_weights != len_center_x or len_center_x != len_center_y:
            raise Exception(f'Arrays widths ({len_widths}), weights ({len_weights}), len_center_x ({len_center_x}), and len_center_y ({len_center_y}) are not of the same length.')

        self.widths = widths
        self.weights = weights
        self.center_x = center_x
        self.center_y = center_y
        self.num_elements = len(widths)

    def compute_values(self, x_values: list, y_values: list):
        super().compute_values(x_values, y_values)

    def compute_log_g(self, x: int, y: int) -> float:
        previous_log_g = self.log_g_of_t[self.index - 1] if self.index > 0 else 0
        return previous_log_g + self.fj

    def compute_dx_log_g(self, x, y) -> float:
        previous_dx = self.dx_log_g_of_t[self.index - 1] if self.index > 0 else 0
        return previous_dx - self.x_bar / self.widths[self.index] * self.fj

    def compute_dy_log_g(self, x, y) -> float:
        previous_dy = self.dy_log_g_of_t[self.index - 1] if self.index > 0 else 0
        return previous_dy - self.y_bar / self.widths[self.index] * self.fj

    def compute_curvature(self, x, y) -> float:
        previous_curvature = self.curvature[self.index - 1] if self.index > 0 else 0
        return previous_curvature + (self.x_bar ** 2 + self.y_bar ** 2 - 2) / self.widths[self.index] * self.fj

    def compute_x_bar(self, x, width_index):
        return (x - self.center_x[width_index]) / self.widths[width_index]

    def compute_y_bar(self, y, width_index):
        return (y - self.center_y[width_index]) / self.widths[width_index]

    def compute_fj(self, x, y, width_index, x_bar, y_bar):
        return self.weights[width_index] * math.e ** (-0.5 * (x_bar ** 2 + y_bar ** 2))

    def display(self, x_values, y_values):
        super().display()
