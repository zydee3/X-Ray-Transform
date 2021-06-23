from .base_metric import BaseMetric


class Euclidean(BaseMetric):
    def __init__(self):
        super().__init__()

    def compute_values(self, x_values: list, y_values: list):
        num_elements = len(x_values)
        self.log_g_of_t += ([0] * num_elements)
        self.dx_log_g_of_t += ([0] * num_elements)
        self.dy_log_g_of_t += ([0] * num_elements)
        self.curvature += ([0] * num_elements)

    def compute_log_g(self, x: int, y: int) -> float:
        return 0

    def compute_dx_log_g(self, x: int, y: int) -> float:
        return 0

    def compute_dy_log_g(self, x: int, y: int) -> float:
        return 0

    def compute_curvature(self, x: int, y: int) -> float:
        return 0

    def display(self, x_values: list, y_values: list):
        super().display()
