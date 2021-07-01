from numba.experimental import jitclass
from numba import njit, types
from numpy import zeros


members = [
    ('log_g_of_t', types.float64[:]),
    ('dx_log_g_of_t', types.float64[:]),
    ('dy_log_g_of_t', types.float64[:]),
    ('curvature', types.float64[:])
]


@njit(fastmath=True, nogil=True)
def parallel_compute_values(num_elements, log_g_of_t, dx_log_g_of_t, dy_log_g_of_t, curvature):
    log_g_of_t = zeros(num_elements)
    dx_log_g_of_t = zeros(num_elements)
    dy_log_g_of_t = zeros(num_elements)
    curvature = zeros(num_elements)


@jitclass(members)
class Euclidean:
    def __init__(self):
        self.log_g_of_t = zeros(0)
        self.dx_log_g_of_t = zeros(0)
        self.dy_log_g_of_t = zeros(0)
        self.curvature = zeros(0)

    def compute_values(self, x_values: list, y_values: list):
        parallel_compute_values(x_values.size, self.log_g_of_t, self.dx_log_g_of_t, self.dy_log_g_of_t, self.curvature)
