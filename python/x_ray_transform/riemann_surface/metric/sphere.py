from numba import types
from numpy import zeros, float_power, square, log, full
from numba.experimental import jitclass
from numba import njit


members = [
    ('radius', types.float64),
    ('radius_squared', types.float64),
    ('log_g_of_t', types.float64[:]),
    ('dx_log_g_of_t', types.float64[:]),
    ('dy_log_g_of_t', types.float64[:]),
    ('curvature', types.float64[:])
]


@njit(fastmath=True, parallel=True, nogil=True)
def parallel_compute_values(x_values, y_values, radius, radius_squared, log_g_of_t, dx_log_g_of_t, dy_log_g_of_t, curvature):
    sum_of_squares = radius_squared + square(x_values) + square(y_values)
    log_g_of_t = log(square(radius_squared) * 4) - log(square(sum_of_squares))
    dx_log_g_of_t = 4 * x_values / sum_of_squares
    dy_log_g_of_t = 4 * y_values / sum_of_squares
    curvature = full(x_values.size, (-1 / radius_squared))


@jitclass(members)
class Sphere:
    def __init__(self, radius):
        self.radius = radius
        self.radius_squared = square(radius)
        self.log_g_of_t = zeros(0)
        self.dx_log_g_of_t = zeros(0)
        self.dy_log_g_of_t = zeros(0)
        self.curvature = zeros(0)

    def compute_values(self, x_values: types.float64, y_values: types.float64):
        parallel_compute_values(x_values, y_values, self.radius, self.radius_squared, self.log_g_of_t, self.dx_log_g_of_t, self.dy_log_g_of_t, self.curvature)
