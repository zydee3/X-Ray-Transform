from numba.experimental import jitclass
from numba import njit, types
from numpy import e, power, square
from numpy import zeros


members = [
    ('log_g_of_t', types.float64[:]),
    ('dx_log_g_of_t', types.float64[:]),
    ('dy_log_g_of_t', types.float64[:]),
    ('curvature', types.float64[:]),
    ('widths', types.float64[:]),
    ('weights', types.float64[:]),
    ('center_x', types.float64[:]),
    ('center_y', types.float64[:])
]


@njit(fastmath=True, parallel=True, nogil=True)
def parallel_compute_values(x_values, y_values, log_g_of_t, dx_log_g_of_t, dy_log_g_of_t, curvature, widths, weights, center_x, center_y):
    log_g_of_t = zeros(x_values.size)
    dx_log_g_of_t = zeros(x_values.size)
    dy_log_g_of_t = zeros(x_values.size)
    curvature = zeros(x_values.size)

    x_bar = x_values - center_x / widths
    y_bar = y_values - center_y / widths
    fj = weights * power(e, -0.5 * square(x_bar) + square(y_bar))

    for index in range(x_values.size):
        previous_log = 0 if index == 0 else log_g_of_t[index - 1]
        previous_dx = 0 if index == 0 else dx_log_g_of_t[index - 1]
        previous_dy = 0 if index == 0 else dy_log_g_of_t[index - 1]
        previous_curvature = 0 if index == 0 else curvature[index - 1]

        precompute = widths[index] * fj[index]
        log_g_of_t[index] = previous_log + fj[index]
        dx_log_g_of_t[index] = previous_dx - x_bar[index] / precompute
        dy_log_g_of_t[index] = previous_dy - y_bar[index] / precompute
        curvature[index] = previous_curvature + (square(x_bar[index] + y_bar[index] - 2))


@jitclass(members)
class Gaussian:
    def __init__(self, widths: types.float64, weights: types.float64, center_x: types.float64, center_y: types.float64):
        len_widths = widths.size
        len_weights = weights.size
        len_center_x = center_x.size
        len_center_y = center_y.size

        if len_widths != len_weights or len_weights != len_center_x or len_center_x != len_center_y:
            raise Exception('len_widths != len_weights or len_weights != len_center_x or len_center_x != len_center_y')

        self.log_g_of_t = zeros(0)
        self.dx_log_g_of_t = zeros(0)
        self.dy_log_g_of_t = zeros(0)
        self.curvature = zeros(0)

        self.widths = widths
        self.weights = weights
        self.center_x = center_x
        self.center_y = center_y

    def compute_values(self, x_values: types.float64, y_values: types.float64):
        parallel_compute_values(x_values, y_values, self.log_g_of_t, self.dx_log_g_of_t, self.dy_log_g_of_t, self.curvature, self.widths, self.weights, self.center_x, self.center_y)
