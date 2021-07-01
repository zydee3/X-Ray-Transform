from numba import types
from numpy import zeros, square, e, power
from numba.experimental import jitclass
from numba import njit


members = [
    ('coefficients', types.float64[:]),
    ('log_g_of_t', types.float64[:]),
    ('dx_log_g_of_t', types.float64[:]),
    ('dy_log_g_of_t', types.float64[:]),
    ('curvature', types.float64[:])
]


@njit(fastmath=True, parallel=True, nogil=True)
def parallel_compute_values(x_values, y_values, coefficients, log_g_of_t, dx_log_g_of_t, dy_log_g_of_t, curvature):
    coefficient_0_x = coefficients[0] * x_values

    log_g_of_t = square(coefficient_0_x) + (coefficients[1] * x_values * y_values) + square(coefficients[2] * y_values) \
                      + (coefficients[3] * x_values) + (coefficients[4] * y_values) + coefficients[5]
    dx_log_g_of_t = (2 * coefficients[0] * x_values) + (coefficients[1] * y_values) + coefficients[3]
    dy_log_g_of_t = (coefficients[1] * x_values) + (2 * coefficients[2] * y_values) + coefficients[4]
    curvature = -1 * power(e, log_g_of_t) * (coefficients[0] + coefficients[1])


@jitclass(members)
class Polynomial:
    def __init__(self, coefficients):
        self.coefficients = coefficients
        self.log_g_of_t = zeros(0)
        self.dx_log_g_of_t = zeros(0)
        self.dy_log_g_of_t = zeros(0)
        self.curvature = zeros(0)

        if coefficients.size != 6:
            raise Exception("Number of coefficients is not 6.")

    def compute_values(self, x_values: types.float64, y_values: types.float64):
        parallel_compute_values(x_values, y_values, self.coefficients, self.log_g_of_t, self.dx_log_g_of_t, self.dy_log_g_of_t, self.curvature)

