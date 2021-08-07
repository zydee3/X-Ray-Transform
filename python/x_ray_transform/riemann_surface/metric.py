from numba import njit, types
from numba.experimental import jitclass
from numpy import zeros, array, power, square, log, full, e

from x_ray_transform.constants import type_metric_constant_curvature, type_metric_gaussian, type_metric_polynomial

#merge euclidean, hyperbolic, sphere => constant_curvature


@njit(fastmath=True, parallel=True, nogil=True)
def parallel_compute_values_constant_curvature(radius_squared, x_values, y_values):
    difference_of_squares = radius_squared - square(x_values) - square(y_values)
    log_g_of_t = log(square(radius_squared) * 4) - log(square(difference_of_squares))
    dx_log_g_of_t = 4 * x_values / difference_of_squares
    dy_log_g_of_t = 4 * y_values / difference_of_squares
    curvature = full(x_values.size, (-1 / radius_squared))

    return log_g_of_t, dx_log_g_of_t, dy_log_g_of_t, curvature


@njit(fastmath=True, parallel=True, nogil=True)
def parallel_compute_values_gaussian(x_values, y_values, widths, weights, center_x, center_y):
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

    return log_g_of_t, dx_log_g_of_t, dy_log_g_of_t, curvature


@njit(fastmath=True, parallel=True, nogil=True)
def parallel_compute_values_polynomial(x_values, y_values, coefficients):
    coefficient_0_x = coefficients[0] * x_values

    log_g_of_t = square(coefficient_0_x) + (coefficients[1] * x_values * y_values) + square(coefficients[2] * y_values) \
                      + (coefficients[3] * x_values) + (coefficients[4] * y_values) + coefficients[5]
    dx_log_g_of_t = (2 * coefficients[0] * x_values) + (coefficients[1] * y_values) + coefficients[3]
    dy_log_g_of_t = (coefficients[1] * x_values) + (2 * coefficients[2] * y_values) + coefficients[4]
    curvature = -1 * power(e, log_g_of_t) * (coefficients[0] + coefficients[1])

    return log_g_of_t, dx_log_g_of_t, dy_log_g_of_t, curvature





members = [
    ('metric_type', types.int8),
    ('log_g_of_t', types.float64[:]),
    ('dx_log_g_of_t', types.float64[:]),
    ('dy_log_g_of_t', types.float64[:]),
    ('curvature', types.float64[:]),
    ('radius', types.float64),
    ('radius_squared', types.float64),
    ('widths', types.float64[:]),
    ('weights', types.float64[:]),
    ('center_x', types.float64[:]),
    ('center_y', types.float64[:]),
    ('coefficients', types.float64[:])
]


@jitclass(members)
class Metric:
    def __init__(self, metric_type):
        self.metric_type = metric_type
        self.radius = 0
        self.radius_squared = 0
        self.log_g_of_t = zeros(0)
        self.dx_log_g_of_t = zeros(0)
        self.dy_log_g_of_t = zeros(0)
        self.curvature = zeros(0)
        self.widths = zeros(0)
        self.weights = zeros(0)
        self.center_x = zeros(0)
        self.center_y = zeros(0)
        self.coefficients = zeros(0)

    def set_widths(self, widths):
        self.widths = widths

    def set_bumps(self, widths, weights, center_x, center_y):
        self.widths = widths
        self.weights = weights
        self.center_x = center_x
        self.center_y = center_y

    def set_radius(self, radius):
        self.radius = radius
        self.radius_squared = square(radius)

    def set_coefficients(self, coefficients):
        self.coefficients = coefficients

    def compute_values(self, x_values: array, y_values: array):
        result = (zeros(0), zeros(0), zeros(0), zeros(0))

        if self.metric_type == type_metric_constant_curvature:
            result = parallel_compute_values_hyperbolic(self.radius_squared, x_values, y_values)
        if self.metric_type == type_euclidean:
            result = parallel_compute_values_euclidean(x_values.size)
        if self.metric_type == type_polynomial:
            result = parallel_compute_values_polynomial(x_values, y_values, self.coefficients)

        self.log_g_of_t = result[0]
        self.dx_log_g_of_t = result[1]
        self.dy_log_g_of_t = result[2]
        self.curvature = result[3]
