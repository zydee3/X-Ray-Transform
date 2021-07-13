import numba
from numpy import zeros, array, empty, power, square, log, full, e
from numba.experimental import jitclass
from numba import njit, types


type_euclidean = 0
type_gaussian = 1
type_hyperbolic = 2
type_polynomial = 3
type_sphere = 4


@njit(fastmath=True, nogil=True)
def parallel_compute_values_euclidean(num_elements, log_g_of_t, dx_log_g_of_t, dy_log_g_of_t, curvature):
    log_g_of_t = zeros(num_elements)
    dx_log_g_of_t = zeros(num_elements)
    dy_log_g_of_t = zeros(num_elements)
    curvature = zeros(num_elements)


@njit(fastmath=True, parallel=True, nogil=True)
def parallel_compute_values_gaussian(x_values, y_values, log_g_of_t, dx_log_g_of_t, dy_log_g_of_t, curvature, widths, weights, center_x, center_y):
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


@njit(fastmath=True, parallel=True, nogil=True)
def parallel_compute_values_hyperbolic(self, x_values, y_values):
    difference_of_squares = self.radius_squared + square(x_values) + square(y_values)
    self.log_g_of_t = log(square(self.radius_squared) * 4) - log(square(difference_of_squares))
    self.dx_log_g_of_t = 4 * x_values / difference_of_squares
    self.dy_log_g_of_t = 4 * y_values / difference_of_squares
    self.curvature = full(x_values.size, (-1 / self.radius_squared))


@njit(fastmath=True, parallel=True, nogil=True)
def parallel_compute_values_polynomial(x_values, y_values, coefficients, log_g_of_t, dx_log_g_of_t, dy_log_g_of_t, curvature):
    coefficient_0_x = coefficients[0] * x_values

    log_g_of_t = square(coefficient_0_x) + (coefficients[1] * x_values * y_values) + square(coefficients[2] * y_values) \
                      + (coefficients[3] * x_values) + (coefficients[4] * y_values) + coefficients[5]
    dx_log_g_of_t = (2 * coefficients[0] * x_values) + (coefficients[1] * y_values) + coefficients[3]
    dy_log_g_of_t = (coefficients[1] * x_values) + (2 * coefficients[2] * y_values) + coefficients[4]
    curvature = -1 * power(e, log_g_of_t) * (coefficients[0] + coefficients[1])


@njit(fastmath=True, parallel=True, nogil=True)
def parallel_compute_values_sphere(self, x_values, y_values):
    sum_of_squares = self.radius_squared + square(x_values) + square(y_values)
    self.log_g_of_t = log(square(self.radius_squared) * 4) - log(square(sum_of_squares))
    self.dx_log_g_of_t = 4 * x_values / sum_of_squares
    self.dy_log_g_of_t = 4 * y_values / sum_of_squares
    self.curvature = full(x_values.size, (-1 / self.radius_squared))


members = [
    ('metric_type', types.int32),
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
        if self.metric_type == type_euclidean:
            parallel_compute_values_euclidean(x_values.size, self.log_g_of_t, self.dx_log_g_of_t, self.dy_log_g_of_t, self.curvature)
        if self.metric_type == type_gaussian:
            parallel_compute_values_gaussian(x_values, y_values, self.log_g_of_t, self.dx_log_g_of_t, self.dy_log_g_of_t, self.curvature, self.widths, self.weights, self.center_x, self.center_y)
        if self.metric_type == type_hyperbolic:
            parallel_compute_values_hyperbolic(self, x_values, y_values)
        if self.metric_type == type_polynomial:
            parallel_compute_values_polynomial(x_values, y_values, self.coefficients, self.log_g_of_t, self.dx_log_g_of_t, self.dy_log_g_of_t, self.curvature)
        if self.metric_type == type_sphere:
            parallel_compute_values_sphere(self, x_values, y_values)



