import numba
from numba.experimental.jitclass.base import JitClassType
from numba.experimental import jitclass
from numba import njit, types, typeof, prange
from numpy import e, sin, cos, power, exp, square, floor_divide, zeros, ones, any, where, append

from .metric import Metric
from .domain import Domain


type_forward_euler = 0
type_backward_euler = 1
type_runge_kutta_4 = 2


@njit(fastmath=True, parallel=True, nogil=True)
def parallel_compute_geodesic(domain_is_inside, domain_min_radius, timeout, step_size, x_values, y_values, theta_values):
    min_radius_squared = square(domain_min_radius)
    num_max_steps = floor_divide(timeout, step_size)
    num_geodesics = x_values.size

    x = zeros(num_geodesics * num_max_steps)
    y = zeros(num_geodesics * num_max_steps)
    theta = zeros(num_geodesics * num_max_steps)

    for i in prange(num_geodesics):
        x[i] = x_values[i]
        y[i] = y_values[i]
        theta[i] = theta_values[i]

    step = 1

    # todo: optimize this
    inside_points = ones(num_geodesics)
    while any(inside_points) and step < num_max_steps:
        inside_points_indexes = zeros(0)
        outside_points_indexes = zeros(0)

        for i in prange(num_geodesics):
            if inside_points[i] == 1:
                append(inside_points_indexes, inside_points[i])
            else:
                append(outside_points_indexes, inside_points[i])



@njit(fastmath=True, parallel=True, nogil=True)
def parallel_compute_geo_step_forward_euler(x_values, y_values, theta_values, step_size, log_g_of_t, dx_log_g_of_t, dy_log_g_of_t):
    cos_theta = cos(theta_values)
    sin_theta = sin(theta_values)
    hh = exp(-0.5 * log_g_of_t) * step_size
    x = x_values + hh * cos_theta
    y = y_values + hh * sin_theta
    theta = theta_values + 0.5 * hh * (cos_theta * dy_log_g_of_t - sin_theta * dx_log_g_of_t)

    return x, y, theta


@njit(fastmath=True, parallel=True, nogil=True)
def parallel_compute_geo_step_backwards_euler(x_values, y_values, theta_values, step_size, log_g_of_t, dx_log_g_of_t, dy_log_g_of_t):
    cos_theta = cos(theta_values)
    sin_theta = sin(theta_values)

    precomputed_exp = exp(-0.5 * log_g_of_t)

    k1x = precomputed_exp * cos_theta
    k1y = precomputed_exp * sin_theta
    k1th = 0.5 * precomputed_exp * (cos_theta * dy_log_g_of_t - sin_theta * dx_log_g_of_t)

    x = x_values + step_size * k1x
    y = y_values + step_size * k1y
    theta = theta_values + step_size * k1th

    return x, y, theta


@njit(fastmath=True, parallel=True, nogil=True)
def parallel_compute_geo_step_runge_kutta(metric, x_values, y_values, theta_values, step_size, log_g_of_t, dx_log_g_of_t, dy_log_g_of_t):
    # first slope
    cos_theta = cos(theta_values)
    sin_theta = sin(theta_values)

    precomputed_exp = exp(-0.5 * log_g_of_t)

    k1x = precomputed_exp * cos_theta
    k1y = precomputed_exp * sin_theta
    k1th = 0.5 * precomputed_exp * (cos_theta * dy_log_g_of_t - sin_theta * dx_log_g_of_t)

    x = x_values + (step_size / 2) * k1x
    y = y_values + (step_size / 2) * k1y
    theta = theta_values + (step_size / 2) * k1th

    # second slope
    metric.compute_values(x, y)
    cos_theta = cos(theta)
    sin_theta = sin(theta)

    precomputed_exp = exp(-0.5 * log_g_of_t)

    k2x = precomputed_exp * cos_theta
    k2y = precomputed_exp * sin_theta
    k2th = 0.5 * precomputed_exp * (cos_theta * dy_log_g_of_t - sin_theta * dx_log_g_of_t)

    x = x_values + (step_size / 2) * k2x
    y = y_values + (step_size / 2) * k2y
    theta = theta_values + (step_size / 2) * k2th

    # third slope
    metric.compute_values(x, y)
    cos_theta = cos(theta)
    sin_theta = sin(theta)

    precomputed_exp = exp(-0.5 * log_g_of_t)

    k3x = precomputed_exp * cos_theta
    k3y = precomputed_exp * sin_theta
    k3th = 0.5 * precomputed_exp * (cos_theta * dy_log_g_of_t - sin_theta * dx_log_g_of_t)

    x = x_values + step_size * k3x
    y = y_values + step_size * k3y
    theta = theta_values + step_size * k3th

    # fourth slope
    metric.compute_values(x, y)
    cos_theta = cos(theta)
    sin_theta = sin(theta)

    precomputed_exp = exp(-0.5 * log_g_of_t)

    k4x = precomputed_exp * cos_theta
    k4y = precomputed_exp * sin_theta
    k4th = 0.5 * precomputed_exp * (cos_theta * dy_log_g_of_t - sin_theta * dx_log_g_of_t)

    x = x_values + (step_size / 6) * (k1x + k2x * 2 + k3x * 3 + k4x)
    y = y_values + (step_size / 6) * (k1y + k2y * 2 + k3y * 3 + k4y)
    theta = theta_values + (step_size / 6) * (k1th + k2th * 2 + k3th * 3 + k4th)

    return x, y, theta


members = [
    ('domain', Domain.class_type.instance_type),
    ('metric', Metric.class_type.instance_type),
    ('step_type', types.int8),
    ('step_size', types.double),
    ('timeout', types.double)
]


@jitclass(members)
class Surface:
    def __init__(self, domain, metric, step_type, step_size, timeout):
        self.domain = domain
        self.metric = metric
        self.step_type = step_type
        self.step_size = step_size
        self.timeout = timeout

    def geo_step(self, x_values, y_values, theta_values):
        self.metric.compute_values(x_values, y_values)
        if self.step_type == type_forward_euler:
            parallel_compute_geo_step_forward_euler(x_values, y_values, theta_values, self.step_size, self.metric.log_g_of_t, self.metric.dx_log_g_of_t, self.metric.dy_log_g_of_t)
        if self.step_type == type_backward_euler:
            prediction = parallel_compute_geo_step_backwards_euler(x_values, y_values, theta_values, self.step_size, self.metric.log_g_of_t, self.metric.dx_log_g_of_t, self.metric.dy_log_g_of_t)
            return parallel_compute_geo_step_backwards_euler(prediction[0], prediction[1], prediction[2], self.step_size, self.metric.log_g_of_t, self.metric.dx_log_g_of_t, self.metric.dy_log_g_of_t)
        if self.step_type == type_runge_kutta_4:
            return parallel_compute_geo_step_runge_kutta(self.metric, x_values, y_values, theta_values, self.step_size, self.metric.log_g_of_t, self.metric.dx_log_g_of_t, self.metric.dy_log_g_of_t)
