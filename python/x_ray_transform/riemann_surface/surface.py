from numba import types, njit
from numba.experimental import jitclass
from numpy import cos, sin, square, floor_divide, zeros, NAN
from x_ray_transform.constants import type_step_forward_euler, type_step_backward_euler, type_step_runge_kutta_4
from .domain import Domain
from x_ray_transform.riemann_surface.meta.surface.geostep import parallel_compute_geo_step_backwards_euler
from x_ray_transform.riemann_surface.meta.surface.geostep import parallel_compute_geo_step_forward_euler
from x_ray_transform.riemann_surface.meta.surface.geostep import parallel_compute_geo_step_runge_kutta
from .metric import Metric


"""
    RiemannSurface
    -Something to plot geodesics (plotGeo)
    -XrayI0,I1
    -I0star (inverts I0)
    -geodesic (to be used in plotGeo)
    -geodesicJacobiAB (and maybe some stuff to find conjugate points idk monard didnt ask abt it)
    (geodesic and geodesicJacobiAB are doing a weird parameterization thing dont worry abt that)
    -StoBeta/BetatoS (StoBeta converts between arclength from 0-Beta to Beta)
    -geodesicFoot or End (projects geodesics onto the boundary as a Beta,Alpha, used in I0star)
"""


@njit(fastmath=True, parallel=True, nogil=True)
def to_x_y_theta(domain_func_boundary, domain_x, domain_y, domain_theta, beta, theta):
    radius = domain_func_boundary(beta - domain_theta)
    x_values = cos(beta) * radius + domain_x
    y_values = sin(beta) * radius + domain_y

    return radius, x_values, y_values


@njit(fastmath=True, nogil=True)
def to_alpha_beta(min_radius, surface_timeout, surface_step_size, domain_x, domain_y, x_values, y_values, theta):
    """
    test

    :param min_radius: dsdsdd
    :param surface_timeout:
    :param surface_step_size:
    :param domain_x:
    :param domain_y:
    :param x_values:
    :param y_values:
    :param theta:
    :return:
    """

    radius = square(min_radius)
    num_max_steps = floor_divide(surface_timeout, surface_step_size)

    num_values = x_values.size
    beta = NAN(num_values)
    x = zeros(num_values)
    y = zeros(num_values)
    theta = zeros(num_values)


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
        to_alpha_beta
        self.metric.compute_values(x_values, y_values)
        if self.step_type == type_step_forward_euler:
            return parallel_compute_geo_step_forward_euler(x_values, y_values, theta_values, self.step_size, self.metric.log_g_of_t, self.metric.dx_log_g_of_t, self.metric.dy_log_g_of_t)
        if self.step_type == type_step_backward_euler:
            prediction = parallel_compute_geo_step_backwards_euler(x_values, y_values, theta_values, self.step_size, self.metric.log_g_of_t, self.metric.dx_log_g_of_t, self.metric.dy_log_g_of_t)
            return parallel_compute_geo_step_backwards_euler(prediction[0], prediction[1], prediction[2], self.step_size, self.metric.log_g_of_t, self.metric.dx_log_g_of_t, self.metric.dy_log_g_of_t)
        if self.step_type == type_step_runge_kutta_4:
            return parallel_compute_geo_step_runge_kutta(self.metric, x_values, y_values, theta_values, self.step_size, self.metric.log_g_of_t, self.metric.dx_log_g_of_t, self.metric.dy_log_g_of_t)
