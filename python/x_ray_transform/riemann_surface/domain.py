from numba.experimental import jitclass
from numba import njit, types
from numpy import sin, cos, power, square, pi, arctan2, zeros, empty, array, sqrt, linspace, amin, amax, array, any, where

from x_ray_transform import type_domain_circle, type_domain_cosine, type_domain_ellipse

# region compute_values


@njit(fastmath=True, parallel=True, nogil=True)
def parallel_compute_alpha_normal(theta, boundary, dtheta_boundary):
    ct = cos(theta)
    st = sin(theta)
    return 0.5 * pi - arctan2(boundary * ct + dtheta_boundary * st, boundary * st - dtheta_boundary * ct)


@njit(parallel=True, nogil=True)
def parallel_compute_values_circle(radius, theta, compute_alpha_normal):
    boundary = empty(theta.size)
    boundary.fill(radius)
    if not compute_alpha_normal:
        return boundary, zeros(0)

    return boundary, parallel_compute_alpha_normal(theta, boundary, 0)


@njit(parallel=True, nogil=True)
def parallel_compute_values_cosine(cycles, radius, amplitude, theta, compute_alpha_normal):
    ct = cycles * theta
    cos_ct = cos(ct)
    boundary = radius + amplitude * cos_ct

    if not compute_alpha_normal:
        return boundary, zeros(0)

    dtheta_boundary = -1 * cycles * amplitude * sin(ct)
    # ddtheta_boundary = -1 * square(cycles) * amplitude * cos_ct

    return boundary, parallel_compute_alpha_normal(theta, boundary, dtheta_boundary)


@njit(fastmath=True, parallel=True)
def parallel_compute_values_ellipse(minor_radius, major_radius, theta_offset, theta, compute_alpha_normal):
    # save redundant computations
    mm = minor_radius * major_radius

    # computing boundary
    delta_theta = theta - theta_offset
    distance = sqrt(major_radius * square(cos(delta_theta)) + minor_radius * square(sin(delta_theta)))
    boundary = mm / distance

    if not compute_alpha_normal:
        return boundary, zeros(0)

    # computing dtheta and ddtheta boundary
    dtheta_boundary = power(boundary, 3) * sin(2 * delta_theta) * (square(major_radius) - square(minor_radius ** 2)) / 2 / mm * mm
    # ddtheta_boundary = square(3 * dtheta_boundary) / boundary + power(boundary, 3) * cos(2 * delta_theta) * (major_radius - minor_radius) / mm_squared

    return boundary, parallel_compute_alpha_normal(theta, boundary, dtheta_boundary)


# endregion


# region internals


@njit(fastmath=True, parallel=True, nogil=True)
def get_bounding_box(theta_offset, radius_samples):
    # radius_samples = compute_values(linspace(0, pi * 2, 1000) - theta_offset)
    euclidean_samples_x = cos(linspace(0, pi * 2, 1000)) * radius_samples
    euclidean_samples_y = sin(linspace(0, pi * 2, 1000)) * radius_samples

    min_x = amin(euclidean_samples_x)
    max_x = amax(euclidean_samples_x)
    min_y = amin(euclidean_samples_y)
    max_y = amax(euclidean_samples_y)

    return min_x, max_x, min_y, max_y


@njit(fastmath=True, parallel=True, nogil=True)
def is_inside(x_offset, y_offset, theta_offset, compute_boundary, x, y, min_radius_squared):
    x_values = x - x_offset
    y_values = y - y_offset
    xy_squared = square(x_values) + square(y_values)
    inside_points = xy_squared <= min_radius_squared

    if not any(inside_points):
        outside_index = where(not inside_points)
        radius = compute_boundary(arctan2(y_values(outside_index), x_values(outside_index)) - theta_offset)
        inside_points[outside_index] = xy_squared(outside_index) <= square(radius)

    return inside_points


# endregion


members = [
    ('domain_type', types.int8),
    ('radius', types.double),
    ('amplitude', types.double),
    ('cycles', types.double),
    ('minor_radius', types.double),
    ('major_radius', types.double),
    ('theta_offset', types.double),
    ('x_offset', types.double),
    ('y_offset', types.double),
    ('origin', types.double)
]


@jitclass(members)
class Domain:
    def __init__(self, domain_type):
        self.domain_type = domain_type
        self.radius = 0
        self.amplitude = 0
        self.cycles = 0
        self.minor_radius = 0
        self.major_radius = 0
        self.theta_offset = 0
        self.x_offset = 0
        self.y_offset = 0
        self.origin = 0

    def compute_values(self, theta, compute_alpha_normal=False):
        if self.domain_type == type_domain_circle:
            return parallel_compute_values_circle(self.radius, theta, compute_alpha_normal)
        if self.domain_type == type_domain_cosine:
            return parallel_compute_values_cosine(self.cycles, self.radius, self.amplitude, theta, compute_alpha_normal)
        if self.domain_type == type_domain_ellipse:
            return parallel_compute_values_ellipse(self.minor_radius, self.major_radius, self.theta_offset, theta, compute_alpha_normal)

    def get_bounding_box(self):
        return get_bounding_box(self.theta_offset)

    def is_inside(self, x_values, y_values, min_radius_squared):
        return is_inside(self.x_offset, self.y_offset, self.theta_offset, self.compute_boundary, x_values, y_values, min_radius_squared)
