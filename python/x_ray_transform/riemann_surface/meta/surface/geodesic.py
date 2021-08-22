from numpy import square, floor_divide, zeros, ones, append
from numba import prange, njit


@njit(fastmath=True, parallel=True, nogil=True)
def compute_geodesic(domain_is_inside, domain_min_radius, timeout, step_size, x_values, y_values, theta_values):
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

    # todo this is terribly slow, revisit
