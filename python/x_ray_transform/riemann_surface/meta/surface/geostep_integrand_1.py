from numba import njit
from numpy import sin, cos, exp

from x_ray_transform.constants import type_step_forward_euler, type_step_backward_euler, type_step_runge_kutta_4


@njit(fastmath=True, parallel=True, nogil=True)
def compute_geo_step_I1(compute_values, integrand_u, integrand_v, step_type, step_size, x_values, y_values, theta_values, u):
    """

    :param compute_values:
    :param integrand_u:
    :param integrand_v:
    :param step_type:
    :param step_size:
    :param x_values:
    :param y_values:
    :param theta_values:
    :param u:
    :return:
    """

    half_step_size = step_size / 2
    values = compute_values(x_values, y_values)
    cos_theta = cos(theta_values)
    sin_theta = sin(theta_values)

    # FORWARDS EULER ===================================================================================================
    if step_type == type_step_forward_euler:
        step_magnitude = exp(-0.5 * values[0]) * step_size
        x_values = x_values + step_magnitude * cos_theta
        y_values = y_values + step_magnitude * sin_theta
        theta_values = theta_values + 0.5 * step_magnitude * (cos_theta * values[2] - sin_theta * values[1])
        u = u + ((integrand_u(x_values, y_values) * cos_theta) + integrand_v(x_values, y_values) * sin_theta) * step_magnitude

        return x_values, y_values, theta_values, u

    # BACKWARDS EULER ==================================================================================================
    elif step_type == type_step_backward_euler:
        # predictor ----------------------------------------------------------------------------------------------------
        exponential_log_of_g_1 = exp(-0.5 * values[0])
        x = exponential_log_of_g_1 * cos_theta
        y = exponential_log_of_g_1 * sin_theta
        theta = 0.5 * exponential_log_of_g_1 * (cos_theta * values[2] - sin_theta * values[1])

        # corrector
        values = compute_values(x_values + step_size * x, y_values + step_size * y)
        new_theta = theta_values + step_size * theta
        cos_theta = cos(new_theta)
        sin_theta = sin(new_theta)

        exponential_log_of_g_2 = exp(-0.5 * values[0])

        u = u + (
                (integrand_u(x_values, y_values) * cos_theta + (integrand_v(x_values, y_values) * sin_theta)) * exponential_log_of_g_1
               + (integrand_u(x_values + step_size * x, y_values + step_size * y) * cos_theta + integrand_v(x_values + step_size * x, y_values + step_size * y) * sin_theta) * exponential_log_of_g_2
        )
        x = x_values + half_step_size * (x + exponential_log_of_g_2 * cos_theta)
        y = y_values + half_step_size * (y + exponential_log_of_g_2 * sin_theta)
        theta = theta_values + half_step_size * (theta + 0.5 * exponential_log_of_g_2 * (cos_theta * values[2] - sin_theta * values[1]))

        return x, y, theta, u

    # RUNGE KUTTA 4 ====================================================================================================
    elif step_type == type_step_runge_kutta_4:
        # first slope --------------------------------------------------------------------------------------------------
        exponential_log_of_g_1 = exp(-0.5 * values[0])
        k1_x = exponential_log_of_g_1 * cos_theta
        k1_y = exponential_log_of_g_1 * sin_theta
        k1_theta = 0.5 * exponential_log_of_g_1 * (cos_theta * values[2] - sin_theta * values[1])

        # second slope -------------------------------------------------------------------------------------------------
        values = compute_values(x_values + half_step_size * k1_x, y_values + half_step_size * k1_y)
        new_theta = theta_values + half_step_size * k1_theta
        cos_theta = cos(new_theta)
        sin_theta = sin(new_theta)

        exponential_log_of_g_2 = exp(-0.5 * values[0])
        k2_x = exponential_log_of_g_2 * cos_theta
        k2_y = exponential_log_of_g_2 * sin_theta
        k2_theta = 0.5 * exponential_log_of_g_2 * (cos_theta * values[2] - sin_theta * values[1])

        # third slope --------------------------------------------------------------------------------------------------
        values = compute_values(x_values + half_step_size * k2_x, y_values + half_step_size * k2_y)
        new_theta = theta_values + half_step_size * k2_theta
        cos_theta = cos(new_theta)
        sin_theta = sin(new_theta)

        exponential_log_of_g_3 = exp(-0.5 * values[0])
        x = exponential_log_of_g_3 * cos_theta
        y = exponential_log_of_g_3 * sin_theta
        theta = 0.5 * exponential_log_of_g_3 * (cos_theta * values[2] - sin_theta * values[1])

        # fourth slope -------------------------------------------------------------------------------------------------
        values = compute_values(x_values + step_size * x, y_values + step_size * y)
        new_theta = theta_values + step_size * theta
        cos_theta = cos(new_theta)
        sin_theta = sin(new_theta)

        exponential_log_of_g_4 = exp(-0.5 * values[0])

        u = u + ((integrand_u(x_values, y_values) * cos_theta + integrand_v(x_values, y_values) * sin_theta) * exponential_log_of_g_1
                 + 2 * (integrand_u(x_values + half_step_size * k1_x, y_values + half_step_size * k1_y) * cos_theta + integrand_v(x_values + half_step_size * k1_x, y_values + half_step_size * k1_y) * sin_theta) * exponential_log_of_g_2
                 + 2 * (integrand_u(x_values + half_step_size * k2_x, y_values + half_step_size * k2_y) * cos_theta + integrand_v(x_values + half_step_size * k2_x, y_values + half_step_size * k2_y) * sin_theta) * exponential_log_of_g_3
                 + (integrand_u(x_values + step_size * x, y_values + step_size * y) * cos_theta + integrand_v(x_values + step_size * x, y_values + step_size * y) * sin_theta) * exponential_log_of_g_4)

        partial_step_size = half_step_size / 3
        x = x_values + partial_step_size * (k1_x + 2 * (k2_x + x) + exponential_log_of_g_4 * cos_theta)
        y = y_values + partial_step_size * (k1_y + 2 * (k2_y + y) + exponential_log_of_g_4 * sin_theta)
        theta = theta_values + partial_step_size * k1_theta + 2 * (k2_theta + theta) + (0.5 * exponential_log_of_g_4 * (cos_theta * values[2] - sin_theta * values[1]))

        return x, y, theta, u
