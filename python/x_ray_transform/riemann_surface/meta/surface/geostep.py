from numba import njit
from numpy import sin, cos, exp

from x_ray_transform.constants import type_step_forward_euler, type_step_backward_euler, type_step_runge_kutta_4


@njit(fastmath=True, parallel=True, nogil=True)
def parallel_compute_geo_step(step_type, metric, step_size, x_values, y_values, theta_values):
    """
    Steps a geodesic forward in accordance to the metric and stepping method

    :param step_type: of type x_ray_transform.constants.type_step_generic
    :param metric: of type x_ray_transform.riemann_surface.metric
    :param step_size:
    :param x_values:
    :param y_values:
    :param theta_values:
    :return:
    """

    # FORWARDS EULER ===================================================================================================
    if step_type == type_step_forward_euler:
        cos_theta = cos(theta_values)
        sin_theta = sin(theta_values)
        hh = exp(-0.5 * metric.log_g_of_t) * step_size
        x = x_values + hh * cos_theta
        y = y_values + hh * sin_theta
        theta = theta_values + 0.5 * hh * (cos_theta * metric.dy_log_g_of_t - sin_theta * metric.dx_log_g_of_t)

        return x, y, theta

    # BACKWARDS EULER ==================================================================================================
    elif step_type == type_step_backward_euler:
        cos_theta = cos(theta_values)
        sin_theta = sin(theta_values)

        precomputed_exp = exp(-0.5 * metric.log_g_of_t)

        k1_x = precomputed_exp * cos_theta
        k1_y = precomputed_exp * sin_theta
        k1_theta = 0.5 * precomputed_exp * (cos_theta * metric.dy_log_g_of_t - sin_theta * metric.dx_log_g_of_t)

        x = x_values + step_size * k1_x
        y = y_values + step_size * k1_y
        theta = theta_values + step_size * k1_theta

        return x, y, theta

    # RUNGE KUTTA 4 ====================================================================================================
    elif step_type == type_step_runge_kutta_4:
        # first slope --------------------------------------------------------------------------------------------------
        cos_theta = cos(theta_values)
        sin_theta = sin(theta_values)

        precomputed_exp = exp(-0.5 * metric.log_g_of_t)

        k1_x = precomputed_exp * cos_theta
        k1_y = precomputed_exp * sin_theta
        k1_theta = 0.5 * precomputed_exp * (cos_theta * metric.dy_log_g_of_t - sin_theta * metric.dx_log_g_of_t)

        x_values = x_values + (step_size / 2) * k1_x
        y_values = y_values + (step_size / 2) * k1_y
        theta_values = theta_values + (step_size / 2) * k1_theta

        # second slope -------------------------------------------------------------------------------------------------
        metric.compute_values(x_values, y_values)
        cos_theta = cos(theta_values)
        sin_theta = sin(theta_values)

        precomputed_exp = exp(-0.5 * metric.log_g_of_t)

        k2_x = precomputed_exp * cos_theta
        k2_y = precomputed_exp * sin_theta
        k2_theta = 0.5 * precomputed_exp * (cos_theta * metric.dy_log_g_of_t - sin_theta * metric.dx_log_g_of_t)

        x_values = x_values + (step_size / 2) * k2_x
        y_values = y_values + (step_size / 2) * k2_y
        theta_values = theta_values + (step_size / 2) * k2_theta

        # third slope --------------------------------------------------------------------------------------------------
        metric.compute_values(x_values, y_values)
        cos_theta = cos(theta_values)
        sin_theta = sin(theta_values)

        precomputed_exp = exp(-0.5 * metric.log_g_of_t)

        k3_x = precomputed_exp * cos_theta
        k3_y = precomputed_exp * sin_theta
        k3_theta = 0.5 * precomputed_exp * (cos_theta * metric.dy_log_g_of_t - sin_theta * metric.dx_log_g_of_t)

        x_values = x_values + step_size * k3_x
        y_values = y_values + step_size * k3_y
        theta_values = theta_values + step_size * k3_theta

        # fourth slope -------------------------------------------------------------------------------------------------
        metric.compute_values(x_values, y_values)
        cos_theta = cos(theta_values)
        sin_theta = sin(theta_values)

        precomputed_exp = exp(-0.5 * metric.log_g_of_t)

        k4_x = precomputed_exp * cos_theta
        k4_y = precomputed_exp * sin_theta
        k4_theta = 0.5 * precomputed_exp * (cos_theta * metric.dy_log_g_of_t - sin_theta * metric.dx_log_g_of_t)

        x_values = x_values + (step_size / 6) * (k1_x + k2_x * 2 + k3_x * 3 + k4_x)
        y_values = y_values + (step_size / 6) * (k1_y + k2_y * 2 + k3_y * 3 + k4_y)
        theta_values = theta_values + (step_size / 6) * (k1_theta + k2_theta * 2 + k3_theta * 3 + k4_theta)

        return x_values, y_values, theta_values
