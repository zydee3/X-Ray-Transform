from numba import njit
from numpy import sin, cos, exp


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
