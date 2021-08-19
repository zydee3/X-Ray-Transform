from numba import njit
from numpy import sin, cos, exp


@njit(fastmath=True, parallel=True, nogil=True)
def compute_jacobi_0(metric_compute_values, step_size, step_type, x_values, y_values, theta_value, b_value, b_dot_value):
    half_step_size = step_size / 2
    values = metric_compute_values(x_values, y_values)
    cos_theta = cos(theta_value)
    sin_theta = sin(theta_value)

    # FORWARDS EULER ===================================================================================================
    if step_type == type_step_forward_euler:
        exponential_log_of_g = exp(-0.5 * values[0]) * step_size
        x = x_values + exponential_log_of_g * cos_theta
        y = y_values + exponential_log_of_g * sin_theta
        theta = theta_value + 0.5 * exponential_log_of_g * (cos_theta * values[2] - sin_theta * values[1])
        b = b_value + step_size * b_dot_value
        b_dot = b_dot_value - step_size * values[3] * b_value

        return x, y, theta, b, b_dot

    # BACKWARDS EULER ==================================================================================================
    elif step_type == type_step_backward_euler:
        # predictor ----------------------------------------------------------------------------------------------------
        exponential_log_of_g = exp(-0.5 * values[0])
        x = exponential_log_of_g * cos_theta
        y = exponential_log_of_g * sin_theta
        theta = 0.5 * exponential_log_of_g * (cos_theta * values[2] - sin_theta * values[1])
        b_dot_1 = -values[3] * b_value

        # corrector ----------------------------------------------------------------------------------------------------
        values = metric_compute_values(x_values + step_size * x, y_values + step_size * y)
        cos_theta = cos(theta_value + step_size * theta)
        sin_theta = sin(theta_value + step_size * theta)

        new_value = -0.5 * values[0]
        x = x_values + half_step_size * (x + exp(new_value) * cos_theta)
        y = y_values + half_step_size * (y + exp(new_value) * sin_theta)
        theta = theta_value + half_step_size * (theta + 0.5 * exp(-0.5 * values[0]) * (cos_theta * values[2] - sin_theta * values[1]))
        b = b_value + half_step_size * (b_dot_1 + values[3] * (step_size * b_dot_1 - b_value))

        return x, y, theta, b, b_dot

    # RUNGE KUTTA 4 ====================================================================================================
    elif step_type == type_step_runge_kutta_4:
        # first slope --------------------------------------------------------------------------------------------------
        exponential_log_of_g = exp(-0.5 * values[0])
        k1_x = exponential_log_of_g * cos_theta
        k1_y = exponential_log_of_g * sin_theta
        k1_theta = 0.5 * exponential_log_of_g * (cos_theta * values[2] - sin_theta * values[1])
        k1_b_dot = -values[3] * b_value

        # second slope -------------------------------------------------------------------------------------------------

        # third slope --------------------------------------------------------------------------------------------------

        # fourth slope -------------------------------------------------------------------------------------------------

