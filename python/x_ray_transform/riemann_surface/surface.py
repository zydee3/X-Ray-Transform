from numba import types
from numba.experimental import jitclass

from x_ray_transform.constants import type_step_forward_euler, type_step_backward_euler, type_step_runge_kutta_4
from .domain import Domain
from .meta.geo_stepper import parallel_compute_geo_step_backwards_euler
from .meta.geo_stepper import parallel_compute_geo_step_forward_euler
from .meta.geo_stepper import parallel_compute_geo_step_runge_kutta
from .metric import Metric


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
        if self.step_type == type_step_forward_euler:
            return parallel_compute_geo_step_forward_euler(x_values, y_values, theta_values, self.step_size, self.metric.log_g_of_t, self.metric.dx_log_g_of_t, self.metric.dy_log_g_of_t)
        if self.step_type == type_step_backward_euler:
            prediction = parallel_compute_geo_step_backwards_euler(x_values, y_values, theta_values, self.step_size, self.metric.log_g_of_t, self.metric.dx_log_g_of_t, self.metric.dy_log_g_of_t)
            return parallel_compute_geo_step_backwards_euler(prediction[0], prediction[1], prediction[2], self.step_size, self.metric.log_g_of_t, self.metric.dx_log_g_of_t, self.metric.dy_log_g_of_t)
        if self.step_type == type_step_runge_kutta_4:
            return parallel_compute_geo_step_runge_kutta(self.metric, x_values, y_values, theta_values, self.step_size, self.metric.log_g_of_t, self.metric.dx_log_g_of_t, self.metric.dy_log_g_of_t)
