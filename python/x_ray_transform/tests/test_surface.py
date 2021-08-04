from x_ray_transform.riemann_surface import Domain, Metric, Surface
from x_ray_transform.meta import time_geo_step
from numba import njit
from numpy import random

from x_ray_transform import type_domain_circle, type_domain_cosine, type_domain_ellipse, type_metric_gaussian, type_metric_polynomial, type_metric_constant_curvature


@njit(fastmath=True, parallel=True)
def test_surface():
    bound = 10
    num_elements = 65536

    domain = Domain(type_domain_circle)
    domain.radius = 2

    metric = Metric(type_metric_constant_curvature)
    metric.radius = 4
    metric.radius_squared = 16

    surface = Surface(domain, metric, 0, 0.01, 20)

    x_values = random.uniform(1, bound, num_elements)
    y_values = random.uniform(1, bound, num_elements)
    theta_values = random.uniform(1, bound, num_elements)

    time_geo_step(surface, x_values, y_values, theta_values, "circle", "hyperbolic")
