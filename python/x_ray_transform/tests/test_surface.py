from x_ray_transform.riemann_surface import Domain, Metric, Surface
from x_ray_transform.meta import time_geo_step
from numba import njit
from numpy import random

type_euclidean = 0
type_gaussian = 1
type_hyperbolic = 2
type_polynomial = 3
type_sphere = 4

type_circle = 0
type_cosine = 1
type_ellipse = 2


@njit(fastmath=True, parallel=True, nogil=True)
def test_surface():
    bound = 10
    num_elements = 65536

    domain = Domain(type_circle)
    domain.radius = 2

    metric = Metric(type_hyperbolic)
    metric.radius = 4
    metric.radius_squared = 16

    surface = Surface(domain, metric, 0, 0.01, 20)

    x_values = random.uniform(1, bound, num_elements)
    y_values = random.uniform(1, bound, num_elements)
    theta_values = random.uniform(1, bound, num_elements)

    time_geo_step(surface, x_values, y_values, theta_values)
