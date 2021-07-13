from x_ray_transform.riemann_surface import Metric
from x_ray_transform.meta import time_computation

from numpy import array, random, power, copy, zeros
from numba import njit, objmode

type_euclidean = 0
type_gaussian = 1
type_hyperbolic = 2
type_polynomial = 3
type_sphere = 4


@njit(fastmath=True, parallel=True, nogil=True)
def test_euclidean(num_elements=100000):
    metric = Metric(type_euclidean)
    time_computation(entity_name="Euclidean", entity=metric, arg1=zeros(num_elements), arg2=zeros(num_elements))


@njit(fastmath=True, parallel=True, nogil=True)
def test_gaussian(bound=2, num_elements=100000):
    x_values = random.uniform(1, bound, num_elements)
    y_values = random.uniform(1, bound, num_elements)
    widths = random.uniform(1, bound, num_elements)
    weights = random.uniform(1, bound, num_elements)
    center_x = random.uniform(1, bound, num_elements)
    center_y = random.uniform(1, bound, num_elements)
    metric = Metric(type_gaussian)
    metric.set_bumps(widths, weights, center_x, center_y)
    time_computation(entity_name="Gaussian", entity=metric, arg1=x_values, arg2=y_values)


@njit(fastmath=True, parallel=True, nogil=True)
def test_hyperbolic(radius=2, bound=1000, num_elements=100000):
    x_values = random.uniform(1, bound, num_elements)
    y_values = random.uniform(1, bound, num_elements)
    metric = Metric(type_hyperbolic)
    metric.set_radius(radius)
    time_computation(entity_name="Hyperbolic", entity=metric, arg1=x_values, arg2=y_values)


@njit(fastmath=True, parallel=True, nogil=True)
def test_sphere(radius=2, bound=1000, num_elements=100000):
    x_values = random.uniform(1, bound, num_elements)
    y_values = random.uniform(1, bound, num_elements)
    metric = Metric(type_sphere)
    metric.set_radius(radius)
    time_computation(entity_name="Sphere", entity=metric, arg1=x_values, arg2=y_values)


@njit(fastmath=True, parallel=True, nogil=True)
def test_polynomial(coefficients=array([0.0, 0.0, 0.0, 0.0, 0.0, 1.0]), bound=1000, num_elements=100000):
    x_values = random.uniform(1, bound, num_elements)
    y_values = random.uniform(1, bound, num_elements)

    metric = Metric(type_polynomial)
    metric.set_coefficients(coefficients)

    time_computation(metric, x_values, y_values, "Polynomial")


@njit(fastmath=True, nogil=True)
def test_metrics():
    num_elements_override = 100000                          # 100,000 elements

    test_gaussian(num_elements=5)                           # 0.0005448 seconds
    test_euclidean(num_elements=num_elements_override)      # 0.0009681 seconds
    test_hyperbolic(num_elements=num_elements_override)     # 0.0017244 seconds
    test_sphere(num_elements=num_elements_override)         # 0.0008538 seconds
    test_polynomial(num_elements=num_elements_override)     # 0.0004439 seconds
