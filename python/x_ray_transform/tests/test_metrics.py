import time

from x_ray_transform import Euclidean, Gaussian, Hyperbolic, Polynomial, Sphere, print_time

from timeit import default_timer as timer
from numpy import array, random, power, copy, zeros
from numba import njit, objmode


@njit(fastmath=True, nogil=True)
def time_computation(metric, x_values, y_values):
    with objmode(start='f8'):
        start = time.perf_counter()

    metric.compute_values(x_values, y_values)

    with objmode():
        end = time.perf_counter()
        num_spaces = 3
        length_name = len(metric.__class__.__name__)
        num_tabs = num_spaces - (length_name - (length_name % 4))/4
        tabs = ""
        for i in range(int(num_tabs)):
            tabs += '\t'

        print('Time Elapsed ({}):{} {}'.format(metric.__class__.__name__, tabs, end - start))


@njit(fastmath=True, parallel=True, nogil=True)
def test_euclidean(num_elements=100000):
    metric = Euclidean()
    time_computation(metric, zeros(num_elements), zeros(num_elements))


@njit(fastmath=True, parallel=True, nogil=True)
def test_gaussian(bound=2, num_elements=100000):
    x_values = random.uniform(1, bound, num_elements)
    y_values = random.uniform(1, bound, num_elements)
    widths = random.uniform(1, bound, num_elements)
    weights = random.uniform(1, bound, num_elements)
    center_x = random.uniform(1, bound, num_elements)
    center_y = random.uniform(1, bound, num_elements)
    metric = Gaussian(widths, weights, center_x, center_y)
    time_computation(metric, x_values, y_values)


@njit(fastmath=True, parallel=True, nogil=True)
def test_hyperbolic(radius=2, bound=1000, num_elements=100000):
    x_values = random.uniform(1, bound, num_elements)
    y_values = random.uniform(1, bound, num_elements)
    metric = Hyperbolic(radius)
    time_computation(metric, x_values, y_values)


@njit(fastmath=True, parallel=True, nogil=True)
def test_sphere(radius=2, bound=1000, num_elements=100000):
    x_values = random.uniform(1, bound, num_elements)
    y_values = random.uniform(1, bound, num_elements)
    metric = Sphere(radius)
    time_computation(metric, x_values, y_values)


@njit(fastmath=True, parallel=True, nogil=True)
def test_polynomial(coefficients=array([0.0, 0.0, 0.0, 0.0, 0.0, 1.0]), bound=1000, num_elements=100000):
    x_values = random.uniform(1, bound, num_elements)
    y_values = random.uniform(1, bound, num_elements)
    metric = Polynomial(coefficients)
    time_computation(metric, x_values, y_values)


def test_metrics():
    num_elements_override = 100000                          # 100,000 elements

    test_gaussian(num_elements=5)                           # 0.0005448 seconds
    test_euclidean(num_elements=num_elements_override)      # 0.0009681 seconds
    test_hyperbolic(num_elements=num_elements_override)     # 0.0017244 seconds
    test_sphere(num_elements=num_elements_override)         # 0.0008538 seconds
    test_polynomial(num_elements=num_elements_override)     # 0.0004439 seconds
