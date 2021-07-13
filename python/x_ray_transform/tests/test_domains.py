from x_ray_transform.riemann_surface import Domain
from x_ray_transform.meta import time_computation
from numpy import array, random, power, copy, zeros
from numba import njit


type_circle = 0
type_cosine = 1
type_ellipse = 2


@njit(fastmath=True, parallel=True, nogil=True)
def test_circle(bound=10, num_elements=100000):
    domain = Domain(type_circle)
    domain.radius = random.randint(bound)
    theta = random.uniform(1, bound, num_elements)
    time_computation(domain, theta, False, "Circle")


@njit(fastmath=True, parallel=True, nogil=True)
def test_cosine(bound=10, num_elements=100000):
    domain = Domain(type_cosine)
    domain.radius = random.randint(bound)
    domain.amplitude = random.randint(bound)
    domain.cycles = random.randint(bound)
    domain.theta_offset = random.randint(bound)
    theta = random.uniform(1, bound, num_elements)
    time_computation(domain, theta, False, "Cosine")


@njit(fastmath=True, parallel=True, nogil=True)
def test_ellipse(bound=10, num_elements=100000):
    domain = Domain(type_ellipse)
    domain.theta_offset = random.randint(bound)
    domain.minor_radius = random.randint(bound)
    domain.major_radius = random.randint(bound)
    theta = random.uniform(1, bound, num_elements)
    time_computation(domain, theta, False, "Ellipse")


@njit(fastmath=True, nogil=True)
def test_domains():
    num_elements_override = 100000

    test_circle(num_elements=num_elements_override)
    test_cosine(num_elements=num_elements_override)
    test_ellipse(num_elements=num_elements_override)
