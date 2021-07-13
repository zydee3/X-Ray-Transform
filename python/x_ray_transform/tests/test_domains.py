from x_ray_transform import Domain
from x_ray_transform.meta import time_computation
from numpy import array, random, power, copy, zeros
from numba import njit


type_circle = 0
type_cosine = 1
type_ellipse = 2


@njit(fastmath=True, parallel=True, nogil=True)
def test_ellipse(bound=10, num_elements=100000):
    theta_offset = random.randint(bound)
    minor_radius = random.randint(bound)
    major_radius = random.randint(bound)
    theta = random.uniform(1, bound, num_elements)
    domain = Domain(type_ellipse)
    time_computation(domain, theta, False, "Ellipse")


@njit(fastmath=True, nogil=True)
def test_domains():
    num_elements_override = 100000

    test_ellipse(num_elements=num_elements_override)
