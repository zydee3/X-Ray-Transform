from numba.experimental import jitclass
from numba import njit, types
from enum import IntEnum


members = [
    ('domain_id', types.int32),
    ('metric_id', types.int32)
]


class Identifier(IntEnum):
    circle_id = 0,
    cosine_id = 1,
    ellipse_id = 2,
    euclidean_id = 3,
    gaussian_id = 4,
    hyperbolic_id = 5,
    polynomial_id = 6,
    sphere_id = 7


@jitclass(members)
class Surface:
    def __init__(self):
        self.domain_id = -1
        self.metric_id = -1
