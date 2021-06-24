from enum import Enum


class MetricType(Enum):
    euclidean = 0,
    gaussian = 1,
    hyperbolic = 2,
    polynomial = 3,
    sphere = 4
