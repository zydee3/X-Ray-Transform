from x_ray_transform.riemann_surface import Metric

from numpy import random, square
from x_ray_transform.constants import type_metric_constant_curvature, type_metric_gaussian, type_metric_polynomial


def generate_metric_constant_curvature(radius=0):
    metric = Metric(type_metric_constant_curvature)
    metric.radius = radius
    metric.radius_squared = square(radius)
    return metric


def generate_metric_gaussian(widths, weights, center_x, center_y):
    metric = Metric(type_metric_gaussian)
    metric.widths = widths
    metric.weights = weights
    metric.center_x = center_x
    metric.center_y = center_y
    return metric


def generate_metric_polynomials(coefficients):
    metric = Metric(type_metric_polynomial)
    metric.coefficients = coefficients
    return metric
