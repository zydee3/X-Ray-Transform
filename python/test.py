from x_ray_transformation import Surface, Circle, Cosine, Ellipse, Hyperbolic, Euclidean, Gaussian, Polynomial, Sphere
from random import randrange
from timeit import default_timer as timer
import traceback


def test_domains():
    pass


def test_metrics(num_elements: int = 1, x_values=None, y_values=None, lower_bound: int = 0, upper_bound: int = 1000) -> None:
    if x_values is None or y_values is None:
        x_values = [randrange(lower_bound, upper_bound, 1) for _ in range(num_elements)]
        y_values = [randrange(lower_bound, upper_bound, 1) for _ in range(num_elements)]

    metrics = [#Euclidean(),
               #Gaussian([0.25, 0.25], [0.45, -0.45], [0.3, -0.3], [0.3, -0.3]),
               #Hyperbolic(2),
               #Polynomial([-1.1, 0.0, -1.1, 0, 0, 3]),
               Sphere(4)]

    try:
        for metric in metrics:
            start = timer()
            domain = Circle(1)
            surface = Surface(domain, metric)
            surface.metric.compute_values(x_values, y_values)

            # print(surface.metric.__str__())
            end = timer()
            print(end - start)
    except Exception as exception:
        print("Exception Raised: ", exception)
        print(traceback.format_exc())

