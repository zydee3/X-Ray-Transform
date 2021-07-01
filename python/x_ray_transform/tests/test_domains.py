from x_ray_transform import Circle, Cosine, Ellipse
from timeit import default_timer as timer

default_time_message = "Total Time Elapsed"


def print_time(message=default_time_message, time=0):
    print(message + ': {0:.4g} Seconds'.format(time))


def unit_test_domains(domain):
    pass


def test_domains():
    start = timer()
    circle = Circle(radius=2)
    cosine = Cosine(radius=2, amplitude=1, cycles=1)
    ellipse = Ellipse(minor_radius=1, major_radius=2, theta_offset=3)
    print_time(time=timer() - start)
