import time

from numba import objmode, njit

default_time_message = "Total Time Elapsed"


@njit(fastmath=True)
def time_computation(entity=None, arg1=None, arg2=None, entity_name=""):
    with objmode(start='f8'):
        start = time.perf_counter()

    entity.compute_values(arg1, arg2)

    # just nice formatting
    with objmode():
        end = time.perf_counter()
        num_spaces = 3
        length_name = len(entity_name)
        num_tabs = num_spaces - (length_name - (length_name % 4))/4
        tabs = ""
        for i in range(int(num_tabs)):
            tabs += '\t'

        print('Time Elapsed ({} - {}):{} {}'.format("entity_type", entity_name, tabs, end - start))


@njit(fastmath=True)
def time_geo_step(surface, x_values, y_values, theta_values, domain_name="", metric_name=""):
    with objmode(start='f8'):
        start = time.perf_counter()

    surface.geo_step(x_values, y_values, theta_values)

    with objmode():
        end = time.perf_counter()
        print('Time Elapsed (Surface - Domain: {}, Metric: {}): {}'.format(domain_name, metric_name, end - start))
