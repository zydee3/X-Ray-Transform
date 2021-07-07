import time
from numba import objmode, njit

default_time_message = "Total Time Elapsed"


@njit(fastmath=True)
def time_computation(entity=None, arg1=None, arg2=None, entity_name=""):
    with objmode(start='f8'):
        start = time.perf_counter()

    entity.compute_values(arg1, arg2)
    # entity_type = entity.__class__.__name__

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
