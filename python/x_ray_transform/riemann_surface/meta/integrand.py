@njit(fastmath=True, parallel=True, nogil=True)
def parallel_compute_integrand(domain_min_radius, timeout, step_size, alpha, beta, integrand):
    min_radius_squared = square(domain_min_radius)
    num_max_steps = floor_divide(timeout, step_size)
