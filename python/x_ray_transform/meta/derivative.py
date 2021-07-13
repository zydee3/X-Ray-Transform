from numba import njit
from numpy import power, float_power


@njit(fastmath=True)
def ddx(f, x: int, epsilon: float = float_power(2, -12)) -> float:
    return (8 * (f(x + epsilon) - f(x - epsilon)) +
            f(x - 2 * epsilon) - f(x + 2 * epsilon)) / (12 * epsilon)


@njit(fastmath=True)
def dddx(f, x: int, epsilon: float = float_power(2, -8)):
    return ((f(x - epsilon * 3) + f(x + epsilon * 3)) / 90 +
            3 * (f(x + epsilon) + f(x - epsilon) -
                 (f(x + epsilon * 2) + f(x - epsilon * 2)) / 10) / 2 -
            49 * f(x) / 18) / power(epsilon, 2)


@njit(fastmath=True)
def ddx(f, x: int, y: int, epsilon: float = float_power(2, -12)):
    return (8 * (f(x + epsilon, y) - f(x - epsilon, y)) +
            f(x - 2 * epsilon, y) - f(x + 2 * epsilon, y)) / (12 * epsilon)


@njit(fastmath=True)
def dddx(f, x: int, y: int, epsilon: float = float_power(2, -8)):
    return ((f(x - epsilon * 3, y) + f(x + epsilon * 3, y)) / 90 +
            3 * (f(x + epsilon, y) + f(x - epsilon, y) -
            (f(x + epsilon * 2, y) + f(x - epsilon * 2, y)) / 10) / 2 -
            49 * f(x, y) / 18) / power(epsilon, 2)


@njit(fastmath=True)
def ddy(f, x: int, y: int, epsilon: float = float_power(2, -12)):
    return (8 * (f(x + epsilon, y) - f(x - epsilon, y)) +
            f(x - 2 * epsilon, y) - f(x + 2 * epsilon, y)) / (12 * epsilon)


@njit(fastmath=True)
def dddy(f, x: int, y: int, epsilon: float = float_power(2, -8)):
    return ((f(x, y - epsilon * 3) + f(x, y + epsilon * 3)) / 90 +
            3 * (f(x, y + epsilon) + f(x, y - epsilon) -
            (f(x, y + epsilon * 2) + f(x, y - epsilon * 2)) / 10) / 2 -
            49 * f(x, y) / 18) / power(epsilon, 2)
