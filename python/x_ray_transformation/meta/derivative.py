from typing import Callable


def ddx(f: Callable, x: int, epsilon: float = (2 ** -12)) -> float:
    return (8 * (f(x + epsilon) - f(x - epsilon)) +
            f(x - 2 * epsilon) - f(x + 2 * epsilon)) / (12 * epsilon)


def dddx(f: Callable, x: int, epsilon: float = (2 ** -8)):
    return ((f(x - epsilon * 3) + f(x + epsilon * 3)) / 90 +
            3 * (f(x + epsilon) + f(x - epsilon) -
                 (f(x + epsilon * 2) + f(x - epsilon * 2)) / 10) / 2 -
            49 * f(x) / 18) / (epsilon ** 2)


def ddx(f: Callable, x: int, y: int, epsilon: float = (2 ** -12)):
    return (8 * (f(x + epsilon, y) - f(x - epsilon, y)) +
            f(x - 2 * epsilon, y) - f(x + 2 * epsilon, y)) / (12 * epsilon)


def dddx(f: Callable, x: int , y: int, epsilon: float = (2 ** -8)):
    return ((f(x - epsilon * 3, y) + f(x + epsilon * 3, y)) / 90 +
            3 * (f(x + epsilon, y) + f(x - epsilon, y) -
            (f(x + epsilon * 2, y) + f(x - epsilon * 2, y)) / 10) / 2 -
            49 * f(x, y) / 18) / (epsilon ** 2)


def ddy(f: Callable, x: int, y: int, epsilon: float = (2 ** -12)):
    return (8 * (f(x + epsilon, y) - f(x - epsilon, y)) +
            f(x - 2 * epsilon, y) - f(x + 2 * epsilon, y)) / (12 * epsilon)


def dddy(f: Callable, x: int , y: int, epsilon: float = (2 ** -8)):
    return ((f(x, y - epsilon * 3) + f(x, y + epsilon * 3)) / 90 +
            3 * (f(x, y + epsilon) + f(x, y - epsilon) -
            (f(x, y + epsilon * 2) + f(x, y - epsilon * 2)) / 10) / 2 -
            49 * f(x, y) / 18) / (epsilon ** 2)
