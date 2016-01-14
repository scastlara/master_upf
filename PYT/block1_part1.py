#!/usr/bin/python3
# encoding: utf-8
import numpy as np


def get_sphere_volume(radius):
    """Calculate the volume of a sphere

    >>> get_sphere_volume(5)
    523.5987755982989"""
    assert radius > 0
    return 4 / 3 * np.pi * radius ** 3


def recursive_factorial(n):
    """Calculates recursively the factorial

    >>> recursive_factorial(5)
    120"""
    assert n >= 0
    if n > 0:
        return n * recursive_factorial(n - 1)
    elif n == 0:
        return 1
    else:
        raise Exception("Error")


def factorial(n):
    """Calculates the factorial in a non recursive form

    >>> factorial(5)
    120"""
    assert n >= 0
    f = 1
    while(n > 1):
        f *= n
        n -= 1
    return f


def count_up(n, odd=False):
    """
    Counts up from 0 to whatever positive number you want.
    Odd is neither even nor odd.

    >>> count_up(5)
    0
    1
    2
    3
    4
    5
    >>> count_up(5, odd=True)
    0
    1
    3
    5
    """
    assert n > 1
    print(0)
    if odd:
        for i in range(1, n + 1, 2):
            print(i)
    else:
        for i in range(1, n + 1):
            print(i)


def get_final_price(price, discount_percentage=10):
    """Return the final price after applying the discount percentage
    >>> get_final_price(100, 10)
    90.0"""
    assert price > 0
    assert discount_percentage <= 100
    assert discount_percentage >= 1
    return price * (100 - discount_percentage) / 100


if __name__ == '__main__':
    import doctest
    doctest.testmod()
