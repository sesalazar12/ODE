#!/usr/bin/env python
import numpy as np



def euler(x0,t,f):

    """
    Example:

    Args:

    Returns

    """

    h = t[1] - t[0]
    x = np.zeros(t.size)
    x[0] = x0
    for i in range(t.size - 1):
        x[i+1] = x[i] + (h * f(x[i],t[i])
    return x


def rk2(x0,t,f):
    """
    Example:

    Args:

    Returns:

    """

    h = t[1] - t[0]
    x = np.zeros(t.size)
    x[0] = x0
    for i in range(t.size-1):
        k1 = h * f(x[i],t[i])
        k2 = h * f(x[i],t[i])
        x[i+1] = x[i] + k2
    return x

