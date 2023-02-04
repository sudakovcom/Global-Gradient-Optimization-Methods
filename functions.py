import math
from interval import imath
from interval import fpu
from interval import interval


def F_0(x):
    return -20 * math.exp(-0.2 * math.sqrt(0.5 * (x[0] * x[0] + x[1] * x[1]))) - math.exp(
        0.5 * (math.cos(2 * math.pi * x[0]) + math.cos(2 * math.pi * x[1]))) + math.exp(1) + 20


def F_0_i(x):
    return -20 * imath.exp(-0.2 * imath.sqrt(0.5 * (x[0] * x[0] + x[1] * x[1]))) - imath.exp(
        0.5 * (imath.cos(2 * math.pi * x[0]) + imath.cos(2 * math.pi * x[1]))) + math.exp(1) + 20


def F_1(x):
    return (x[0] ** 2) / 2 - (x[1] ** 2) / 4 + 3


def F_1_i(x):
    return (x[0] ** 2) / 2 - (x[1] ** 2) / 4 + 3


def F_2(x):
    return (1.5 - x[0] + x[0] * x[1]) ** 2 + (2.25 - x[0] + x[0] * x[1] ** 2) ** 2 + (2.625 - x[0] + x[0] * x[1] ** 3) ** 2


def F_2_i(x):
    return (1.5 - x[0] + x[0] * x[1]) ** 2 + (2.25 - x[0] + x[0] * x[1] ** 2) ** 2 + (2.625 - x[0] + x[0] * x[1] ** 3) ** 2


def F_3(x):
    return (x[0] + 2 * x[1] - 7) ** 2 + (2 * x[0] + x[1] - 5) ** 2


def F_3_i(x):
    return (x[0] + 2 * x[1] - 7) ** 2 + (2 * x[0] + x[1] - 5) ** 2


def F_4(x):
    return 0.26 * (x[0] ** 2 + x[1] ** 2) - 0.48 * x[0] * x[1]


def F_4_i(x):
    return 0.26 * (x[0] ** 2 + x[1] ** 2) - 0.48 * x[0] * x[1]


Functions_1 = [F_0, F_0_i, F_1, F_1_i, F_2, F_2_i, F_3, F_3_i, F_4, F_4_i]
