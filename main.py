import numpy as np
import math
from interval import imath
from interval import fpu
from interval import interval


def left(interval):
    return fpu.max(interval)[0]


def right(interval):
    return fpu.max(interval)[1]


def mid(interval):
    return (fpu.max(interval)[0] + fpu.max(interval)[1]) / 2


def norm(p_1, p_2):
    r = 0
    for i in range(len(p_1)):
        r += (p_1[i] - p_2[i]) ** 2
    return math.sqrt(r)

def GoldenRatio(F, index, a, b, p, e_d, e_f):  #f(function), i(index of direction),
    # a(left border), b(right border), p(current point), e_d(error of d), e_f(error of f)
    phi = (1 + math.sqrt(5)) / 2  #constant of golden ratio
    x_1 = b - (b - a) / phi
    x_2 = a + (b - a) / phi
    p_1 = p.copy()  #current point
    p_2 = p.copy()  #current point
    p_1[index] = x_1
    p_2[index] = x_2
    f_1 = F(p_1)  #value in 1-st point
    f_2 = F(p_2)  #value in 2-nd point
    while (b - a > e_d) | (abs(f_1 - f_2) > e_f):  #termination criteria
        if f_1 <= f_2:
            b = x_2
            x_2 = x_1
            x_1 = b - (b - a) / phi

            p_1[index] = x_1
            p_2[index] = x_2

            f_2 = f_1
            f_1 = F(p_1)
        else:
            a = x_1
            x_1 = x_2
            x_2 = a + (b - a) / phi

            p_1[index] = x_1
            p_2[index] = x_2

            f_1 = f_2
            f_2 = F(p_2)

    best_point = []
    for i in range(len(p)):
        best_point.append((p_1[i] + p_2[i]) / 2)

    return best_point  #point of extremum with error e_d



def FastSearch(F, D, p, e_d, e_f, method):  #F(function), D(set), p(start point), e(error)
    dimension = len(p)

    while True:
        p_0 = p
        for index in range(dimension):
            p = method(F, index, D[index][0], D[index][1], p, e_d, e_f)
        if norm(p_0, p) < e_d:
            break
        if F(p_0) - F(p) < e_f:
            break

    best_point = p
    best_value = F(p)

    return best_point, best_value


def F(x):
    return (x[0] ** 2) / 2 - (x[1] ** 2) / 4 + 3


p, v = FastSearch(F, [[-10, 10], [-10, 10]], [1, 1], 0.001, 0.001, GoldenRatio)  #point of minimum
print(p)
print(v)
