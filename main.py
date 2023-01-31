import numpy as np
from numpy import linalg as LA
import math


def F(x):
    return (x[0] ** 2) / 2 - (x[1] ** 2) / 4 + 3


def Gradient(p, f, step):
    gradient = np.array([None] * len(p))
    for i in range(len(p)):
        p_ = p.copy()
        p_[i] += step
        gradient[i] = (f(p_) - f(p)) / step
    return gradient


def r(p_1, p_2):
    r_ = 0
    for i in range(len(p_1)):
        r_ += (p_1[i] - p_2[i]) ** 2
    return math.sqrt(r_)


def argmin1(f, gradient, D, p, e):  # f(function), i(index of direction), a(left border), b(right border), p(current point), e(error)
    phi = (1 + math.sqrt(5)) / 2  # constant of golden ratio
    borders = Borders(p, gradient, D)
    a = borders[0]
    b = borders[1]
    p_1 = b - (b - a) / phi
    p_2 = a + (b - a) / phi
    f_1 = f(p_1)  # value in 1-st point
    f_2 = f(p_2)  # value in 2-nd point
    while np.linalg.norm(b - a) > e:  # termination criteria
        if f_1 <= f_2:
            b = p_2
            p_2 = p_1
            p_1 = b - (b - a) / phi

            f_2 = f_1
            f_1 = f(p_1)
        else:
            a = p_1
            p_1 = p_2
            p_2 = a + (b - a) / phi

            f_1 = f_2
            f_2 = f(p_2)
    mid = p_1
    for i in range(len(mid)):
        mid[i] = (p_1[i] + p_2[i]) / 2
    return mid  # point of extremum with error e


def fast_search1(f, D, p, e, e_n, method):  # f(function), D(set), p(start point), e(error)
    p_0 = p
    gradient = Gradient(p, f, 0.000001)
    while np.linalg.norm(gradient) > 0.01:
        p_0 = p
        p = argmin1(f, gradient, D, p, e)
        gradient = Gradient(p, f, 0.000001)
        # cont = False
        # for i in range(len(p)):
        #     if abs(p[i] - D[i][0]) > 0.001 and abs(p[i] - D[i][1]) > 0.001 and abs(gradient[i]) > 0.001:
        #         cont = True
        # if not cont:
        #     break
        if np.linalg.norm(p - p_0) < 0.00001:
            break
    return p


def Borders(p, gradient, D):
    max_t = math.inf
    min_t = -math.inf
    dir = gradient / np.linalg.norm(gradient)
    for i in range(len(p)):
        if dir[i] > 0:
            max_t = min((D[i][1] - p[i]) / dir[i], max_t)
            min_t = max((D[i][0] - p[i]) / dir[i], min_t)
        elif gradient[i] < 0:
            min_t = max((D[i][1] - p[i]) / dir[i], min_t)
            max_t = min((D[i][0] - p[i]) / dir[i], max_t)

    borders = [p + dir * min_t, p + dir * max_t]
    return borders


p = fast_search1(F, [[-10, 10], [-10, 10]], [1, 1], 0.01, 0.001, argmin1)  # point of minimum
y = F(p)  # minimum value
print(p)
print(y)
