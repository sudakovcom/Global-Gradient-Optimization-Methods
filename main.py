import numpy as np
import math
from interval import imath
from interval import fpu
from interval import interval
from functions import Functions


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
def GoldenRatio(func_index, index, a, b, p, e_d, e_f):  # f(function), i(index of direction),
    # a(left border), b(right border), p(current point), e_d(error of d), e_f(error of f)
    F = Functions[func_index * 2]
    phi = (1 + math.sqrt(5)) / 2  # constant of golden ratio
    x_1 = b - (b - a) / phi
    x_2 = a + (b - a) / phi
    p_1 = p.copy()  # current point
    p_2 = p.copy()  # current point
    p_1[index] = x_1
    p_2[index] = x_2
    f_1 = F(p_1)  # value in 1-st point
    f_2 = F(p_2)  # value in 2-nd point
    while (b - a > e_d) | (abs(f_1 - f_2) > e_f):  # termination criteria
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

    return best_point  # point of extremum with error e_d


def MooreSkelboe(func_index, index, a, b, p, e_d, e_f):  # F(function), i(index of direction),
    # a(left border), b(right border), p(current point), e_d(error of d), e_f(error of f)
    F = Functions[func_index * 2 + 1]
    interval_d = []
    for i in range(len(p)):
        interval_d.append(interval[p[i], p[i]])
    interval_d[index] = interval[a, b]

    interval_f = F(interval_d)
    set_of_intervals = [[interval_d, interval_f]]
    U = right(interval_f)
    w_f = right(interval_f) - left(interval_f)
    w_d = right(interval_d[index]) - left(interval_d[index])
    best_interval = set_of_intervals[0]
    while (w_d > e_d) | (w_f > e_f):
        set_of_intervals.pop(0)
        mid_p = mid(best_interval[0][index])
        interval_1 = best_interval[0].copy()
        interval_2 = best_interval[0].copy()
        interval_1[index] = interval[left(best_interval[0][index]), mid_p]
        interval_1_f = F(interval_1)
        interval_2[index] = interval[mid_p, right(best_interval[0][index])]
        interval_2_f = F(interval_2)
        U = min(U, right(interval_1_f))
        U = min(U, right(interval_2_f))

        for i in range(len(set_of_intervals)):
            if U < left(set_of_intervals[i][1]):
                set_of_intervals = set_of_intervals[:i]
                break

        set_of_intervals.append([interval_1, interval_1_f])
        set_of_intervals.append([interval_2, interval_2_f])

        set_of_intervals.sort(key=lambda item: left(item[1]))
        best_interval = set_of_intervals[0]
        w_f = right(best_interval[1]) - left(best_interval[1])
        w_d = right(best_interval[0][index]) - left(best_interval[0][index])
        print(w_f)

    best_point = []
    for i in range(len(p)):
        best_point.append(mid(best_interval[0][i]))

    return best_point


def FastSearch(func_index, D, p, e_d, e_f, method):  # F(function), D(set), p(start point), e(error)
    while True:
        p_0 = p
        for index in range(len(p)):
            p = method(func_index, index, D[index][0], D[index][1], p, e_d, e_f)
        if (norm(p_0, p) < e_d) & (Functions[func_index * 2](p_0) - Functions[func_index * 2](p) < e_f):
            break

    best_point = p
    best_value = Functions[func_index * 2](p)

    return best_point, best_value


def print_result(func_num, D, p, e_d, e_f):
    p_1, v_1 = FastSearch(func_num, D, p, e_d, e_f, GoldenRatio)  # point of minimum
    p_2, v_2 = FastSearch(func_num, D, p, e_d, e_f, MooreSkelboe)  # point of minimum
    print("min of Golden Ratio:", p_1, v_1)
    print("min of Moore-Skelboe:", p_2, v_2)


def Gradient(D, p, f, step):
    gradient = np.array([None] * len(p))
    for i in range(len(p)):
        p_ = p.copy()
        p_[i] += step
        gradient[i] = (f(p_) - f(p)) / step
        if (abs(p[i] - D[i][0]) < step) & (p[i] < 0):
            p[i] = 0
        if (abs(p[i] - D[i][1]) < step) & (p[i] > 0):
            p[i] = 0
    return gradient


def FastSearchG(func_index, D, p, e_d, e_f, method):  # F(function), D(set), p(start point), e(error)
    F = Functions[func_index * 2]
    gradient = Gradient(D, p, F, e_d)
    while np.linalg.norm(gradient) > e_f:
        p_0 = p
        p = method(func_index, gradient, D, p, e_d, e_f)
        gradient = Gradient(D, p, F, e_d)
        if norm(p, p_0) < e_d:
            break

    best_point = p
    best_value = Functions[func_index * 2](p)

    return best_point, best_value


def Borders(p, gradient, D):
    max_t = math.inf
    min_t = -math.inf
    dir = gradient / np.linalg.norm(gradient)
    for i in range(len(p)):
        if dir[i] > 0:
            max_t = min((D[i][1] - p[i]) / dir[i], max_t)
            min_t = max((D[i][0] - p[i]) / dir[i], min_t)
        elif dir[i] < 0:
            min_t = max((D[i][1] - p[i]) / dir[i], min_t)
            max_t = min((D[i][0] - p[i]) / dir[i], max_t)

    borders = [p + dir * min_t, p + dir * max_t]
    return borders


def GoldenRatioG(func_index, gradient, D, p, e_d, e_f):  # f(function), i(index of direction),
    # a(left border), b(right border), p(current point), e_d(error of d), e_f(error of f)
    F = Functions[func_index * 2]
    phi = (1 + math.sqrt(5)) / 2  # constant of golden ratio
    a = Borders(p, gradient, D)[0]
    b = Borders(p, gradient, D)[1]
    x_1 = b - (b - a) / phi
    x_2 = a + (b - a) / phi
    p_1 = x_1
    p_2 = x_2
    f_1 = F(p_1)  # value in 1-st point
    f_2 = F(p_2)  # value in 2-nd point
    while (norm(b, a) > e_d) | (abs(f_1 - f_2) > e_f):  # termination criteria
        if f_1 <= f_2:
            b = x_2
            x_2 = x_1
            x_1 = b - (b - a) / phi

            p_1 = x_1
            p_2 = x_2

            f_2 = f_1
            f_1 = F(p_1)
        else:
            a = x_1
            x_1 = x_2
            x_2 = a + (b - a) / phi

            p_1 = x_1
            p_2 = x_2

            f_1 = f_2
            f_2 = F(p_2)

    best_point = []
    for i in range(len(p)):
        best_point.append((p_1[i] + p_2[i]) / 2)

    return best_point  # point of extremum with error e_d



def MooreSkelboeG(func_index, gradient, D, p, e_d, e_f):  # F(function), i(index of direction),
    # a(left border), b(right border), p(current point), e_d(error of d), e_f(error of f)
    F = Functions[func_index * 2 + 1]
    interval_d = []
    for i in range(len(p)):
        interval_d.append(interval[p[i], p[i]])
    interval_d[index] = interval[a, b]

    interval_f = F(interval_d)
    set_of_intervals = [[interval_d, interval_f]]
    U = right(interval_f)
    w_f = right(interval_f) - left(interval_f)
    w_d = right(interval_d[index]) - left(interval_d[index])
    best_interval = set_of_intervals[0]
    while (w_d > e_d) | (w_f > e_f):
        set_of_intervals.pop(0)
        mid_p = mid(best_interval[0][index])
        interval_1 = best_interval[0].copy()
        interval_2 = best_interval[0].copy()
        interval_1[index] = interval[left(best_interval[0][index]), mid_p]
        interval_1_f = F(interval_1)
        interval_2[index] = interval[mid_p, right(best_interval[0][index])]
        interval_2_f = F(interval_2)
        U = min(U, right(interval_1_f))
        U = min(U, right(interval_2_f))

        for i in range(len(set_of_intervals)):
            if U < left(set_of_intervals[i][1]):
                set_of_intervals = set_of_intervals[:i]
                break

        set_of_intervals.append([interval_1, interval_1_f])
        set_of_intervals.append([interval_2, interval_2_f])

        set_of_intervals.sort(key=lambda item: left(item[1]))
        best_interval = set_of_intervals[0]
        w_f = right(best_interval[1]) - left(best_interval[1])
        w_d = right(best_interval[0][index]) - left(best_interval[0][index])
        print(w_f)

    best_point = []
    for i in range(len(p)):
        best_point.append(mid(best_interval[0][i]))

    return best_point


print(FastSearchG(3, [[-10, 10], [-10, 10]], [1, 1], 0.01, 0.01, GoldenRatioG))








