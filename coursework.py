import numpy as np
from interval import imath
from interval import fpu
from interval import interval
from functions import Functions


# возвращает левую границу интервала
def left(interval):
    return fpu.max(interval)[0]


# возвращает правую границу интервала
def right(interval):
    return fpu.max(interval)[1]


# возвращает середину интервала
def mid(interval):
    return (fpu.max(interval)[0] + fpu.max(interval)[1]) / 2


# возвращает евклидово расстояние между 2 векторами
def dist(p_1, p_2):
    return np.sqrt(np.sum(np.square(np.array(p_1) - np.array(p_2))))


def GoldenRatio(func_index, index, a, b, p, e_d,
                e_f):  # f(function), i(index of direction),
    # a(left border), b(right border), p(current point), e_d(error of d), e_f(error of f)
    F = Functions[func_index * 2]
    phi = (1 + np.sqrt(5)) / 2  # constant of golden ratio
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


def MooreSkelboe(func_index, index, a, b, p, e_d,
                 e_f):  # func_index(number of function in list of functions), index(index of direction),
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

        val_1 = left(interval_1_f)
        val_2 = left(interval_2_f)

        if (len(set_of_intervals) == 0) or (val_1 > left(
                set_of_intervals[-1][1])):
            set_of_intervals.append([interval_1, interval_1_f])
        else:
            l = 0
            r = len(set_of_intervals) - 1
            while l < r:
                m = int((l + r) / 2)
                if left(set_of_intervals[m][1]) > val_1:
                    r = m
                else:
                    l = m + 1
            set_of_intervals.insert(l, [interval_1, interval_1_f])

        if (len(set_of_intervals) == 0) or (
                val_2 > left(set_of_intervals[-1][1])):
            set_of_intervals.append([interval_2, interval_2_f])
        else:
            l = 0
            r = len(set_of_intervals) - 1
            while l < r:
                m = int((l + r) / 2)
                if left(set_of_intervals[m][1]) > val_2:
                    r = m
                else:
                    l = m + 1
            set_of_intervals.insert(l, [interval_2, interval_2_f])

        best_interval = set_of_intervals[0]
        w_f = right(best_interval[1]) - left(best_interval[1])
        w_d = right(best_interval[0][index]) - left(best_interval[0][index])

    best_point = []
    for i in range(len(p)):
        best_point.append(mid(best_interval[0][i]))

    return best_point


def Gradient(D, p, f, step):  # градиетн с отражением
    gradient = np.array([0] * len(p))
    p_ = p.copy()
    for i in range(len(p)):
        p_[i] += step
        gradient[i] = (f(p_) - f(p)) / step
        if (abs(p[i] - D[i][0]) < step) & (gradient[i] < 0):
            gradient[i] = 0
        if (abs(p[i] - D[i][1]) < step) & (gradient[i] > 0):
            gradient[i] = 0
        p_[i] -= step
    return gradient


def FastSearch(func_index, D, p, e_d, e_f,
               method):  # F(function), D(set), p(start point), e(error)
    while True:
        p_0 = p
        for index in range(len(p)):
            p = method(func_index, index, D[index][0], D[index][1], p, e_d,
                       e_f)
        if (dist(p_0, p) < e_d) & (
                Functions[func_index * 2](p_0) - Functions[func_index * 2](
            p) < e_f):
            break

    best_point = p
    best_value = Functions[func_index * 2](p)
    return best_point, best_value


def FastSearchG(func_index, D, p, e_d, e_f,
                method):  # F(function), D(set), p(start point), e(error)
    F = Functions[func_index * 2]
    gradient = Gradient(D, p, F, e_d)
    while np.linalg.norm(gradient) > e_f:
        p_0 = p
        p = method(func_index, gradient, D, p, e_d, e_f)
        gradient = Gradient(D, p, F, e_d)
        if dist(p, p_0) < e_d:
            break

    best_point = p
    best_value = Functions[func_index * 2](p)
    return best_point, best_value


def Borders(p, gradient, D):
    max_t = np.inf
    min_t = -np.inf
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


def GoldenRatioG(func_index, gradient, D, p, e_d,
                 e_f):  # f(function), i(index of direction),
    # a(left border), b(right border), p(current point), e_d(error of d), e_f(error of f)
    F = Functions[func_index * 2]
    phi = (1 + np.sqrt(5)) / 2  # constant of golden ratio
    a = Borders(p, gradient, D)[0]
    b = Borders(p, gradient, D)[1]
    x_1 = b - (b - a) / phi
    x_2 = a + (b - a) / phi
    p_1 = x_1
    p_2 = x_2
    f_1 = F(p_1)  # value in 1-st point
    f_2 = F(p_2)  # value in 2-nd point
    while (dist(b, a) > e_d) | (abs(f_1 - f_2) > e_f):  # termination criteria
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


def MooreSkelboeG(func_index, gradient, D, p, e_d,
                  e_f):  # p(current point), e_d(error of d), e_f(error of f)
    F = Functions[func_index * 2 + 1]
    interval_d = [None] * len(p)

    a = Borders(p, gradient, D)[0]
    b = Borders(p, gradient, D)[1]

    for i in range(len(p)):
        interval_d[i] = interval[a[i], b[i]]

    interval_f = F(interval_d)
    set_of_intervals = [[interval_d, interval_f]]
    U = right(interval_f)
    w_f = right(interval_f) - left(interval_f)
    for index in range(len(p)):
        w_d = np.square(right(interval_d[index]) - left(interval_d[index]))
    w_d = np.sqrt(w_d)
    best_interval = set_of_intervals[0]
    while (w_d > e_d) | (w_f > e_f):
        set_of_intervals.pop(0)
        mid_p = [None] * len(p)
        for index in range(len(p)):
            mid_p[index] = mid(best_interval[0][index])

        interval_1 = best_interval[0].copy()
        interval_2 = best_interval[0].copy()
        for index in range(len(p)):
            interval_1[index] = interval[left(best_interval[0][index]), mid_p[index]]
            interval_2[index] = interval[mid_p[index], right(best_interval[0][index])]

        interval_1_f = F(interval_1)
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

        for index in range(len(p)):
            w_d = np.square(right(best_interval[0][index]) - left(best_interval[0][index]))
        w_d = np.sqrt(w_d)
        # print(w_f)

    best_point = []
    for i in range(len(p)):
        best_point.append(mid(best_interval[0][i]))

    return best_point


# print(FastSearchG(3, [[-10, 10], [-10, 10]], [1, 1], 0.01, 0.01, GoldenRatioG))

n = 2
l = -500
r = 500
D = [[l, r]] * n
p = [2] * n
print(FastSearchG(18, D, p, 0.001, 0.001, MooreSkelboeG)[1])
print(FastSearchG(18, D, p, 0.001, 0.001, GoldenRatioG)[1])

print("finished")
