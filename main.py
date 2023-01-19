from sympy import *
from sympy.abc import x
from interval import imath
from interval import fpu
from interval import interval


def F(x):
    return (6 * imath.exp(x) / x + imath.sin(x) * x - imath.cos(x) * x * x) / x / x


def F2(set):
    return (set[0] * set[0]) + (set[0] * set[1]) - (set[1] * set[1]) + set[0]


def F1(set):
    return (6 * imath.exp(set[0]) / set[0] + imath.sin(set[0]) * set[0] - imath.cos(set[0]) * set[0] * set[0]) / set[0] / set[0]


def Multidimensional_Moore_Skelboe(e, set_d):
    interval_e = F2(set_d)
    list_of_intervals = [[set_d, interval_e]]

    U = fpu.max(interval_e)[1]
    w = fpu.max(interval_e)[1] - fpu.max(interval_e)[0]
    best_set = list_of_intervals[0]

    while w > e:
        list_of_intervals.pop(0)
        set_1 = best_set[0].copy()
        set_2 = best_set[0].copy()
        sep_index = 0
        max_len = fpu.max(best_set[0][0])[1] - fpu.max(best_set[0][0])[0]
        for i in range(len(best_set[0])):
            if max_len < fpu.max(best_set[0][i])[1] - fpu.max(best_set[0][i])[0]:
                max_len = fpu.max(best_set[0][i])[1] - fpu.max(best_set[0][i])[0]
                sep_index = i

        mid = (fpu.max(best_set[0][sep_index])[0] + fpu.max(best_set[0][sep_index])[1]) / 2

        set_1[sep_index] = interval[fpu.max(set_1[sep_index])[0], mid]
        set_2[sep_index] = interval[mid, fpu.max(set_2[sep_index])[1]]

        interval_1e = F2(set_1)
        interval_2e = F2(set_2)
        U = min(U, fpu.max(interval_1e)[1])
        U = min(U, fpu.max(interval_2e)[1])
        list_of_intervals.append([set_1, interval_1e])
        list_of_intervals.append([set_2, interval_2e])
        for el in list_of_intervals:
            if U < fpu.max(el[1])[0]:
                list_of_intervals.remove(el)
        list_of_intervals.sort(key=lambda item: fpu.max(item[1])[0])
        best_set = list_of_intervals[0]
        w = fpu.max(best_set[0][1])[1] - fpu.max(best_set[0][1])[0]
    return list_of_intervals[0]


def Moore_Skelboe(e, interval_d):
    interval_e = F(interval_d)
    list_of_intervals = [[interval_d, interval_e]]
    U = fpu.max(interval_e)[1]
    w = len(interval_e)
    best_interval = list_of_intervals[0]
    while w > e:
        list_of_intervals.pop(0)
        mid = (fpu.max(best_interval[0])[0] + fpu.max(best_interval[0])[1]) / 2
        interval_1 = interval[fpu.max(best_interval[0])[0], mid]
        interval_1e = F(interval_1)
        interval_2 = interval[mid, fpu.max(best_interval[0])[1]]
        interval_2e = F(interval_2)
        U = min(U, fpu.max(interval_1e)[1])
        U = min(U, fpu.max(interval_2e)[1])
        list_of_intervals.append([interval_1, interval_1e])
        list_of_intervals.append([interval_2, interval_2e])
        for el in list_of_intervals:
            if U < fpu.max(el[1])[0]:
                list_of_intervals.remove(el)
        list_of_intervals.sort(key=lambda item: fpu.max(item[1])[0])
        best_interval = list_of_intervals[0]
        w = fpu.max(best_interval[0])[1] - fpu.max(best_interval[0])[0]
    return list_of_intervals[0][1]


# def Moore_Skelboe(f, index, D, p):
#     e = 0.001
#     interval_d = D.copy()
#     for i in range(len(p)):
#         interval_d[i] = interval[p[i], p[i]]
#     interval_d[index] = D[index]
#
#     interval_e = f(interval_d)
#     list_of_intervals = [[interval_d, interval_e]]
#     U = fpu.max(interval_e)[1]
#     w = len(interval_e)
#     best_interval = list_of_intervals[0]
#     while w > e:
#         list_of_intervals.pop(0)
#         mid = (fpu.max(best_interval[0][index])[0] + fpu.max(best_interval[0][index])[1]) / 2
#         interval_1 = best_interval[0].copy()
#         interval_2 = best_interval[0].copy()
#         interval_1[index] = interval[fpu.max(best_interval[0][index])[0], mid]
#         interval_1e = f(interval_1)
#         interval_2[index] = interval[mid, fpu.max(best_interval[0][index])[1]]
#         interval_2e = f(interval_2)
#         U = min(U, fpu.max(interval_1e)[1])
#         U = min(U, fpu.max(interval_2e)[1])
#         list_of_intervals.append([interval_1, interval_1e])
#         list_of_intervals.append([interval_2, interval_2e])
#         for el in list_of_intervals:
#             if U < fpu.max(el[1])[0]:
#                 list_of_intervals.remove(el)
#         list_of_intervals.sort(key=lambda item: fpu.max(item[1])[0])
#         best_interval = list_of_intervals[0]
#         w = fpu.max(best_interval[0][index])[1] - fpu.max(best_interval[0][index])[0]
#     return fpu.max(list_of_intervals[0][1])[0]

a = interval[-1, 1]
b = interval[-1, 1]
print(Multidimensional_Moore_Skelboe(0.0001, [a, b]))