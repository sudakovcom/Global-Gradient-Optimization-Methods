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


def Multidimensional_Moore_Skelboe(func_index, D, e_d, e_f):
    F = Functions[func_index * 2 + 1]
    interval_f = F(D)
    list_of_intervals = [[D, interval_f]]

    U = right(interval_f)
    w_f = right(interval_f) - left(interval_f)
    best_set = list_of_intervals[0]

    w_d = right(best_set[0][0]) - left(best_set[0][0])
    for i in range(len(best_set[0])):
        if w_d < right(best_set[0][i]) - left(best_set[0][i]):
            w_d = right(best_set[0][i]) - left(best_set[0][i])

    while (w_d > e_d) | (w_f > e_f):
        list_of_intervals.pop(0)
        set_1 = best_set[0].copy()
        set_2 = best_set[0].copy()
        sep_index = 0
        w_d = right(best_set[0][0]) - left(best_set[0][0])
        for i in range(len(best_set[0])):
            if w_d < right(best_set[0][i]) - left(best_set[0][i]):
                w_d = right(best_set[0][i]) - left(best_set[0][i])
                sep_index = i

        mid = (left(best_set[0][sep_index]) + right(
            best_set[0][sep_index])) / 2

        set_1[sep_index] = interval[left(set_1[sep_index]), mid]
        set_2[sep_index] = interval[mid, right(set_2[sep_index])]

        interval_1f = F(set_1)
        interval_2f = F(set_2)
        U = min(U, right(interval_1f))
        U = min(U, right(interval_2f))
        list_of_intervals.append([set_1, interval_1f])
        list_of_intervals.append([set_2, interval_2f])
        for el in list_of_intervals:
            if U < left(el[1]):
                list_of_intervals.remove(el)
        list_of_intervals.sort(key=lambda item: left(item[1]))
        best_set = list_of_intervals[0]
        w_f = right(best_set[1]) - left(best_set[1])

    min_value = (right(best_set[1]) + left(best_set[1])) / 2
    return min_value


n = 1
l = -500
r = 500
D = [interval[l, r]] * n

print(Multidimensional_Moore_Skelboe(18, D, 0.001, 0.001))

print("finished")
