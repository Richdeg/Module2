#библиотека для параметризации функции и построения графиков
from sympy import *

#для решения СЛАУ
from Project1 import Gaussa
import numpy

#для формирования таблицы
from prettytable import PrettyTable

x = symbols('x')
y = (x - 1) ** 6 #исходная функция
ab = [1, 2] #отрезок из крайних точек
K = 3 #степень многочлена
N = 4 #количество интервалов

def approksim(y, ab, N, K):
    #получение массивов опорных точек и значений функции
    global xi, yi
    xi = []
    yi = []
    a = ab[0]
    b = ab[1]
    h = (b - a) / N  # шаг
    while a < b:
        xi.append(a)
        yi.append(y.evalf(subs={'x': a}))
        a += h
    if round(a - h, 5) != b:
        xi.append(b)
        yi.append(y.evalf(subs={'x': b}))

    #поиск констант
    A = numpy.zeros([K + 1, K + 1])
    for i in range(K + 1):
        for j in range(K + 1):
            for t in xi: A[i, j] += t ** (2 * K - i - j)
    a = numpy.zeros([K + 1])
    for i in range(K + 1):
        for j in range(len(xi)): a[i] += yi[j] * xi[j] ** (K - i)
    const = []
    for j in Gaussa.Metod_Gaussa(A, a): const.append(float(j))

    #составление уравнения
    f = 0
    for i in range(len(const)):
        f += const[i] * x ** (len(const) - 1 - i)

    return f

#графики
p = plot(y, (x, ab[0], ab[1]), line_color='r', show=False)
p.extend(plot(approksim(y, ab, N, K), (x, ab[0], ab[1]), show=False))
p.show()

#значения линейной функции, полученной МНК, в тех же точках
f1 = []
for i in xi: f1.append(approksim(y, ab, N, 1).evalf(subs={'x': i}))
#отклонения (невязки в точках)
e1 = []
for i in range(len(xi)): e1.append(abs(f1[i] - yi[i]))

#значения квадратичной функции, полученной МНК, в тех же точках
f2 = []
for i in xi: f2.append(approksim(y, ab, N, 2).evalf(subs={'x': i}))
#отклонения (невязки в точках)
e2 = []
for i in range(len(xi)): e2.append(abs(f2[i] - yi[i]))

#значения кубической функции, полученной МНК, в тех же точках
f3 = []
for i in xi: f3.append(approksim(y, ab, N, 3).evalf(subs={'x': i}))
#отклонения (невязки в точках)
e3 = []
for i in range(len(xi)): e3.append(abs(f3[i] - yi[i]))

#таблица
th = ['x', 'f(x)', 'f1', 'e1', 'f2', 'e2', 'f3', 'e3']
table = PrettyTable(th)
data = [[] for i in range(N + 1)]
for i in range(N + 1):
    data[i].append(xi[i])
    data[i].append(yi[i])
    data[i].append(f1[i])
    data[i].append(e1[i])
    data[i].append(f2[i])
    data[i].append(e2[i])
    data[i].append(f3[i])
    data[i].append(e3[i])

#построчное заполнение таблицы
for i in range(N + 1):
    dat = data[:]
    while dat:
        table.add_row(dat[i][:len(th)])
        dat = dat[i][len(th):]
print(table)

#Суммы квадратов отклонений (невязок)
e21 = 0
for i in e1: e21 += i ** 2
e22 = 0
for i in e2: e22 += i ** 2
e23 = 0
for i in e3: e23 += i ** 2
print(e21)
print(e22)
print(e23)
