#библиотека для параметризации функции и построения графиков
from sympy import *

#для решения СЛАУ
from Project1 import Gaussa
import numpy

#для решения трехдиагональных СЛАУ
from Project1 import Progonka

#для формирования таблицы
from prettytable import PrettyTable

def lin_interp(y, ab, N):
    #получение массивов точек и значений функции
    global xi, yi
    xi = []
    yi = []
    a = ab[0]
    b = ab[1]
    h = (b - a) / N #шаг
    while a < b:
        xi.append(a)
        yi.append(y.evalf(subs={'x':a}))
        a += h
    if round(a - h, 5) != b:
        xi.append(b)
        yi.append(y.evalf(subs={'x': b}))

    #вычисление уравнений кусков линейного сплайна с привязкой к своим отрезкам
    yn = {}
    for i in range(N):
        ai = (yi[i + 1] - yi[i]) / (xi[i + 1] - xi[i])
        bi = yi[i] - ai * xi[i]
        yn[xi[i], xi[i + 1]] = ai * x + bi

    #построение сплайна
    #p = plot(y, (x, ab[0], ab[1]), line_color = 'r', show = False)
    #for i, j in yn:
        #p.extend(plot(yn[i, j], (x, i, j), show = False))
    #p.show()
    return yn

def left_square_interp(y, ab, N):
    #получение массивов точек и значений функции
    xi = []
    yi = []
    a = ab[0]
    b = ab[1]
    h = (b - a) / N #шаг
    while a < b:
        xi.append(a)
        yi.append(y.evalf(subs={'x':a}))
        a += h
    if round(a - h, 5) != b:
        xi.append(b)
        yi.append(y.evalf(subs={'x': b}))

    #ищем константы уравнений
    const = [[] for i in range(N)]
    for i in range(N):
        A = numpy.array([[1, xi[i], xi[i] ** 2], [1, xi[i + 1], xi[i + 1] ** 2], [0, 1, 2 * xi[i]]])
        if i == 0: a = numpy.array([yi[i], yi[i + 1], 0])
        else: a = numpy.array([yi[i], yi[i + 1], const[i - 1][1] + 2 * const[i - 1][2] * xi[i]])
        for j in Gaussa.Metod_Gaussa(A, a): const[i].append(float(j))

    #вычисление уравнений кусков определённого слева квадратичного сплайна с привязкой к своим отрезкам
    yn = {}
    for i in range(N):
        yn[xi[i], xi[i + 1]] = const[i][0] + const[i][1] * x + const[i][2] * x ** 2

    #построение сплайна

    return yn

def right_square_interp(y, ab, N):
    #получение массивов точек и значений функции
    xi = []
    yi = []
    a = ab[0]
    b = ab[1]
    h = (b - a) / N #шаг
    while a < b:
        xi.append(a)
        yi.append(y.evalf(subs={'x':a}))
        a += h
    if round(a - h, 5) != b:
        xi.append(b)
        yi.append(y.evalf(subs={'x': b}))

    #ищем константы уравнений
    const = [[] for i in range(N)]
    for i in range(N, 0, -1):
        A = numpy.array([[1, xi[i], xi[i] ** 2], [1, xi[i - 1], xi[i - 1] ** 2], [0, 1, 2 * xi[i]]])
        if i == N: a = numpy.array([yi[i], yi[i - 1], 0])
        else: a = numpy.array([yi[i], yi[i - 1], const[N - 1 - i][1] + 2 * const[N - 1 - i][2] * xi[i]])
        for j in Gaussa.Metod_Gaussa(A, a): const[N - i].append(float(j))

    #переворачиваем значения в списке
    const_new = []
    for i in reversed(const): const_new.append(i)

    #вычисление уравнений кусков определённого справа квадратичного сплайна с привязкой к своим отрезкам
    yn = {}
    for i in range(N):
        yn[xi[i], xi[i + 1]] = const_new[i][0] + const_new[i][1] * x + const_new[i][2] * x ** 2

    return yn

def cube_interp(y, ab, N):
    #получение массивов точек и значений функции
    xi = []
    yi = []
    a = ab[0]
    b = ab[1]
    h = (b - a) / N #шаг
    while a < b:
        xi.append(a)
        yi.append(y.evalf(subs={'x':a}))
        a += h
    if round(a - h, 5) != b:
        xi.append(b)
        yi.append(y.evalf(subs={'x': b}))

    #ищем константы ci уравнений
    A = numpy.zeros([N - 1, N - 1])
    A[0, 0] = 2 * (xi[2] - xi[0])
    A[0, 1] = xi[2] - xi[1]
    A[N - 2, N - 3] = xi[N - 1] - xi[N - 2]
    A[N - 2, N - 2] = 2 * (xi[N] - xi[N - 2])
    for i in range(2, N - 1):
        A[i - 1, i - 2] = xi[i] - xi[i - 1]
        A[i - 1, i - 1] = 2 * (xi[i + 1] - xi[i - 1])
        A[i - 1, i] = xi[i + 1] - xi[i]

    a = numpy.zeros([N - 1])
    for i in range(1, N):
        a[i - 1] = 6 * ((yi[i + 1] - yi[i]) / (xi[i + 1] - xi[i]) - (yi[i] - yi[i - 1]) / (xi[i] - xi[i - 1]))

    #с помощью метода прогонки находим константы ci
    ci = [0]
    for j in Progonka.Metod_Progonki(A, a):
        ci.append(float(j))
    ci.append(0)

    #находим константы di
    di = [0]
    for i in range(1, N + 1):
        di.append((ci[i] - ci[i - 1]) / (xi[i] - xi[i - 1]))

    #находим константы bi
    bi = [0]
    for i in range(1, N + 1):
        hi = xi[i] - xi[i - 1]
        bi.append(ci[i] * hi / 2 - di[i] * hi ** 2 / 6 + (yi[i] - yi[i - 1]) / hi)

    #находим константы ai
    ai = []
    for i in yi: ai.append(i)

    #вычисление уравнений кусков кубического сплайна с привязкой к своим отрезкам
    yn = {}
    for i in range(1, N + 1):
        yn[xi[i - 1], xi[i]] = ai[i] + bi[i] * (x - xi[i]) + ci[i] * (x - xi[i]) ** 2 / 2 + di[i] * (x - xi[i]) ** 3 / 6

    return yn

#запуск программы
x = symbols('x')
y = (x - 1) ** 6 #исходная функция
ab = [1, 5] #отрезок из крайних точек

N = 10 #количество интервалов

print(lin_interp(y, ab, N))
print(left_square_interp(y, ab, N))
print(right_square_interp(y, ab, N))
print(cube_interp(y, ab, N))

#таблицы аргументов и значений функции
th = ['x', 'f(x)']
table = PrettyTable(th)
data = [[] for i in range(N + 1)]
for i in range(N + 1):
    data[i].append(xi[i])
    data[i].append(yi[i])

#построчное заполнение таблицы
for i in range(N + 1):
    dat = data[:]
    while dat:
        table.add_row(dat[i][:len(th)])
        dat = dat[i][len(th):]
print(table)

N = 50 #количество интервалов

print(lin_interp(y, ab, N))
print(left_square_interp(y, ab, N))
print(right_square_interp(y, ab, N))
print(cube_interp(y, ab, N))

#таблицы аргументов и значений функции
th = ['x', 'f(x)']
table = PrettyTable(th)
data = [[] for i in range(N + 1)]
for i in range(N + 1):
    data[i].append(xi[i])
    data[i].append(yi[i])

#построчное заполнение таблицы
for i in range(N + 1):
    dat = data[:]
    while dat:
        table.add_row(dat[i][:len(th)])
        dat = dat[i][len(th):]
print(table)

#Поочередное построение сплайнов
p = plot(y, (x, ab[0], ab[1]), show=False)
for i, j in lin_interp(y, ab, 10):
    p.extend(plot(lin_interp(y, ab, 10)[i, j], (x, i, j), line_color='r', show=False))
for i, j in lin_interp(y, ab, 50):
    p.extend(plot(lin_interp(y, ab, 50)[i, j], (x, i, j), line_color='g', show=False))
p.show()

p = plot(y, (x, ab[0], ab[1]), show=False)
for i, j in left_square_interp(y, ab, 10):
    p.extend(plot(left_square_interp(y, ab, 10)[i, j], (x, i, j), line_color='r', show=False))
for i, j in left_square_interp(y, ab, 50):
    p.extend(plot(left_square_interp(y, ab, 50)[i, j], (x, i, j), line_color='g', show=False))
p.show()

p = plot(y, (x, ab[0], ab[1]), show=False)
for i, j in right_square_interp(y, ab, 10):
    p.extend(plot(right_square_interp(y, ab, 10)[i, j], (x, i, j), line_color='r', show=False))
for i, j in right_square_interp(y, ab, 50):
    p.extend(plot(right_square_interp(y, ab, 50)[i, j], (x, i, j), line_color='g', show=False))
p.show()

p = plot(y, (x, ab[0], ab[1]), show=False)
for i, j in cube_interp(y, ab, 10):
    p.extend(plot(cube_interp(y, ab, 10)[i, j], (x, i, j), line_color='r', show=False))
for i, j in cube_interp(y, ab, 50):
    p.extend(plot(cube_interp(y, ab, 50)[i, j], (x, i, j), line_color='g', show=False))
p.show()

#Построение всех сплайнов на одном графике
N = 10
p = plot(y, (x, ab[0], ab[1]), show=False)

for i, j in lin_interp(y, ab, N):
    p.extend(plot(lin_interp(y, ab, N)[i, j], (x, i, j), line_color='r', show=False))

for i, j in left_square_interp(y, ab, N):
    p.extend(plot(left_square_interp(y, ab, N)[i, j], (x, i, j), line_color='g', show=False))

for i, j in right_square_interp(y, ab, N):
    p.extend(plot(right_square_interp(y, ab, N)[i, j], (x, i, j), line_color='b', show=False))

for i, j in cube_interp(y, ab, N):
    p.extend(plot(cube_interp(y, ab, N)[i, j], (x, i, j), line_color='y', show=False))
p.show()
