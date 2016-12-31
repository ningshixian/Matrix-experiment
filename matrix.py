#!/usr/bin/python
# -*- coding: utf-8 -*-
# coding=utf-8
import numpy as np
import matplotlib.pyplot as plt
import math
import matplotlib.pyplot as plt


def cal_fanshu(x, type):
    if x != None:
        sum = 0
        if type == 1:
            for i in range(len(x)):
                y = abs(x[i])
                sum += y
            return sum
        elif type == 2:
            for i in range(len(x)):
                y = x[i] * x[i]
                sum += y
            return np.power(sum, 0.5)  # 数组元素求n次方
        elif type == 0:  # 0代表无穷
            max = x[0]
            for i in range(len(x)):
                if max < x[i]:
                    max = x[i]
            return max
        else:
            return '范数类型只可是1或2或0'
    else:
        return '输入不可为空'


def cal_2():
    left = -(10 ** (-15))
    right = (10 ** (-15))
    # x = np.random.randint(left, right, size=10)  # 制定生成随机数范围和数组大小
    x = np.asarray(np.linspace(left, right, 100))  # 均匀分布
    print x
    # f = (np.log(1 + x)) / x
    f1 = []
    for i in range(len(x)):
        if x[i] == 0:
            f1.append(1)
        else:
            f1.append((np.log(1 + x[i])) / x[i])

    f2 = []
    for i in range(len(x)):
        d = 1 + x[i]
        if d == 1:
            f2.append(1)
        else:
            f2.append((np.log(d) / (d - 1)))

    plt.figure()  # 实例化作图变量
    plt.xlabel('x')  # x轴文本
    plt.ylabel('f')  # y轴文本
    plt.axis([left, right, -1, 2])  # x轴范围-12到12，y轴范围-1到1
    plt.grid(True)  # 是否绘制网格线
    plt.plot(x, f1, 'g-', label="fun1")  # 颜色green，形式为线条
    plt.plot(x, f2, 'r-', label="fun2")  # 颜色red
    plt.legend()  # 绘制图例
    plt.show()  # 展示图像


def cal_3(A, x):  # 多项式的系数以及给定点
    coe = A[0]
    for i in range(1, len(A)):
        coe = coe * x + A[i]
    return coe


def cal_4(A):
    (rSize, cSize) = A.shape
    L = np.eye(rSize, cSize)
    U = np.zeros([rSize, cSize])
    for r in range(rSize):
        U[0][r] = A[0][r]  # U的第一行赋值
        if r > 0:
            L[r][0] = A[r][0] / U[0][0]  # L的第一列赋值
    for i in range(1, rSize):
        for j in range(i, rSize):
            sum = 0
            for k in range(0, i - 1):
                sum = sum + L[i][k] * U[k][j]
            U[i][j] = A[i][j] - sum
        if i < rSize:
            for j in range(i + 1, rSize):
                sum = 0
                for k in range(i - 1):
                    sum = sum + L[j][k] * U[k][i]
                L[j][i] = (A[j][i] - sum) / U[i][i]
    return [L, U]


def generateA(n):
    A = np.zeros((n, n))
    for i in range(n):
        for j in range(n):
            x = i + 1
            y = j + 1
            A[i][j] = 1.0 / (x + y - 1.0)
    return np.asarray(A)


def cal_5(A):
    (rSize, cSize) = A.shape
    n = rSize
    L = np.zeros((n, n))
    print L
    L[0][0] = math.sqrt(A[0][0])
    for i in range(1, n):
        for k in range(i):
            sum = 0
            for x in range(i):
                sum += L[i][x] * L[k][x]
            L[i][k] = (A[i][k] - sum) / L[k][k]

        sum = 0
        for k in range(i):
            sum += L[i][k] * L[i][k]
        L[i][i] = np.sqrt(A[i][i] - sum)

    return L, L.T


def cal_6_1(x):
    shape = np.shape(x)
    n = max(shape)

    e1 = np.zeros(n)
    e1[1] = 1

    u = x + cal_fanshu(x, 2) * e1
    u = u / np.sqrt(sum(u * u))
    return u


def cal_6_2(w):
    shape = np.shape(w)
    n = max(shape)

    H = np.eye(n) - 2 * ((w * w.T) / sum(w * w))
    return H


def cal_6_3(A):
    u = cal_6_1(A[:, 0])
    HxA = A - 2 * (u * u.T) * A
    return HxA


def cal_7_jacobi(A, b, epoch=10):
    n = np.shape(A)
    n = n[0]
    D = np.eye(n) * A
    L = -np.tril(A, -1)
    U = -np.triu(A, 1)

    B = np.linalg.inv(D) * (L + U)
    f = np.linalg.inv(D) * b

    x = np.array([0, 0, 0]).T
    oldx = x

    for i in range(epoch):
        x = B * x + f
        error = abs(x - oldx)
        print error
        oldx = x
    return x


def cal_7_GS(A, b, epoch=10):
    n = np.shape(A)
    n = n[0]
    D = np.eye(n) * A
    L = -np.tril(A, -1)
    U = -np.triu(A, 1)

    B = np.linalg.inv(D - L) * U
    f = np.linalg.inv(D - L) * b

    x = np.array([0, 0, 0]).T
    oldx = x

    for i in range(epoch):
        x = B * x + f
        error = abs(x - oldx)
        print error
        oldx = x
    return x


def cal_8(x, la=1, epoch=10):
    x = x
    y = x ^ 3 + 2 * x ^ 2 + 10 * x - 100
    y_derivative = 3 * x ^ 2 + 4 * x + 10
    x = x - la * y / y_derivative

    error = []
    for i in range(epoch):
        oldx = x
        y = x ^ 3 + 2 * x ^ 2 + 10 * x - 100
        y_derivative = 3 * x ^ 2 + 4 * x + 10
        x = x - la * y / y_derivative

        error.append(abs(x - oldx))
    return x, error


def func(x):
    y = float(x ** 3 + 2 * x ** 2 + 10 * x - 100)
    return y


def cal_8_2(x0, x1, epoch=10):
    old0 = x0
    old1 = x1
    derivative = (func(old1) - func(old0)) / float(old1 - old0)
    x = old1 - func(old1) / derivative

    error = []
    for i in range(epoch):
        old0 = old1
        old1 = x
        derivative = (func(old1) - func(old0)) / (old1 - old0)
        x = old1 - func(old1) / derivative

        error.append(abs(x - old1))
    return x, error


def func9(x):
    y = np.exp(x) * np.cos(x) + 2
    return y


def cal_9(lowerBound, upperBound, epoch):
    f0 = func9(lowerBound)
    f1 = func9(upperBound)
    if f0 * f1 > 0:
        root = lowerBound
    x = lowerBound
    for i in range(epoch):
        old = x
        x = (lowerBound + upperBound) / 2.0
        f = func9(x)
        f0 = func9(lowerBound)
        f1 = func9(upperBound)
        if f * f0 <= 0:
            upperBound = x
        if f * f1 <= 0:
            lowerBound = x
    print x


def func_1(x):
    y = np.sin(np.pi * x)
    return y


def cal_10(low, high, n):
    x = np.linspace(low, high, n)
    y = func_1(x)
    previusAverageDeviation = y
    currentAverageDeviation = y
    polCoefficients = np.zeros(n).T

    for i in range(n):
        polCoefficients[i] = previusAverageDeviation[0]
        for j in range(n - i):
            currentAverageDeviation[j] = (previusAverageDeviation[j + 1] - previusAverageDeviation[j]) / (
                x[j + i] - x[j])
        previusAverageDeviation = currentAverageDeviation
    return polCoefficients


def func_11(x):
    y = 1 / (1 + x ** 2)
    return y


def cal_11_lagrange(interpolations, x):
    n = max(np.shape(interpolations))
    y = func_11(interpolations)
    px = np.zeros(np.shape(x))
    for i in range(n):
        numerator = 1.0
        denominator = 1.0
        for j in range(i - 1):
            numerator = numerator * (x - interpolations[j])
            denominator = denominator * (interpolations[i] - interpolations[j])
        for j in range(i + 1, n):
            numerator = numerator * (x - interpolations[j])
            denominator = denominator * (interpolations[i] - interpolations[j])
        px += y[i] * (numerator / denominator)
    return px


def f(x):
    y = np.exp(3 * x) * np.cos(np.pi * x)
    return y


def complexSimpson(lowerBound, upperBound, n):
    x = np.linspace(lowerBound, upperBound, n)
    integration = 0.0
    integration += f(lowerBound)
    integration += f(upperBound)
    for i in range(1, n):
        integration += 2 * f(x[i])
    for i in range(n - 1):
        integration += 4 * f((x[i + 1] + x[i]) / 2.0)
    integration = integration * (upperBound - lowerBound) / (6 * n)
    return integration


def complexT(lowerBound, upperBound, n):
    x = np.linspace(lowerBound, upperBound, n + 1)
    integration = 0.0
    integration += f(lowerBound)
    integration += f(upperBound)

    for i in range(1, n):
        integration += 2 * f(x[i])
    integration = integration * (upperBound - lowerBound) / (2 * n)
    return integration


def f13(x):
    y = x ** 2
    return y


def cal_13_gaussChebyshev(n):
    integration = 0.0
    for i in range(n):
        integration += (np.pi / (n + 1)) * f13(np.cos((2 * i + 1) * np.pi / (2 * n + 2)))
    return integration


def f14(t, y):
    yt = (t * y - y * 2) / ((t + 1) ** 2)
    return yt


def cal_14_euler(startPoint, initialValue, endPoint, stepLength):
    stepCount = (endPoint - startPoint) / stepLength
    u = np.zeros(stepCount)
    u[0] = initialValue

    currentPoint = startPoint

    for i in range(int(stepCount)-1):
        u[i + 1] = u[i] + stepLength * f14(currentPoint, u[i])
        currentPoint += stepLength
    return u


def cal_14_improvedEuler(startPoint, initialValue, endPoint, stepLength):
    stepCount = (endPoint - startPoint) / stepLength
    u = np.zeros(stepCount)
    u[1] = initialValue

    currentPoint = startPoint

    for i in range(int(stepCount)-1):
        k1 = f14(currentPoint, u[i])
        k2 = f14(currentPoint + stepLength, u[i] + stepLength * k1)
        u[i + 1] = u[i] + 0.5 * stepLength * (k1 + k2)
        currentPoint += stepLength
    return u


def cal_14_RungeKutta(startPoint, initialValue, endPoint, stepLength):
    stepCount = (endPoint - startPoint) / stepLength
    u = np.zeros(stepCount)
    u[1] = initialValue

    currentPoint = startPoint

    for i in range(int(stepCount) - 1):
        k1 = f14(currentPoint, u[i])
        k2 = f14(currentPoint + 0.5 * stepLength, u[i] + 0.5 * stepLength * k1)
        k3 = f14(currentPoint + 0.5 * stepLength, u[i] + 0.5 * stepLength * k2)
        k4 = f14(currentPoint + stepLength, u[i] + stepLength * k3)
        u[i + 1] = u[i] + stepLength / 6. * (k1 + 2. * k2 + 2 * k3 + k4)
        currentPoint += stepLength
    return u


if __name__ == "__main__":
    result = cal_14_RungeKutta(0, 2, 1, 0.1)
    x = np.linspace(0, 1, 10)
    plt.figure()
    plt.plot(x, result, 'r')
    plt.grid(True)
    plt.xlabel('x')
    plt.ylabel('y')
    plt.show()







# x = np.linspace(-5, 5, 100)
# for i in range(4, 10):
# interpolations = np.linspace(-5, 5, i)
# plt.figure()
# plt.plot(x, cal_11_lagrange(interpolations, x), 'r')
# plt.grid(True)
#     plt.xlabel('x')
#     plt.ylabel('y')
#     plt.show()



