# import math
# from enum import Enum
#
# '''
# Exercise 1 -
# '''
#
# # Global Vars #
# eps = 0.0001
# V_na = 55  # mV
# V_k = -72  # mV
# V_l = -49.3  # mV
# g_na = 0.080  # mS
# g_k = 0.030  # mS
# g_l = 0.0003  # mS
# C_m = 1  # uF
# dt = 0.0005  # mSec
# simulation_time = 100  # mSec
#
#
# def power_10(y):
#     return math.pow(10, y)
#
#
# def calc_alpha_m(V_tag):
#     numerator = -0.1 * (35 + V_tag)
#     denum = -1 + math.exp(-(35 + V_tag) / 10)
#     return numerator / denum
#
#
# def calc_alpha_n(V_tag):
#     numerator = -0.01 * (50 + V_tag)
#     denum = -1 + math.exp(-(50 + V_tag) / 10)
#     return numerator / denum
#
# def calc_alpha_h(V_tag):
#     return 0.07*math.exp(-(60+V_tag) / 20)
#
#
# def calc_beta_m(V_tag):
#     return 4 * math.exp(-(60+V_tag) / 18)
#
#
# def calc_beta_n(V_tag):
#     return 0.125 * math.exp(-(60+V_tag) / 80)
#
#
# def calc_beta_h(V_tag):
#     numerator = 1
#     denum = 1 + math.exp(-(30 + V_tag) / 10)
#     return numerator / denum
#
#
# def calc_x_dot(x, alpha, beta):
#     return (alpha * (1 - x)) - (beta * x)
#
#
# def calc_V_dot(m, h, n, V):
#     first_item = -Im
#     second_item = g_na * math.pow(m, 3) * h * (V - V_na)
#     third_item = g_k * math.pow(n, 4) * (V - V_k)
#     fourth_item = g_l * (V - V_l)
#     return (first_item + second_item + third_item + fourth_item) / (-C_m)
#
#
#
# def euler_func_V(V, m, h, n):
#     return V + dt * calc_V_dot(m=m, V=V, n=n, h=h)
#
#
# def euler_func(x, x_dot):
#     return x + (dt*x_dot)
#
#
#
# Im = 0
# time = 0.0
# V = -60
# m = 0
# n = 0
# h = 0
#
# while time < simulation_time:
#     m_dot = calc_x_dot(V, alpha=calc_alpha_m(V), beta=calc_beta_m(V))
#     n_dot = calc_x_dot(V, alpha=calc_alpha_n(V), beta=calc_beta_n(V))
#     h_dot = calc_x_dot(V, alpha=calc_alpha_h(V), beta=calc_beta_h(V))
#     m = euler_func(x=m, x_dot=m_dot)
#     n = euler_func(x=n, x_dot=n_dot)
#     h = euler_func(x=h, x_dot=h_dot)
#     V = euler_func_V(V=V, m=m, n=n, h=h)
#     time += dt
#     print("time = " + str(time)[:6] + " V = " + str(V)[:6] + " m = " + str(m)[:6] + " n = " + str(n)[:6] + " h = " + str(h)[:6])
import numpy as np
import math


def alphaM(V):
    return (-0.1 * (V + 35)) / (np.exp(- (V + 35)/10) - 1)


def betaM(V):
    return 4 * np.exp(-(V + 60) / 18)


def alphaH(V):
    return 0.07 * np.exp(-(V + 60) / 20)


def betaH(V):
    return 1 / (np.exp(-(V + 30)/10) + 1)


def alphaN(V):
    return ((- 0.01) * (V + 50)) / (np.exp(-(V + 50)/10) - 1)


def betaN(V):
    return 0.125 * np.exp(-(V + 60) / 80)


def HH(I0, T0):
    dt = 0.001
    T = math.ceil(T0 / dt)  # [ms]
    gNa0 = 80  # [mS/cm^2]
    ENa = 55  # [mV]
    gK0 = 30  # [mS/cm^2]
    EK = -72  # [mV]
    gL0 = 0.3  # [mS/cm^2]
    EL = -49.3  # [mV]

    t = np.arange(0, T) * dt
    V = np.zeros([T, 1])
    m = np.zeros([T, 1])
    h = np.zeros([T, 1])
    n = np.zeros([T, 1])

    V[0] = -60.0
    m[0] = 0.1
    h[0] = 0.6
    n[0] = 0.3

    for i in range(0, T - 1):
        V[i + 1] = V[i] + dt * (
                    gNa0 * m[i] ** 3 * h[i] * (V[i] - ENa) +
                    gK0 * n[i] ** 4 * (V[i] -EK) +
                    gL0 * (V[i] - EL) + I0)
        m[i + 1] = m[i] + dt * (alphaM(V[i]) * (1 - m[i]) - betaM(V[i]) * m[i])
        h[i + 1] = h[i] + dt * (alphaH(V[i]) * (1 - h[i]) - betaH(V[i]) * h[i])
        n[i + 1] = n[i] + dt * (alphaN(V[i]) * (1 - n[i]) - betaN(V[i]) * n[i])
    print("V = " + str(V[T - 1]))
    print("m = " + str(m[T - 1]))
    print("h = " + str(h[T - 1]))
    print("n = " + str(n[T - 1]))
    print("t = " + str(t[T - 1]))


HH(I0=0, T0=1)
