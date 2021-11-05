import math
import matplotlib.pyplot as plt

'''
xercise 1 -
'''

# Global Vars #
eps = 0.0001
V_na = 55  # mV
V_k = -72  # mV
V_l = -49.3  # mV
g_na = 80  # mS
g_k = 30  # mS
g_l = 0.3  # mS
C_m = 1  # uF
dt = 0.0005  # mSec


def calc_alpha_m(V_tag):
    numerator = -0.1 * (35 + V_tag)
    denum = -1 + math.exp(-(35 + V_tag) / 10)
    return numerator / denum


def calc_alpha_n(V_tag):
    numerator = -0.01 * (50 + V_tag)
    denum = -1 + math.exp(-(50 + V_tag) / 10)
    return numerator / denum


def calc_alpha_h(V_tag):
    return 0.07 * math.exp(-(60 + V_tag) / 20)


def calc_beta_m(V_tag):
    return 4 * math.exp(-(60 + V_tag) / 18)


def calc_beta_n(V_tag):
    return 0.125 * math.exp(-(60 + V_tag) / 80)


def calc_beta_h(V_tag):
    numerator = 1
    denum = 1 + math.exp(-(30 + V_tag) / 10)
    return numerator / denum


def calc_x_dot(x, alpha, beta):
    return (alpha * (1 - x)) - (beta * x)


def calc_V_dot(m, h, n, V, Im):
    first_item = -Im
    second_item = g_na * math.pow(m, 3) * h * (V - V_na)
    third_item = g_k * math.pow(n, 4) * (V - V_k)
    fourth_item = g_l * (V - V_l)
    return (first_item + second_item + third_item + fourth_item) / (-C_m)


def euler_func_V(V, m, h, n,Im):
    return V + dt * calc_V_dot(m=m, V=V, n=n, h=h,Im=Im)


def euler_func(x, x_dot):
    return x + (dt * x_dot)


def Hodgkings_Huksley(Im, simulation_time, M0=0, N0=0, H0=0, pulse_time = 0, pulse_length = 0):
    time = 0.0
    V = -60
    m = M0
    n = N0
    h = H0
    while time < simulation_time:
        if pulse_time < time < pulse_time + pulse_length:
            I = Im
        else:
            I = 0
        m_dot = calc_x_dot(m, alpha=calc_alpha_m(V), beta=calc_beta_m(V))
        n_dot = calc_x_dot(n, alpha=calc_alpha_n(V), beta=calc_beta_n(V))
        h_dot = calc_x_dot(h, alpha=calc_alpha_h(V), beta=calc_beta_h(V))
        m = euler_func(x=m, x_dot=m_dot)
        n = euler_func(x=n, x_dot=n_dot)
        h = euler_func(x=h, x_dot=h_dot)
        V = euler_func_V(V=V, m=m, n=n, h=h,Im=I)
        time += dt
    return V

simulation_times = [x * 0.1 for x in range(170, 340)]

# 0. Resting
end_value = []
for i in simulation_times:
    end_value.append(float(Hodgkings_Huksley(Im=0, simulation_time=i, M0=0.05, N0=0.34, H0=0.54)))


fig = plt.figure()

plt.plot(simulation_times, end_value)
plt.show()

end_value = []
for i in simulation_times:
    end_value.append(float(Hodgkings_Huksley(Im=0, simulation_time=i)))

fig = plt.figure()

plt.plot(simulation_times, end_value)
plt.show()


# 1. Resting

end_value = []
for i in simulation_times:
    end_value.append(float(Hodgkings_Huksley(Im=30, simulation_time=i, M0=0.05, N0=0.34, H0=0.54, pulse_time=25, pulse_length=0.001)))

fig = plt.figure()

plt.plot(simulation_times, end_value)
plt.show()