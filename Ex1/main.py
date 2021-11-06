import datetime
import math
import matplotlib.pyplot as plt
import datetime as dt
'''
Exercise 1 - Neourons Modeling
Martha Efraim - 
Bar Weitzman - 
Yehonatan Oliker - 206035883
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
RESOLUTION = 0.1
plot_index = 0

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


def Hodgkings_Huksley(Im, simulation_time, M0=0, N0=0, H0=0, pulse_time = [0], pulse_length = 0):
    time = 0.0
    V = -60
    m = M0
    n = N0
    h = H0
    while time < simulation_time:
        for item in pulse_time:
            if item < time < item + pulse_length:
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
    return [float(V), float(m), float(n), float(h)]




# 0. Resting
# Run time - 0-27 seconds, gaps of 0.1 seconds
simulation_times = [x * RESOLUTION for x in range(10, 270)]
V_values = []
current_results = []
for i in simulation_times:
    current_results = Hodgkings_Huksley(Im=0, simulation_time=i, M0=0.05, N0=0.34, H0=0.54)
    V_values.append(current_results[0])


fig = plt.figure()
plt.title("Section 0.A - M_init = 0.05, N_init = 0.34, H_init = 0.54")
plt.plot(simulation_times, V_values)
plot_index +=1
plt.savefig(str(plot_index) + '.png')

plt.show()


V_values = []
current_results = []
for i in simulation_times:
    current_results = Hodgkings_Huksley(Im=0, simulation_time=i)
    V_values.append(current_results[0])

fig = plt.figure()
plt.title("Section 0.B - M_init = N_init = H_init = 0")
plt.plot(simulation_times, V_values)
plot_index +=1
plt.savefig(str(plot_index) + '.png')

plt.show()


# 1. Active potential
# Run time - 1-50 seconds, gaps of 0.1 seconds (stabilization around 23 seconds)
simulation_times = [x * RESOLUTION for x in range(10, 500)]
V_values = []
h_values = []
n_values = []
m_values = []

current_results = []
for i in simulation_times:
    current_results = Hodgkings_Huksley(Im=30000, simulation_time=i, M0=0.05, N0=0.34, H0=0.54, pulse_time=[25], pulse_length=0.001)
    V_values.append(current_results[0])
    m_values.append(current_results[1])
    n_values.append(current_results[2])
    h_values.append(current_results[3])

# A - Parameters
fig = plt.figure()
plt.title("Section 1.A - m(t)")
plt.plot(simulation_times, m_values)
plot_index +=1
plt.savefig(str(plot_index) + '.png')

plt.show()

fig = plt.figure()
plt.title("Section 1.A - n(t)")
plt.plot(simulation_times, n_values)
plot_index +=1
plt.savefig(str(plot_index) + '.png')

plt.show()

fig = plt.figure()
plt.title("Section 1.A - h(t)")
plt.plot(simulation_times, h_values)
plot_index +=1
plt.savefig(str(plot_index) + '.png')

plt.show()

# B - Membrane voltage
fig = plt.figure()
plt.title("Section 1.B - Membrane voltage")
plt.plot(simulation_times, V_values)
plot_index +=1
plt.savefig(str(plot_index) + '.png')

plt.show()

# C - G calculations:

G_na_values = []
G_k_values = []
for index, item in enumerate(m_values):
    G_na_values.append(g_na*h_values[index]*m_values[index]**3)
    G_k_values.append(g_k * n_values[index] ** 4)


fig = plt.figure()
plt.title("Section 1.C - g_Na(t)")
plt.plot(simulation_times, G_na_values)
plot_index +=1
plt.savefig(str(plot_index) + '.png')

plt.show()

fig = plt.figure()
plt.title("Section 1.C - g_K(t)")
plt.plot(simulation_times, G_k_values)
plot_index +=1
plt.savefig(str(plot_index) + '.png')

plt.show()



# D - I (current) calculations:

I_na_values = []
I_k_values = []
I_total_values = []

for index, item in enumerate(m_values):
    I_na_values.append(G_na_values[index]*(V_values[index]-V_na))
    I_k_values.append(G_k_values[index] * (V_values[index] - V_na))
    if index > 0:
        I_total_values.append(I_k_values[index] +
                              I_na_values[index] +
                              (g_l*(V_values[index]-V_l)) +
                              (C_m*(V_values[index]-V_values[index-1])/RESOLUTION))

I_total_values_adjusted = [I_total_values[0]]
I_total_values_adjusted.extend(I_total_values)


fig = plt.figure()
plt.title("Section 1.D - i_Na(t)")
plt.plot(simulation_times, I_na_values)
plot_index +=1
plt.savefig(str(plot_index) + '.png')

plt.show()

fig = plt.figure()
plt.title("Section 1.D - i_K(t)")
plt.plot(simulation_times, I_k_values)
plot_index +=1
plt.savefig(str(plot_index) + '.png')

plt.show()

fig = plt.figure()
plt.title("Section 1.D - i_total(t)")
plt.plot(simulation_times, I_total_values_adjusted)
plot_index +=1
plt.savefig(str(plot_index) + '.png')

plt.show()


# 2. Threshold potential finding
# Run time - 20-30 seconds, gaps of 0.1 seconds - only to recognize threshold value (stabilization around 23 seconds)
simulation_times = [x * 0.1 for x in range(200, 300)]

I_optional_values = range(8000, 30000, 1000) ## NEED TO CHANGE 8000 TO 1000
I_max = 0
found_I_max = False
for i_value in I_optional_values:
    current_results = []
    V_values = []
    for i in simulation_times:
        current_results = Hodgkings_Huksley(Im=i_value, simulation_time=i, M0=0.05, N0=0.34, H0=0.54, pulse_time=[25], pulse_length=0.001)
        V_values.append(current_results[0])
        # Membrane voltage
    fig = plt.figure()
    plt.title("Section 2.A - Find Threshold - 1000 gaps, Im = " + str(i_value))
    plt.plot(simulation_times, V_values)
    plot_index +=1
    plt.savefig(str(plot_index) + '.png')

    plt.show()
    for item in V_values:
        if item > 0:
            I_max = i_value
            found_I_max = True
            break
    if found_I_max:
        break

found_I_max = False
I_optional_values = range(I_max-900, I_max, 100)
for i_value in I_optional_values:
    current_results = []
    V_values = []
    for i in simulation_times:
        current_results = Hodgkings_Huksley(Im=i_value, simulation_time=i, M0=0.05, N0=0.34, H0=0.54, pulse_time=[25], pulse_length=0.001)
        V_values.append(current_results[0])
        # Membrane voltage
    fig = plt.figure()
    plt.title("Section 2.A - Find Threshold - 100 gaps, Im = " + str(i_value))
    plt.plot(simulation_times, V_values)
    plot_index +=1
    plt.savefig(str(plot_index) + '.png')

    plt.show()
    for item in V_values:
        if item > 0:
            I_max = i_value
            found_I_max = True
            break
    if found_I_max:
        break


# According to the 30 graphs printed above, we have decided to run on Im = I_max - 100
# Which is the highest value which didn't spike.
# So, repeating graphs just as in section 1:
I_th = I_max - 100


# Run time - 1-50 seconds, gaps of 0.1 seconds (stabilization around 23 seconds)
simulation_times = [x * RESOLUTION for x in range(10, 500)]
V_values = []
h_values = []
n_values = []
m_values = []

current_results = []
for i in simulation_times:
    current_results = Hodgkings_Huksley(Im=I_th, simulation_time=i, M0=0.05, N0=0.34, H0=0.54, pulse_time=[25], pulse_length=0.001)
    V_values.append(current_results[0])
    m_values.append(current_results[1])
    n_values.append(current_results[2])
    h_values.append(current_results[3])

# A - Parameters
fig = plt.figure()
plt.title("Section 2.B.A - m(t)")
plt.plot(simulation_times, m_values)
plot_index +=1
plt.savefig(str(plot_index) + '.png')

plt.show()

fig = plt.figure()
plt.title("Section 2.B.A - n(t)")
plt.plot(simulation_times, n_values)
plot_index +=1
plt.savefig(str(plot_index) + '.png')

plt.show()

fig = plt.figure()
plt.title("Section 2.B.A - h(t)")
plt.plot(simulation_times, h_values)
plot_index +=1
plt.savefig(str(plot_index) + '.png')

plt.show()

# B - Membrane voltage
fig = plt.figure()
plt.title("Section 2.B.B - Membrane voltage")
plt.plot(simulation_times, V_values)
plot_index +=1
plt.savefig(str(plot_index) + '.png')

plt.show()

# C - G calculations:

G_na_values = []
G_k_values = []
for index, item in enumerate(m_values):
    G_na_values.append(g_na*h_values[index]*m_values[index]**3)
    G_k_values.append(g_k * n_values[index] ** 4)


fig = plt.figure()
plt.title("Section 2.B.C - g_Na(t)")
plt.plot(simulation_times, G_na_values)
plot_index +=1
plt.savefig(str(plot_index) + '.png')

plt.show()

fig = plt.figure()
plt.title("Section 2.B.C - g_K(t)")
plt.plot(simulation_times, G_k_values)
plot_index += 1
plt.savefig(str(plot_index) + '.png')

plt.show()



# D - I (current) calculations:

I_na_values = []
I_k_values = []
I_total_values = []

for index, item in enumerate(m_values):
    I_na_values.append(G_na_values[index]*(V_values[index]-V_na))
    I_k_values.append(G_k_values[index] * (V_values[index] - V_na))
    if index > 0:
        I_total_values.append(I_k_values[index] +
                              I_na_values[index] +
                              (g_l*(V_values[index]-V_l)) +
                              (C_m*(V_values[index]-V_values[index-1])/RESOLUTION))

I_total_values_adjusted = [I_total_values[0]]
I_total_values_adjusted.extend(I_total_values)


fig = plt.figure()
plt.title("Section 2.B.D - i_Na(t)")
plt.plot(simulation_times, I_na_values)
plot_index +=1
plt.savefig(str(plot_index) + '.png')

plt.show()

fig = plt.figure()
plt.title("Section 2.B.D - i_K(t)")
plt.plot(simulation_times, I_k_values)
plot_index +=1
plt.savefig(str(plot_index) + '.png')

plt.show()

fig = plt.figure()
plt.title("Section 2.B.D - i_total(t)")
plt.plot(simulation_times, I_total_values_adjusted)
plot_index +=1
plt.savefig(str(plot_index) + '.png')

plt.show()




# 3. Different excitations
# Run time - 1-50 seconds, gaps of 0.1 seconds (stabilization around 23 seconds)
simulation_times = [x * RESOLUTION for x in range(10, 500)]
V_values = []
V_values_double = []

# A. Very big excitation (twice as high as original)
current_results_original = []
current_results_double = []
for i in simulation_times:
    current_results_original = Hodgkings_Huksley(Im=30000, simulation_time=i, M0=0.05, N0=0.34, H0=0.54, pulse_time=[25], pulse_length=0.001)
    V_values.append(current_results_original[0])
    current_results_double = Hodgkings_Huksley(Im=30000, simulation_time=i, M0=0.05, N0=0.34, H0=0.54, pulse_time=[25], pulse_length=0.001)
    V_values_double.append(current_results_double[0])

fig = plt.figure()
plt.title("Section 3.A - Membrane voltage, original and doubled")
plt.plot(simulation_times, V_values)
plt.plot(simulation_times, V_values_double)
plot_index +=1
plt.savefig(str(plot_index) + '.png')

plt.show()


# Run time - 1-50 seconds, gaps of 0.1 seconds (stabilization around 23 seconds)
simulation_times = [x * RESOLUTION for x in range(10, 500)]
V_values = []

# B. 2 pulses lower than threshold value

for i in simulation_times:
    current_results = Hodgkings_Huksley(Im=I_th, simulation_time=i, M0=0.05, N0=0.34, H0=0.54, pulse_time=[25,25.1], pulse_length=0.001)
    V_values.append(current_results[0])


fig = plt.figure()
plt.title("Section 3.B - Double almost threshold pulses")
plt.plot(simulation_times, V_values)
plot_index +=1
plt.savefig(str(plot_index) + '.png')

plt.show()
