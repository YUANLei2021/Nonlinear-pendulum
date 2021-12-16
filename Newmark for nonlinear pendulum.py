import numpy as np
import matplotlib.pyplot as plt

# ------------------------ input  section -------------------
t_total = np.pi * 4  # total time to calculate
dt = 0.001  # delta time

mass = 1.0  # mass of pendulum
length = 1.0  # length of pendulum
g = 1.0  # gravity acceleration

x_0 = 0.0  # initial displacement
v_0 = 1.999999238456499  # initial velocity
a_0 = 0.0  # initial acceleration

print("period of linear system: ", 2*np.pi/np.sqrt(g/length))

# ---------------------- calculate section -------------------
t = 0
x_list = [x_0]  # record u
t_list = [0]  # record t
delta = 0.5  # parameter of new-mark method
gamma = 0.5  # parameter

z0 = 1 / dt ** 2 / delta
z1 = gamma / dt / delta
z2 = 1 / dt / delta
z3 = 1 / 2 / delta - 1
z4 = gamma / delta - 1
z5 = dt / 2 * (gamma / delta - 2)
z6 = dt * (1 - gamma)
z7 = gamma * dt

y1 = mass*length**2
y2 = mass*g*length

x_t = x_0
v_t = v_0
a_t = a_0
while t < t_total:
    p = y1*z0+y2
    q = y1*(z0 * x_t + z2 * v_t + z3 * a_t)
    x_t2 = q / p
    a_t2 = z0 * (x_t2 - x_t) - z2 * v_t - z3 * a_t
    v_t2 = v_t + z6 * a_t + z7 * a_t2

    x_t = x_t2
    v_t = v_t2
    a_t = a_t2
    x_list.append(x_t2)
    t_list.append(t)
    t += dt

# ------------ find the local min and local max points -------------
for i in range(1, len(t_list)-1):

    t_i = t_list[i]
    x_0 = x_list[i-1]
    x_1 = x_list[i]
    x_2 = x_list[i+1]
    if x_1 < x_0 and x_1 < x_2:
        print(t_i, x_1, 'min')
    if x_1 > x_0 and x_1 > x_2:
        print(t_i, x_1, 'max')

# ------------------   plot figure -----------------------------
plt.figure()
plt.plot(t_list, x_list)
plt.show()
