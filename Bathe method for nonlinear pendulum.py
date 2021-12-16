import numpy as np
import matplotlib.pyplot as plt

# -------------------input section------------------------
# the total time to calculate and delta time
t_total = 100
dt = 0.001

# mass and length of the pendulum
mass = 1.0
length = 1.0
# gravitational acceleration
g = 1.0

# initial conditions
x_0 = 0.0
v_0 = 1.999999238456499
a_0 = 0.0

# threshold of newton raphson method
dx_threshold = 0.1**8

# ----------------------------------calculate section------------------------
print(f'The period linear pendulum：{(2*np.pi/np.sqrt(g/length)):_.4f}')

t = 0
u_list = [x_0]  # record u
t_list = [0]   # record t

x_t = x_0
v_t = v_0
a_t = a_0

z0 = 16 / dt**2 * mass * length**2
z1 = mass * g * length
z2 = mass * length**2
z3 = 8 / dt * mass * length**2
z4 = 9 / dt**2 * mass * length**2
z5 = 12 / dt**2 * mass * length**2
z6 = 3 / dt**2 * mass * length**2
z7 = 4 / dt * mass * length**2
z8 = 1 / dt * mass * length**2
while t < t_total:
    x_2t = x_t
    dx = 1.0
    while dx >= dx_threshold:  # first sub-step
        dx = -(z0*(x_2t-x_t)-z3*v_t-z2 * a_t + z1*np.sin(x_2t))/(z0+z1*np.cos(x_2t))  # newton-raphson method
        x_2t += dx
    v_2t = (x_2t-x_t)/(dt/4) - v_t
    dx = 1.0
    x_3t = x_2t
    while dx >= dx_threshold:  # second sub-step
        dx = -(z8*v_t - z7 * v_2t + z6*x_t - z5 * x_2t + z4 * x_3t+z1*np.sin(x_3t))/(z4 + z1*np.cos(x_3t))
        x_3t += dx
    v_3t = x_t/dt - x_2t * 4 / dt + x_3t * 3 / dt
    a_3t = v_t/dt - v_2t * 4 / dt + v_3t * 3 / dt
    x_t = x_3t
    v_t = v_3t
    a_t = a_3t
    t = np.round(t+dt, 5)
    u_list.append(x_3t)
    t_list.append(t)


#  -------------------------find local min/max points -----------------------
for i in range(1, len(t_list)-1):

    t_i = t_list[i]
    u_0 = u_list[i-1]
    u_1 = u_list[i]
    u_2 = u_list[i+1]
    if u_1 < u_0 and u_1 < u_2:
        print(t_i, u_1, 'local min')
    if u_1 > u_0 and u_1 > u_2:
        print(t_i, u_1, 'local max')

#  -----------------------plot figure------------------
plt.figure()
plt.plot(t_list, u_list)
plt.xlabel('t')
plt.ylabel('u')
plt.xlim(0, t_total)
plt.show()
