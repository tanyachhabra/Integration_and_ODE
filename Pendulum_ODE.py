############For simple pendulum##########
import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import odeint


# Constants
g = 9.81  # gravitational acceleration (m/s^2)
L = 1.0   # pendulum length (m)
m = 1.0   # pendulum mass (kg)

# Initial conditions
theta0 = np.pi/4  # initial angle (radians)
omega0 = 0.0      # initial angular velocity (rad/s)

# Time array
t = np.linspace(0, 10, 1000)
dt = t[1] - t[0]

# ODEs for the simple pendulum
def derivatives(state, t):
    theta, omega = state
    return [omega, -g/L * np.sin(theta)]

# Euler method
def euler_step(state, dt):
    theta, omega = state
    theta_new = theta + omega * dt
    omega_new = omega - g/L * np.sin(theta) * dt
    return [theta_new, omega_new]

# RK4 method
def rk4_step(state, dt):
    k1 = derivatives(state, 0)
    k2 = derivatives([state[i] + 0.5*dt*k1[i] for i in range(2)], 0.5*dt)
    k3 = derivatives([state[i] + 0.5*dt*k2[i] for i in range(2)], 0.5*dt)
    k4 = derivatives([state[i] + dt*k3[i] for i in range(2)], dt)
    return [state[i] + (dt/6)*(k1[i] + 2*k2[i] + 2*k3[i] + k4[i]) for i in range(2)]

# Energy calculation
def calculate_energy(state):
    theta, omega = state
    kinetic = 0.5 * m * (L * omega)**2
    potential = m * g * L * (1 - np.cos(theta))
    return kinetic + potential


# Solve using odeint
odeint_solution = odeint(derivatives, [theta0, omega0], t)

# Solve using Euler method
euler_solution = [[theta0, omega0]]
for _ in range(1, len(t)):
    euler_solution.append(euler_step(euler_solution[-1], dt))
euler_solution = np.array(euler_solution)

# Solve using RK4 method
rk4_solution = [[theta0, omega0]]
for _ in range(1, len(t)):
    rk4_solution.append(rk4_step(rk4_solution[-1], dt))
rk4_solution = np.array(rk4_solution)



# Calculate energies
odeint_energy = [calculate_energy(state) for state in odeint_solution]
euler_energy = [calculate_energy(state) for state in euler_solution]
rk4_energy = [calculate_energy(state) for state in rk4_solution]

# Plotting
plt.figure(figsize=(12, 10))

# Theta plot
plt.figure()
plt.plot(t, odeint_solution[:, 0], 'r.', markersize=2, label='odeint')
plt.plot(t, euler_solution[:, 0], 'b-', label='Euler')
plt.plot(t, rk4_solution[:, 0], 'g-', label='RK4')
plt.title('Angle (θ) vs Time')
plt.xlabel('Time (s)')
plt.ylabel('Angle (rad)')
plt.legend()
plt.grid(True)
plt.savefig("theta0.png")

# Energy plot
plt.figure()
plt.plot(t, odeint_energy, 'r.', markersize=2, label='odeint')
plt.plot(t, euler_energy, 'b-', label='Euler')
plt.plot(t, rk4_energy, 'g-', label='RK4')
plt.title('Total Energy vs Time')
plt.xlabel('Time (s)')
plt.ylabel('Energy (J)')
plt.legend()
plt.grid(True)
plt.savefig("Energy.png")
plt.show()


euler_error = np.abs(euler_solution[:, 0] - odeint_solution[:, 0])
rk4_error = np.abs(rk4_solution[:, 0] - odeint_solution[:, 0])

# Plotting errors
plt.figure()
plt.plot(t, euler_error, 'b-', label='Euler Error')
plt.plot(t, rk4_error, 'g-', label='RK4 Error')
plt.title('Error vs Time')
plt.xlabel('Time (s)')
plt.ylabel('Error (rad)')
plt.legend()
plt.grid(True)
plt.savefig("error.png")
plt.show()

###############for double pendulum #############
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation

# Constants and initial conditions
fps = 60
steps_per_frame = 32
line_time = 0.25

m_1, l_1 = 0.8, 0.5  # mass and length of pendulum 1
m_2, l_2 = 0.8, 0.5  # mass and length of pendulum 2
g = 9.81  # acceleration due to gravity
friction_coefficient = 0.1  # friction

# Initial angles (radians) and angular velocities (rad/s)
theta_1_0 = np.pi / 4  # Starting angle for pendulum 1 (45 degrees)
theta_2_0 = -np.pi / 4  # Starting angle for pendulum 2 (-45 degrees)
theta_dot_1_0 = 5.0  # Initial angular velocity for pendulum 1 (positive for upward motion)
theta_dot_2_0 = 0  # Initial angular velocity for pendulum 2

u_vector = [theta_1_0, theta_2_0, theta_dot_1_0, theta_dot_2_0]
delta_t = 1 / fps / steps_per_frame
u_vector_time_snapshots = []

def symplectic_euler_step():
    theta_1, theta_2, theta_dot_1, theta_dot_2 = u_vector

    # Calculate angular accelerations using the equations of motion
    theta_double_dot_1 = (-m_2 * l_1 * theta_dot_1**2 * np.sin(theta_1 - theta_2) * np.cos(theta_1 - theta_2) +
                          m_2 * g * np.sin(theta_2) * np.cos(theta_1 - theta_2) -
                          m_2 * l_2 * theta_dot_2**2 * np.sin(theta_1 - theta_2) -
                          (m_1 + m_2) * g * np.sin(theta_1)) / \
                         ((m_1 + m_2) * l_1 - m_2 * l_1 * np.cos(theta_1 - theta_2)**2) - \
                         friction_coefficient * np.sign(theta_dot_1)

    theta_double_dot_2 = (m_2 * l_2 * theta_dot_2**2 * np.sin(theta_1 - theta_2) * np.cos(theta_1 - theta_2) +
                          (m_1 + m_2) * g * np.sin(theta_1) * np.cos(theta_1 - theta_2) +
                          l_1 * theta_dot_1**2 * np.sin(theta_1 - theta_2) * (m_1 + m_2) -
                          g * np.sin(theta_2) * (m_1 + m_2)) / \
                         (l_2 * (m_1 + m_2) - m_2 * l_2 * np.cos(theta_1 - theta_2)**2) - \
                         friction_coefficient * np.sign(theta_dot_2)

    u_vector[2] += theta_double_dot_1 * delta_t
    u_vector[3] += theta_double_dot_2 * delta_t
    u_vector[0] += u_vector[2] * delta_t
    u_vector[1] += u_vector[3] * delta_t

    u_vector_time_snapshots.append(u_vector.copy())

radius_of_graph_axes = max(l_1, l_2) * 2.1
fig, ax = plt.subplots(figsize=(8, 8))

# Setting up the plot for the pendulum motion
ax.set_xlim(-radius_of_graph_axes, radius_of_graph_axes)
ax.set_ylim(-radius_of_graph_axes, radius_of_graph_axes)
ax.set_aspect('equal')
ax.set_facecolor((1, 1, 1))
ax.set_title('Double Pendulum Animation')
line0, = ax.plot([], [], lw=0.5, color='blue')
line1, = ax.plot([], [], lw=0.5, color='red')
bar0, = ax.plot([], [], lw=3, color='black')
bar1, = ax.plot([], [], lw=3, color='black')

def get_axis_coordinates(particle, axis):
    if axis == 0:
        return (l_1 * np.sin(u_vector_time_snapshots[-1][0]) if particle == 0
                else l_1 * np.sin(u_vector_time_snapshots[-1][0]) + l_2 * np.sin(u_vector_time_snapshots[-1][1]))
    else:
        return (-l_1 * np.cos(u_vector_time_snapshots[-1][0]) if particle == 0
                else -l_1 * np.cos(u_vector_time_snapshots[-1][0]) - l_2 * np.cos(u_vector_time_snapshots[-1][1]))

def animate(i):
    for _ in range(steps_per_frame):
        symplectic_euler_step()

    line0.set_data(get_axis_coordinates(0, 0), get_axis_coordinates(0, 1))
    line1.set_data(get_axis_coordinates(1, 0), get_axis_coordinates(1, 1))

    bar0.set_data([0, l_1 * np.sin(u_vector[0])], [0, -l_1 * np.cos(u_vector[0])])
    bar1.set_data([l_1 * np.sin(u_vector[0]), l_1 * np.sin(u_vector[0]) + l_2 * np.sin(u_vector[1])],
                  [-l_1 * np.cos(u_vector[0]), -l_1 * np.cos(u_vector[0]) - l_2 * np.cos(u_vector[1])])

ani = animation.FuncAnimation(fig, animate, frames=500, interval=1000.0 / fps)
plt.show()

# Save the animation
ani.save('double_pendulum_animation.gif', writer='pillow', fps=10)
