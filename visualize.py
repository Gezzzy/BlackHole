
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation, PillowWriter

data = np.loadtxt("output.txt", delimiter=":")

lam = data[:, 0]
x = data[:, 1]
y = data[:, 2]
z = data[:, 3]
vr = data[:, 4]
vtheta = data[:, 5]
vphi = data[:, 6]

step = 10
x = x[::step]
y = y[::step]
z = z[::step]

fig = plt.figure()
ax = fig.add_subplot(projection='3d')

max_range = np.max([np.max(np.abs(x)), np.max(np.abs(y)), np.max(np.abs(z))])
padding = 0.1 * max_range
ax.set_xlim(-max_range-padding, max_range+padding)
ax.set_ylim(-max_range-padding, max_range+padding)
ax.set_zlim(-max_range-padding, max_range+padding)

ax.set_xlabel("X")
ax.set_ylabel("Y")
ax.set_zlabel("Z")

# Draw black hole
r_s = 2 * 1 * 1

phi_s = np.linspace(0, 2*np.pi, 30)
theta_s = np.linspace(0, np.pi, 30)
phi_s, theta_s = np.meshgrid(phi_s, theta_s)

X_s = r_s * np.sin(theta_s) * np.cos(phi_s)
Y_s = r_s * np.sin(theta_s) * np.sin(phi_s)
Z_s = r_s * np.cos(theta_s)

ax.plot_surface(X_s, Y_s, Z_s, color='k', alpha=1.0)

line, = ax.plot([], [], [], lw=1, color='blue')
point, = ax.plot([], [], [], marker='o', color='red')

def init():
    line.set_data([], [])
    line.set_3d_properties([])
    point.set_data([], [])
    point.set_3d_properties([])
    return line, point

def update(frame):
    print("Frame: ", frame, "/", len(x))
    print("\033[A\033[A")

    line.set_data(x[:frame], y[:frame])
    line.set_3d_properties(z[:frame])

    point.set_data(x[frame:frame+1], y[frame:frame+1])
    point.set_3d_properties(z[frame:frame+1])

    # Rotate camera for better visualization
    ax.view_init(elev=30, azim=frame * 0.3)

    return line, point

ani = FuncAnimation(
    fig,
    update,
    frames=len(x),
    init_func=init,
    interval=20,
    blit=True
)

ani.save("orbit.gif", writer=PillowWriter(fps=30))

plt.show()
