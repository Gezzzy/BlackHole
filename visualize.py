
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation, PillowWriter

# Load data
data = np.loadtxt("output.txt", delimiter=":")

lam = data[:, 0]
x = data[:, 1]
y = data[:, 2]
z = data[:, 3]

# Speed-up (frame skipping)
step = 10
x = x[::step]
y = y[::step]
z = z[::step]

# Create figure
fig = plt.figure()
ax = fig.add_subplot(projection='3d')

# Fixed scale: find the max extent
max_range = np.max([np.max(np.abs(x)), np.max(np.abs(y)), np.max(np.abs(z))])

ax.set_xlim(-max_range, max_range)
ax.set_ylim(-max_range, max_range)
ax.set_zlim(-max_range, max_range)

ax.set_xlabel("X")
ax.set_ylabel("Y")
ax.set_zlabel("Z")

# Central mass (black hole)
ax.scatter(0, 0, 0, s=100, color='k')

# Orbit line + moving point
line, = ax.plot([], [], [], lw=1, color='blue')
point, = ax.plot([], [], [], marker='o', color='red')

def init():
    line.set_data([], [])
    line.set_3d_properties([])
    point.set_data([], [])
    point.set_3d_properties([])
    return line, point

def update(frame):
    line.set_data(x[:frame], y[:frame])
    line.set_3d_properties(z[:frame])

    point.set_data(x[frame:frame+1], y[frame:frame+1])
    point.set_3d_properties(z[frame:frame+1])

    # Rotate camera for better visualization
    ax.view_init(elev=30, azim=frame * 0.3)

    return line, point

# Create animation
ani = FuncAnimation(
    fig,
    update,
    frames=len(x),
    init_func=init,
    interval=20,
    blit=True
)

# Save as GIF
ani.save("orbit.gif", writer=PillowWriter(fps=30))

plt.show()
