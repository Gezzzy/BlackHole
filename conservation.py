import numpy as np
import matplotlib.pyplot as plt

# Constants (must match your C code)
G = 1.0
M = 1.0

# Load data: lambda:x:y:z:vr:vtheta:vphi
data = np.loadtxt("output.txt", delimiter=":")

lam = data[:, 0]
x = data[:, 1]
y = data[:, 2]
z = data[:, 3]
vr = data[:, 4]
vtheta = data[:, 5]
vphi = data[:, 6]

# Reconstruct spherical coordinates
r = np.sqrt(x**2 + y**2 + z**2)
theta = np.arccos(z / r)

# Precompute trig
sin_theta = np.sin(theta)

# Metric factor
f = 1.0 - (2.0 * G * M) / r

# --- Angular Momentum ---
L = r**2 * (sin_theta**2) * vphi

# --- Energy ---
E = np.sqrt(
    f * (
            1.0
            + (vr**2) / f
            + r**2 * (vtheta**2 + (sin_theta**2) * vphi**2)
        )
)

# Initial values
L0 = L[0]
E0 = E[0]

# Relative errors
L_err = (L - L0) / L0
E_err = (E - E0) / E0

# =========================
# PLOTS
# =========================

# Angular momentum
plt.figure()
plt.plot(lam, L)
plt.xlabel("Proper time (λ)")
plt.ylabel("L")
plt.title("Angular Momentum Conservation")
plt.grid()

# Energy
plt.figure()
plt.plot(lam, E)
plt.xlabel("Proper time (λ)")
plt.ylabel("E")
plt.title("Energy Conservation")
plt.grid()

# Relative errors
plt.figure()
plt.plot(lam, L_err, label="L error")
plt.plot(lam, E_err, label="E error")
plt.xlabel("Proper time (λ)")
plt.ylabel("Relative error")
plt.title("Conservation Errors")
plt.legend()
plt.grid()

plt.show()
