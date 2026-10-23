
import numpy as np
import matplotlib.pyplot as plt

# Parameters copied from authors' MATLAB code
q = 0.5002
p = 1 - q
c = 1 / 24
r = q / p
N = 1000

# Calculate extinction time using recurrence
Delta = np.zeros(N)
T = np.zeros(N)

# Initial condition
Delta[0] = -(c / p) * np.log(1 - 1 / r)

# Recurrence equation
for i in range(1, N):
    Delta[i] = r * Delta[i - 1] - c / (p * i)

# Cumulative extinction time
T = np.cumsum(Delta)

# Plot extinction duration
plt.figure(figsize=(7, 5))

plt.plot(
    range(1, N + 1),
    np.log10(T),
    color="black",
    label="Recurrence"
)

plt.xlabel("Starting population size (n)")
plt.ylabel("Log10 extinction time (days)")
plt.title(f"Time to extinction (q/p = {r:.6f})")

plt.legend()
plt.grid(True, alpha=0.4)
plt.tight_layout()
plt.show()
