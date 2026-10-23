

import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import quad

# Parameters from the authors' MATLAB code
q = 0.5002
p = 1 - q
c = 1 / 24
r = q / p
N = 1000

# METHOD 1: Recurrence
Delta = np.zeros(N)
T = np.zeros(N)

Delta[0] = -(c / p) * np.log(1 - 1 / r)

for i in range(1, N):
    Delta[i] = r * Delta[i - 1] - c / (p * i)

T = np.cumsum(Delta)

# METHOD 2: Analytical expression
T_an = np.zeros(N)

for k in range(2, N + 1):
    term1, error = quad(
        lambda t: (1 - t**(k-1)) / (1-t),
        0, 1/r
    )

    m = k - 1
    term2 = (
        np.log(m) + 0.57721
        + 1/(2*m) - 1/(12*m**3)
    )

    term3 = (r**k - 1) * np.log(r/(r-1))

    T_an[k-1] = (
        -c/(p*(r-1))
        * (r**k * term1 - term2 - term3)
    )

# Plot both methods
population = np.arange(1, N + 1)

plt.figure(figsize=(7, 5))

plt.plot(
    population,
    np.log10(T),
    "k--",
    label="Recurrence"
)

plt.plot(
    population[1:],
    np.log10(T_an[1:]),
    "r-",
    label="Analytical"
)

plt.xlabel("Starting population size (n)")
plt.ylabel("Log10 extinction time (days)")
plt.title(f"Extinction duration (q/p = {r:.6f})")

plt.legend()
plt.grid(True, alpha=0.4)
plt.tight_layout()
plt.show()

