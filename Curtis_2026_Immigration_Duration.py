
import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import quad
from scipy.special import digamma

# Parameters from authors' MATLAB code
q = 0.5002
p = 1 - q
c = 1 / 24
r = q / p
N = 1000

# METHOD 1: Recurrence
Delta = np.zeros(N)
M = np.zeros(N)

# Starting at 1, time to reach 1 is zero
Delta[0] = 0

for i in range(1, N):
    Delta[i] = r * Delta[i-1] + c / (p*i)

M = np.cumsum(Delta)

# METHOD 2: Analytical
M_an = np.zeros(N)

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

    M_an[k-1] = (
        c / (p*(r-1))
        * (r**k * term1 - term2)
    )

# Plot both methods
population = np.arange(1, N + 1)

plt.figure(figsize=(7, 5))

plt.plot(
    population[1:],
    np.log10(M[1:]),
    "k--",
    label="Recurrence"
)

plt.plot(
    population[1:],
    np.log10(M_an[1:]),
    "r-",
    label="Analytical"
)

plt.xlabel("Target population size (n)")
plt.ylabel("Log10 time (days)")
plt.title(f"Time to grow from 1 to n (q/p = {r:.6f})")

plt.legend()
plt.grid(True, alpha=0.4)
plt.tight_layout()
plt.show()
