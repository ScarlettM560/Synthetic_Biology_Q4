
import numpy as np
import matplotlib.pyplot as plt
from scipy.special import digamma

# Parameters
c = 1 / 24
N = 1000
capacities = [10**4, 10**16]

# Exact harmonic numbers using digamma
def H(m):
    if m == 0:
        return 0.0
    return digamma(m + 1) + np.euler_gamma

# Original MATLAB harmonic approximation
def H_approx(m):
    return (
        np.log(m) + 0.57721
        + 1/(2*m) - 1/(12*m**3)
    )

def original_recurrence(capacity):
    delta = np.zeros(N)
    delta[0] = 1/N + 2*c*H_approx(capacity-1)

    for i in range(1, N):
        delta[i] = delta[i-1] - 2*c/i

    return np.cumsum(delta)

def corrected_recurrence(capacity):
    delta = np.zeros(N)
    delta[0] = 1/capacity + 2*c*H(capacity-1)

    for i in range(1, N):
        delta[i] = delta[i-1] - 2*c/i

    return np.cumsum(delta)

def analytical(capacity):
    population = np.arange(1, N+1)
    return (
        2*c*population *
        (H(capacity-1) - np.array(
            [H(i-1) for i in population]
        ))
        + population/capacity
        + 2*c*(population-1)
    )

plt.figure(figsize=(8, 5))

for capacity, color in zip(capacities, ["blue", "red"]):
    original = original_recurrence(capacity)
    corrected = corrected_recurrence(capacity)
    analytic = analytical(capacity)

    population = np.arange(1, N+1)

    plt.plot(
        population, np.log10(original),
        color=color, linestyle=":",
        label=f"Original {capacity:.0e}"
    )

    plt.plot(
        population, np.log10(corrected),
        color=color, linestyle="--",
        label=f"Corrected {capacity:.0e}"
    )

    plt.plot(
        population, np.log10(analytic),
        color=color, linestyle="-",
        label=f"Analytical {capacity:.0e}"
    )

    print(f"\nCapacity: {capacity:.0e}")
    for i in [2, 5, 10, 100, 1000]:
        print(
            f"n={i}: original difference = "
            f"{abs(original[i-1]-analytic[i-1]):.6f}, "
            f"corrected difference = "
            f"{abs(corrected[i-1]-analytic[i-1]):.6f}"
        )

plt.xlabel("Starting population size")
plt.ylabel("Log10 extinction time (days)")
plt.title("Validation of equal birth/death calculations")
plt.legend()
plt.grid(True, alpha=0.4)
plt.tight_layout()
plt.show()
