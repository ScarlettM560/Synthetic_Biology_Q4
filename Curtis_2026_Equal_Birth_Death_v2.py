
import numpy as np
import matplotlib.pyplot as plt

# Parameters from the authors' MATLAB code
p = q = 0.5
c = 1 / 24
N = 1000
capacities = [10**4, 10**16]

def harmonic_approx(m):
    return (
        np.log(m) + 0.57721
        + 1/(2*m) - 1/(12*m**3)
    )

def recurrence(capacity):
    delta = np.zeros(N)
    H = harmonic_approx(capacity - 1)

    # Initial condition from MATLAB
    delta[0] = 1/N + 2*c*H

    for i in range(1, N):
        delta[i] = delta[i-1] - c/(p*i)

    return np.cumsum(delta)

def analytical(capacity):
    T = np.zeros(N)
    H_capacity = harmonic_approx(capacity - 1)

    for k in range(2, N+1):
        H_i = harmonic_approx(k-1)

        T[k-1] = (
            2*c*k*(H_capacity - H_i)
            + k/capacity
            + 2*c*(H_i - 1)
        )

    return T

population = np.arange(1, N+1)

plt.figure(figsize=(8, 5))

for capacity, color in zip(capacities, ["blue", "red"]):
    T_rec = recurrence(capacity)
    T_an = analytical(capacity)

    plt.plot(
        population,
        np.log10(T_rec),
        color=color,
        linestyle="--",
        label=f"Recurrence n={capacity:.0e}"
    )

    plt.plot(
        population[1:],
        np.log10(T_an[1:]),
        color=color,
        linestyle="-",
        label=f"Analytical n={capacity:.0e}"
    )

    # Numerical comparison for small populations
    print(f"\nCarrying capacity = {capacity:.0e}")

    for k in [2, 5, 10, 100, 1000]:
        difference = abs(T_rec[k-1] - T_an[k-1])
        print(f"Population {k}: difference = {difference:.6f} days")

plt.xlabel("Starting population size")
plt.ylabel("Log10 extinction time (days)")
plt.title("Equal birth/death: analytical vs recurrence")

plt.legend()
plt.grid(True, alpha=0.4)
plt.tight_layout()
plt.show()
