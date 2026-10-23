
import numpy as np
import matplotlib.pyplot as plt

# Parameters from the authors' MATLAB code
p = 0.5
q = 0.5
r = q / p
c = 1 / 24
N = 1000

def unbiased_extinction(N, capacity):
    Delta = np.zeros(N)

    # Harmonic-number approximation from MATLAB
    m = capacity - 1
    H = (
        np.log(m) + 0.57721
        + 1 / (2*m)
        - 1 / (12*m**3)
    )

    # Initial condition
    Delta[0] = 1 / N + 2*c*H

    # Recurrence equation
    for i in range(1, N):
        Delta[i] = (
            r * Delta[i-1]
            - c / (p*i)
        )

    return np.cumsum(Delta)

# Calculate extinction durations
T1 = unbiased_extinction(N, 10**16)
T2 = unbiased_extinction(N, 10**4)

population = np.arange(1, N+1)

plt.figure(figsize=(7, 5))

plt.plot(
    population,
    np.log10(T1),
    color="red",
    label="Capacity = 10^16"
)

plt.plot(
    population,
    np.log10(T2),
    color="blue",
    label="Capacity = 10^4"
)

plt.xlabel("Starting population size (n)")
plt.ylabel("Log10 extinction time (days)")
plt.title("Equal birth and death probabilities (p=q=0.5)")

plt.legend()
plt.grid(True, alpha=0.4)
plt.tight_layout()
plt.show()
