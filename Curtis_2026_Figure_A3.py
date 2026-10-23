
import numpy as np
import matplotlib.pyplot as plt

# Parameters from Curtis et al. (2026)
q_base = 0.5002
c = 1 / 24
N = 1000

# Reproducible random number generator
seed = 42
rng = np.random.default_rng(seed)

# Noise range from Figure A3 caption
noise = 0.0002

def extinction_time(q, N):
    p = 1 - q
    r = q / p

    Delta = np.zeros(N)

    Delta[0] = -(c / p) * np.log(1 - 1/r)

    for i in range(1, N):
        Delta[i] = r * Delta[i-1] - c / (p*i)

    return np.cumsum(Delta)

# Baseline calculation
T_base = extinction_time(q_base, N)

plt.figure(figsize=(7, 5))

# Generate 20 random noise trajectories
for run in range(20):
    q_noise = q_base + rng.uniform(-noise, noise)

    # Keep q above 0.5 for this extinction formula
    if q_noise <= 0.5:
        continue

    T_noise = extinction_time(q_noise, N)

    plt.plot(
        np.arange(1, N+1),
        np.log10(T_noise),
        color="red",
        alpha=0.15
    )

# Plot baseline on top
plt.plot(
    np.arange(1, N+1),
    np.log10(T_base),
    color="black",
    linewidth=2.5,
    label="Baseline"
)

plt.xlabel("n")
plt.ylabel("Time log10(days)")
plt.title("Time to extinction q=0.5002")

plt.xlim(0, 1000)
plt.ylim(-0.5, 2.5)

plt.legend()
plt.tight_layout()
plt.show()

print("Random seed:", seed)
