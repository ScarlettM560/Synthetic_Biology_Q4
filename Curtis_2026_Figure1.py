
import numpy as np
import matplotlib.pyplot as plt

# Starting and target population sizes
i = 1
N = 10000

# Ratio of death probability to birth probability
ratios = np.linspace(0.5, 1.5, 1000)

# Calculate probability of reaching target population
probabilities = []

for r in ratios:
    if np.isclose(r, 1.0):
        P = i / N
    else:
        P = (1 - r**i) / (1 - r**N)

    probabilities.append(P)

# Plot results
plt.plot(ratios, probabilities)

plt.xlabel("Death / Birth probability (q/p)")
plt.ylabel("Probability of reaching target population")
plt.title("Immigration probability - Curtis et al. (2026)")

plt.grid(True)
plt.tight_layout()
plt.show()

