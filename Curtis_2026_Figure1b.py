
import numpy as np
import matplotlib.pyplot as plt

# Starting and target populations
i = 1
N = 10000

# Death-to-birth probability ratios
ratios = np.linspace(0.92, 1.04, 1000)

probabilities = []

for r in ratios:
    if r == 1.0:
        P = i / N
    elif r < 1.0:
        P = (1 - r**i) / (1 - r**N)
    else:
        P = (r - 1) * r**(-N) / (1 - r**(-N))

    probabilities.append(P)

# Create Figure 1b
plt.figure(figsize=(7, 5))

plt.plot(ratios, probabilities, color="black")

plt.xlabel("Ratio of deaths to births (q/p)")
plt.ylabel("Probability of reaching 10,000")
plt.title("Curtis et al. (2026) - Figure 1b")

plt.ylim(0, 0.09)
plt.grid(True, alpha=0.4)

plt.tight_layout()
plt.show()
