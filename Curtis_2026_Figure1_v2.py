
import numpy as np
import matplotlib.pyplot as plt

# Figure 1a parameters
i = 1
N = 10000

# Ratios focused around equal birth and death
ratios = np.linspace(0.91, 1.04, 1000)

log_probabilities = []

for r in ratios:
    if r == 1.0:
        log_p = np.log10(i / N)
    elif r < 1.0:
        p = (1 - r**i) / (1 - r**N)
        log_p = np.log10(p)
    else:
        # Numerically stable logarithmic calculation
        log_p = (
            np.log10(r - 1)
            - N * np.log10(r)
            - np.log10(1 - r**(-N))
        )

    log_probabilities.append(log_p)

plt.figure(figsize=(7, 5))
plt.plot(ratios, log_probabilities, color="black")

plt.axvline(1.0, color="black", linestyle=":")
plt.xlabel("Ratio of deaths to births (q/p)")
plt.ylabel("Log10 probability of reaching 10,000")
plt.title("Curtis et al. (2026) - Figure 1a")

plt.ylim(-180, 10)
plt.tight_layout()
plt.show()

