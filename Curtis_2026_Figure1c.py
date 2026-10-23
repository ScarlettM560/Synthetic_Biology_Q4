
import numpy as np
import matplotlib.pyplot as plt

# Figure 1c: starting with one bacterium
i = 1

# Target population sizes: 10^2 to 10^8
log_N = np.linspace(2, 8, 100)
N_values = 10**log_N

# Death-to-birth ratios shown in the paper
ratios = [1.000100, 1.000075, 1.000050, 1.000025]

plt.figure(figsize=(7, 5))

for r in ratios:
    # Stable calculation of log10 establishment probability
    log_P = (
        np.log10(r - 1)
        - N_values * np.log10(r)
        - np.log10(-np.expm1(-N_values * np.log(r)))
    )

    # Only display values in the published vertical range
    log_P = np.where(log_P >= -50, log_P, np.nan)

    plt.plot(
        log_N, log_P,
        marker="o",
        markersize=3,
        linestyle="None",
        label=f"{r:.6f}"
    )

plt.xlabel("Log10 target population size (N)")
plt.ylabel("Log10 establishment probability")
plt.title("Curtis et al. (2026) - Figure 1c")

plt.xlim(1, 8.5)
plt.ylim(-50, 0)
plt.legend(title="q/p")
plt.grid(True, alpha=0.4)

plt.tight_layout()
plt.show()
