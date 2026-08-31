import csv
from pathlib import Path
import matplotlib.pyplot as plt

DATA = Path(__file__).resolve().parents[1] / "first_1000_prime_square_plus_one.csv"
OUTPUT = Path(__file__).resolve().parent / "p2_plus_1_curve_first_1000_primes.png"

p_values = []
y_values = []

with DATA.open("r", encoding="utf-8-sig", newline="") as f:
    reader = csv.DictReader(f)
    for row in reader:
        p_values.append(int(row["p"]))
        y_values.append(int(row["p^2 + 1"]))

fig, ax = plt.subplots(figsize=(10, 6))
ax.plot(p_values, y_values, linewidth=1.2)
ax.scatter(p_values, y_values, s=5)
ax.set_xlabel("p (prime)")
ax.set_ylabel(r"$p^2 + 1$")
ax.set_title(r"$p^2 + 1$ for the first 1000 primes")
ax.grid(True, alpha=0.25)
fig.tight_layout()
fig.savefig(OUTPUT, dpi=180)
print(OUTPUT)
