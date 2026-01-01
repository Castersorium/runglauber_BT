import math
import matplotlib.pyplot as plt

# Beam energies (GeV)
energies = [17.3, 62.4, 200.0]
mp = 0.938

# Compute y_beam
ybeam = [math.acosh(e / (2 * mp)) for e in energies]

# Data: delta y for each alpha
data = {
    1: [4.5199, 4.81795, 5.76754],
    2: [2.80001, 2.91977, 3.39811],
    3: [2.19617, 2.27954, 2.59785],
    4: [1.90186, 1.9551, 2.19879],
    5: [1.71851, 1.76704, 1.96459],
}
# BRAHMS 62.4 GeV
plt.errorbar(
    4.20+0.08, 2.01,
    yerr=[[0.14 + 0.12], [0.14 + 0.12]],  # 合并统计+系统误差
    fmt='s',
    capsize=3,
    label='BRAHMS 62.4 GeV'
)

# BRAHMS 200 GeV（不对称误差）
plt.errorbar(
    5.36+0.08, 2.05,
    yerr=[[0.60], [0.40]],
    fmt='^',
    capsize=3,
    label='BRAHMS 200 GeV'
)

# NA49 17.3 GeV（无误差）
plt.errorbar(
    2.90+0.08, 1.76,
    marker='D',
    label='NA49 17.3 GeV'
)

# Plot points
for alpha, dy in data.items():
    plt.scatter(ybeam, dy, label=f"alpha = {alpha}")

plt.xlim(left=0)
plt.xlim(right=10)
plt.ylim(bottom=0)

# Plot reference line delta y = 0.58 * ybeam
x_line = [0.0] + ybeam
y_line = [0.0] + [0.58 * y for y in ybeam]
plt.plot(x_line, y_line, label=r"$\delta y = 0.58\, y_{\mathrm{beam}}$")


plt.xlabel(r"$y_{\mathrm{beam}}$")
plt.ylabel(r"$rapidity loss <\delta y>$")
plt.legend()
plt.tight_layout()
plt.show()
