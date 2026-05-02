import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from scipy.stats import linregress as linreg

path = lambda T: f"p1_T{str(T)}.csv"

T_vals = [1473, 1373, 1273, 973]

N_VALS = [1, 2, 3] # Potential N values
# N_VALS = [3] # 3 was definitely linear, 2 felt a little more nonliner sometimes

K_vals = []
R = 8.3145 # J / mol K

for T in T_vals:
    tpath = path(T)
    df = pd.read_csv(tpath)

    P = df["P(Pa)"]
    Ctheta = df["ug adsorbed/ g sample"]

    for n in N_VALS:
        x = P**(-1/n)
        y = 1/Ctheta
        plt.plot(x, y, label=f"N={n}")

    n = 3
    x = P**(-1/n)
    y = 1/Ctheta
    # plt.plot(x, y, label=f"N={n}")

    res = linreg(x, y)
    m = res.slope
    b = res.intercept
    C = 1/b
    K = (C*m)**(1/n)
    K_vals.append(K)

    plt.xlabel("1/(P^(1/N))")
    plt.ylabel(f"1/Cθ, for T = {T} K")
    plt.legend()
    plt.show()

print(K_vals)
print(T_vals)
y = np.log(np.array(K_vals))
x = 1 / np.array(T_vals)
res = linreg(x, y)
m = res.slope


plt.plot(x, y, label=f"Isotherm observations")
plt.xlabel("1/T")
plt.ylabel("ln(K)")

plt.legend()
plt.show()
dH = -m*R
print(f"ΔH⁰: {dH} J / mol")

    