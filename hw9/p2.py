import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from scipy.stats import linregress as linreg

path = lambda beta: f"p2_beta{str(beta).replace(".", "p")}.csv"

beta_vals = [0.1, 1, 10] # K / s

R = 8.3145 # J / mol K

init_covs = []
Tp_vals = []


for beta in beta_vals:
    bpath = path(beta)
    df = pd.read_csv(bpath)
    T = df["T (K)"] # K, Seem to be uniformly spaced!
    t = T/beta
    t -= t[0]
    h = 1/beta # increments of 1k to increments in seconds
    r = df["Release Rate (1/s)"]

    init_coverage = sum(r)*h
    theta = np.array([
        init_coverage - sum(r[:i+1])*h for i in range(len(r))
    ])
    init_covs.append(init_coverage)
    print(f"Initial coverage: {init_coverage}")
    # plt.plot(T/beta, r, label="Release Rate over time")
    # plt.xlabel("t (s)")
    # plt.ylabel("rrate (1/s)")
    # plt.legend()
    # plt.show()

    # Coverage over time
    # plt.plot(t, theta, label=f"β = {beta}")

    imax = r.argmax()
    Tp = T[imax]
    Tp_vals.append(Tp)

#     N_EST = 2
#     N = N_EST
#     y = np.log(r / (theta**N))
#     x = 1 / T
#     plt.plot(x, y, label=f"β = {beta}")
# plt.xlabel("t (s)")
# plt.ylabel("θ (-)")
# plt.legend()
# plt.show()

x = np.log(np.array(Tp_vals) / np.array(beta_vals))
y = 1 / np.array(Tp_vals)

res = linreg(x, y)
m = res.slope
E_A = m * R # K * J/molK = J / mol
print(f"Activation energy: {E_A} J / mol")
plt.scatter(x, y, label="Peaks")
plt.xlabel("1/Tp")
plt.ylabel("ln(Tp/β)")
plt.legend()
plt.show()

