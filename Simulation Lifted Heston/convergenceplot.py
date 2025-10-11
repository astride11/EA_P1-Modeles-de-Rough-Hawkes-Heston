import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import os
os.chdir(r"C:\Users\achle\EA_P1-Modeles-de-Rough-Hawkes-Heston\Simulation Lifted Heston")

df_S = pd.read_csv("S_compare_N.csv")
df_V = pd.read_csv("V_compare_N.csv")

fig, axes = plt.subplots(1, 2, figsize=(14, 6))

# Sous-figure pour S
for col in df_S.columns[1:]:
    axes[0].plot(df_S['time'], df_S[col], label=col)

# Calcul et affichage de la norme sup entre courbes consécutives pour S
normes_S = []
for i in range(1, len(df_S.columns)-1):
    c1 = df_S.columns[i]
    c2 = df_S.columns[i+1]
    dist = np.max(np.abs(df_S[c1] - df_S[c2]))
    normes_S.append(dist)
    axes[0].text(0.98, 0.95 - 0.05*i, f"max|{c1}-{c2}|={dist:.3g}", transform=axes[0].transAxes, ha='right', va='top', fontsize=10, color='black')

axes[0].set_title("Trajectoires S(t) pour différents N")
axes[0].set_xlabel("Temps t")
axes[0].set_ylabel("S(t)")
axes[0].legend()
axes[0].grid(True)

# Sous-figure pour V
for col in df_V.columns[1:]:
    axes[1].plot(df_V['time'], df_V[col], label=col)

# Calcul et affichage de la norme sup entre courbes consécutives pour V
normes_V = []
for i in range(1, len(df_V.columns)-1):
    c1 = df_V.columns[i]
    c2 = df_V.columns[i+1]
    dist = np.max(np.abs(df_V[c1] - df_V[c2]))
    normes_V.append(dist)
    axes[1].text(0.98, 0.95 - 0.05*i, f"max|{c1}-{c2}|={dist:.3g}", transform=axes[1].transAxes, ha='right', va='top', fontsize=10, color='black')

axes[1].set_title("Trajectoires V(t) pour différents N")
axes[1].set_xlabel("Temps t")
axes[1].set_ylabel("V(t)")
axes[1].legend()
axes[1].grid(True)

plt.tight_layout()
plt.savefig("convergence_N.png", dpi=300)
plt.show()