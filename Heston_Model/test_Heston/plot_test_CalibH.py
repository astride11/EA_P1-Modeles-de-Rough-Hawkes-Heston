import pandas as pd
import matplotlib.pyplot as plt
import math

def plot_smiles(datafile, outname):
    # Lecture du CSV exporté depuis ton code C++
    df = pd.read_csv(datafile)

    # Maturités uniques triées
    maturities = sorted(df["T"].unique())
    n = len(maturities)

    # 2 subplots par ligne
    cols = 2
    rows = math.ceil(n / cols)

    fig, axes = plt.subplots(rows, cols, figsize=(6*cols, 4*rows), sharey=True)
    axes = axes.flatten()  # on met tout dans une seule liste pour itérer facilement

    for i, T in enumerate(maturities):
        ax = axes[i]
        sub = df[df["T"] == T]

        ax.plot(sub["K"], sub["iv_market"], "o-", label="Vol vraie", linewidth=2)
        ax.plot(sub["K"], sub["iv_model"], "s--", label="Vol calibrée", linewidth=2)

        ax.set_title(f"Maturité T = {T:.2f} ans", fontsize=11)
        ax.set_xlabel("Strike K")
        ax.grid(True, linestyle="--", alpha=0.6)

    # Pour les figures inutilisées (si nb de maturités impair)
    for j in range(i + 1, len(axes)):
        fig.delaxes(axes[j])

    axes[0].set_ylabel("Volatilité implicite")
    axes[0].legend(loc="best", fontsize=9)
 
    
    plt.suptitle("Comparaison des smiles de volatilité implicite (Heston)", fontsize=14, weight="bold")
    plt.tight_layout(rect=[0, 0, 1, 0.95])
    plt.savefig(f"Heston_Model/test_Heston/{outname}", dpi=150)



if __name__ == "__main__":
    plot_smiles("Heston_Model/test_Heston/smile_data_calib.csv", "smile_data_calib.png")
    # plot_smiles("smile_data_pert.csv", "smile_data_pert.png")