import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import os


plt.style.use('seaborn-v0_8')
fig, axes = plt.subplots(2, 3, figsize=(15, 10))
fig.suptitle('Impact des Paramètres Heston sur le Smile de Volatilité', fontsize=16)

# Chemin vers le dossier contenant les CSV
data_dir = "Heston_Model/Impacts_Heston"

# Lecture des données avec le bon chemin
kappa_df = pd.read_csv(f'{data_dir}/kappa_impact.csv')
theta_df = pd.read_csv(f'{data_dir}/theta_impact.csv')
sigma_df = pd.read_csv(f'{data_dir}/sigma_impact.csv')
rho_df = pd.read_csv(f'{data_dir}/rho_impact.csv')
v0_df = pd.read_csv(f'{data_dir}/v0_impact.csv')
S0 = 100  # Prix spot

# 1. Impact de Kappa
axes[0,0].plot(kappa_df['Strike'], kappa_df['Kappa_0.5'], label='κ=0.5 ', linewidth=2)
axes[0,0].plot(kappa_df['Strike'], kappa_df['Kappa_1.0'], label='κ=1.0', linewidth=2)
axes[0,0].plot(kappa_df['Strike'], kappa_df['Kappa_2.0'], label='κ=2.0 ', linewidth=2)
axes[0,0].set_title('Impact de Kappa\n(Vitesse de retour)')
axes[0,0].set_xlabel('Strike')
axes[0,0].set_ylabel('Volatilité Implicite')
axes[0,0].legend()
axes[0,0].grid(True, alpha=0.3)

# 2. Impact de Theta
axes[0,1].plot(theta_df['Strike'], theta_df['Theta_0.01'], label='θ=0.01 ', linewidth=2)
axes[0,1].plot(theta_df['Strike'], theta_df['Theta_0.04'], label='θ=0.04', linewidth=2)
axes[0,1].plot(theta_df['Strike'], theta_df['Theta_0.09'], label='θ=0.09 ', linewidth=2)
axes[0,1].set_title('Impact de Theta\n(Variance long terme)')
axes[0,1].set_xlabel('Strike')
axes[0,1].legend()
axes[0,1].grid(True, alpha=0.3)

# 3. Impact de Sigma
axes[0,2].plot(sigma_df['Strike'], sigma_df['Sigma_0.1'], label='σ=0.1 ', linewidth=2)
axes[0,2].plot(sigma_df['Strike'], sigma_df['Sigma_0.3'], label='σ=0.3 ', linewidth=2)
axes[0,2].plot(sigma_df['Strike'], sigma_df['Sigma_0.5'], label='σ=0.5 ', linewidth=2)
axes[0,2].set_title('Impact de Sigma\n(Volatilité de la volatilité)')
axes[0,2].set_xlabel('Strike')
axes[0,2].legend()
axes[0,2].grid(True, alpha=0.3)

# 4. Impact de Rho
axes[1,0].plot(rho_df['Strike'], rho_df['Rho_-0.9'], label='ρ=0.5 ', linewidth=2)
axes[1,0].plot(rho_df['Strike'], rho_df['Rho_-0.5'], label='ρ=-0.5 ', linewidth=2)
axes[1,0].plot(rho_df['Strike'], rho_df['Rho_0.0'], label='ρ=0.0 ', linewidth=2)
axes[1,0].set_title('Impact de Rho\n(Corrélation spot/vol)')
axes[1,0].set_xlabel('Strike')
axes[1,0].set_ylabel('Volatilité Implicite')
axes[1,0].legend()
axes[1,0].grid(True, alpha=0.3)

# 5. Impact de V0
axes[1,1].plot(v0_df['Strike'], v0_df['V0_0.01'], label='v₀=0.01 ', linewidth=2)
axes[1,1].plot(v0_df['Strike'], v0_df['V0_0.04'], label='v₀=0.04 ', linewidth=2)
axes[1,1].plot(v0_df['Strike'], v0_df['V0_0.09'], label='v₀=0.09 ', linewidth=2)
axes[1,1].set_title('Impact de V0\n(Variance initiale)')
axes[1,1].set_xlabel('Strike')
axes[1,1].legend()
axes[1,1].grid(True, alpha=0.3)

# Suppression du subplot vide
fig.delaxes(axes[1,2])

plt.tight_layout()
script_dir = os.path.dirname(os.path.abspath(__file__))
plt.savefig(os.path.join(script_dir, 'heston_parameters_impact.png'), dpi=300, bbox_inches='tight')
plt.show()