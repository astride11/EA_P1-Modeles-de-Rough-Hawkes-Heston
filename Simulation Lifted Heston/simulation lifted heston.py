import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import os
os.chdir(r"C:\Users\achle\EA_P1-Modeles-de-Rough-Hawkes-Heston\Simulation Lifted Heston")


# Liste des fichiers de simulation
files = [
    'simulation_N1.csv',
    'simulation_N20.csv',
    'simulation_N50.csv',
    'simulation_N100.csv',
    'simulation_N200.csv',
    'simulation_N500.csv'
]

# Valeurs de N correspondantes
N_values = [1, 20, 50, 100, 200, 500]

# Paramètres du modèle (à ajuster selon vos paramètres réels)
params = {
    'kappa': 0.3,
    'theta': 0.05,
    'sigma': 0.1,
    'rho': -0.5,
    'V0': 0.05,
    'H': 0.1,
    'S0': 100.0,
    'r': 0.01,
    'T': 1.0,
    'N_steps': 1000
}

# Créer une figure avec 2 lignes et 3 colonnes
fig, axes = plt.subplots(2, 3, figsize=(18, 10))
fig.suptitle('Lifted Heston Simulations - Variance Process', fontsize=16, fontweight='bold')

# Tracer les variances
for idx, (file, N) in enumerate(zip(files, N_values)):
    row = idx // 3
    col = idx % 3
    ax = axes[row, col]
    
    # Lire le fichier CSV
    df = pd.read_csv(file)
    
    # Calculer dt
    dt = params['T'] / params['N_steps']
    time = np.arange(len(df)) * dt
    
    # Tracer la variance
    ax.plot(time, df['V'], linewidth=1.5, color='blue')
    ax.set_xlabel('Time (years)', fontsize=10)
    ax.set_ylabel('Variance V(t)', fontsize=10)
    ax.set_title(f'N = {N}', fontsize=12, fontweight='bold')
    ax.grid(True, alpha=0.3)
    
    # Ajouter un texte avec les paramètres
    textstr = f'H={params["H"]}\nκ={params["kappa"]}\nθ={params["theta"]}\nσ={params["sigma"]}\nρ={params["rho"]}'
    ax.text(0.02, 0.98, textstr, transform=ax.transAxes, fontsize=8,
            verticalalignment='top', bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.5))

plt.tight_layout()
plt.savefig('variance_comparison.png', dpi=300, bbox_inches='tight')
plt.show()

# Créer une figure pour les prix
fig, axes = plt.subplots(2, 3, figsize=(18, 10))
fig.suptitle('Lifted Heston Simulations - Stock Price Process', fontsize=16, fontweight='bold')

# Tracer les prix
for idx, (file, N) in enumerate(zip(files, N_values)):
    row = idx // 3
    col = idx % 3
    ax = axes[row, col]
    
    # Lire le fichier CSV
    df = pd.read_csv(file)
    
    # Calculer dt
    dt = params['T'] / params['N_steps']
    time = np.arange(len(df)) * dt
    
    # Tracer le prix
    ax.plot(time, df['S'], linewidth=1.5, color='green')
    ax.set_xlabel('Time (years)', fontsize=10)
    ax.set_ylabel('Stock Price S(t)', fontsize=10)
    ax.set_title(f'N = {N}', fontsize=12, fontweight='bold')
    ax.grid(True, alpha=0.3)
    ax.axhline(y=params['S0'], color='r', linestyle='--', alpha=0.5, label=f'S₀={params["S0"]}')
    ax.legend(fontsize=8)
    
    # Ajouter un texte avec les paramètres
    textstr = f'H={params["H"]}\nκ={params["kappa"]}\nθ={params["theta"]}\nσ={params["sigma"]}\nρ={params["rho"]}'
    ax.text(0.02, 0.98, textstr, transform=ax.transAxes, fontsize=8,
            verticalalignment='top', bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.5))

plt.tight_layout()
plt.savefig('price_comparison.png', dpi=300, bbox_inches='tight')
plt.show()

print("Graphiques sauvegardés : variance_comparison.png et price_comparison.png")