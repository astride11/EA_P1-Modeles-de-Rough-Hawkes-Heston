import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import os
os.chdir(r"C:\Users\achle\EA_P1-Modeles-de-Rough-Hawkes-Heston\Simulation Lifted Heston")
df = pd.read_csv('vix_paths.csv')

# Tracé des trajectoires de VIX pour chaque path
fig, ax = plt.subplots()
for key, group in df.groupby('path'):
    ax.plot(group['t'], group['VIX'], label=f'Path {int(key)}')
ax.set_xlabel('t')
ax.set_ylabel('VIX')
ax.set_title('Simulation des trajectoires VIX par path')
ax.legend()
plt.savefig('VIX simulation.png', dpi=300, bbox_inches='tight')
plt.show()