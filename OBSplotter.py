import numpy as np
import math
import matplotlib.pyplot as plt 
import os

# ==========================================
# 1. CONFIGURACIÓN GENERAL
# ==========================================
suffix = "3D_NVT_SPHERE_N100_rho0.250RAD_v0.0_pert_period"

obs_path = f"results/obs_{suffix}.dat"
tray_path = f"results/tray_{suffix}.dat"

output_path_4panel = f"results/plot_energy_{suffix}.png"
output_path_temp = f"results/plot_temp_{suffix}.png"
output_path_tray = f"results/plot_tray_{suffix}.png"

# ==========================================
# 2. LECTURA DE OBSERVABLES Y CÁLCULOS
# ==========================================
try:
    t_obs, K, U, E, T = np.loadtxt(obs_path, delimiter=" ", skiprows=1, unpack=True)
except FileNotFoundError:
    print(f"Error: No se encontró el archivo {obs_path}")
    exit()

# Cálculos de fluctuaciones y promedios
E_0 = E[0]
delta_E_inst = 100.0 * (E - E_0) / np.abs(E_0) if E_0 != 0 else np.zeros_like(E)

n_steps = np.arange(1, len(E) + 1)
K_avg, U_avg, E_avg = np.cumsum(K)/n_steps, np.cumsum(U)/n_steps, np.cumsum(E)/n_steps
T_avg = np.cumsum(T)/n_steps # Promedio acumulado de temperatura
delta_E_avg = 100.0 * (E_avg - E_0) / np.abs(E_0) if E_0 != 0 else np.zeros_like(E)

# ==========================================
# 3. FIGURA 1: PANEL 2x2 (ANÁLISIS DE ENERGÍA)
# ==========================================
fig1, axs = plt.subplots(2, 2, figsize=(10, 7), gridspec_kw={'wspace': 0.15, 'hspace': 0.15})
fig1.suptitle(f'Energy Analysis & Fluctuations\n[{suffix}]', fontsize=14, fontweight='bold')

# E(t)
axs[0, 0].plot(t_obs, K, 'k-', alpha=0.7, label='$E_{kin}$')
axs[0, 0].plot(t_obs, U, 'r-', alpha=0.7, label='$E_{pot}$')
axs[0, 0].plot(t_obs, E, 'g-', linewidth=2, label='$E_{tot}$')
axs[0, 0].set_ylabel('E(t)'); axs[0, 0].grid(True); axs[0, 0].legend()

# ΔE(t)
axs[0, 1].plot(t_obs, delta_E_inst, 'g-', linewidth=1)
axs[0, 1].axhline(0, color='black', lw=1); axs[0, 1].set_ylabel(r'$\Delta E(t)$ (%)')
axs[0, 1].yaxis.tick_right(); axs[0, 1].yaxis.set_label_position("right"); axs[0, 1].grid(True)

# <E(t)>
axs[1, 0].plot(t_obs, K_avg, 'k--', label=r'$\langle E_{kin} \rangle$')
axs[1, 0].plot(t_obs, U_avg, 'r--', label=r'$\langle E_{pot} \rangle$')
axs[1, 0].plot(t_obs, E_avg, 'g-', linewidth=2, label=r'$\langle E_{tot} \rangle$')
axs[1, 0].set_ylabel(r'$\langle E(t) \rangle$'); axs[1, 0].set_xlabel('Time'); axs[1, 0].grid(True); axs[1, 0].legend()

# Δ<E(t)>
axs[1, 1].plot(t_obs, delta_E_avg, 'g-', linewidth=1.5)
axs[1, 1].axhline(0, color='black', lw=1); axs[1, 1].set_ylabel(r'$\Delta \langle E(t) \rangle$ (%)')
axs[1, 1].set_xlabel('Time'); axs[1, 1].yaxis.tick_right(); axs[1, 1].yaxis.set_label_position("right"); axs[1, 1].grid(True)

fig1.savefig(output_path_4panel, dpi=300, bbox_inches='tight')

# ==========================================
# 4. FIGURA 2: TEMPERATURA (Solo para 3D NVT)
# ==========================================
if "3D_" in suffix and "NVT" in suffix:
    fig2, ax_t = plt.subplots(figsize=(8, 5))
    ax_t.plot(t_obs, T, color='purple', alpha=0.3, label='Instantaneous T')
    ax_t.plot(t_obs, T_avg, color='darkviolet', linewidth=2, label='Accumulated Average $\langle T \rangle$')
    
    ax_t.set_title(f'Temperature Evolution (NVT Thermostat)\n{suffix}', fontweight='bold')
    ax_t.set_xlabel('Time')
    ax_t.set_ylabel('Temperature')
    ax_t.grid(True, linestyle='--', alpha=0.7)
    ax_t.legend()
    
    fig2.savefig(output_path_temp, dpi=300, bbox_inches='tight')
    print(f"Gráfica de temperatura guardada en: {output_path_temp}")

# ==========================================
# 5. FIGURA 3: TRAYECTORIAS (Solo para 1D y 2D)
# ==========================================
if "3D_" not in suffix:
    # (Mantiene la lógica de lectura de trayectorias que ya teníamos para 1D/2D)
    try:
        t_tray = []
        particle_positions = {}  
        with open(tray_path, 'r') as f:
            lines = f.readlines()
        i = 0
        while i < len(lines):
            line = lines[i].strip()
            if not line: i += 1; continue
            N_part = int(line); i += 1
            curr_t = float(lines[i].strip().split('=')[1]); t_tray.append(curr_t); i += 1
            for _ in range(N_part):
                p = lines[i].split()
                p_id = int(p[0])
                if p_id not in particle_positions: particle_positions[p_id] = {'x': [], 'y': []}
                particle_positions[p_id]['x'].append(float(p[1]))
                particle_positions[p_id]['y'].append(float(p[2]))
                i += 1
        
        fig3, ax3 = plt.subplots(figsize=(7, 4))
        if "2D_" in suffix:
            for p_id, pos in particle_positions.items():
                r = [math.sqrt(x**2 + y**2) for x, y in zip(pos['x'], pos['y'])]
                ax3.plot(t_tray, r, '.', markersize=1.5)
            ax3.set_ylabel('$r_i$')
        else:
            for p_id, pos in particle_positions.items():
                ax3.plot(t_tray, pos['x'])
            ax3.set_ylabel('Position x')
            
        ax3.set_title(f'Trajectories/Distances - {suffix}', fontweight='bold')
        ax3.set_xlabel('Time'); ax3.grid(True)
        fig3.savefig(output_path_tray, dpi=300, bbox_inches='tight')
    except FileNotFoundError:
        pass

plt.show()








#os.system("cd compound && convert -delay 10 -loop 0 *.png compound.gif && xdg-open compound.gif")