import numpy as np
import math
import matplotlib.pyplot as plt 
import re
import os
from pathlib import Path
import concurrent.futures 
import time

# Usamos un backend no interactivo por defecto en los workers para que el multiprocesamiento 
# no crashee intentando abrir múltiples ventanas al mismo tiempo.
plt.switch_backend('agg')

# ==========================================
# 1. FUNCIONES DE GRAFICADO (WORKERS PARALELOS)
# ==========================================

def plot_benchmark():
    results_dir = Path("results")
    input_path = results_dir / "benchmark_CT_vs_N.dat"
    output_path = results_dir / "benchmark_analysis_dual.pdf" 

    try:
        data = np.loadtxt(input_path, comments="#")
    except FileNotFoundError:
        return f"[Aviso] No se encontró {input_path}. Saltando Benchmark."

    if len(data) < 2:
        return "[Aviso] No hay suficientes datos en el benchmark para graficar."

    mid = len(data) // 2
    data_NVE, data_NVT = data[:mid], data[mid:]
    
    if len(data_NVE) == 0 or len(data_NVT) == 0:
        return "[Aviso] Formato de benchmark incompleto. Saltando."

    N_NVE, CT_NVE = data_NVE[:, 0], data_NVE[:, 1]
    N_NVT, CT_NVT = data_NVT[:, 0], data_NVT[:, 1]

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(13, 5))

    # RAW
    ax1.plot(N_NVE, CT_NVE, marker='o', linestyle='-', color='blue', linewidth=2, label='NVE')
    ax1.plot(N_NVT, CT_NVT, marker='s', linestyle='-', color='red', linewidth=2, label='NVT')
    ax1.set_xlabel('Number of Particles (N)')
    ax1.set_ylabel('Computing Time (s)')
    ax1.set_title('Raw Performance: N vs CT', fontweight='bold')
    ax1.legend(); ax1.grid(True)

    # LOG-LOG
    log_N_NVE, log_CT_NVE = np.log10(N_NVE), np.log10(CT_NVE)
    log_N_NVT, log_CT_NVT = np.log10(N_NVT), np.log10(CT_NVT)

    m_NVE, b_NVE = np.polyfit(log_N_NVE, log_CT_NVE, 1)
    m_NVT, b_NVT = np.polyfit(log_N_NVT, log_CT_NVT, 1)

    ax2.plot(log_N_NVE, log_CT_NVE, 'bo', markersize=8)
    ax2.plot(log_N_NVT, log_CT_NVT, 'rs', markersize=8)
    ax2.plot(log_N_NVE, (m_NVE * log_N_NVE) + b_NVE, 'b--', lw=2, label=f'NVE O(N^{m_NVE:.2f})')
    ax2.plot(log_N_NVT, (m_NVT * log_N_NVT) + b_NVT, 'r--', lw=2, label=f'NVT O(N^{m_NVT:.2f})')

    ax2.set_xlabel('log10(N)')
    ax2.set_ylabel('log10(CT)')
    ax2.set_title('Linearization: Log-Log Plot', fontweight='bold')
    ax2.legend(); ax2.grid(True)

    plt.tight_layout()
    fig.savefig(output_path, format='pdf', bbox_inches='tight') 
    plt.close(fig)
    return f"¡Exito! Benchmark guardado en {output_path.name}"

def plot_rdf(target_ensemble="NVE"):
    results_dir = Path("results")
    output_path = results_dir / f"plot_rdf_comparison_{target_ensemble}.pdf" 
    file_list = list(results_dir.glob(f"rdf_3D_{target_ensemble}_*.dat"))

    if not file_list:
        return f"[Aviso] No hay archivos RDF para el ensamble {target_ensemble}."

    data_files = []
    for filepath in file_list:
        match = re.search(r'rho([0-9]+\.[0-9]+)', filepath.name)
        if match:
            data_files.append((float(match.group(1)), filepath))
    
    data_files.sort(key=lambda x: x[0])

    fig, ax = plt.subplots(figsize=(8, 5))
    ax.set_title(f'Radial Distribution Function $g(r)$\nEnsemble {target_ensemble}', fontsize=14, fontweight='bold')
    colors = plt.cm.plasma(np.linspace(0, 0.85, len(data_files)))

    max_r, max_g = 10, 3
    for i, (density, filepath) in enumerate(data_files):
        try:
            r_rdf, g_r = np.loadtxt(filepath, skiprows=1, unpack=True)
            ax.plot(r_rdf, g_r, '-', color=colors[i], linewidth=2, label=rf'$\rho = {density:.3f}$')
            max_r, max_g = max(max_r, max(r_rdf)), max(max_g, max(g_r))
        except Exception:
            pass

    ax.axhline(1.0, color='gray', linestyle='--', alpha=0.7, label='Ideal Gas')
    ax.set_xlabel('Distance $r$', fontsize=12)
    ax.set_ylabel('$g(r)$', fontsize=12)
    ax.set_xlim(0, max_r); ax.set_ylim(0, max_g * 1.1)
    ax.grid(True, linestyle=':', alpha=0.6)
    ax.legend(title="Density", loc='lower right', framealpha=0.9)

    plt.tight_layout()
    fig.savefig(output_path, format='pdf', bbox_inches='tight') 
    plt.close(fig)
    return f"¡Exito! RDF ({target_ensemble}) guardado en {output_path.name}"

def plot_observables(suffix):
    results_dir = Path("results")
    obs_path = results_dir / f"obs_{suffix}.dat"
    tray_path = results_dir / f"tray_{suffix}.dat"
    
    out_4panel = results_dir / f"plot_energy_{suffix}.pdf" 
    out_temp = results_dir / f"plot_temp_{suffix}.pdf"     
    out_tray = results_dir / f"plot_tray_{suffix}.pdf"     

    msgs = []
    
    # 1. ENERGÍA
    try:
        t_obs, K, U, E, T = np.loadtxt(obs_path, delimiter=" ", skiprows=1, unpack=True)
        # TODO Futuro: Aquí añadiremos la Presión cuando la incluyas en el motor (Virial)
        
        E_0 = E[0]
        delta_E_inst = 100.0 * (E - E_0) / np.abs(E_0) if E_0 != 0 else np.zeros_like(E)
        n_steps = np.arange(1, len(E) + 1)
        K_avg, U_avg, E_avg, T_avg = np.cumsum(K)/n_steps, np.cumsum(U)/n_steps, np.cumsum(E)/n_steps, np.cumsum(T)/n_steps
        delta_E_avg = 100.0 * (E_avg - E_0) / np.abs(E_0) if E_0 != 0 else np.zeros_like(E)

        fig1, axs = plt.subplots(2, 2, figsize=(12, 8), gridspec_kw={'wspace': 0.25, 'hspace': 0.35})
        fig1.suptitle(f'Energy Analysis [{suffix}]', fontweight='bold', fontsize=14)
        
        # [0,0] Energía Instantánea
        axs[0,0].plot(t_obs, K, 'k-', alpha=0.7, label=r'$E_{kin}$')
        axs[0,0].plot(t_obs, U, 'r-', alpha=0.7, label=r'$E_{pot}$')
        axs[0,0].plot(t_obs, E, 'g-', lw=2, label=r'$E_{tot}$')
        axs[0,0].set_title('Instantaneous Energy')
        axs[0,0].set_xlabel('Time ($t$)'); axs[0,0].set_ylabel('Energy')
        axs[0,0].grid(True); axs[0,0].legend()
        
        # [0,1] Fluctuación de Energía Instantánea (CON FÓRMULA)
        axs[0,1].plot(t_obs, delta_E_inst, 'g-', lw=1)
        axs[0,1].axhline(0, color='k', lw=1)
        axs[0,1].set_title("Instantaneous Fluctuation\n" + r"$\Delta E = \frac{E(t) - E_0}{|E_0|} \times 100$")
        axs[0,1].set_xlabel('Time ($t$)'); axs[0,1].set_ylabel(r'$\Delta E$ (%)')
        axs[0,1].yaxis.tick_right(); axs[0,1].yaxis.set_label_position("right"); axs[0,1].grid(True)
        
        # [1,0] Energías Acumuladas Promedio
        axs[1,0].plot(t_obs, K_avg, 'k--', label=r'$\langle E_{kin} \rangle$')
        axs[1,0].plot(t_obs, U_avg, 'r--', label=r'$\langle E_{pot} \rangle$')
        axs[1,0].plot(t_obs, E_avg, 'g-', lw=2, label=r'$\langle E_{tot} \rangle$')
        axs[1,0].set_title('Cumulative Average Energy')
        axs[1,0].set_xlabel('Time ($t$)'); axs[1,0].set_ylabel(r'$\langle E \rangle$')
        axs[1,0].grid(True); axs[1,0].legend()
        
        # [1,1] Fluctuación de Energía Acumulada Promedio (CON FÓRMULA)
        axs[1,1].plot(t_obs, delta_E_avg, 'g-', lw=1.5)
        axs[1,1].axhline(0, color='k', lw=1)
        axs[1,1].set_title("Average Fluctuation\n" + r"$\langle \Delta E \rangle = \frac{\langle E(t) \rangle - E_0}{|E_0|} \times 100$")
        axs[1,1].set_xlabel('Time ($t$)'); axs[1,1].set_ylabel(r'$\langle \Delta E \rangle$ (%)')
        axs[1,1].yaxis.tick_right(); axs[1,1].yaxis.set_label_position("right"); axs[1,1].grid(True)
        
        fig1.savefig(out_4panel, format='pdf', bbox_inches='tight') 
        plt.close(fig1)
        msgs.append(f"Energía -> {out_4panel.name}")

        # TEMPERATURA (Si es NVT)
        if "NVT" in suffix:
            fig2, ax_t = plt.subplots(figsize=(8, 5))
            ax_t.plot(t_obs, T, color='purple', alpha=0.3, label='Instantaneous T')
            ax_t.plot(t_obs, T_avg, color='darkviolet', lw=2, label=r'$\langle T \rangle$')
            ax_t.set_title(f'Temperature Evolution [{suffix}]', fontweight='bold')
            ax_t.set_xlabel('Time ($t$)')
            ax_t.set_ylabel('Temperature ($T$)')
            ax_t.grid(True); ax_t.legend()
            fig2.savefig(out_temp, format='pdf', bbox_inches='tight') 
            plt.close(fig2)
            msgs.append(f"Temperatura -> {out_temp.name}")

    except Exception as e:
        msgs.append(f"Error procesando energía/temp para {suffix}: {e}")

    # 2. TRAYECTORIAS (Si no es 3D)
    if "3D_" not in suffix and tray_path.exists():
        try:
            t_tray, particle_positions = [], {}
            with open(tray_path, 'r') as f:
                lines = f.readlines()
            i = 0
            while i < len(lines):
                line = lines[i].strip()
                if not line: i += 1; continue
                N_part = int(line); i += 1
                curr_t = float(lines[i].strip().split('=')[1]); t_tray.append(curr_t); i += 1
                for _ in range(N_part):
                    p = lines[i].split(); p_id = int(p[0])
                    if p_id not in particle_positions: particle_positions[p_id] = {'x': [], 'y': []}
                    particle_positions[p_id]['x'].append(float(p[1]))
                    particle_positions[p_id]['y'].append(float(p[2]))
                    i += 1
            
            fig3, ax3 = plt.subplots(figsize=(7, 4))
            if "2D_" in suffix:
                for p_id, pos in particle_positions.items():
                    r = [math.sqrt(x**2 + y**2) for x, y in zip(pos['x'], pos['y'])]
                    ax3.plot(t_tray, r, '.', markersize=1.5)
                ax3.set_ylabel('Radius ($r$)')
            else:
                for p_id, pos in particle_positions.items():
                    ax3.plot(t_tray, pos['x'])
                ax3.set_ylabel('Position ($x$)')
                
            ax3.set_title(f'Trajectories [{suffix}]', fontweight='bold')
            ax3.set_xlabel('Time ($t$)')
            ax3.grid(True)
            
            fig3.savefig(out_tray, format='pdf', bbox_inches='tight') 
            plt.close(fig3)
            msgs.append(f"Trayectoria -> {out_tray.name}")
        except Exception:
            pass

    return "\n".join(msgs)

# ==========================================
# 2. INTERFAZ Y MANEJO DEL PROGRAMA PRINCIPAL
# ==========================================
def get_latest_suffix():
    results_dir = Path("results")
    obs_files = list(results_dir.glob("obs_*.dat"))
    if not obs_files:
        return None
    latest_file = max(obs_files, key=os.path.getmtime)
    return latest_file.name[4:-4] 

def main():
    print("=" * 50)
    print(" 📊 STATISTICAL MDS ENGINE - PLOTTER MASTER ")
    print("=" * 50)
    print("1. Default Rápido (Observables + RDF de la última simulación)")
    print("2. Graficar TODAS las corridas (Multiprocesamiento)")
    print("3. Analizar Benchmark (N vs Tiempo)")
    print("4. Seleccionar simulación manualmente")
    print("=" * 50)
    
    choice = input("Elige una opción (1-4): ").strip()
    start_time = time.time()

    with concurrent.futures.ProcessPoolExecutor() as executor:
        futures = []

        if choice == '1':
            latest = get_latest_suffix()
            if not latest:
                print("No hay datos en results/")
                return
            print(f"\n[+] Lanzando procesos para la simulación: {latest}")
            futures.append(executor.submit(plot_observables, latest))
            
            ens = "NVT" if "NVT" in latest else "NVE"
            futures.append(executor.submit(plot_rdf, ens))

        elif choice == '2':
            results_dir = Path("results")
            obs_files = list(results_dir.glob("obs_*.dat"))
            print(f"\n[+] Se encontraron {len(obs_files)} archivos de observables. ¡Enviando a los núcleos!")
            for obs in obs_files:
                suffix = obs.name[4:-4]
                futures.append(executor.submit(plot_observables, suffix))
            
            futures.append(executor.submit(plot_rdf, "NVE"))
            futures.append(executor.submit(plot_rdf, "NVT"))
            futures.append(executor.submit(plot_benchmark))

        elif choice == '3':
            futures.append(executor.submit(plot_benchmark))

        elif choice == '4':
            results_dir = Path("results")
            obs_files = list(results_dir.glob("obs_*.dat"))
            if not obs_files:
                print("No hay datos.")
                return
            print("\nSimulaciones disponibles:")
            for i, obs in enumerate(obs_files):
                print(f"{i+1}. {obs.name[4:-4]}")
            sel = int(input("\nIngresa el número: ").strip()) - 1
            if 0 <= sel < len(obs_files):
                futures.append(executor.submit(plot_observables, obs_files[sel].name[4:-4]))
            else:
                print("Selección inválida.")
                return

        print("\nGenerando PDFs en paralelo... Por favor espera...")
        for f in concurrent.futures.as_completed(futures):
            print(f.result())

    end_time = time.time()
    print(f"\n✅ ¡Todo listo! Tiempo total de ploteo: {end_time - start_time:.2f} segundos.")

if __name__ == '__main__':
    main()