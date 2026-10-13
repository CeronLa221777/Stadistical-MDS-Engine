import numpy as np
import math
import matplotlib.pyplot as plt 
import re
import os
from pathlib import Path
import concurrent.futures 
import time
from collections import defaultdict
from matplotlib.lines import Line2D

# Use a non-interactive backend by default in workers so multiprocessing
# doesn't crash trying to open multiple windows simultaneously.
plt.switch_backend('agg')

# ==========================================
# 1. PLOTTING FUNCTIONS (PARALLEL WORKERS)
# ==========================================

def plot_benchmark():
    results_dir = Path("results")
    input_path = results_dir / "benchmark_CT_vs_N.dat"
    output_path = results_dir / "benchmark_analysis_dual.pdf" 

    try:
        # CORRECCIÓN: usecols=(0,1) fuerza a leer solo N y Computing_Time
        data = np.loadtxt(input_path, comments="#", usecols=(0, 1))
    except FileNotFoundError:
        return f"[Warning] {input_path} not found. Skipping Benchmark."

    if len(data) < 2:
        return "[Warning] Not enough data in benchmark to plot."

    mid = len(data) // 2
    data_NVE, data_NVT = data[:mid], data[mid:]
    
    if len(data_NVE) == 0 or len(data_NVT) == 0:
        return "[Warning] Incomplete benchmark format. Skipping."

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
    
    import warnings
    with warnings.catch_warnings():
        warnings.simplefilter('ignore', np.RankWarning)
        m_NVE, b_NVE = np.polyfit(log_N_NVE, log_CT_NVE, 1)
        m_NVT, b_NVT = np.polyfit(log_N_NVT, log_CT_NVT, 1)

    ax2.plot(log_N_NVE, log_CT_NVE, 'bo', markersize=8)
    ax2.plot(log_N_NVT, log_CT_NVT, 'rs', markersize=8)
    ax2.plot(log_N_NVE, (m_NVE * log_N_NVE) + b_NVE, 'b--', lw=2, label=f'NVE O(N^{m_NVE:.2f})')
    ax2.plot(log_N_NVT, (m_NVT * log_N_NVT) + b_NVT, 'r--', lw=2, label=f'NVT O(N^{m_NVT:.2f})')
    ax2.set_xlabel('log10(N)'); ax2.set_ylabel('log10(CT)')
    ax2.set_title('Linearization: Log-Log Plot', fontweight='bold')
    ax2.legend(); ax2.grid(True)

    plt.tight_layout()
    fig.savefig(output_path, format='pdf', bbox_inches='tight') 
    plt.close(fig)
    return f"Success! Benchmark saved to {output_path.name}"

def plot_rdf(target_ensemble="NVE"):
    results_dir = Path("results")
    file_list = list(results_dir.glob(f"rdf_3D_{target_ensemble}_*.dat"))

    if not file_list:
        return f"[Warning] No RDF files found for {target_ensemble}."

    if target_ensemble == "NVE":
        # NVE: Single Plot
        output_path = results_dir / f"plot_rdf_comparison_{target_ensemble}.pdf" 
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
        ax.set_xlabel('Distance $r$'); ax.set_ylabel('$g(r)$')
        ax.set_xlim(0, max_r); ax.set_ylim(0, max_g * 1.1)
        ax.grid(True, linestyle=':', alpha=0.6); ax.legend(title="Density", loc='upper right', bbox_to_anchor=(1.25, 1))
        plt.tight_layout(); fig.savefig(output_path, format='pdf', bbox_inches='tight') 
        plt.close(fig)
        return f"Success! RDF ({target_ensemble}) saved to {output_path.name}"

    else:
        # NVT: Group by Isotherm (Temperature) to prevent overlapping mess
        files_by_T = defaultdict(list)
        for filepath in file_list:
            match_rho = re.search(r'rho([0-9]+\.[0-9]+)', filepath.name)
            match_T = re.search(r'_T([0-9]+\.[0-9]+)', filepath.name)
            if match_rho and match_T:
                rho = float(match_rho.group(1))
                T_val = float(match_T.group(1))
                files_by_T[T_val].append((rho, filepath))
        
        msgs = []
        for T_val, data_files in files_by_T.items():
            data_files.sort(key=lambda x: x[0])
            fig, ax = plt.subplots(figsize=(8, 5))
            ax.set_title(f'Radial Distribution Function $g(r)$\nEnsemble {target_ensemble} - Isotherm T = {T_val}', fontsize=14, fontweight='bold')
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
            ax.set_xlabel('Distance $r$'); ax.set_ylabel('$g(r)$')
            ax.set_xlim(0, max_r); ax.set_ylim(0, max_g * 1.1)
            ax.grid(True, linestyle=':', alpha=0.6)
            ax.legend(title="Density", loc='upper right', bbox_to_anchor=(1.25, 1), fontsize='small')
            
            out_path = results_dir / f"plot_rdf_comparison_{target_ensemble}_T{T_val}.pdf"
            plt.tight_layout(); fig.savefig(out_path, format='pdf', bbox_inches='tight'); plt.close(fig)
            msgs.append(f"Success! RDF NVT (T={T_val}) saved to {out_path.name}")
        return "\n".join(msgs)

def plot_observables(suffix):
    results_dir = Path("results")
    obs_path = results_dir / f"obs_{suffix}.dat"
    out_4panel = results_dir / f"plot_energy_{suffix}.pdf" 
    msgs = []
    
    try:
        t_obs, K, U, E, T, P = np.loadtxt(obs_path, delimiter=" ", skiprows=1, unpack=True)
        
        E_0 = E[0]
        delta_E_inst = 100.0 * (E - E_0) / np.abs(E_0) if E_0 != 0 else np.zeros_like(E)
        n_steps = np.arange(1, len(E) + 1)
        K_avg, U_avg, E_avg, T_avg = np.cumsum(K)/n_steps, np.cumsum(U)/n_steps, np.cumsum(E)/n_steps, np.cumsum(T)/n_steps
        P_avg = np.cumsum(P)/n_steps

        fig1, axs = plt.subplots(2, 2, figsize=(12, 8), gridspec_kw={'wspace': 0.25, 'hspace': 0.35})
        fig1.suptitle(f'Dynamic Observables [{suffix}]', fontweight='bold', fontsize=14)
        
        axs[0,0].plot(t_obs, K, 'k-', alpha=0.7, label=r'$E_{kin}$')
        axs[0,0].plot(t_obs, U, 'r-', alpha=0.7, label=r'$E_{pot}$')
        axs[0,0].plot(t_obs, E, 'g-', lw=2, label=r'$E_{tot}$')
        axs[0,0].set_title('Instantaneous Energy')
        axs[0,0].set_xlabel('Time ($t$)'); axs[0,0].set_ylabel('Energy')
        axs[0,0].grid(True); axs[0,0].legend()
        
        axs[0,1].plot(t_obs, P, color='orange', alpha=0.4, label='Inst. Pressure')
        axs[0,1].plot(t_obs, P_avg, color='darkorange', lw=2, label=r'$\langle P \rangle$')
        axs[0,1].set_title("Pressure Evolution")
        axs[0,1].set_xlabel('Time ($t$)'); axs[0,1].set_ylabel('Pressure ($P$)')
        axs[0,1].grid(True); axs[0,1].legend()
        
        axs[1,0].plot(t_obs, K_avg, 'k--', label=r'$\langle E_{kin} \rangle$')
        axs[1,0].plot(t_obs, U_avg, 'r--', label=r'$\langle E_{pot} \rangle$')
        axs[1,0].plot(t_obs, E_avg, 'g-', lw=2, label=r'$\langle E_{tot} \rangle$')
        axs[1,0].set_title('Cumulative Average Energy')
        axs[1,0].set_xlabel('Time ($t$)'); axs[1,0].set_ylabel(r'$\langle E \rangle$')
        axs[1,0].grid(True); axs[1,0].legend()
        
        delta_E_avg = 100.0 * (E_avg - E_0) / np.abs(E_0) if E_0 != 0 else np.zeros_like(E)
        axs[1,1].plot(t_obs, delta_E_avg, 'g-', lw=1.5)
        axs[1,1].axhline(0, color='k', lw=1)
        axs[1,1].set_title("Average Energy Fluctuation\n" + r"$\langle \Delta E \rangle = \frac{\langle E(t) \rangle - E_0}{|E_0|} \times 100$")
        axs[1,1].set_xlabel('Time ($t$)'); axs[1,1].set_ylabel(r'$\langle \Delta E \rangle$ (%)')
        axs[1,1].yaxis.tick_right(); axs[1,1].yaxis.set_label_position("right"); axs[1,1].grid(True)
        
        fig1.savefig(out_4panel, format='pdf', bbox_inches='tight') 
        plt.close(fig1)
        msgs.append(f"Dynamics -> {out_4panel.name}")

    except Exception as e:
        pass # Handle silently in batch mode

    return "\n".join(msgs)

# ==========================================
# THERMODYNAMIC STATE ANALYSIS & G(R) COMPARISON
# ==========================================
def plot_thermodynamics(target_ensemble="NVT"):
    results_dir = Path("results")
    obs_files = list(results_dir.glob(f"obs_3D_{target_ensemble}_*.dat"))
    
    if not obs_files:
        return f"[Warning] Not enough suite data to plot {target_ensemble} thermodynamics."

    records = []
    U_RCUT = 4.0 * (2.5**-12 - 2.5**-6)

    for obs_path in obs_files:
        match_rho = re.search(r'rho([0-9]+\.[0-9]+)', obs_path.name)
        if not match_rho: continue
        rho = float(match_rho.group(1))
        
        match_T = re.search(r'_T([0-9]+\.[0-9]+)', obs_path.name)
        if not match_T:
            match_T = re.search(r'_Tinit([0-9]+\.[0-9]+)', obs_path.name)
        T_target = float(match_T.group(1)) if match_T else 1.0

        match_N = re.search(r'N(\d+)_', obs_path.name)
        N_part = int(match_N.group(1)) if match_N else 300
        
        try:
            # 1. Dynamic Averages
            t, K, U, E, T_inst, P = np.loadtxt(obs_path, delimiter=" ", skiprows=1, unpack=True)
            cut = len(t) // 5
            if len(t) <= cut: continue
            
            K_eq, U_eq, T_eq, P_eq = K[cut:], U[cut:], T_inst[cut:], P[cut:]
            avg_T = np.mean(T_eq)
            
            if target_ensemble == "NVT":
                Cv = (np.mean(K_eq) / avg_T) + (N_part * np.var(U_eq) / (avg_T**2))
            else:
                var_K = np.var(K_eq)
                den = 1.0 - (N_part * var_K / (np.mean(K_eq) * avg_T))
                Cv = (np.mean(K_eq) / avg_T) / den if den > 0.01 else np.nan

            # 2. Static Integrals from g(r)
            rdf_name = obs_path.name.replace("obs_", "rdf_")
            rdf_path = results_dir / rdf_name
            
            P_gr_val, U_gr_val = np.nan, np.nan
            if rdf_path.exists():
                r, g_r = np.loadtxt(rdf_path, skiprows=1, unpack=True)
                mask = (r > 0.1) & (r <= 2.5) 
                r_m, g_m = r[mask], g_r[mask]
                
                u_r = 4.0 * (r_m**-12 - r_m**-6) - U_RCUT
                rf_r = 24.0 * (2.0 * r_m**-12 - r_m**-6)
                
                U_gr_val = 2.0 * np.pi * rho * np.trapz((r_m**2) * u_r * g_m, r_m)
                P_gr_val = (rho * avg_T) + ( (2.0/3.0) * np.pi * (rho**2) * np.trapz((r_m**2) * rf_r * g_m, r_m) )

            records.append({
                'rho': rho,
                'T_target': T_target,
                'P_dyn': np.mean(P_eq), 'U_dyn': np.mean(U_eq), 'Cv_dyn': Cv,
                'P_gr': P_gr_val, 'U_gr': U_gr_val
            })
            
        except Exception as e:
            pass

    if not records:
        return "[Warning] Could not extract averages."

    # --- ERROR CALCULATION (Normalized Mean Absolute Error) ---
    P_dyn_arr = np.array([r['P_dyn'] for r in records if not np.isnan(r['P_gr'])])
    P_gr_arr = np.array([r['P_gr'] for r in records if not np.isnan(r['P_gr'])])
    err_P = (np.mean(np.abs(P_dyn_arr - P_gr_arr)) / (np.max(P_dyn_arr) - np.min(P_dyn_arr))) * 100 if len(P_dyn_arr) > 1 else 0.0

    U_dyn_arr = np.array([r['U_dyn'] for r in records if not np.isnan(r['U_gr'])])
    U_gr_arr = np.array([r['U_gr'] for r in records if not np.isnan(r['U_gr'])])
    err_U = (np.mean(np.abs(U_dyn_arr - U_gr_arr)) / (np.max(U_dyn_arr) - np.min(U_dyn_arr))) * 100 if len(U_dyn_arr) > 1 else 0.0

    print(f"      [📊] {target_ensemble} Method Comparison Deviation -> P: {err_P:.2f}%, U/N: {err_U:.2f}%")
    # --------------------------------------------------------

    if target_ensemble == "NVE":
        records.sort(key=lambda x: x['rho'])
        densities = [r['rho'] for r in records]
        
        fig, (ax1, ax2, ax3) = plt.subplots(1, 3, figsize=(15, 5))
        fig.suptitle('Thermodynamic Equation of State - NVE Ensemble', fontsize=16, fontweight='bold')

        ax1.plot(densities, [r['P_dyn'] for r in records], 'o-', color='darkorange', lw=2, label=r'Dynamic $\langle P \rangle$ (Virial)')
        ax1.plot(densities, [r['P_gr'] for r in records], 's--', color='purple', alpha=0.7, label=r'Static $P$ from $g(r)$')
        ax1.set_xlabel(r'Density ($\rho$)'); ax1.set_ylabel('Pressure ($P$)')
        ax1.set_title(fr'Equation of State: P vs $\rho$'+'\n'+fr'(Mean Dev: {err_P:.2f}%)')
        ax1.grid(True); ax1.legend()

        ax2.plot(densities, [r['U_dyn'] for r in records], 'o-', color='crimson', lw=2, label=r'Dynamic $\langle U/N \rangle$')
        ax2.plot(densities, [r['U_gr'] for r in records], 's--', color='navy', alpha=0.7, label=r'Static $U/N$ from $g(r)$')
        ax2.set_xlabel(r'Density ($\rho$)'); ax2.set_ylabel('Internal Energy ($U/N$)')
        ax2.set_title(fr'Internal Energy vs $\rho$'+'\n'+fr'(Mean Dev: {err_U:.2f}%)')
        ax2.grid(True); ax2.legend()

        ax3.plot(densities, [r['Cv_dyn'] for r in records], 'o-', color='teal', lw=2, label=f'$C_v$ (NVE)')
        ax3.set_xlabel(r'Density ($\rho$)'); ax3.set_ylabel('Specific Heat ($C_v$)')
        ax3.set_title(r'Specific Heat vs $\rho$'); ax3.grid(True); ax3.legend()

        out_eos = results_dir / f"plot_thermo_EOS_NVE.pdf"
        plt.tight_layout(); fig.savefig(out_eos, format='pdf', bbox_inches='tight'); plt.close(fig)
        return f"Success! NVE Thermodynamics saved to {out_eos.name}"

    elif target_ensemble == "NVT":
        unique_Ts = sorted(list(set([r['T_target'] for r in records])))
        unique_rhos = sorted(list(set([r['rho'] for r in records])))

        # ----------------------------------------------------
        # 1. ISOTHERMS (Plotting vs Density for each Temperature)
        # ----------------------------------------------------
        fig1, axs1 = plt.subplots(1, 3, figsize=(16, 5))
        fig1.suptitle(f'Isothermal Equations of State (NVT Ensemble)\nGlobal Method Deviation -> P: {err_P:.2f}%, U/N: {err_U:.2f}%', fontsize=14, fontweight='bold')
        colors_T = plt.cm.jet(np.linspace(0, 1, len(unique_Ts)))

        for i, T_val in enumerate(unique_Ts):
            sub = sorted([r for r in records if r['T_target'] == T_val], key=lambda x: x['rho'])
            if not sub: continue
            r_vals = [s['rho'] for s in sub]
            
            # Gráfico Comparativo: Dinámico (Círculo Continuo) vs g(r) (X Punteada)
            axs1[0].plot(r_vals, [s['P_dyn'] for s in sub], 'o-', color=colors_T[i], markersize=4)
            axs1[0].plot(r_vals, [s['P_gr'] for s in sub], 'x--', color=colors_T[i], alpha=0.7)
            
            axs1[1].plot(r_vals, [s['U_dyn'] for s in sub], 'o-', color=colors_T[i], markersize=4)
            axs1[1].plot(r_vals, [s['U_gr'] for s in sub], 'x--', color=colors_T[i], alpha=0.7)
            
            axs1[2].plot(r_vals, [s['Cv_dyn'] for s in sub], 'o-', color=colors_T[i], label=f'T={T_val:.1f}')

        # Leyenda personalizada para explicar los marcadores
        custom_lines = [Line2D([0], [0], color='k', marker='o', linestyle='-', lw=1.5, label='Dynamic (Virial / Fluct)'),
                        Line2D([0], [0], color='k', marker='x', linestyle='--', lw=1.5, alpha=0.7, label='Static from $g(r)$')]

        axs1[0].set_xlabel(r'Density ($\rho$)'); axs1[0].set_ylabel('Pressure ($P$)'); axs1[0].set_title(r'Isotherms: P vs $\rho$')
        axs1[0].grid(True); axs1[0].legend(handles=custom_lines, loc='upper left', fontsize='small')
        
        axs1[1].set_xlabel(r'Density ($\rho$)'); axs1[1].set_ylabel('Internal Energy ($U/N$)'); axs1[1].set_title(r'Isotherms: U/N vs $\rho$')
        axs1[1].grid(True); axs1[1].legend(handles=custom_lines, loc='upper left', fontsize='small')
        
        axs1[2].set_xlabel(r'Density ($\rho$)'); axs1[2].set_ylabel('Specific Heat ($C_v$)'); axs1[2].set_title(r'Isotherms: $C_v$ vs $\rho$')
        axs1[2].grid(True); axs1[2].legend(fontsize='small', loc='upper right', bbox_to_anchor=(1.25, 1))

        out_iso_t = results_dir / f"plot_thermo_Isotherms_NVT.pdf"
        fig1.tight_layout(); fig1.savefig(out_iso_t, format='pdf', bbox_inches='tight'); plt.close(fig1)

        # ----------------------------------------------------
        # 2. ISOCHORES (Plotting vs Temperature for each Density)
        # ----------------------------------------------------
        fig2, axs2 = plt.subplots(1, 3, figsize=(16, 5))
        fig2.suptitle(f'Isochoric Thermodynamic Dependencies (NVT Ensemble)\nGlobal Method Deviation -> P: {err_P:.2f}%, U/N: {err_U:.2f}%', fontsize=14, fontweight='bold')
        colors_rho = plt.cm.viridis(np.linspace(0, 1, len(unique_rhos)))

        for i, rho_val in enumerate(unique_rhos):
            sub = sorted([r for r in records if r['rho'] == rho_val], key=lambda x: x['T_target'])
            if not sub: continue
            t_vals = [s['T_target'] for s in sub]
            
            axs2[0].plot(t_vals, [s['P_dyn'] for s in sub], 'o-', color=colors_rho[i], markersize=4)
            axs2[0].plot(t_vals, [s['P_gr'] for s in sub], 'x--', color=colors_rho[i], alpha=0.7)
            
            axs2[1].plot(t_vals, [s['U_dyn'] for s in sub], 'o-', color=colors_rho[i], markersize=4)
            axs2[1].plot(t_vals, [s['U_gr'] for s in sub], 'x--', color=colors_rho[i], alpha=0.7)
            
            axs2[2].plot(t_vals, [s['Cv_dyn'] for s in sub], 'o-', color=colors_rho[i], label=rf'$\rho$={rho_val:.2f}')

        axs2[0].set_xlabel(r'Temperature ($T$)'); axs2[0].set_ylabel('Pressure ($P$)'); axs2[0].set_title(r'Isochores: P vs $T$')
        axs2[0].grid(True); axs2[0].legend(handles=custom_lines, loc='upper left', fontsize='small')
        
        axs2[1].set_xlabel(r'Temperature ($T$)'); axs2[1].set_ylabel('Internal Energy ($U/N$)'); axs2[1].set_title(r'Isochores: U/N vs $T$')
        axs2[1].grid(True); axs2[1].legend(handles=custom_lines, loc='upper right', fontsize='small')
        
        axs2[2].set_xlabel(r'Temperature ($T$)'); axs2[2].set_ylabel('Specific Heat ($C_v$)'); axs2[2].set_title(r'Isochores: $C_v$ vs $T$')
        axs2[2].grid(True); axs2[2].legend(fontsize='small', loc='upper right', bbox_to_anchor=(1.25, 1))

        out_iso_r = results_dir / f"plot_thermo_Isochores_NVT.pdf"
        fig2.tight_layout(); fig2.savefig(out_iso_r, format='pdf', bbox_inches='tight'); plt.close(fig2)

        return f"Success! NVT Isotherms saved to {out_iso_t.name}\nSuccess! NVT Isochores saved to {out_iso_r.name}"

# ==========================================
# 3. INTERFACE AND MAIN PROGRAM EXECUTION
# ==========================================
def get_latest_suffix():
    results_dir = Path("results")
    obs_files = list(results_dir.glob("obs_*.dat"))
    if not obs_files: return None
    return max(obs_files, key=os.path.getmtime).name[4:-4] 

def main():
    print("=" * 55)
    print(" 📊 STATISTICAL MDS ENGINE - PLOTTER MASTER ")
    print("=" * 55)
    print("1. Quick Default (Latest simulation)")
    print("2. Plot ALL and generate Thermodynamic Suite Analysis")
    print("3. Thermodynamics Only (P, U/N, Cv vs Density/Temp)")
    print("=" * 55)
    
    choice = input("Choose an option (1-3): ").strip()
    start_time = time.time()

    with concurrent.futures.ProcessPoolExecutor() as executor:
        futures = []

        if choice == '1':
            latest = get_latest_suffix()
            if latest:
                print(f"\n[+] Processing: {latest}")
                futures.append(executor.submit(plot_observables, latest))
                ens = "NVT" if "NVT" in latest else "NVE"
                futures.append(executor.submit(plot_rdf, ens))
            else:
                print("No data found in results/")

        elif choice == '2':
            results_dir = Path("results")
            obs_files = list(results_dir.glob("obs_*.dat"))
            print(f"\n[+] Dispatching {len(obs_files)} files to parallel cores...")
            for obs in obs_files:
                futures.append(executor.submit(plot_observables, obs.name[4:-4]))
            
            futures.append(executor.submit(plot_rdf, "NVE"))
            futures.append(executor.submit(plot_rdf, "NVT"))
            futures.append(executor.submit(plot_thermodynamics, "NVE"))
            futures.append(executor.submit(plot_thermodynamics, "NVT"))
            futures.append(executor.submit(plot_benchmark))

        elif choice == '3':
            print("\n[+] Extracting state averages and calculating g(r) integrals...")
            futures.append(executor.submit(plot_thermodynamics, "NVE"))
            futures.append(executor.submit(plot_thermodynamics, "NVT"))

        print("Generating PDFs... Please wait...")
        for f in concurrent.futures.as_completed(futures):
            print(f.result())

    print(f"\n✅ All done! Total time: {time.time() - start_time:.2f} s.")

if __name__ == '__main__':
    main()