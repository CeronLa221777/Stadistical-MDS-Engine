import numpy as np
import matplotlib.pyplot as plt

# 1. Rutas de los archivos
input_path = "results/benchmark_CT_vs_N.dat"
output_path = "results/benchmark_analysis_dual.png"

# 2. Cargar los datos (ignora el encabezado con '#')
data = np.loadtxt(input_path, comments="#")

# --- SEPARAR LOS DATOS ---
# Asumiendo que hiciste 6 corridas de NVE primero y 6 de NVT después.
# Si los hiciste en otro orden, puedes ajustar los índices [0:6] y [6:12]
data_NVE = data[0:6]
data_NVT = data[6:12]

# Extraer N y Tiempo de Cómputo (CT) para NVE
N_NVE = data_NVE[:, 0]
CT_NVE = data_NVE[:, 1]

# Extraer N y Tiempo de Cómputo (CT) para NVT
N_NVT = data_NVT[:, 0]
CT_NVT = data_NVT[:, 1]

# 3. Crear la figura con dos subgráficas
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(13, 5))

# ==========================================
# GRÁFICA 1: Datos Crudos (Comportamiento Parabólico)
# ==========================================
ax1.plot(N_NVE, CT_NVE, marker='o', linestyle='-', color='blue', linewidth=2, label='NVE (No Thermostat)')
ax1.plot(N_NVT, CT_NVT, marker='s', linestyle='-', color='red', linewidth=2, label='NVT (Andersen Thermostat)')

ax1.set_xlabel('Number of Particles (N)', fontsize=12)
ax1.set_ylabel('Computing Time (s)', fontsize=12)
ax1.set_title('Raw Performance: N vs CT', fontsize=14, fontweight='bold')
ax1.legend()
ax1.grid(True)

# ==========================================
# GRÁFICA 2: Linealización Log-Log (Comportamiento Lineal)
# ==========================================
# Transformación Logarítmica
log_N_NVE, log_CT_NVE = np.log10(N_NVE), np.log10(CT_NVE)
log_N_NVT, log_CT_NVT = np.log10(N_NVT), np.log10(CT_NVT)

# Regresión lineal para NVE
m_NVE, b_NVE = np.polyfit(log_N_NVE, log_CT_NVE, 1)
fit_NVE = (m_NVE * log_N_NVE) + b_NVE

# Regresión lineal para NVT
m_NVT, b_NVT = np.polyfit(log_N_NVT, log_CT_NVT, 1)
fit_NVT = (m_NVT * log_N_NVT) + b_NVT

# Graficar puntos
ax2.plot(log_N_NVE, log_CT_NVE, marker='o', linestyle='', color='blue', markersize=8)
ax2.plot(log_N_NVT, log_CT_NVT, marker='s', linestyle='', color='red', markersize=8)

# Graficar líneas de ajuste
ax2.plot(log_N_NVE, fit_NVE, linestyle='--', color='darkblue', linewidth=2, 
         label=f'NVE Fit: Slope = {m_NVE:.2f}')
ax2.plot(log_N_NVT, fit_NVT, linestyle='--', color='darkred', linewidth=2, 
         label=f'NVT Fit: Slope = {m_NVT:.2f}')

ax2.set_xlabel('log10(N)', fontsize=12)
ax2.set_ylabel('log10(Computing Time)', fontsize=12)
ax2.set_title('Linearization: Log-Log Plot', fontsize=14, fontweight='bold')
ax2.legend()
ax2.grid(True)

# 4. Guardar y mostrar
plt.tight_layout()
plt.savefig(output_path, dpi=300, bbox_inches='tight')
plt.show()

# Imprimir los resultados en la terminal
print("-" * 40)
print("RESULTADOS DEL BENCHMARK")
print("-" * 40)
print(f"Complejidad Ensamble NVE: O(N^{m_NVE:.3f})")
print(f"Complejidad Ensamble NVT: O(N^{m_NVT:.3f})")
print("-" * 40)