import numpy as np
import matplotlib.pyplot as plt
import re
from pathlib import Path # <-- Reemplaza a 'os' y 'glob'

# ==========================================
# 1. CONFIGURACIÓN DEL ANÁLISIS
# ==========================================
# Define qué ensamble quieres comparar ("NVE" o "NVT")
target_ensemble = "NVE"  

# Si usas NVT, define la temperatura objetivo para el título
T_target = 1.0  

# Definir el directorio usando pathlib
data_dir = Path("results")
output_path = data_dir / f"plot_rdf_comparison_{target_ensemble}.png"

# ==========================================
# 2. BÚSQUEDA Y LECTURA DE ARCHIVOS
# ==========================================
# Buscar todos los archivos de RDF usando el glob integrado en pathlib
file_list = list(data_dir.glob(f"rdf_3D_{target_ensemble}_*.dat"))

if not file_list:
    print(f"Error: No se encontraron archivos de correlación (rdf_3D_{target_ensemble}_*.dat) en '{data_dir}'.")
    exit()

# Lista para guardar tuplas (densidad, ruta_del_archivo)
data_files = []

for filepath in file_list:
    # filepath.name reemplaza a os.path.basename()
    filename = filepath.name
    
    # Extraer la densidad (rho) usando expresiones regulares
    match = re.search(r'rho([0-9]+\.[0-9]+)', filename)
    
    if match:
        density = float(match.group(1))
        data_files.append((density, filepath))
    else:
        print(f"[Aviso] No se pudo extraer la densidad del archivo: {filename}")

# Ordenar los archivos de menor a mayor densidad
data_files.sort(key=lambda x: x[0])

# ==========================================
# 3. CREACIÓN DE LA GRÁFICA COMPACTA
# ==========================================
fig, ax = plt.subplots(figsize=(8, 5))

# Generar un título dinámico según el ensamble
if target_ensemble == "NVT":
    title_str = f'Radial Distribution Function $g(r)$\nEnsemble {target_ensemble} ($T = {T_target}$)'
else:
    title_str = f'Radial Distribution Function $g(r)$\nEnsemble {target_ensemble} (Constant Energy)'

ax.set_title(title_str, fontsize=14, fontweight='bold')

# Usar un colormap para que los colores transicionen suavemente según la densidad
colors = plt.cm.plasma(np.linspace(0, 0.85, len(data_files)))

# Iterar sobre los archivos ordenados y graficar
for i, (density, filepath) in enumerate(data_files):
    try:
        # np.loadtxt lee el objeto Path sin problemas
        r_rdf, g_r = np.loadtxt(filepath, skiprows=1, unpack=True)
        ax.plot(r_rdf, g_r, '-', color=colors[i], linewidth=2, label=rf'$\rho = {density:.3f}$')
    except Exception as e:
        print(f"Error leyendo {filepath.name}: {e}")

# Línea base del gas ideal
ax.axhline(1.0, color='gray', linestyle='--', alpha=0.7, label='Ideal Gas ($g(r)=1$)')

# Ajustes estéticos y de layout
ax.set_xlabel('Distance $r$', fontsize=12)
ax.set_ylabel('$g(r)$', fontsize=12)
ax.set_xlim(0, max(r_rdf) if 'r_rdf' in locals() else 10) 
ax.set_ylim(0, max(g_r) * 1.1 if 'g_r' in locals() else 3)
ax.grid(True, linestyle=':', alpha=0.6)

# Leyenda compacta
ax.legend(title="System Density", loc='lower right', fontsize='small', title_fontsize='medium', framealpha=0.9)

# Guardar y mostrar (savefig acepta Path)
plt.tight_layout()
plt.savefig(output_path, dpi=300, bbox_inches='tight')
print(f"¡Gráfica comparativa guardada exitosamente en: {output_path}!")

plt.show()