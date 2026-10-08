# Statistical-MDS-Engine
*(Statistical Molecular Dynamics Simulation Engine)*

> **Acknowledgments & Origin:** This project is a derivative engine heavily inspired by the original **PRISM** *(Parallelization & Research of interacting Research Model)* project. It has been adapted and restructured to offer native and full cross-platform compatibility (Linux and Windows) and to serve as a foundation for new applications in statistical mechanics.

A versatile, from-scratch Molecular Dynamics engine written in C++17. This project simulates particle interactions across 1D, 2D, and 3D spaces, allowing for the exploration of thermodynamic observables, phase transitions, and statistical mechanics principles.

Soon, the engine will expand to include **new advanced statistical analysis applications**, going beyond classical molecular dynamics.

## 🚀 Core Features

* **Cross-Platform Compatibility (Linux & Windows):** The source code has been rewritten using `std::filesystem` in C++ and `pathlib` in Python. This ensures that folder creation, as well as reading and writing heavy data files, works flawlessly regardless of the operating system, resolving classic path conflicts (`/` vs `\`).
* **Multi-Dimensional Support:** Seamlessly run simulations in 1D, 2D, or 3D. The mathematical engine dynamically adapts degrees of freedom and neutralizes unused axes to prevent division-by-zero (NaN) artifacts in boundary calculations.
* **Density-Driven Box Scaling:** The simulation box volume automatically scales based on the specified particle count (N) and target density (rho), ensuring rigorous statistical mechanics conditions.
* **Flexible Initializations:** 
  * **Uniform:** Particles are scattered uniformly across the entire simulation box.
  * **Spherical/Circular:** Particles are uniformly packed into a centered sphere or circle (using *rejection sampling* and proper radial weighting) to study expansion and equilibrium.
* **Smart Output Management:** Automated creation of a `results/` directory. Output files are dynamically named using the simulation parameters (e.g., `obs_3D_NVT_UNI_N1000_rho0.250.dat`) to prevent accidental overwriting and keep data organized.

## ⚛️ Physics & Thermodynamics

* **Ensemble Control:** Toggle effortlessly between Microcanonical (NVE) and Canonical (NVT) ensembles.
* **Andersen Thermostat:** Features a stochastic Andersen Thermostat to simulate coupling to a phantom heat bath, maintaining a constant target temperature via randomized Maxwell-Boltzmann velocity reassignments.
* **Normalized Observables:** Tracks Total Energy, Kinetic Energy, and Potential Energy *per particle* (E/N, K/N, U/N), alongside the instantaneous Temperature (T).
* **Structural Analysis (3D):** The C++ engine natively calculates the **Radial Distribution Function / Pair Correlation Function** $g(r)$ during the simulation using the Minimum Image Convention. To maximize computational efficiency, the $O(N^2)$ histogram sampling is performed strategically every $X$ steps rather than every integration step.
* **Reduced Units:** Operates in Lennard-Jones/reduced units ($k_B = 1, m = 1$), making the engine universally applicable to different atomic species by applying the proper scaling factors post-simulation.

## 🛠️ Compilation & Usage (Cross-Platform)

The project requires a compiler that supports **C++17** and utilizes a universal `Makefile` that automatically detects the operating system.

### 🐧 For Linux / macOS users
```bash
# 1. Compile the code
make

# 2. Run the simulation
./sim3D.x
```

### 🪟 For Windows users (PowerShell)
You need to have the GCC compiler installed via [MSYS2](https://www.msys2.org/) (packages `mingw-w64-ucrt-x86_64-gcc` and `mingw-w64-ucrt-x86_64-make`). Make sure the MSYS2 `bin` folder is added to your system's Environment Variables.

```powershell
# 1. Compile the code (ensure mingw32-make is aliased to 'make')
make

# 2. Run the simulation
.\sim3D.exe
```

## 📊 Data Visualization

The repository includes cross-platform Python scripts (powered by `pathlib`) that utilize `numpy` and `matplotlib` to parse the `.dat` output files and generate production-ready plots.

**Prerequisites:**
```bash
pip install numpy matplotlib
```

* **`OBSplotter.py` (Observables & Trajectories):** Generates a comprehensive $2\times2$ panel analyzing instantaneous energies ($E_{kin}$, $E_{pot}$, $E_{tot}$), accumulated averages ($\langle E \rangle$), and percentage fluctuations ($\Delta E$) to verify algorithmic energy conservation. It also plots temperature evolution (if NVT is detected) and visualizes 1D trajectories or 2D orbital planes.
* **`Radialdisfunc.py` (Structural Comparison):** Scans the `results/` directory for multiple $g(r)$ files belonging to the same ensemble, automatically extracts the density ($\rho$) using Regex, and plots a consolidated, color-graded comparative graph where phase transitions (gas $\rightarrow$ liquid $\rightarrow$ solid) become visually evident.
* **`NvsCTplotter.py` (Performance Benchmarking):** Analyzes the computational efficiency by parsing benchmark data to compare execution times across NVE and NVT ensembles. It automatically performs a log-log linear regression to calculate and print the algorithm's empirical time complexity (e.g., $O(N^2)$) directly to the terminal.

All generated graphs are saved directly to the `results/` folder as high-resolution `.png` images.

---
*Future Roadmap: Upcoming updates will include new statistical mechanics applications and advanced data-driven physics analysis.*