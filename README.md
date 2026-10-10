# Statistical-MDS-Engine
*(Statistical Molecular Dynamics Simulation Engine)*

> **Acknowledgments & Origin:** This project is a derivative engine heavily inspired by the original **PRISM** *(Parallelization & Research of interacting Research Model)* project. It has been adapted and restructured to offer native and full cross-platform compatibility (Linux and Windows) and to serve as a foundation for new applications in statistical mechanics.

A versatile, highly-optimized Molecular Dynamics engine written in C++17. This project simulates particle interactions across 1D, 2D, and 3D spaces, allowing for the exploration of thermodynamic observables, phase transitions, and statistical mechanics principles.

Soon, the engine will expand to include **new advanced statistical analysis applications**, going beyond classical molecular dynamics.

## 🚀 Core Features

- **OpenMP Multi-Core Parallelization:** The engine utilizes `#pragma omp` directives to distribute heavy $O(N^2)$ workloads (force interactions and structural sampling) across all available CPU threads, drastically reducing computation time while avoiding thread-management overhead on lighter $O(N)$ routines.
- **Optimized Disk I/O:** Trajectories and thermodynamic observables are strategically sampled at specific intervals (e.g., every 100 steps) rather than every integration step. This completely eliminates disk-write bottlenecks and prevents massive file sizes without sacrificing statistical independence.
- **Cross-Platform Compatibility (Linux & Windows):** The source code has been rewritten using `std::filesystem` in C++ and `pathlib` in Python. This ensures that folder creation and data handling work flawlessly regardless of the operating system, resolving classic path conflicts (`/` vs `\`).
- **Multi-Dimensional Support:** Seamlessly run simulations in 1D, 2D, or 3D. The mathematical engine dynamically adapts degrees of freedom and neutralizes unused axes to prevent division-by-zero (NaN) artifacts.
- **Smart Output Management:** Automated creation of a `results/` directory. Output files are dynamically named using the simulation parameters (e.g., `obs_3D_NVT_N300_rho0.200_Tinit1.0_period.dat`) to prevent accidental overwriting and keep data organized.

## ⚛️ Physics & Thermodynamics

- **Shifted Lennard-Jones Potential:** Implements a rigorous Lennard-Jones potential truncated at $r_c = 2.5\sigma$. The potential energy is smoothly shifted ($U(r) - U(r_c)$) to neutralize boundary discontinuities, guaranteeing perfect total energy conservation in the microcanonical ensemble.
- **Ensemble Control:** Toggle effortlessly between Microcanonical (NVE) and Canonical (NVT) ensembles.
- **Maxwell-Boltzmann Initialization:** Particles are initialized uniformly in space, but their velocities are assigned from a strict Gaussian distribution based on a target initial temperature, instantly neutralizing center-of-mass drift.
- **Berendsen Thermostat:** Replaces standard stochastic collisions with a continuous velocity-rescaling algorithm to reach the target temperature smoothly, preserving the physical trajectory and temporal correlation of the particles.
- **Normalized Observables:** Tracks Total Energy, Kinetic Energy, and Potential Energy *per particle* (E/N, K/N, U/N), alongside the instantaneous Temperature (T).
- **Structural Analysis (3D):** Natively calculates the **Radial Distribution Function / Pair Correlation Function** $g(r)$ using the Minimum Image Convention. 

## 🛠️ Compilation & Usage (Cross-Platform)

The project requires a compiler that supports **C++17** and **OpenMP**. It utilizes a universal `Makefile` that automatically detects the operating system.

### 🐧 For Linux / macOS users
```bash
# 1. Compile the code (optimized with -O3 and -fopenmp)
make

# 2. Run the simulation
./sim3D.x

# 3. Generate analysis PDFs automatically
make plot
```

### 🪟 For Windows users (PowerShell)
You need to have the GCC compiler installed via [MSYS2](https://www.msys2.org/) (packages `mingw-w64-ucrt-x86_64-gcc` and `mingw-w64-ucrt-x86_64-make`). Make sure the MSYS2 `bin` folder is added to your system's Environment Variables.

```powershell
# 1. Compile the code (ensure mingw32-make is aliased to 'make')
make

# 2. Run the simulation
.\sim3D.exe

# 3. Generate analysis PDFs automatically
make plot
```

## 📊 Data Visualization (Multiprocessing)

The repository relies on a unified, highly efficient master script (`PlotterMaster.py`) powered by `concurrent.futures`, `pathlib`, `numpy`, and `matplotlib` to parse the `.dat` output files and generate production-ready vector graphics.

**Prerequisites:**
```bash
pip install numpy matplotlib
```

By simply running `make plot`, the engine launches an interactive Python menu that allows you to:
- **Graph the latest simulation automatically.**
- **Batch-process ALL simulation files simultaneously:** Utilizing actual CPU Multiprocessing to render dozens of graphs in parallel.
- **Extract Benchmarks:** Automatically analyze `benchmark_CT_vs_N.dat`, performing log-log linear regressions to print the empirical algorithm time complexity (e.g., $O(N^2)$) in the terminal.
- **Render High-Quality Vectors:** All visualizations (Thermodynamic Panels, $g(r)$ Comparisons, and Trajectories) are natively exported as `.pdf` files, ensuring infinite zoom resolution for academic papers and reports.

---
*Future Roadmap: Upcoming updates will include pressure calculation via the Virial Theorem, specific heat ($C_v$) extraction, and advanced data-driven physics analysis.*