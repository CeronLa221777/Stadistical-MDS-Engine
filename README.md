# PRISM: Parallelization & Research of interacting Research Model

A versatile, from-scratch Molecular Dynamics engine written in C++17. This project simulates particle interactions across 1D, 2D, and 3D spaces, allowing for the exploration of thermodynamic observables, phase transitions, and statistical mechanics principles.

## 🚀 Core Features

* **Multi-Dimensional Support:** Seamlessly run simulations in 1D, 2D, or 3D. The mathematical engine dynamically adapts degrees of freedom and neutralizes unused axes to prevent division-by-zero (NaN) artifacts in boundary calculations.
* **Density-Driven Box Scaling:** The simulation box volume automatically scales based on the specified particle count (N) and target density (rho), ensuring rigorous statistical mechanics conditions.
* **Flexible Initializations:** * **Uniform:** Particles are scattered uniformly across the entire simulation box.
  * **Spherical/Circular:** Particles are uniformly packed into a centered sphere or circle (using rejection sampling and proper radial weighting) to study expansion and equilibrium.
* **Smart Output Management:** Automated creation of a `results/` directory. Output files are dynamically named using the simulation parameters (e.g., `obs_3D_NVT_UNI_N1000_rho0.250.dat`) to prevent accidental overwriting and keep data organized.
* **Performance Benchmarking:** Built-in execution timer (using `<chrono>`) that appends N vs. Computing Time to a benchmark file, allowing for empirical testing of the algorithm's O(N^2) time complexity.

## ⚛️ Physics & Thermodynamics

* **Ensemble Control:** Toggle effortlessly between Microcanonical (NVE) and Canonical (NVT) ensembles.
* **Andersen Thermostat:** Features a stochastic Andersen Thermostat to simulate coupling to a phantom heat bath, maintaining a constant target temperature via randomized Maxwell-Boltzmann velocity reassignments.
* **Normalized Observables:** Tracks Total Energy, Kinetic Energy, and Potential Energy *per particle* (E/N, K/N, U/N), alongside the instantaneous Temperature (T).
* **Structural Analysis (3D):** The C++ engine natively calculates the **Radial Distribution Function / Pair Correlation Function** $g(r)$ during the simulation using the Minimum Image Convention. To maximize computational efficiency, the $O(N^2)$ histogram sampling is performed strategically every $X$ steps rather than every integration step.
* **Reduced Units:** Operates in Lennard-Jones/reduced units ($k_B = 1, m = 1$), making the engine universally applicable to different atomic species by applying the proper scaling factors post-simulation.

## 🛠️ Compilation & Usage

This project utilizes the `<filesystem>` library for file management, requiring a compiler that supports **C++17** or higher. This project also utilizes a `Makefile` for streamlined compilation, automatically applying flags (`-O3`) for maximum performance and modern file system management.
**To compile (Linux):**
```bash
make
```

## To run:
```bash
./sim3D.x     # On Linux/macOS
```

## 📊 Data Visualization
The repository includes Python scripts utilizing `numpy` and `matplotlib` to parse the output `.dat` files and generate production-ready plots of the system's thermodynamic evolution.

The visualization scripts automatically detect the ensemble type (NVE vs NVT) from the C++ generated filename and adjust titles and axis limits accordingly. All generated graphs are saved directly to the `results/` folder.

* **`OBSplotter.py` (Observables & Trajectories):** * Generates a comprehensive $2\times2$ panel analyzing instantaneous energies ($E_{kin}$, $E_{pot}$, $E_{tot}$), accumulated averages ($\langle E \rangle$), and percentage fluctuations ($\Delta E$) to verify algorithmic energy conservation.
  * Plots temperature evolution (if NVT is detected).
  * Automatically plots 1D trajectories or 2D orbital planes depending on the system's dimensionality.

* **`Radialdisfunc.py` (Structural Comparison):** * A specialized script that scans the `results/` directory for multiple $g(r)$ files belonging to the same ensemble.
  * Automatically extracts the density ($\rho$) from each filename using Regex.
  * Plots a consolidated, color-graded comparative graph of the Radial Distribution Functions across different densities, making phase transitions (gas $\rightarrow$ liquid $\rightarrow$ solid) visually evident.

* **`NvsCTplotter.py` (Performance Benchmarking):** * Analyzes the computational efficiency of the engine by parsing benchmark data to compare execution times across NVE and NVT ensembles.
  * Generates a dual-panel figure: the first plot visualizes the raw computing time versus the number of particles (demonstrating the expected polynomial growth), while the second plot displays a log-log linearization.
  * Automatically performs a linear regression on the logarithmic data to extract the slope, explicitly calculating and printing the algorithm's empirical time complexity (e.g., O(N^2)) directly to the terminal.

All generated graphs are saved directly to the `results/` folder as high-resolution `.png` images.
