#include <iostream>
#include <sstream>
#include <string>
#include <iomanip>
#include <cmath>
#include <vector>
#include <fstream>
#include <random>
#include <chrono>
#include <filesystem>
#include "verlet.hpp"
#include "observables.hpp"

enum class Dimension {D1, D2, D3};      

int main() {
    // ==========================================
    // 1. USER INPUTS & INTERACTIVE MENU
    // ==========================================
    int base_N;
    std::cout << "Enter the base number of particles (e.g., 300): ";
    std::cin >> base_N;
    if (base_N < 250 || base_N > 550) {
        std::cout << "[Notice] Operating outside the recommended range (250-550).\n";
    }

    char ens_choice;
    std::cout << "Use NVT ensemble with thermostat? (y/n, 'n' defaults to NVE): ";
    std::cin >> ens_choice;
    bool use_thermostat = (ens_choice == 'y' || ens_choice == 'Y');

    char bc_choice;
    std::cout << "Use Periodic Boundary Conditions? (y/n, 'n' defaults to Reflective): ";
    std::cin >> bc_choice;
    bool periodicB = (bc_choice == 'y' || bc_choice == 'Y');
    bool reflectiveB = !periodicB;

    char cv_choice;
    std::cout << "Keep Volume constant and vary N across densities? (y/n): ";
    std::cin >> cv_choice;
    bool constant_volume = (cv_choice == 'y' || cv_choice == 'Y');
    
    // We standardize the fixed volume using the base N at a reference density of 0.7
    double fixed_V = base_N / 0.7; 
    if (constant_volume) {
        std::cout << "  [+] Fixed Volume set to " << fixed_V << "\n";
    }

    bool iso_energetic = false;
    if (!use_thermostat) {
        char iso_choice;
        std::cout << "Initialize NVE with precise Target Total Energy (E/N) instead of T_init? (y/n): ";
        std::cin >> iso_choice;
        iso_energetic = (iso_choice == 'y' || iso_choice == 'Y');
    }

    int run_mode;
    std::cout << "\nSelect Sweep Resolution Mode:\n";
    std::cout << "1) Isochore Focus  (Cv vs T/E)   -> 5 Densities, 15 Temp/Energy points\n";
    std::cout << "2) Isotherm Focus  (Cv vs rho)   -> 15 Densities, 5 Temp/Energy points\n";
    std::cout << "3) Single Run      (1 specific point)\n";
    std::cout << "Choice: ";
    std::cin >> run_mode;

    // === definition of basic system conditions ===
    Dimension sim_dim = Dimension::D3;  
    bool static_initialization = false; 

    int current_dim = 3;
    if (sim_dim == Dimension::D1) current_dim = 1;
    else if (sim_dim == Dimension::D2) current_dim = 2;

    double tau = 0.1;               
    double dt = 0.001;
    int steps = 50000;

    int sample_interval_obs = 100;    
    int sample_interval_traj = 100;   
    int sample_interval_rdf = 100;    

    double d_f = 3.0;
    if(sim_dim == Dimension::D1) d_f = 1.0;
    else if(sim_dim == Dimension::D2) d_f = 2.0;

    namespace fs = std::filesystem;
    fs::path out_dir("results"); 
    if(!fs::exists(out_dir)){
        fs::create_directory(out_dir);
    }

    // ==========================================
    // 2. DYNAMIC DENSITY & TARGET ARRAYS
    // ==========================================
    std::vector<double> densities;
    std::vector<double> targets; 

    if (run_mode == 1) {
        // ISOCHORE FOCUS: 5 Curves (Densities), 15 X-axis points (Temperatures/Energies)
        densities = {0.2, 0.5, 0.8, 1.0, 1.2};
        if (iso_energetic) {
            targets = {-3.0, -2.5, -2.0, -1.5, -1.0, -0.5, 0.0, 0.5, 1.0, 1.5, 2.0, 2.5, 3.0, 3.5, 4.0};
        } else {
            targets = {0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.1, 1.2, 1.3, 1.4, 1.5, 1.7, 1.9, 2.1, 2.4};
        }
    } else if (run_mode == 2) {
        // ISOTHERM FOCUS: 15 X-axis points (Densities), 5 Curves (Temperatures/Energies)
        densities = {0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.1, 1.2, 1.3, 1.4, 1.5};
        if (iso_energetic) {
            targets = {-2.0, -1.0, 0.0, 1.0, 2.0};
        } else {
            targets = {0.8, 1.0, 1.2, 1.5, 2.0};
        }
    } else {
        // SINGLE RUN
        double d_val, t_val;
        std::cout << "Enter exact Density (e.g., 0.7): ";
        std::cin >> d_val;
        std::cout << (iso_energetic ? "Enter exact Target E/N (e.g., -1.5): " : "Enter exact Target T (e.g., 1.0): ");
        std::cin >> t_val;
        densities = {d_val};
        targets = {t_val};
    }

    int total_runs = densities.size() * targets.size();
    std::cout << "\nStarting multi-parameter simulation suite (" << total_runs << " runs)...\n";
    std::cout << "=========================================================\n";

    // ==========================================
    // 3. MAIN SIMULATION LOOPS
    // ==========================================
    for (double current_target : targets) {
        for (double rho : densities) {
            
            double L = 0.0;
            int N = base_N;

            // Geometry calculation based on Constant Volume mode
            if (constant_volume) {
                if (sim_dim == Dimension::D3) {
                    L = std::cbrt(fixed_V);
                } else if (sim_dim == Dimension::D2) {
                    L = std::sqrt(fixed_V);
                } else {
                    L = fixed_V;
                }
                N = std::round(rho * fixed_V);
                if (N < 2) N = 2; // Safeguard for extremely low densities
            } else {
                switch (sim_dim){
                    case Dimension::D1: L = N / rho; break;
                    case Dimension::D2: L = std::sqrt(N/rho); break;
                    case Dimension::D3: L = std::pow(N / rho, 1.0/3.0); break;
                }
            }

            double Lx = L, Ly = L, Lz = L;
            if(sim_dim == Dimension::D1){Ly = 1.0; Lz = 1.0;}
            if(sim_dim == Dimension::D2){ Lz = 1.0;}

            double x_min = -Lx/2.0, x_max = Lx/2.0;
            double y_min = -Ly/2.0, y_max = Ly/2.0;
            double z_min = -Lz/2.0, z_max = Lz/2.0;

            std::cout << "\n-> Simulating " << (iso_energetic ? "E/N" : "T") << " = " << current_target 
                      << " | rho = " << rho << " | N = " << N << " | Box L = " << L << "...\n";

            std::vector<Particle3D> particles(N);
            std::mt19937 gen(42); 

            std::uniform_real_distribution<double> dist_x(-Lx/2.0, Lx/2.0);
            std::uniform_real_distribution<double> dist_y(-Ly/2.0, Ly/2.0);
            std::uniform_real_distribution<double> dist_z(-Lz/2.0, Lz/2.0);

            // ==========================================
            // SMART INITIALIZATION: Random vs Lattice
            // ==========================================
            if (use_thermostat && rho < 0.7) {
                std::cout << "      [+] Low density detected (NVT). Using Random Insertion...\n";
                for(int i = 0; i < N; i++){
                    double base_x, base_y, base_z;
                    switch(sim_dim){
                        case Dimension::D1:
                            do { base_x = dist_x(gen); } 
                            while(tooClose(particles, base_x, 0.0, 0.0, i, 1.0, periodicB, Lx, Ly, Lz));
                            particles[i].x = base_x; particles[i].y = 0.0; particles[i].z = 0.0;
                            break;
                        case Dimension::D2:
                            do { base_x = dist_x(gen); base_y = dist_y(gen); } 
                            while(tooClose(particles, base_x, base_y, 0.0, i, 1.0, periodicB, Lx, Ly, Lz));
                            particles[i].x = base_x; particles[i].y = base_y; particles[i].z = 0.0;
                            break;
                        case Dimension::D3:
                            do { base_x = dist_x(gen); base_y = dist_y(gen); base_z = dist_z(gen); } 
                            while(tooClose(particles, base_x, base_y, base_z, i, 1.0, periodicB, Lx, Ly, Lz));
                            particles[i].x = base_x; particles[i].y = base_y; particles[i].z = base_z;
                            break;
                    }
                }
            } else {
                std::cout << "      [+] Using Crystalline Lattice initialization for structural stability...\n";
                initializeLattice(particles, Lx, Ly, Lz, current_dim);
            }
            
            // ==========================================
            // THERMAL / ENERGETIC VELOCITY SETUP
            // ==========================================
            double T_initial = current_target;         
            double T_target = current_target; // Used only if NVT

            if (iso_energetic && !use_thermostat) {
                // Calculate lattice potential energy to find required kinetic energy
                double U_init = potentialEnergy3D(particles, periodicB, Lx, Ly, Lz);
                double u_init = U_init / N;
                double k_req = current_target - u_init; // K/N = E/N - U/N
                
                if (k_req <= 0.0) {
                    std::cout << "      [!] Target E/N (" << current_target << ") is mathematically impossible "
                              << "from lattice state (U/N = " << u_init << "). Skipping.\n";
                    continue; // Skip this simulation run entirely
                }
                
                // T = (2 * K) / (d_f * N). Since k_req is K/N:
                T_initial = (2.0 * k_req) / d_f;
                std::cout << "      [+] Iso-Energetic Target Matched: U/N=" << u_init 
                          << " -> Set T_init=" << T_initial << " to reach E/N=" << current_target << "\n";
            }

            if (static_initialization) {
                for (auto& p : particles) {
                    p.vx = 0.0; p.vy = 0.0; p.vz = 0.0;
                }
            } else {
                initializeMaxwellBoltzmann(particles, T_initial, current_dim, gen);
            }

            // ==========================================
            // FILE NAMING CONVENTION
            // ==========================================
            std::stringstream ss;
            switch (sim_dim){
                case Dimension::D1: ss << "1D_"; break;
                case Dimension::D2: ss << "2D_"; break;
                case Dimension::D3: ss << "3D_"; break;
            }

            ss << (use_thermostat ? "NVT_" : "NVE_");
            ss << "N" << N << "_rho" << std::fixed << std::setprecision(3) << rho;  
            
            if (static_initialization) {
                ss << "_Tinit0";
            } else if (iso_energetic) {
                ss << "_Einit" << std::fixed << std::setprecision(1) << current_target;
            } else {
                ss << "_T" << std::fixed << std::setprecision(1) << current_target;
            }         
            
            ss << (periodicB ? "_period" : "_box");                               

            std::string suffix = ss.str();
            fs::path traj_filename = out_dir / ("tray_" + suffix + ".dat");
            fs::path obs_filename  = out_dir / ("obs_" + suffix + ".dat");
            fs::path rdf_filename  = out_dir / ("rdf_" + suffix + ".dat"); 

            std::ofstream traj(traj_filename); 
            std::ofstream obs(obs_filename); 

            double dr_rdf = 0.05;                       
            double max_dist_rdf = L / 2.0;              
            int num_bins_rdf = static_cast<int>(max_dist_rdf / dr_rdf);
            std::vector<double> rdf_hist(num_bins_rdf, 0.0);
            int rdf_snapshots = 0;             

            obs << "# t K/N U/N E/N T P\n";

            auto start_time = std::chrono::high_resolution_clock::now();
            double step_Virial = 0.0;

            for(int i = 0; i < steps; i++){
                double t = i * dt;

                step_Virial = velocityVerlet3D(particles, dt, x_min, x_max, y_min, y_max, z_min, z_max, reflectiveB, periodicB, Lx, Ly, Lz);

                if (use_thermostat) {
                    applyBerendsenThermostat(particles, T_target, tau, dt, current_dim);
                }

                if (sim_dim == Dimension::D3 && (i % sample_interval_rdf == 0)) {
                    updateRDF3D(particles, rdf_hist, dr_rdf, periodicB, Lx, Ly, Lz);
                    rdf_snapshots++;
                }

                if (i % sample_interval_traj == 0) {
                    traj << particles.size() <<"\n";
                    traj << "#t = " << t << "\n";
                    for(size_t j = 0; j < particles.size(); j++){
                         traj << j << " " << particles[j].x << " " << particles[j].y << " " << particles[j].z << "\n";
                    }
                }

                if (i % sample_interval_obs == 0) {
                    double K_total = kineticEnergy3D(particles);
                    double U_total = potentialEnergy3D(particles, periodicB, Lx, Ly, Lz);
                    double E_total = K_total + U_total;

                    double P_inst = pressure3D(K_total, step_Virial, Lx, Ly, Lz);

                    double K_norm = K_total / N;
                    double U_norm = U_total / N;
                    double E_norm = E_total / N;
                    double T_inst = (2.0*K_norm) / d_f;

                    obs << t << " " << K_norm << " " << U_norm << " " << E_norm << " " << T_inst << " " << P_inst << "\n";
                }
            }

            auto end_time = std::chrono::high_resolution_clock::now();
            std::chrono::duration<double> elapsed_seconds = end_time - start_time;

            if (sim_dim == Dimension::D3 && rdf_snapshots > 0) {
                normalizeAndSaveRDF3D(rdf_hist, rdf_filename.string(), N, Lx, Ly, Lz, rdf_snapshots, dr_rdf);
            }

            std::cout << "      [+] Done! Computing time: " << elapsed_seconds.count() << " seconds.\n";

            fs::path bench_filename = out_dir / "benchmark_CT_vs_N.dat"; 
            std::ofstream bench_file(bench_filename, std::ios::app); 
            std::ifstream bench_check(bench_filename);
            bench_check.seekg(0, std::ios::end);
            if (bench_check.tellg() == 0) {
                bench_file << "# N Computing_Time(s) Steps Density Target\n";
            }
            
            bench_file << N << " " << elapsed_seconds.count() << " " << steps << " " << rho << " " << current_target << "\n";

            bench_file.close();
            traj.close();
            obs.close();
        } 
    } 

    std::cout << "\nAll simulations completed successfully!\n";
    return 0;
}