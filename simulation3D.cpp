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
    // 1. USER INPUTS (INTERACTIVE MENU)
    // ==========================================
    int N;
    std::cout << "Enter the number of particles (Warning: Recommended between 250 and 550): ";
    std::cin >> N;
    if (N < 250 || N > 550) {
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

    // === definition of basic system conditions ===
    Dimension sim_dim = Dimension::D3;  
    bool static_initialization = false; 

    int current_dim = 3;
    if (sim_dim == Dimension::D1) current_dim = 1;
    else if (sim_dim == Dimension::D2) current_dim = 2;

    double tau = 0.1;               
    double dt = 0.001;
    int steps = 50000;

    // --- SAMPLING FREQUENCIES (OPTIMIZATION) ---
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
    // 2. DENSITY & TEMPERATURE NESTED LOOPS
    // ==========================================
    std::vector<double> densities = {0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.1, 1.2};
    std::vector<double> temperatures;

    if (use_thermostat) {
        // 10 temperaturas para el estudio termodinámico completo en NVT
        temperatures = {0.6, 0.8, 1.0, 1.2, 1.4, 1.6, 1.8, 2.0, 2.2, 2.4};
    } else {
        // NVE por defecto corre a una energía inicial dada por T=1.0
        temperatures = {1.0}; 
    }

    int total_runs = densities.size() * temperatures.size();
    std::cout << "\nStarting multi-parameter simulation suite (" << total_runs << " runs)...\n";
    std::cout << "=========================================================\n";

    for (double current_T : temperatures) {
        for (double rho : densities) {
            
            double T_initial = current_T;         
            double T_target = current_T;

            double L = 0.0;
            switch (sim_dim){
                case Dimension::D1: L = N / rho; break;
                case Dimension::D2: L = std::sqrt(N/rho); break;
                case Dimension::D3: L = std::pow(N / rho, 1.0/3.0); break;
            }

            double Lx = L, Ly = L, Lz = L;
            if(sim_dim == Dimension::D1){Ly = 1.0; Lz = 1.0;}
            if(sim_dim == Dimension::D2){ Lz = 1.0;}

            double x_min = -Lx/2.0, x_max = Lx/2.0;
            double y_min = -Ly/2.0, y_max = Ly/2.0;
            double z_min = -Lz/2.0, z_max = Lz/2.0;

            std::cout << "\n-> Simulating T = " << current_T << " | rho = " << rho << " | Box L = " << L << "...\n";

            std::vector<Particle3D> particles(N);
            std::mt19937 gen(42); 

            std::uniform_real_distribution<double> dist_x(-Lx/2.0, Lx/2.0);
            std::uniform_real_distribution<double> dist_y(-Ly/2.0, Ly/2.0);
            std::uniform_real_distribution<double> dist_z(-Lz/2.0, Lz/2.0);

            // ==========================================
            // SMART INITIALIZATION: Random vs Lattice
            // ==========================================
            if (rho < 0.7) {
                std::cout << "      [+] Low density detected. Using Random Insertion...\n";
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
                std::cout << "      [+] High density (>= 0.7) detected. Using Crystalline Lattice...\n";
                initializeLattice(particles, Lx, Ly, Lz, current_dim);
            }
            
            if (static_initialization) {
                for (auto& p : particles) {
                    p.vx = 0.0; p.vy = 0.0; p.vz = 0.0;
                }
            } else {
                initializeMaxwellBoltzmann(particles, T_initial, current_dim, gen);
            }

            std::stringstream ss;
            switch (sim_dim){
                case Dimension::D1: ss << "1D_"; break;
                case Dimension::D2: ss << "2D_"; break;
                case Dimension::D3: ss << "3D_"; break;
            }

            ss << (use_thermostat ? "NVT_" : "NVE_");
            ss << "N" << N << "_rho" << std::fixed << std::setprecision(3) << rho;  
            
            // Actualizado para reflejar la temperatura correcta en el archivo
            if (static_initialization) {
                ss << "_Tinit0";
            } else {
                ss << "_T" << std::fixed << std::setprecision(1) << current_T;
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
                bench_file << "# N Computing_Time(s) Steps Density Temp\n";
            }
            
            // Añadí la temperatura al benchmark para tener registro
            bench_file << N << " " << elapsed_seconds.count() << " " << steps << " " << rho << " " << current_T << "\n";

            bench_file.close();
            traj.close();
            obs.close();
        } // End of Density Loop
    } // End of Temperature Loop

    std::cout << "\nAll simulations completed successfully!\n";
    return 0;
}