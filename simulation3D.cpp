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
    constexpr double PI = 3.14159265358979323846;
    
    int N = 300;                    
    double rho = 0.20;              
    double T_initial = 1.0;         
    
    // === definition of basic system conditions ===
    Dimension sim_dim = Dimension::D3;  
    bool static_initialization = false; // true: v=0 / false: Maxwell-Boltzmann distribution

    int current_dim = 3;
    if (sim_dim == Dimension::D1) current_dim = 1;
    else if (sim_dim == Dimension::D2) current_dim = 2;

    // Switches to choose boundary conditions and thermostat
    bool periodicB = true;
    bool reflectiveB = false;

    bool use_thermostat = false;     
    double T_target = 1.0;          
    double tau = 0.1;               

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

    std::cout << "Simulating N = " << N << " with density rho = " << rho << std::endl;
    std::cout << "Calculated box size L = " << L << std::endl;

    double dt = 0.001;
    int steps = 50000;

    // --- SAMPLING FREQUENCIES (OPTIMIZATION) ---
    int sample_interval_obs = 100;    // Calculate energy and temperature every 100 steps
    int sample_interval_traj = 100;   // Save positions to disk every 100 steps
    int sample_interval_rdf = 100;    // Take g(r) sample every 100 steps
    
    std::vector<Particle3D> particles(N);
    std::mt19937 gen(42);

    std::uniform_real_distribution<double> dist_x(-Lx/2.0, Lx/2.0);
    std::uniform_real_distribution<double> dist_y(-Ly/2.0, Ly/2.0);
    std::uniform_real_distribution<double> dist_z(-Lz/2.0, Lz/2.0);

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
    
    if (static_initialization) {
        for (auto& p : particles) {
            p.vx = 0.0; p.vy = 0.0; p.vz = 0.0;
        }
    } else {
        initializeMaxwellBoltzmann(particles, T_initial, current_dim, gen);
    }

    namespace fs = std::filesystem;
    fs::path out_dir("results"); 
    if(!fs::exists(out_dir)){
        fs::create_directory(out_dir);
    }

    std::stringstream ss;
    switch (sim_dim){
        case Dimension::D1: ss << "1D_"; break;
        case Dimension::D2: ss << "2D_"; break;
        case Dimension::D3: ss << "3D_"; break;
    }

    ss << (use_thermostat ? "NVT_" : "NVE_");
    ss << "N" << N << "_rho" << std::fixed << std::setprecision(3) << rho;  
    ss << (static_initialization ? "_Tinit0" : ("_Tinit" + std::to_string(T_initial).substr(0,3)));          
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

    double d_f = 3.0;
    if(sim_dim == Dimension::D1) d_f = 1.0;
    else if(sim_dim == Dimension::D2) d_f = 2.0;

    obs << "# t K/N U/N E/N T \n";

    auto start_time = std::chrono::high_resolution_clock::now();

    for(int i = 0; i < steps; i++){
        double t = i * dt;

        // 1. DYNAMICS
        velocityVerlet3D(particles, dt, x_min, x_max, y_min, y_max, z_min, z_max, reflectiveB, periodicB, Lx, Ly, Lz);

        if (use_thermostat) {
            applyBerendsenThermostat(particles, T_target, tau, dt, current_dim);
        }

        // 2. STRUCTURAL SAMPLING: Only every N steps
        if (sim_dim == Dimension::D3 && (i % sample_interval_rdf == 0)) {
            updateRDF3D(particles, rdf_hist, dr_rdf, periodicB, Lx, Ly, Lz);
            rdf_snapshots++;
        }

        // 3. TRAJECTORY SAMPLING (Optimization)
        if (i % sample_interval_traj == 0) {
            traj << particles.size() <<"\n";
            traj << "#t = " << t << "\n";
            for(size_t j = 0; j < particles.size(); j++){
                 traj << j << " " << particles[j].x << " " << particles[j].y << " " << particles[j].z << "\n";
            }
        }

        // 4. THERMODYNAMIC SAMPLING (CUTTING EXECUTION TIME IN HALF)
        if (i % sample_interval_obs == 0) {
            double K_total = kineticEnergy3D(particles);
            
            // Heavy call (O(N^2)) relegated to occur only 1 out of every 100 times
            double U_total = potentialEnergy3D(particles, periodicB, Lx, Ly, Lz);
            double E_total = K_total + U_total;

            double K_norm = K_total / N;
            double U_norm = U_total / N;
            double E_norm = E_total / N;
            double T_inst = (2.0*K_norm) / d_f;

            obs << t << " " << K_norm << " " << U_norm << " " << E_norm << " " << T_inst << "\n";
        }
    }

    auto end_time = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed_seconds = end_time - start_time;

    if (sim_dim == Dimension::D3 && rdf_snapshots > 0) {
        normalizeAndSaveRDF3D(rdf_hist, rdf_filename.string(), N, Lx, Ly, Lz, rdf_snapshots, dr_rdf);
    }

    std::cout << "Simulation finished. Computing time: " << elapsed_seconds.count() << " seconds.\n";

    fs::path bench_filename = out_dir / "benchmark_CT_vs_N.dat"; 
    std::ofstream bench_file(bench_filename, std::ios::app); 

    std::ifstream bench_check(bench_filename);
    bench_check.seekg(0, std::ios::end);
    if (bench_check.tellg() == 0) {
        bench_file << "# N Computing_Time(s) Steps Density\n";
    }
    
    bench_file << N << " " << elapsed_seconds.count() << " " << steps << " " << rho << "\n";

    bench_file.close();
    traj.close();
    obs.close();

    return 0;
}