#include "observables.hpp"
#include <cmath>
#include <fstream>
#include <iostream>
#include <algorithm> 
#include <omp.h> 

// 1. Cutoff Radius Adjustment: r_c = 2.5 sigma -> r_c^2 = 6.25
constexpr double RCUT2 = 6.25; 
constexpr double PI = 3.14159265358979323846;

// 2. Truncation and Shift: Precalculate U(r_c) to subtract
// U(2.5) = 4 * ( (1/2.5)^12 - (1/2.5)^6 ) = -0.016316891136
constexpr double U_RCUT = -0.016316891136;

//---------3D-------------------------------------------
double kineticEnergy3D(const std::vector<Particle3D>& particles)
{
    double K = 0.0;
    int N = particles.size();

    
    for(int i = 0; i < N; i++){
        K += 0.5 *( particles[i].vx * particles[i].vx + 
                    particles[i].vy * particles[i].vy + 
                    particles[i].vz * particles[i].vz ); // mass = 1
    }
    return K;
}

double potentialEnergy3D(const std::vector<Particle3D>& particles,
                         bool usePBounds,
                         double Lx, double Ly, double Lz)
{
    double U = 0.0;
    int N = particles.size();
    double rcut2 = RCUT2;

   
    #pragma omp parallel for schedule(dynamic) reduction(+:U)
    for(int i = 0; i < N; i++){
        for(int j = i + 1; j < N; j++){

            double dx = particles[i].x - particles[j].x;
            double dy = particles[i].y - particles[j].y;
            double dz = particles[i].z - particles[j].z;

            if(usePBounds){
                dx -= Lx * std::round(dx / Lx);
                dy -= Ly * std::round(dy / Ly);
                dz -= Lz * std::round(dz / Lz);
            }

            double r2 = dx*dx + dy*dy + dz*dz;

            if(r2 < rcut2){
                double r2_inv = 1.0 / r2;
                double r6_inv = r2_inv * r2_inv * r2_inv;
                double r12_inv = r6_inv * r6_inv;

                // Original interacting Lennard-Jones potential
                double u_lj = 4.0 * (r12_inv - r6_inv);
                
                // Truncated and shifted u(r) = U(r) - U(r_c)
                U += (u_lj - U_RCUT);
            }
        }
    }

    return U;
}

// ------------------------------------------------------------------
// FUNCTION 1: g(r) Histogram Collector
// ------------------------------------------------------------------
void updateRDF3D(const std::vector<Particle3D>& particles,
                 std::vector<double>& rdf_hist,
                 double dr,
                 bool usePBounds,
                 double Lx, double Ly, double Lz)
{
    int N = particles.size();
    
    // By minimum image convention, measuring beyond L/2 has no physical meaning
    double max_dist = std::min({Lx, Ly, Lz}) / 2.0;

    
    #pragma omp parallel for schedule(dynamic)
    for (int i = 0; i < N; i++) {
        for (int j = i + 1; j < N; j++) {
            
            double dx = particles[i].x - particles[j].x;
            double dy = particles[i].y - particles[j].y;
            double dz = particles[i].z - particles[j].z;

            // minimum image convention
            if (usePBounds) {
                dx -= Lx * std::round(dx / Lx);
                dy -= Ly * std::round(dy / Ly);
                dz -= Lz * std::round(dz / Lz);
            }

            double r = std::sqrt(dx*dx + dy*dy + dz*dz);

            // Only add if it is within the minimum measurable box
            if (r < max_dist) {
                size_t bin = static_cast<size_t>(r / dr);
                if (bin < rdf_hist.size()) {
                    // Vital atomic update to avoid thread collisions
                    #pragma omp atomic
                    rdf_hist[bin] += 2.0; 
                }
            }
        }
    }
}

// ------------------------------------------------------------------
// FUNCTION 2: g(r) Normalizer and Output
// ------------------------------------------------------------------
void normalizeAndSaveRDF3D(const std::vector<double>& rdf_hist,
                           const std::string& filename,
                           int N, double Lx, double Ly, double Lz,
                           int num_snapshots, double dr)
{
    std::ofstream out(filename);
    if (!out.is_open()) {
        std::cerr << "[ERROR] Could not create the file " << filename << std::endl;
        return;
    }

    double V = Lx * Ly * Lz;
    double rho = N / V; // Ideal density

    out << "# r g(r)\n";

    // This loop runs only for the number of 'bins' (e.g., 100 times at the very end). 
    // It is instantaneous on a single core.
    for (size_t bin = 0; bin < rdf_hist.size(); bin++) {
        double r_inner = bin * dr;
        double r_outer = (bin + 1) * dr;
        double r_center = r_inner + 0.5 * dr;

        // Spherical shell volume V = 4/3 * PI * (r_outer^3 - r_inner^3)
        double dV = (4.0 / 3.0) * PI * (r_outer * r_outer * r_outer - r_inner * r_inner * r_inner);
        
        // How many particles would be in this shell if it were an ideal gas without interactions?
        double ideal_particles = rho * dV;

        // Normalize the count using the total snapshots taken, N particles, and the ideal gas
        double g_r = 0.0;
        if (ideal_particles > 0.0) {
            g_r = rdf_hist[bin] / (N * num_snapshots * ideal_particles);
        }

        out << r_center << " " << g_r << "\n";
    }
    
    out.close();
}