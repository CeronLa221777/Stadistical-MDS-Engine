#include "verlet.hpp"
#include <cmath>
#include <random>
#include <omp.h> 

constexpr double RCUT2 = 6.25; 

void computeAccelerations3D(const std::vector<Particle3D>& particles,
                            std::vector<double>& acc_x,
                            std::vector<double>& acc_y,
                            std::vector<double>& acc_z,
                            bool usePBounds,
                            double Lx, double Ly, double Lz)
{   
    int N = particles.size();
    double rcut2 = RCUT2;

    std::fill(acc_x.begin(), acc_x.end(), 0.0);
    std::fill(acc_y.begin(), acc_y.end(), 0.0);
    std::fill(acc_z.begin(), acc_z.end(), 0.0);

    
    #pragma omp parallel for schedule(dynamic)
    for(int i = 0; i < N; i++){
        for(int j = 0; j < N; j++){
            if (i == j) continue; 

            double dx = particles[i].x - particles[j].x;
            double dy = particles[i].y - particles[j].y;
            double dz = particles[i].z - particles[j].z;

            if(usePBounds){
                dx -= Lx * std::round(dx / Lx);
                dy -= Ly * std::round(dy / Ly);
                dz -= Lz * std::round(dz / Lz);
            }

            double r2 = dx*dx + dy*dy + dz*dz;

            if (r2 < rcut2){
                double r2_inv = 1.0 / r2;
                double r6_inv = r2_inv * r2_inv * r2_inv;
                double f_scalar = 24.0 * r6_inv * r2_inv * (2.0 * r6_inv - 1.0);

                acc_x[i] += f_scalar*dx;
                acc_y[i] += f_scalar*dy;
                acc_z[i] += f_scalar*dz;
            }
        }
    }
}

void velocityVerlet3D(std::vector<Particle3D>& particles, double dt,
                      double xmin, double xmax,
                      double ymin, double ymax,
                      double zmin, double zmax,
                      bool useBoundaries, bool usePBounds, 
                      double Lx, double Ly, double Lz)
{
    int N = particles.size();
    std::vector<double> acc_x(N, 0.0);
    std::vector<double> acc_y(N, 0.0);
    std::vector<double> acc_z(N, 0.0);

    computeAccelerations3D(particles, acc_x, acc_y, acc_z, usePBounds, Lx, Ly, Lz);         
    

    for(int i = 0; i < N; i++){
        particles[i].vx += 0.5 * acc_x[i] * dt;                 
        particles[i].vy += 0.5 * acc_y[i] * dt;
        particles[i].vz += 0.5 * acc_z[i] * dt;
        
        particles[i].x += particles[i].vx * dt;                 
        particles[i].y += particles[i].vy * dt;
        particles[i].z += particles[i].vz * dt;
    }

    if (useBoundaries){
        applyReflectiveBC3D(particles, xmin, xmax, ymin, ymax, zmin, zmax);     
    }
    if (usePBounds){
        applyPeriodicBoundary(particles, Lx, Ly, Lz);  
    }

    computeAccelerations3D(particles, acc_x, acc_y, acc_z, usePBounds, Lx, Ly, Lz);         
    

    for(int i = 0; i < N; i++){
        particles[i].vx += 0.5 * acc_x[i] * dt;                 
        particles[i].vy += 0.5 * acc_y[i] * dt;
        particles[i].vz += 0.5 * acc_z[i] * dt;
    }
}

void applyReflectiveBC3D(std::vector<Particle3D>& particles,
                         double xmin, double xmax,
                         double ymin, double ymax,
                         double zmin, double zmax)
{

    for (int i = 0; i < particles.size(); i++) {
        auto& p = particles[i];
        if (p.x < xmin) { p.x = 2.0 * xmin - p.x; p.vx = -p.vx; }
        if (p.x > xmax) { p.x = 2.0 * xmax - p.x; p.vx = -p.vx; }
        if (p.y < ymin) { p.y = 2.0 * ymin - p.y; p.vy = -p.vy; }
        if (p.y > ymax) { p.y = 2.0 * ymax - p.y; p.vy = -p.vy; }
        if (p.z < zmin) { p.z = 2.0 * zmin - p.z; p.vz = -p.vz; }
        if (p.z > zmax) { p.z = 2.0 * zmax - p.z; p.vz = -p.vz; }
    }
}

void applyPeriodicBoundary(std::vector<Particle3D>& particles,
                           double Lx, double Ly, double Lz)
{
    double half_Lx = Lx / 2.0;
    double half_Ly = Ly / 2.0;
    double half_Lz = Lz / 2.0;


    for (int i = 0; i < particles.size(); i++) {
        auto& p = particles[i];
        if (p.x >= half_Lx) p.x -= Lx;
        else if (p.x < -half_Lx) p.x += Lx;

        if (p.y >= half_Ly) p.y -= Ly;
        else if (p.y < -half_Ly) p.y += Ly;

        if (p.z >= half_Lz) p.z -= Lz;
        else if (p.z < -half_Lz) p.z += Lz;
    }
}

bool tooClose(const std::vector<Particle3D>& particles,
              double x, double y, double z,
              int current, double minDist, bool usePBounds,
              double Lx, double Ly, double Lz)
{
    double minDist2 = minDist * minDist;

    for(int j = 0; j < current; j++){
        double dx = x - particles[j].x;
        double dy = y - particles[j].y;
        double dz = z - particles[j].z;

        if(usePBounds){
            dx -= Lx * std::round(dx / Lx);
            dy -= Ly * std::round(dy / Ly);
            dz -= Lz * std::round(dz / Lz);
        }

        double r2 = dx*dx + dy*dy + dz*dz;

        if(r2 < minDist2) return true;
    }
    return false;
}

void initializeMaxwellBoltzmann(std::vector<Particle3D>& particles, 
                                double T_initial, int dim, std::mt19937& gen)
{
    std::normal_distribution<double> dist_vel(0.0, std::sqrt(T_initial));
    
    double sum_vx = 0.0, sum_vy = 0.0, sum_vz = 0.0;
    int N = particles.size();

    for (auto& p : particles) {
        p.vx = (dim >= 1) ? dist_vel(gen) : 0.0;
        p.vy = (dim >= 2) ? dist_vel(gen) : 0.0;
        p.vz = (dim == 3) ? dist_vel(gen) : 0.0;
        
        sum_vx += p.vx;
        sum_vy += p.vy;
        sum_vz += p.vz;
    }

    double vcm_x = sum_vx / N;
    double vcm_y = sum_vy / N;
    double vcm_z = sum_vz / N;


    for (int i = 0; i < N; i++) {
        if (dim >= 1) particles[i].vx -= vcm_x;
        if (dim >= 2) particles[i].vy -= vcm_y;
        if (dim == 3) particles[i].vz -= vcm_z;
    }
}

void applyBerendsenThermostat(std::vector<Particle3D>& particles,
                              double T_target, double tau, double dt, int dim)
{
    int N = particles.size();
    double K = 0.0;
    

    for (int i = 0; i < N; i++) {
        K += 0.5 * (particles[i].vx*particles[i].vx + particles[i].vy*particles[i].vy + particles[i].vz*particles[i].vz);
    }
    
    double T_current = (2.0 * K) / (N * dim);
    if (T_current < 1e-8) return; 

    double ratio = dt / tau;
    double lambda = std::sqrt(1.0 + ratio * ((T_target / T_current) - 1.0));

   
    for (int i = 0; i < N; i++) {
        particles[i].vx *= lambda;
        particles[i].vy *= lambda;
        particles[i].vz *= lambda;
    }
}