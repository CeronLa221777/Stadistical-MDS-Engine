#include "verlet.hpp"
#include <cmath>
#include <random>

constexpr double RCUT2 = 16.0; //radio de cort definido arriba en el codigo



void computeAccelerations3D(const std::vector<Particle3D>& particles,
                            std::vector<double>& acc_x,
                            std::vector<double>& acc_y,
                            std::vector<double>& acc_z,
                            const std::vector<double>& k,
                            bool usePBounds,
                            double Lx, double Ly, double Lz)
{   int N = particles.size();
    double rcut2 = RCUT2;

    double kx = k[0];
    double ky = k[1];
    double kz = k[2];

    std::fill(acc_x.begin(), acc_x.end(), 0.0);
    std::fill(acc_y.begin(), acc_y.end(), 0.0);
    std::fill(acc_z.begin(), acc_z.end(), 0.0);

    // Trampa armónica
    for(int i = 0; i < N; i++){
        acc_x[i] = -kx * particles[i].x;
        acc_y[i] = -ky * particles[i].y;
        acc_z[i] = -kz * particles[i].z;
    }

    // Interacción por pares
    for(int i = 0; i < N; i++){
        for(int j = i + 1; j < N; j++){

            double dx = particles[i].x - particles[j].x;
            double dy = particles[i].y - particles[j].y;
            double dz = particles[i].z - particles[j].z;

            // minimum image convention
            if(usePBounds){
                dx -= Lx * std::round(dx / Lx);
                dy -= Ly * std::round(dy / Ly);
                dz -= Lz * std::round(dz / Lz);
            }

            double r2 = dx*dx + dy*dy + dz*dz;

            if (r2 < rcut2){

                double r4 = r2*r2;
                double r6 = r4*r2;
                double r14 = r6*r6*r2;

                double f_scalar = 12.0 / r14;

                acc_x[i] += f_scalar*dx;
                acc_y[i] += f_scalar*dy;
                acc_z[i] += f_scalar*dz;

                acc_x[j] -= f_scalar*dx;
                acc_y[j] -= f_scalar*dy;
                acc_z[j] -= f_scalar*dz;
            }
        }
    }
}

//implementacion velocity verlet para las diferentes dimensiones
void velocityVerlet3D(std::vector<Particle3D>& particles, double dt,
                        const std::vector<double>& k, double xmin, double xmax,
                                  double ymin, double ymax,
                                  double zmin, double zmax,
                                bool useBoundaries, bool usePBounds,double Lx, double Ly,double Lz)
{
    int N = particles.size();
    std::vector<double> acc_x(N, 0.0);
    std::vector<double> acc_y(N, 0.0);
    std::vector<double> acc_z(N, 0.0);

    computeAccelerations3D(particles, acc_x, acc_y, acc_z, k, usePBounds,Lx,Ly,Lz);         //Calculamos aceleración en configuración "inicial"
    for(int i = 0; i < N; i++){
        //velocidades en x,y,z
        particles[i].vx += 0.5 * acc_x[i] * dt;                 //Actualizamos velocidades medio paso
        particles[i].vy += 0.5 * acc_y[i] * dt;
        particles[i].vz += 0.5 * acc_z[i] * dt;
        //Movimiento en x,y,z 
        particles[i].x += particles[i].vx * dt;                 //Corregimos la posición paso entero
        particles[i].y += particles[i].vy * dt;
        particles[i].z += particles[i].vz * dt;

    }

    if (useBoundaries){
    applyReflectiveBC3D(particles, xmin, xmax, ymin, ymax, zmin, zmax);     //Checkeamos condicion: agregar condiciones de frontera
    }
    if (usePBounds){
      applyPeriodicBoundary(particles,Lx, Ly,Lz);  //Function periodic boundaries 
    }
    



    computeAccelerations3D(particles, acc_x, acc_y, acc_z, k, usePBounds,Lx,Ly,Lz);         //Calculamos la aceleración con la nueva configuración de posiciones
    for(int i = 0; i < N; i++){
        particles[i].vx += 0.5 * acc_x[i] * dt;                 //Nueva velocidad de la partícula
        particles[i].vy += 0.5 * acc_y[i] * dt;
        particles[i].vz += 0.5 * acc_z[i] * dt;
    }
}


void applyReflectiveBC3D(std::vector<Particle3D>& particles,
                         double xmin, double xmax,
                         double ymin, double ymax,
                         double zmin, double zmax)
{
    for (auto& p : particles) {

        // eje x
        if (p.x < xmin) {
            p.x = 2.0 * xmin - p.x;
            p.vx = -p.vx;
        }
        if (p.x > xmax) {
            p.x = 2.0 * xmax - p.x;
            p.vx = -p.vx;
        }

        // eje y
        if (p.y < ymin) {
            p.y = 2.0 * ymin - p.y;
            p.vy = -p.vy;
        }
        if (p.y > ymax) {
            p.y = 2.0 * ymax - p.y;
            p.vy = -p.vy;
        }

        // eje z
        if (p.z < zmin) {
            p.z = 2.0 * zmin - p.z;
            p.vz = -p.vz;
        }
        if (p.z > zmax) {
            p.z = 2.0 * zmax - p.z;
            p.vz = -p.vz;
        }

    }
}




void applyPeriodicBoundary(std::vector<Particle3D>& particles,
                           double Lx, double Ly, double Lz)
{
    for (auto& p : particles) {

        if (p.x >= Lx)
            p.x -= Lx;
        else if (p.x < 0.0)
            p.x += Lx;

        if (p.y >= Ly)
            p.y -= Ly;
        else if (p.y < 0.0)
            p.y += Ly;

        if (p.z >= Lz)
            p.z -= Lz;
        else if (p.z < 0.0)
            p.z += Lz;
    }
}
//Distribución condiciones iniciales


bool tooClose(const std::vector<Particle3D>& particles,
              double x, double y, double z,
              int current,
              double minDist,
              bool usePBounds,
              double Lx, double Ly, double Lz){
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

        if(r2 < minDist2)
            return true;
    }

    return false;
}

// Termostato de Andersen
void applyAndersenThermostat(std::vector<Particle3D>& particles,
                             double T_target,
                             double nu,
                             double dt,
                             int dim,
                             std::mt19937& gen){
    // Distribución uniforme para la probabilidad de colisión
    std::uniform_real_distribution<double> dist_prob(0.0, 1.0);
    
    // Distribución normal (Maxwell-Boltzmann) para las nuevas velocidades
    std::normal_distribution<double> dist_vel(0.0, std::sqrt(T_target));
    
    double collision_prob = nu * dt;
    
    for(auto& p : particles){
        // ¿Choca esta partícula con el baño térmico?
        if(dist_prob(gen) < collision_prob){
            
            // Asignar nuevas velocidades dependiendo de la dimensionalidad
            if(dim >= 1) p.vx = dist_vel(gen);
            if(dim >= 2) p.vy = dist_vel(gen);
            if(dim == 3) p.vz = dist_vel(gen);
        }
    }
}
