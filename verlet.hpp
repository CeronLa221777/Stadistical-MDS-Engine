#ifndef VERLET_HPP
#define VERLET_HPP

#include <vector>
#include <random>

struct Particle3D {
    double x, y, z;
    double vx, vy, vz;
};

void computeAccelerations3D(const std::vector<Particle3D>& particles,
                            std::vector<double>& acc_x,
                            std::vector<double>& acc_y,
                            std::vector<double>& acc_z,
                            bool usePBounds,
                            double Lx, double Ly, double Lz);

void velocityVerlet3D(std::vector<Particle3D> &particles,
                      double dt,
                      double xmin, double xmax,
                      double ymin, double ymax,
                      double zmin, double zmax,
                      bool useBoundaries, bool usePBounds, 
                      double Lx, double Ly, double Lz);

void applyReflectiveBC3D(std::vector<Particle3D>& particles,
                         double xmin, double xmax,
                         double ymin, double ymax,
                         double zmin, double zmax);

void applyPeriodicBoundary(std::vector<Particle3D>& particles,
                           double Lx, double Ly, double Lz);

bool tooClose(const std::vector<Particle3D>& particles,
              double x, double y, double z,
              int current, double minDist, bool usePBounds,
              double Lx, double Ly, double Lz);

void initializeMaxwellBoltzmann(std::vector<Particle3D>& particles, 
                                double T_initial, int dim, std::mt19937& gen);

void applyBerendsenThermostat(std::vector<Particle3D>& particles,
                              double T_target, double tau, double dt, int dim);

#endif