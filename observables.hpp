#ifndef OBSERVABLES_HPP
#define OBSERVABLES_HPP

#include <vector>
#include <string> 
#include "verlet.hpp"



double kineticEnergy3D(const std::vector<Particle3D>& particles);
double potentialEnergy3D(const std::vector<Particle3D>& particles,
                         bool usePBounds,
                         double Lx, double Ly, double Lz);

// Ecuacion de estado problema 4 numerico
double potentialEnergy3D(const std::vector<Particle3D>& particles,
                         bool usePBounds,
                         double Lx, double Ly, double Lz);

// ---- Presión por el teorema del virial (punto 4) ----
// W = sum_{i<j} f(r_ij) . r_ij   (virial de pares, truncado en r_c)
double virial3D(const std::vector<Particle3D>& particles,
                bool usePBounds,
                double Lx, double Ly, double Lz);

// P = rho*kB*T_inst + W/(3V),  con kB*T_inst = 2K/(3N)  =>  P = (2K + W)/(3V)
double pressure3D(const std::vector<Particle3D>& particles,
                  bool usePBounds,
                  double Lx, double Ly, double Lz);

// Nuevas funciones para el g(r)


// Nuevas funciones para el g(r)
void updateRDF3D(const std::vector<Particle3D>& particles,
                 std::vector<double>& rdf_hist,
                 double dr,
                 bool usePBounds,
                 double Lx, double Ly, double Lz);

void normalizeAndSaveRDF3D(const std::vector<double>& rdf_hist,
                           const std::string& filename,
                           int N, double Lx, double Ly, double Lz,
                           int num_snapshots, double dr);

#endif