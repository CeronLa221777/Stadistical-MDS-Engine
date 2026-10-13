#ifndef OBSERVABLES_HPP
#define OBSERVABLES_HPP

#include <vector>
#include <string> 
#include "verlet.hpp"

double kineticEnergy3D(const std::vector<Particle3D>& particles);
double potentialEnergy3D(const std::vector<Particle3D>& particles,
                         bool usePBounds,
                         double Lx, double Ly, double Lz);

// ------------------------------------------------------------------
// Equation of state (Numerical calculation)
// Instantaneous pressure: P = (2K + W) / (3V)
// (The Virial W is now calculated in verlet.cpp)
// ------------------------------------------------------------------
double pressure3D(double K, double W, double Lx, double Ly, double Lz);

// ------------------------------------------------------------------
// Radial Distribution Function g(r) functions
// ------------------------------------------------------------------
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