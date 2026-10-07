#ifndef OBSERVABLES_HPP
#define OBSERVABLES_HPP

#include <vector>
#include <string> // Necesario para los nombres de archivo
#include "verlet.hpp"

//aquí declararemos las variables de los observables a calcular

double kineticEnergy3D(const std::vector<Particle3D>& particles);
double potentialEnergy3D(const std::vector<Particle3D>& particles,
                         const std::vector<double>& k,
                         bool usePBounds,
                         double Lx, double Ly, double Lz);

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