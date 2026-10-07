#ifndef VERLET_HPP
#define VERLET_HPP

#include <vector>
#include <random>



struct Particle3D {
    double x, y, z;
    double vx, vy, vz;
};

//definicion de la función de aceleracion usadas en el resto del proyecto


void computeAccelerations3D(const std::vector<Particle3D>& particles,
                            std::vector<double>& acc_x,
                            std::vector<double>& acc_y,
                            std::vector<double>& acc_z,
                            const std::vector<double>& k,
                            bool usePBounds,
                            double Lx, double Ly, double Lz);

//definicion de función para el algoritmo velocity verlet
void velocityVerlet3D(std::vector<Particle3D> &particles,
                        double dt,
                        const std::vector<double>& k,
                        double xmin, double xmax,
                        double ymin, double ymax,
                        double zmin, double zmax,
                        bool useBoundaries, bool usePBounds,double Lx, double Ly,double Lz);

//funciones pertinenetes a las condiciones de frontera periodicas
void applyReflectiveBC3D(std::vector<Particle3D>& particles,
                         double xmin, double xmax,
                         double ymin, double ymax,
                         double zmin, double zmax);

//Condiciones periódicas
void applyPeriodicBoundary(std::vector<Particle3D>& particles,
                           double Lx, double Ly, double Lz);
//std::vec..bla bla bla es el sistema de N partículas, & se pasa por ref no se copia memoria Estos es la interfaz del integrador

//Ciudar las condiciones iniciales para que no se sobrepongan las partículas
bool tooClose(const std::vector<Particle3D>& particles,
              double x, double y, double z,
              int current,
              double minDist,
              bool usePBounds,
              double Lx, double Ly, double Lz);

// Termostato de Andersen
void applyAndersenThermostat(std::vector<Particle3D>& particles,
                             double T_target,
                             double nu,
                             double dt,
                             int dim,
                             std::mt19937& gen);


#endif