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


enum class Dimension {D1, D2, D3};      // clase para las dimensiones posibles del experimento (1D, 2D, 3D)
enum class Placement {Sphere, Uniform}; // Clase para las posiciones posibles de las particulas (Esfera al rededor del orignes, distribucion uniforme)

int main() {
    constexpr double PI = 3.14159265358979323846;
    int N = 100;                     // número de partículas
    double rho = 0.20;              // densidad del sistema
    double v_initial = 0.0;         // velocidad inicial
    
    // === definicion condiciones basicas del sistema ===
    Dimension sim_dim = Dimension::D3;  // escogemos la dimension del sistema poniendo D1, D2, D3
    Placement placement = Placement::Sphere; //escogemos si las particulas se distribuyen de forma uniforme o en una "esfera"

    // Switches para elegir condiones del sistema
    bool use_rotation = false;       // rotación 2D o 3D
    bool perturbation = true;        // perturbación en posiciones
    bool periodicB =true;
    bool reflectiveB= false;

    //switches para encender termostato de andersen
    bool use_thermostat = false;     // true: ensamble NVT / false: ensamble NVE
    double T_target = 1.0;          // Temperatura objetivo del baño térmico
    double nu = 10.0;               // Frecuencia de colisión con el baño (probabilidad)

    // =========================================================================
    // STEP 1: CREACION DE LA CAJA DINAMICA DEPENDIENTE DE N Y RHO
    // =========================================================================

    // === calculo dinamico de la caja (L) dependiendo de la dimension, N y rho
    double L = 0.0;
    switch (sim_dim){
        case Dimension::D1:
            L = N / rho;               //densidad lineal: rho = N / L
            break;
        case Dimension::D2:
            L = std::sqrt(N/rho);      //densidad superficial: rho = N / l^2
            break;
        case Dimension::D3:
            L = std::pow(N / rho, 1.0/3.0);  //densidad volumetrica: = N / L^3
            break;
        }

    //asiganmos dimensiones a la caja y determinamos que la esfera ocupa el 80% de ella
    double radius = 0.4 * L;
    double Lx = L, Ly = L, Lz = L;

    //en casos de 1D o 2D reducimos las dimensiones de la caja por coordenada
    if(sim_dim == Dimension::D1){Ly = 1.0; Lz = 1.0;}
    if(sim_dim == Dimension::D2){ Lz = 1.0;}

    //caja centrada en el origen
    double x_min = -Lx/2.0, x_max = Lx/2.0;
    double y_min = -Ly/2.0, y_max = Ly/2.0;
    double z_min = -Lz/2.0, z_max = Lz/2.0;

    std::cout << "Simulando N = " << N << " con densidad rho = " << rho << std::endl;
    std::cout << "Tamano de caja calculado L = " << L << std::endl;

    // Paso temporal
    double dt = 0.001;
    int steps = 50000;

    // Acoplamientos
    std::vector<double> k_harmonic = {0.0, 0.0, 0.0};  // z=0

    // Vector de partículas
    std::vector<Particle3D> particles(N);

    // Perturbación aleatoria
    double noise_amp = 0.2;
    std::mt19937 gen(42);
    std::uniform_real_distribution<double> dist(-noise_amp, noise_amp);

    // =========================================================================
    // STEP 2: INICIALIZACIÓN DE PARTÍCULAS
    // =========================================================================

    // Distribuciones para generar puntos aleatorios
    std::uniform_real_distribution<double> dist_01(0.0, 1.0);
    std::uniform_real_distribution<double> dist_theta(0.0, 2.0 * PI);
    
    for(int i = 0; i < N; i++){
        switch(sim_dim){

            case Dimension::D1:{
                // === Inicialización 1D ===
                //pone las particulas al rededor de un radio o del L de la caja
                double base_x;
                std::uniform_real_distribution<double> dist1D_box(-Lx/2.0, Lx/2.0);
                std::uniform_real_distribution<double> dist1D_sphere(-radius, radius);

                do{
                    if (placement == Placement::Uniform) {
                        base_x = dist1D_box(gen);
                    } else {
                        // En 1D, una "esfera" es simplemente un segmento de línea centrado
                        base_x = dist1D_sphere(gen);
                    }
                    
                    if (perturbation) base_x += dist(gen);
                    
                }while(tooClose(particles, base_x, 0.0, 0.0, i, 1.0, periodicB, Lx, Ly, Lz));

                particles[i].x = base_x;
                particles[i].y = 0.0;
                particles[i].z = 0.0;

                // Velocidades en 1D
                particles[i].vx = (i % 2 == 0) ? v_initial : -v_initial; // Mitad a la izq, mitad a la der
                particles[i].vy = 0.0;
                particles[i].vz = 0.0;
                break;
            }

            case Dimension::D2:{
                // === Inicialización 2D ===
                //pone las particulas uniformemente en la caja, o en un circula definido por L
                double base_x, base_y;
                std::uniform_real_distribution<double> dist2D_x(-Lx/2.0, Lx/2.0);
                std::uniform_real_distribution<double> dist2D_y(-Ly/2.0, Ly/2.0);

                do {
                    if(placement == Placement::Uniform){
                        base_x = dist2D_x(gen);
                        base_y = dist2D_y(gen);
                    }else{
                        // En 2D, distribución uniforme dentro de un círculo usando polares
                        double r = radius * std::sqrt(dist_01(gen));
                        double theta = dist_theta(gen);
                        base_x = r * std::cos(theta);
                        base_y = r * std::sin(theta);
                    }

                    if(perturbation){
                        base_x += dist(gen);
                        base_y += dist(gen);
                    }

                }while(tooClose(particles, base_x, base_y, 0.0, i, 1.0, periodicB, Lx, Ly, Lz));

                particles[i].x = base_x;
                particles[i].y = base_y;
                particles[i].z = 0.0;

                // Velocidades 2D
                double r_planar = std::sqrt(base_x*base_x + base_y*base_y);
                if (r_planar > 1e-12) {
                    if (use_rotation) {
                        particles[i].vx = -v_initial * base_y / r_planar;
                        particles[i].vy =  v_initial * base_x / r_planar;
                    } else {
                        particles[i].vx = v_initial * base_x / r_planar;
                        particles[i].vy = v_initial * base_y / r_planar;
                    }
                }else{
                    particles[i].vx = particles[i].vy = 0.0;
                }
                particles[i].vz = 0.0;
                break;
            }

            case Dimension::D3:{
                // === Inicialización 3D ===
                //en la caja o dentro de una esfera de radio dependiente de L
                double base_x, base_y, base_z;
                std::uniform_real_distribution<double> dist3D_x(-Lx/2.0, Lx/2.0);
                std::uniform_real_distribution<double> dist3D_y(-Ly/2.0, Ly/2.0);
                std::uniform_real_distribution<double> dist3D_z(-Lz/2.0, Lz/2.0);
                std::uniform_real_distribution<double> dist3D_sphere(-radius, radius);

                do{
                    if(placement == Placement::Uniform){
                        base_x = dist3D_x(gen);
                        base_y = dist3D_y(gen);
                        base_z = dist3D_z(gen);
                    }else{
                        // Rejection sampling para llenar la esfera uniformemente
                        do{
                            base_x = dist3D_sphere(gen);
                            base_y = dist3D_sphere(gen);
                            base_z = dist3D_sphere(gen);
                        }while(base_x*base_x + base_y*base_y + base_z*base_z > radius*radius);
                    }

                    if(perturbation){
                        base_x += dist(gen);
                        base_y += dist(gen);
                        base_z += dist(gen);
                    }

                }while(tooClose(particles, base_x, base_y, base_z, i, 1.0, periodicB, Lx, Ly, Lz));

                particles[i].x = base_x;
                particles[i].y = base_y;
                particles[i].z = base_z;

                // Velocidades 3D
                double mag = std::sqrt(base_x*base_x + base_y*base_y + base_z*base_z);
                if (mag > 1e-12) {
                    if(use_rotation){
                        // Rotación cilíndrica alrededor del eje Z
                        double r_planar = std::sqrt(base_x*base_x + base_y*base_y);
                        if(r_planar > 1e-12){
                            particles[i].vx = -v_initial * base_y / r_planar;
                            particles[i].vy =  v_initial * base_x / r_planar;
                            particles[i].vz = 0.0;
                        } else {
                            particles[i].vx = particles[i].vy = particles[i].vz = 0.0;
                        }
                    } else {
                        // Expansión radial en 3D
                        particles[i].vx = v_initial * base_x / mag;
                        particles[i].vy = v_initial * base_y / mag;
                        particles[i].vz = v_initial * base_z / mag;
                    }
                } else {
                    particles[i].vx = particles[i].vy = particles[i].vz = 0.0;
                }
                break;
            }
        }
    }


// =========================================================================
    // STEP 3: NOMBRADO Y GUARDADO DE ARCHIVOS
    // =========================================================================

    // crear la carpeta "results" si no existe
    namespace fs = std::filesystem;
    if(!fs::exists("results")){
        fs::create_directory("results");
    }

    // generar los nombres automáticos para los archivos de datos
    std::stringstream ss;

    // Dimensionalidad
    switch (sim_dim){
        case Dimension::D1: ss << "1D_"; break;
        case Dimension::D2: ss << "2D_"; break;
        case Dimension::D3: ss << "3D_"; break;
    }

    ss << (use_thermostat ? "NVT_" : "NVE_");

    // tipo de distribución
    ss << (placement == Placement::Uniform ? "UNI_" : "SPHERE_");

    // Parámetros físicos
    ss << "N" << N << "_rho" << std::fixed << std::setprecision(3) << rho;  
    ss << (use_rotation ? "ROT" : "RAD");                                   
    ss << "_v" << std::fixed << std::setprecision(1) << v_initial;          
    ss << (perturbation ? "_pert" : "_clean");                              
    ss << (periodicB ? "_period" : "_box");                               

    // Construir nombres de archivo finales
    std::string suffix = ss.str();
    std::string traj_filename = "results/tray_" + suffix + ".dat";
    std::string obs_filename  = "results/obs_" + suffix + ".dat";
    std::string rdf_filename  = "results/rdf_" + suffix + ".dat"; 

    // Guardando los archivos
    std::ofstream traj(traj_filename); 
    std::ofstream obs(obs_filename); 

    // =========================================================================
    // STEP 4: CALCULO DE TRAYECTORIAS Y OBSERVABLES
    // =========================================================================

    // Configuración para la Función de Distribución Radial (RDF)
    double dr_rdf = 0.05;                       // Resolución del histograma
    double max_dist_rdf = L / 2.0;              // Límite por Imagen Mínima
    int num_bins_rdf = static_cast<int>(max_dist_rdf / dr_rdf);
    std::vector<double> rdf_hist(num_bins_rdf, 0.0);
    int rdf_snapshots = 0;
    int sample_interval_rdf = 50;               // Muestrear cada 50 pasos para eficiencia

    // determinar los grados de libertad del sistema para calculo de la temperatura
    double d_f = 3.0;
    if(sim_dim == Dimension::D1) d_f = 1.0;
    else if(sim_dim == Dimension::D2) d_f = 2.0;

    // encabezado para el archivo de observables
    obs << "# t K/N U/N E/N T \n";

    // === INICIAR EL CRONOMETRO ===
    auto start_time = std::chrono::high_resolution_clock::now();

    for(int i = 0; i < steps; i++){

        double t = i * dt;

        velocityVerlet3D(particles, dt, k_harmonic, x_min, x_max, y_min, y_max, z_min, z_max, reflectiveB, periodicB, Lx, Ly, Lz);

        // aplicar el termostato si el sistema esta en NVT
        if (use_thermostat) {
            int current_dim = 3;
            if (sim_dim == Dimension::D1) current_dim = 1;
            else if (sim_dim == Dimension::D2) current_dim = 2;
            
            applyAndersenThermostat(particles, T_target, nu, dt, current_dim, gen);
        }

        // --- Muestreo de la RDF (Solo si es 3D y cada N pasos) ---
        if (sim_dim == Dimension::D3 && i % sample_interval_rdf == 0) {
            updateRDF3D(particles, rdf_hist, dr_rdf, periodicB, Lx, Ly, Lz);
            rdf_snapshots++;
        }

        // Trayectoria partículas
        traj << particles.size() <<"\n";
        traj << "#t = " << t << "\n";

        for(size_t j = 0; j < particles.size(); j++){
             traj << j << " "
                  << particles[j].x << " "
                  << particles[j].y << " "
                  << particles[j].z << "\n";
        }

        // archivo observables
        double K_total = kineticEnergy3D(particles);
        double U_total = potentialEnergy3D(particles, k_harmonic, periodicB, Lx, Ly, Lz);
        double E_total = K_total + U_total;

        // normalizar dividiendo por N
        double K_norm = K_total / N;
        double U_norm = U_total / N;
        double E_norm = E_total / N;

        double T_inst = (2.0*K_norm) / d_f;

        obs << t << " " << K_norm << " " << U_norm << " " << E_norm << " " << T_inst << "\n";
    }

    // --- FIN DEL CRONÓMETRO ---
    auto end_time = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed_seconds = end_time - start_time;

    // --- GUARDAR RDF (Solo si es 3D) ---
    if (sim_dim == Dimension::D3 && rdf_snapshots > 0) {
        normalizeAndSaveRDF3D(rdf_hist, rdf_filename, N, Lx, Ly, Lz, rdf_snapshots, dr_rdf);
    }

    std::cout << "Simulacion terminada. Tiempo de computo: " 
              << elapsed_seconds.count() << " segundos.\n";

    // Guardar los datos de N y Tiempo en un archivo de benchmark
    std::string bench_filename = "results/benchmark_CT_vs_N.dat";
        
    // std::ios::app abre el archivo en modo "append" (añadir al final)
    // Si el archivo no existe, lo crea.
    std::ofstream bench_file(bench_filename, std::ios::app); 

        // Si el archivo está vacío (recién creado), le ponemos un encabezado
    std::ifstream bench_check(bench_filename);
    bench_check.seekg(0, std::ios::end);
    if (bench_check.tellg() == 0) {
        bench_file << "# N Tiempo_Computo(s) Pasos Densidad\n";
    }
    
    // Escribir los datos de esta ejecución
    bench_file << N << " " 
               << elapsed_seconds.count() << " " 
               << steps << " " 
               << rho << "\n";

    bench_file.close();
    traj.close();
    obs.close();

    return 0;
}