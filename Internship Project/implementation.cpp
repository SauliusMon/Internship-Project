#include <iostream>
#include "implementation.h" 
#include "calc.h"
#include <vector>
#include <cmath>
#include <fstream>


const double eps_0 = 8.854e-12;         // Vacuum dielectric permittivity (F/m)
const double eps = 3.0;                 // Material dielectric permittivity
const double d = 1e-6;                  // Thickness d = 1 micrometer (m)
const double S = 4e-6;                  // Area S = 2x2 mm^2 (m^2)
const double U = 3.0;                   // Voltage U = 3 V (or 6 V)
const double mu_n = 1e-7;               // Electron mobility (m^2/Vs)
const double mu_p = 5e-12;               // Hole mobility (m^2/Vs)


const int Nx = 1000;                     // Number of space layers
const double dx = d / Nx;               // Space step (dx)


const double dt = 1e-10;     //d*d / (mu_p * U * Nx * 1000);                // Picosecond order step
const int steps = 10000000;                                  // Number of simulation steps


std::vector<double> n(Nx, 0.0);
std::vector<double> p(Nx, 0.0);
std::vector<double> E(Nx, U / d);
std::vector<double> E_prev = E;


const double k_B = 1.380649e-23;                            // Bolcman constant
const double e_charge = 1.602176634e-19;                    // Electric charge constant
const double T = 300;                                       // Temperature (Kelvin)


const double D_n = (k_B * T * mu_n) / e_charge;             // Diffusion constant electrons (m^2/s)
const double D_p = (k_B * T * mu_p) / e_charge;             // Diffusion constant holes (m^2/s)


const double gamma = (e_charge * (mu_n + mu_p)) / (eps * eps_0);    // Gamma coefficient (probability that electron and hole collide)


const double alpha = 20e6;
const double L = 0.1;
const double C = (eps_0 * eps * S) / L;


void mainProgram() {
    //double N0 = 1e17;                    // Generated couples density

    //for (int i = 0; i < Nx; ++i) {
    //    n[i] = N0 * std::exp(-alpha * (i * dx));
    //    p[i] = n[i];
    //}


    double N0 = (L * C * U * alpha) / ((1.0 - std::exp(-alpha * d)) * e_charge * S);

    for (int i = 0; i < Nx; ++i) {
        n[i] = N0 * std::exp(-alpha * ( i * dx));
        p[i] = n[i];
    }


    std::ofstream outfile("C:\\Users\\Sauls\\source\\repos\\Internship Project\\Internship Project\\tof_rezultatai.csv");
    outfile << "Laikas_s,Srove_A\n";
    outfile.flush();

    
    double last_I_total = 0;
    int same_I_count = 0;

    std::vector<double> n_old = n;
    std::vector<double> p_old = p;


    for (int t_step = 0; t_step < steps; ++t_step) {

        std::vector<double> n_driftUpdate = n;
        std::vector<double> p_driftUpdate = p;

        // Calculate drift
        drift(n, p, n_driftUpdate, p_driftUpdate, Nx, mu_n, mu_p, dt, dx, E);
        
        //n = n_driftUpdate;
        //p = p_driftUpdate;



        std::vector<double> n_diffusionUpdate = n;
        std::vector<double> p_diffusionUpdate = p;

        // Calculate diffusion
        diffusion(n, p, n_diffusionUpdate, p_diffusionUpdate, D_n, D_p, dx, dt, T, mu_n, mu_p, k_B, e_charge, Nx);

        //n = n_diffusionUpdate;
        //p = p_diffusionUpdate;



        std::vector<double> n_recombinationUpdate = n;
        std::vector<double> p_recombinationUpdate = p;

        // Calculate recombination
        //recombination(n, p, n_recombinationUpdate, p_recombinationUpdate, gamma, dt, Nx);

        //n = n_recombinationUpdate;
        //p = p_recombinationUpdate;


        // Updating values
        for (int i = 0; i < Nx; ++i) {
            double delta_n_drift = n_driftUpdate[i] - n[i];
            double delta_n_diff = n_diffusionUpdate[i] - n[i];
            double delta_n_rec = n_recombinationUpdate[i] - n[i];

            double delta_p_drift = p_driftUpdate[i] - p[i];
            double delta_p_diff = p_diffusionUpdate[i] - p[i];
            double delta_p_rec = p_recombinationUpdate[i] - p[i];

            n[i] = n[i] + delta_n_drift + delta_n_diff + delta_n_rec;
            p[i] = p[i] + delta_p_drift + delta_p_diff + delta_p_rec;

            if (delta_n_rec != 0.0) {
                std::cout << "";
            }

            if (n[i] < 0.0) {
                std::string error_msg = "Numerical instability: Electron concentration went negative at cell";
                throw std::runtime_error(error_msg);
            }
            if (p[i] < 0.0) {
                std::string error_msg = "Numerical instability: Hole concentration went negative at cell";
                throw std::runtime_error(error_msg);
            }
        }


        // Puasono equation: dE/dx = e(n-p)/(eps*eps_0)
        double integral_E = E[0] * dx;
        for (int i = 1; i < Nx; ++i) {
            double rho = e_charge * (p[i] - n[i]);
            E[i] = E[i - 1] + rho / (eps * eps_0) * dx;
            integral_E += E[i] * dx;
        }

        // Voltage correction: U = integral(E dx)
        double correction = (U - integral_E) / d;
        for (int i = 0; i < Nx; ++i) {
            E[i] += correction;
        }


        // Calculating total current
        double I_total = 0.0;
        for (int i = 0; i < Nx; ++i) {
            // Dreifo srovės tankis: j_dr = e * E * (n*mu_n + p*mu_p) TIEK KIEK ISEINA -> n[i]
            //double j_dr = e_charge * E[i] * ((n_driftUpdate[i] - n_old[i]) * mu_n + (p_driftUpdate[i] - p_old[i]) * mu_p);
            //double j_dr = e_charge * E[i] * (n[i] * mu_n + p[i] * mu_p);
            double j_dr = e_charge * ((n_old[i] - n_driftUpdate[i]) + (p_old[i] - p_driftUpdate[i])) * dx / dt;

            // Poslinkio srovės tankis: j_sl = eps*eps_0 * dE/dt
            double j_sl = eps * eps_0 * (E[i] - E_prev[i]) / dt;

            // Pilnutinė srovė
            I_total += (j_dr + j_sl) * dx / (S * Nx);
        }

        if (t_step * dt > 6e-6) {
            break;
        }

        E_prev = E;


        if (t_step % 100 == 0) {
            std::ofstream appendfile("C:\\Users\\Sauls\\source\\repos\\Internship Project\\Internship Project\\tof_rezultatai.csv", std::ios::app);
            appendfile << t_step * dt << "," << I_total << "\n";
            appendfile.close();
            outfile.flush();

            //std::cout << t_step * dt << "," << I_total << "\n";
            //std::cout.flush();

            if (std::abs(I_total - last_I_total) < 1e-15) {
                same_I_count++;
            }
            else {
                same_I_count = 1;
                last_I_total = I_total;
            }

            if (same_I_count >= 10) {
                break;
            }
        }
    }

    outfile.close();
    std::cout << "Calculations done";

    system("\"C:\\Users\\Sauls\\AppData\\Local\\Programs\\Python\\Python312\\python.exe\" plotData.py");
}