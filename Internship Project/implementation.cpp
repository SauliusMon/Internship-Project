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


const double dt = 1e-11;                    //d*d / (mu_p * U * Nx * 1000);                // Picosecond order step
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
const double L = 0.01;
const double C = (eps_0 * eps * S) / d;


void mainProgram() {
    system("start cmd /k \"C:\\Users\\Sauls\\AppData\\Local\\Programs\\Python\\Python312\\python.exe plotData.py\"");


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


    for (int t_step = 0; t_step < steps; ++t_step) {

        std::vector<double> n_current = n;
        std::vector<double> p_current = p;



        // Calculate drift
        std::vector<double> n_deltaDrift(Nx, 0.0);
        std::vector<double> p_deltaDrift(Nx, 0.0);

        drift(n, p, n_current, p_current, Nx, mu_n, mu_p, dt, dx, E, n_deltaDrift, p_deltaDrift);


        // Calculate diffusion
        //diffusion(n, p, n_current, p_current, D_n, D_p, dx, dt, T, mu_n, mu_p, k_B, e_charge, Nx);


        // Calculate recombination
        //recombination(n, p, n_current, p_current, gamma, dt, Nx);



        // Updating values
        for (int i = 0; i < Nx; ++i) {
            n[i] = n_current[i];
            p[i] = p_current[i];


            if (n[i] < 0.0) {
                std::string error_msg = "Numerical instability: Electron concentration went negative at cell";
                throw std::runtime_error(error_msg);
            }
            if (p[i] < 0.0) {
                std::string error_msg = "Numerical instability: Hole concentration went negative at cell";
                throw std::runtime_error(error_msg);
            }
        }


        // Puason equation: dE/dx = e(n-p)/(eps*eps_0)
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


        if (t_step % 100 == 0) {
            // Calculating total current
            double I_total = 0.0;
            for (int i = 0; i < Nx; ++i) {
                
                // Drift current density: j_dr = e * E * (n*mu_n + p*mu_p) TIEK KIEK ISEINA -> n[i]
                double j_dr = e_charge * ((n_deltaDrift[i]) + (p_deltaDrift[i])) * dx / dt;


                // Diffusion current density: j_diff = e * (D_n * dn/dx - D_p * dp/dx):
                int i_next = (i + 1 < Nx) ? i + 1 : i;
                int i_prev = (i - 1 >= 0) ? i - 1 : i;
                double denominator = (i_next != i_prev) ? (2.0 * dx) : dx;

                double dn_dx = (n[i_next] - n[i_prev]) / denominator;
                double dp_dx = (p[i_next] - p[i_prev]) / denominator;

                double j_diff = e_charge * (D_n * dn_dx - D_p * dp_dx);
                

                // Displacement current density: j_sl = eps*eps_0 * dE/dt - integralas 0
                double j_sl = eps * eps_0 * (E[i] - E_prev[i]) / dt;


                // Full current
                I_total += (j_dr) * S / Nx;
            }

            if (t_step * dt > 6e-6) {
                break;
            }

            // Writing TIME-&&-CURRENT data to text file
            std::ofstream appendfile("C:\\Users\\Sauls\\source\\repos\\Internship Project\\Internship Project\\tof_rezultatai.csv", std::ios::app);
            appendfile << t_step * dt << "," << I_total << "\n";
            appendfile.close();
            outfile.flush();


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
        E_prev = E;
    }

    outfile.close();
    std::cout << "Calculations done";

}