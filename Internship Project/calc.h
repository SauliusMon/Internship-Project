#pragma once
#include <vector>


void drift(const std::vector<double>& n, const std::vector<double>& p, std::vector<double>& n_driftUpdate, std::vector<double>& p_driftUpdate,
    int Nx, double mu_n, double mu_p, double dt, double dx, const std::vector<double>& E);

void diffusion(const std::vector<double>& n, const std::vector<double>& p, std::vector<double>& n_diffusionUpdate, std::vector<double>& p_diffusionUpdate,
    double D_n, double D_p,
    double dx, double dt, double T, double mu_n, double mu_p, double k_B, double e_charge,
    int Nx);

void recombination(const std::vector<double>& n, const std::vector<double>& p, std::vector<double>& n_recombinationUpdate, std::vector<double>& p_recombinationUpdate,
    double gamma, double dt, int Nx);