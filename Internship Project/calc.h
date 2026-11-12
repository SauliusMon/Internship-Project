#pragma once
#include <vector>


void drift(const std::vector<double>& n, const std::vector<double>& p, std::vector<double>& n_current, std::vector<double>& p_current,
    int Nx, double mu_n, double mu_p, double dt, double dx, std::vector<double>& E_change,
    std::vector<double>& n_deltaDrift, std::vector<double>& p_deltaDrift);


void diffusion(const std::vector<double>& n, const std::vector<double>& p, std::vector<double>& n_current, std::vector<double>& p_current,
    double D_n, double D_p,
    double dx, double dt, double T, double mu_n, double mu_p, double k_B, double e_charge,
    int Nx);

void recombination(const std::vector<double>& n, const std::vector<double>& p, std::vector<double>& n_current, std::vector<double>& p_current,
    double gamma, double dt, int Nx);