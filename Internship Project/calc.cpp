#include <vector>
#include <string>
#include <stdexcept>
#include <iostream>


void drift(const std::vector<double>& n, const std::vector<double>& p, std::vector<double>& n_current, std::vector<double>& p_current,
    int Nx, double mu_n, double mu_p, double dt, double dx, std::vector<double>& E_change,
    std::vector<double>& n_deltaDrift, std::vector<double>& p_deltaDrift)
{
    //Absolute vertes kitoje vietoje
    n_deltaDrift[0] = n[0] * (mu_n * std::abs(E_change[0])) * dt / dx;
    p_deltaDrift[0] = p[0] * (mu_p * std::abs(E_change[0])) * dt / dx;

    if (E_change[0] > 0) {
        n_current[0] -= n_deltaDrift[0];
        n_current[1] += n_deltaDrift[0];

        p_current[0] -= p_deltaDrift[0];
    }
    else if (E_change[0] < 0) {
        n_deltaDrift[0] = 0;

        p_current[0] -= p_deltaDrift[0];
        p_current[1] += p_deltaDrift[0];
    }

    int last = Nx - 1;

    n_deltaDrift[last] = n[last] * (mu_n * std::abs(E_change[last])) * dt / dx;
    p_deltaDrift[last] = p[last] * (mu_p * std::abs(E_change[last])) * dt / dx;

    if (E_change[last] > 0) {
        p_current[last] -= p_deltaDrift[last];
        p_current[last - 1] += p_deltaDrift[last];

        n_current[last] -= n_deltaDrift[last];
    }
    else if (E_change[last] < 0) {
        p_deltaDrift[last] = 0;

        n_current[last] -= n_deltaDrift[last];
        n_current[last - 1] += n_deltaDrift[last];
    }


    for (int i = 1; i < Nx - 1; ++i) {

        n_deltaDrift[i] = n[i] * (mu_n * std::abs(E_change[i])) * dt / dx;

        if (E_change[i] > 0) {
            n_current[i] -= n_deltaDrift[i];
            n_current[i + 1] += n_deltaDrift[i];
        }
        else if (E_change[i] < 0) {
            n_current[i] -= n_deltaDrift[i];
            n_current[i - 1] += n_deltaDrift[i];
        }

        p_deltaDrift[i] = p[i] * (mu_p * std::abs(E_change[i])) * dt / dx;

        if (E_change[i] < 0) {
            p_current[i] -= p_deltaDrift[i];
            p_current[i + 1] += p_deltaDrift[i];
        }
        else if (E_change[i] > 0) {
            p_current[i] -= p_deltaDrift[i];
            p_current[i - 1] += p_deltaDrift[i];
        }
    }
}

/*
DERRIVED FORMULAS

Q = e * V{volume} * n{concentration} = e * S * dx * n;
I = dQ/dt; 
I = j{current density} * S;
(0)

(Provided formula)
j = e * D * (dn/dx); 
(1)

(0)==>>
I = e * S * dx * (dn/dt);
(0)==>>
j * S = e * S * dx * (dn/dt);
==>>
j = e * dx * (dn/dt);
(1)==>>
e * D * (dn/dx) = e * dx * (dn/dt);
==>>
D * (dn {difference between neighboring cells} /dx) * (1/dx) * dt = dn;

*/

void diffusion(const std::vector<double>& n, const std::vector<double>& p, std::vector<double>& n_current, std::vector<double>& p_current,
    double D_n, double D_p,
    double dx, double dt, double T, double mu_n, double mu_p, double k_B, double e_charge,
    int Nx) 
{

    // Assuuming constant dx & dt during the working time:
    double alpha_n = (D_n * dt) / (dx * dx);
    double alpha_p = (D_p * dt) / (dx * dx);

    
    // Solving the boundary of i = 0 separately:
    // Flow of Electrons due to diffusion TO the side of the current
    double flowC_n = alpha_n * (n[0] - n[1]);
    // Flow of Electrons due to diffusion TO the opposite side of the current
    double flowOppC_n = alpha_n * (n[0] - 0.0);

    // In a case of positive / negative numbers - the flow either adds or subtracts
    n_current[0] = n_current[0] - flowC_n - flowOppC_n;


    // Same thing here for the holes (can't pass through left side)
    double flowC_p = alpha_p * (p[0] - p[1]);
    p_current[0] = p_current[0] - flowC_p;


    // Solving the boundary of i = Nx separately:
    int last = Nx - 1;

    flowC_n = alpha_n * (n[last] - n[last - 1]);
    n_current[last] = n_current[last] - flowC_n;

    flowC_p = alpha_p * (p[last] - p[last - 1]);
    double flowOppC_p = alpha_p * (p[last] - 0.0);
    p_current[last] = p_current[last] - flowC_p - flowOppC_p;



    for (int i = 1; i < Nx - 1; ++i) {

        // Same thing here for the electrons
        flowC_n = alpha_n * (n[i] - n[i + 1]);
        flowOppC_n = alpha_n * (n[i] - n[i - 1]);
        n_current[i] = n_current[i] - flowC_n - flowOppC_n;
        
        // Same thing here for the holes
        double flowC_p = alpha_p * (p[i] - p[i + 1]);
        double flowOppC_p = alpha_p * (p[i] - p[i - 1]);
        p_current[i] = p_current[i] - flowC_p - flowOppC_p;

        if (n_current[i] < 0.0) {
            std::string error_msg = "Numerical instability: Electron concentration went negative at cell";
            throw std::runtime_error(error_msg);
        }
        if (p_current[i] < 0.0) {
            std::string error_msg = "Numerical instability: Hole concentration went negative at cell";
            throw std::runtime_error(error_msg);
        }
    }
}



/*
FORMULAS USED

gamma = e(mu_n + mu_p) / (eps_0 * eps) (1);
n_updated = gamma * n * p * dt;

*/

void recombination(const std::vector<double>& n, const std::vector<double>& p, std::vector<double>& n_current, std::vector<double>& p_current,
    double gamma, double dt, int Nx) {
    for (int i = 0; i < Nx; ++i) {
        double recombined =  1 * gamma * n[i] * p[i] * dt;
        n_current[i] = std::max(0.0, n_current[i] - recombined);
        p_current[i] = std::max(0.0, p_current[i] - recombined);
    }
}