#include <vector>
#include <string>
#include <stdexcept>
#include <iostream>


void drift(const std::vector<double>& n, const std::vector<double>& p, std::vector<double>& n_driftUpdate, std::vector<double>& p_driftUpdate, 
    int Nx, double mu_n, double mu_p, double dt, double dx, const std::vector<double>& E) 
{
    std::vector<double> pushed_n(Nx, 0.0);
    std::vector<double> pushed_p(Nx, 0.0);

    pushed_n[0] = n[0] * (mu_n * std::abs(E[0])) * dt / dx;
    pushed_p[0] = p[0] * (mu_p * std::abs(E[0])) * dt / dx;

    if (E[0] > 0) {
        n_driftUpdate[0] -= pushed_n[0];
        n_driftUpdate[1] += pushed_n[0];

        p_driftUpdate[0] -= pushed_p[0];
    }
    else if (E[0] < 0) {
        n_driftUpdate[0] -= pushed_n[0];

        p_driftUpdate[0] -= pushed_p[0];
        p_driftUpdate[1] += pushed_p[0];
    }

    int last = Nx - 1;

    pushed_n[last] = n[last] * (mu_n * std::abs(E[last])) * dt / dx;
    pushed_p[last] = p[last] * (mu_p * std::abs(E[last])) * dt / dx;

    if (E[last] > 0) {
        p_driftUpdate[last] -= pushed_p[last];
        p_driftUpdate[last - 1] += pushed_p[last];

        n_driftUpdate[last] -= pushed_n[last];
    }
    else if (E[last] < 0) {
        p_driftUpdate[last] -= pushed_p[last];

        n_driftUpdate[last] -= pushed_n[last];
        n_driftUpdate[last - 1] += pushed_n[last];
    }


    for (int i = 1; i < Nx - 1; ++i) {

        pushed_n[i] = n[i] * (mu_n * std::abs(E[i])) * dt / dx;

        if (E[i] > 0) {
            n_driftUpdate[i] -= pushed_n[i];
            n_driftUpdate[i + 1] += pushed_n[i];
        }
        else if (E[i] < 0) {
            n_driftUpdate[i] -= pushed_n[i];
            n_driftUpdate[i - 1] += pushed_n[i];
        }

        pushed_p[i] = p[i] * (mu_p * std::abs(E[i])) * dt / dx;

        if (E[i] < 0) {
            p_driftUpdate[i] -= pushed_p[i];
            p_driftUpdate[i + 1] += pushed_p[i];
        }
        else if (E[i] > 0) {
            p_driftUpdate[i] -= pushed_p[i];
            p_driftUpdate[i - 1] += pushed_p[i];
        }

        if (pushed_n[i] <= 0 || pushed_p[i] <= 0) {
            std::cout << "";
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

void diffusion(const std::vector<double>& n, const std::vector<double>& p, std::vector<double>& n_diffusionUpdate, std::vector<double>& p_diffusionUpdate,
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
        n_diffusionUpdate[0] = n[0] - flowC_n - flowOppC_n;


        // Same thing here for the holes (can't pass through left side)
        double flowC_p = alpha_p * (p[0] - p[1]);
        p_diffusionUpdate[0] = p[0] - flowC_p;


    // Solving the boundary of i = Nx separately:
        
        int last = Nx - 1;

        flowC_n = alpha_n * (n[last] - n[last - 1]);
        n_diffusionUpdate[last] = n[last] - flowC_n;

        flowC_p = alpha_p * (p[last] - p[last - 1]);
        double flowOppC_p = alpha_p * (p[last] - 0.0);
        p_diffusionUpdate[last] = p[last] - flowC_p - flowOppC_p;



    for (int i = 1; i < Nx - 1; ++i) {

        // Same thing here for the electrons
        flowC_n = alpha_n * (n[i] - n[i + 1]);
        flowOppC_n = alpha_n * (n[i] - n[i - 1]);
        n_diffusionUpdate[i] = n[i] - flowC_n - flowOppC_n;
        
        // Same thing here for the holes
        double flowC_p = alpha_p * (p[i] - p[i + 1]);
        double flowOppC_p = alpha_p * (p[i] - p[i - 1]);
        p_diffusionUpdate[i] = p[i] - flowC_p - flowOppC_p;

        if (n_diffusionUpdate[i] < 0.0) {
            std::string error_msg = "Numerical instability: Electron concentration went negative at cell";
            throw std::runtime_error(error_msg);
        }
        if (p_diffusionUpdate[i] < 0.0) {
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

void recombination(const std::vector<double>& n, const std::vector<double>& p, std::vector<double>& n_recombinationUpdate, std::vector<double>& p_recombinationUpdate, 
    double gamma, double dt, int Nx) {
    for (int i = 0; i < Nx; ++i) {
        double recombined =  5000 * gamma * n[i] * p[i] * dt;
        n_recombinationUpdate[i] = std::max(0.0, n[i] - recombined);
        p_recombinationUpdate[i] = std::max(0.0, p[i] - recombined);
    }
}