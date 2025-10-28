#include "Heston.hpp"
#include <complex>
#include <vector>
#include <cmath>
#include <algorithm> 
#include <numeric>
#include <stdexcept>

using cdouble = std::complex<double>;

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif
// ===== Constructeur par défaut =====
Heston::Heston() 
    /*
    Nous initialisons les attributs avec des valeurs par defaut.
    Par defaut ce constructeur cree un Call avec les parametres suivants :
    t0 = 0.0
    T = 1.0 
    S0 = 100.0              
    K = 100.0
    r = 0.05 (5%)
    d = 0.0
    market_price = 0.0
    vega = 0.0
    kappa = 0.0
    theta = 0.0
    sigma = 0.0
    rho = 0.0
    v0 = 0.0
    */
    : Option(), kappa(0.0), theta(0.0), sigma(0.0), rho(0.0), v0(0.0) {}

// ===== Constructeur paramétré =====
Heston::Heston(int phi_, double S0_, double K_, double T_, double t0_,
               double r_, double q_, double market_price_,
               double kappa_, double theta_, double sigma_, double rho_, double v0_,
               double vega_)
    /*
    Nous initialisons les attributs avec les valeurs passees en parametres.
    */
    : Option(phi_, S0_, K_, T_, t0_, r_, q_, market_price_, vega_),
      kappa(kappa_), theta(theta_), sigma(sigma_), rho(rho_), v0(v0_) {}


// Fonction pour calculer le prix du call européen selon le modèle de Heston
double Heston::prix_call() const {
    double tau = T - t0;
    if (tau <= 0) {
        return std::max(0.0, S0 - K); // Valeur intrinsèque pour un call
    }
    double p1 = heston_probability(1, tau);
    double p2 = heston_probability(2, tau);
    double call_price = S0 * exp(-q * (T - t0)) * p1 - K * exp(-r * tau) * p2; 

    return call_price;
}

// Fonction caractéristique de Heston - CONFORME MATLAB
cdouble Heston::heston_characteristic(double u, int j, double tau) const {
    double a = kappa * theta;
    double u_param = (j == 1) ? 0.5 : -0.5;
    double b = (j == 1) ? (kappa - rho * sigma) : kappa;

    cdouble i(0, 1);
    cdouble u_i = u * i;

    // Calcul de d - conforme MATLAB
    cdouble d = sqrt(pow(rho * sigma * u_i - b, 2) - sigma * sigma * (2.0 * u_param * u_i - u * u));

    // Calcul de g - conforme MATLAB
    cdouble g = (b - rho * sigma * u_i + d) / (b - rho * sigma * u_i - d);

    // Calcul de C et D - conforme MATLAB
    cdouble exp_dt = exp(d * tau);
    cdouble C = (r - q) * u_i * tau + (a / (sigma * sigma)) * 
                ((b - rho * sigma * u_i + d) * tau - 2.0 * log((1.0 - g * exp_dt) / (1.0 - g)));
    
    cdouble D = ((b - rho * sigma * u_i + d) / (sigma * sigma)) * 
                ((1.0 - exp_dt) / (1.0 - g * exp_dt));

    // Fonction caractéristique
    cdouble phi = exp(C + D * v0 + u_i * log(S0));

    return phi; 
}

// Fonction pour calculer la probabilité P_j
double Heston::heston_probability(int j, double tau) const {
    auto integrand = [&](double u) -> double {
        if (u < 1e-10) u = 1e-10;

        cdouble char_func = heston_characteristic(u, j, tau);
        cdouble numerator = exp(cdouble(0, -1) * u * log(K)) * char_func;
        cdouble denominator = cdouble(0, 1) * u;

        return (numerator / denominator).real();
    };

    // Intégration numérique avec Simpson
    const int N = 1000;
    const double U_max = 100.0;
    double h = U_max / N;
    double integral = 0.0;

    for (int i = 0; i <= N; ++i) {
        double u = i * h;
        double weight = (i == 0 || i == N) ? 1 : (i % 2 == 0) ? 2 : 4;
        integral += weight * integrand(u);
    }
    integral *= h / 3.0;

    return 0.5 + integral / M_PI;
}

// Fonction prix() avec parité call-put
double Heston::prix() const {
    if (phi != 1 && phi != -1) {
        throw std::invalid_argument("Le paramètre phi doit être +1 (Call) ou -1 (Put).");
    }
    
    if (phi == 1) {
        return prix_call();
    } else {
        double call_price = prix_call();
        double put_price = call_price - S0 * exp(-q * (T - t0)) + K * exp(-r * (T - t0)); // ✅ Parité call-put standard
        return put_price;
    }
}
Heston::~Heston() {
    // Destructeur
    // Pas de ressources dynamiques à libérer dans cette classe de base
}

double Heston::getKappa() const {
    return kappa;
}   

double Heston::getTheta() const {
    return theta;
}
double Heston::getSigma() const {
    return sigma;
}
double Heston::getRho() const {
    return rho;
}
double Heston::getV0() const {
    return v0;
}
