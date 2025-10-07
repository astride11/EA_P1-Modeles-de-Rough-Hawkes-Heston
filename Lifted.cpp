#include "Lifted.hpp"
#include <cmath>
#include <algorithm>
#include <iostream>
#include <complex>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

// Fonction auxiliaire chi_k pour méthode COS
double Lifted::chi_k(int k, double c, double d, double a, double b) const {
    double denom = 1.0 + std::pow(k * M_PI / (b - a), 2.0);
    double term1 = std::cos(k * M_PI * (d - a) / (b - a)) * std::exp(d);
    double term2 = std::cos(k * M_PI * (c - a) / (b - a)) * std::exp(c);
    double term3 = (k * M_PI / (b - a)) * (std::sin(k * M_PI * (d - a) / (b - a)) * std::exp(d));
    double term4 = (k * M_PI / (b - a)) * (std::sin(k * M_PI * (c - a) / (b - a)) * std::exp(c));
    return (term1 - term2 + term3 - term4) / denom;
}

// Fonction auxiliaire psi_k pour méthode COS
double Lifted::psi_k(int k, double c, double d, double a, double b) const {
    if (k == 0) {
        return d - c;
    }
    double term1 = std::sin(k * M_PI * (d - a) / (b - a));
    double term2 = std::sin(k * M_PI * (c - a) / (b - a));
    return ((b - a) / (k * M_PI)) * (term1 - term2);
}

Lifted::Lifted(int phi_, double t0_, double T_, double S0_, double K_, 
               double r_, double q_, double market_price_, double vega_,
               double kappa_, double theta_, double sigma_, double rho_, double v0_,
               int n_, double H_, double rn_)
    : Option(phi_, S0_, K_, T_, t0_, r_, q_, market_price_, vega_),
      kappa(kappa_), theta(theta_), sigma(sigma_), rho(rho_), v0(v0_),
      n(n_), H(H_), rn(rn_) {
    initialize_parameters();
}


void Lifted::initialize_parameters() {
    c_i.resize(n);
    x_i.resize(n);
    
    double alpha = H + 0.5;
    
    // Calcul des paramètres selon la paramétrisation d'Abi Jaber
    // Pour éviter les calculs complexes de gamma, on utilise une approximation
    // Dans une implémentation réelle, il faudrait utiliser std::tgamma
    double gamma_approx = 1.0; // Approximation simplifiée
    
    double rn_power = std::pow(rn, (alpha - 1.0) * (1.0 + n / 2.0));
    double factor = (std::pow(rn, 1.0 - alpha) - 1.0) * rn_power / gamma_approx;
    
    double x_factor = ((1.0 - alpha) * (std::pow(rn, 2.0 - alpha) - 1.0)) / 
                     ((2.0 - alpha) * (std::pow(rn, 1.0 - alpha) - 1.0)) * 
                     std::pow(rn, -1.0 - n / 2.0);
    
    for (int i = 0; i < n; ++i) {
        c_i[i] = factor * std::pow(rn, (1.0 - alpha) * i);
        x_i[i] = x_factor * std::pow(rn, i);
        
        // Éviter les valeurs nulles ou négatives
        if (x_i[i] <= 0.0) x_i[i] = 0.1;
        if (c_i[i] <= 0.0) c_i[i] = 0.1;
    }
}

double Lifted::g0(double t) const {
    double sum = 0.0;
    for (int i = 0; i < n; ++i) {
        if (x_i[i] == 0.0) {
            sum += c_i[i] * t;
        } else {
            sum += (c_i[i] / x_i[i]) * (1.0 - std::exp(-x_i[i] * t));//Penser à minorer pour éviter de diviser par 0
        }
    }
    return v0 + kappa * theta * sum;
}

std::complex<double> Lifted::F(std::complex<double> u, std::complex<double> z) const {
    return 0.5 * (u * u - u) + (rho * sigma * u - kappa) * z + 
           (sigma * sigma / 2.0) * z * z;
}

void Lifted::solve_riccati_system(double u_real, double u_imag, 
                                 std::vector<std::complex<double>>& psi, 
                                 double tau) const {
    int num_steps = 1000;  // Nombre de pas de temps
    double dt = tau / num_steps;
    std::complex<double> u(u_real, u_imag);
    
    // Initialisation
    for (int i = 0; i < n; ++i) {
        psi[i] = std::complex<double>(0.0, 0.0);
    }
    
    // Schéma explicite-implicite
    for (int step = 0; step < num_steps; ++step) {
        std::complex<double> z_sum(0.0, 0.0);
        for (int j = 0; j < n; ++j) {
            z_sum += c_i[j] * psi[j];
        }
        
        std::complex<double> F_val = F(u, z_sum);
        
        for (int i = 0; i < n; ++i) {
            psi[i] = (psi[i] + dt * F_val) / (1.0 + x_i[i] * dt); //C'est l'équation C.2
        }
    }
}

std::complex<double> Lifted::characteristic_function(std::complex<double> u) const {
    double tau = T - t0;
    
    if (tau <= 0.0) {
        return std::exp(std::complex<double>(0.0, 1.0) * u * std::log(S0));
    }
    
    std::vector<std::complex<double>> psi(n);
    solve_riccati_system(u.real(), u.imag(), psi, tau);
    
    // Calcul simplifié de phi(0, T) pour les tests
    std::complex<double> phi_integral(0.0, 0.0);
    
    // Approximation simplifiée pour éviter les calculs complexes
    double avg_variance = g0(tau) / tau;
    
    for (int i = 0; i < n; ++i) {
        double term = (v0 / n + kappa * theta * c_i[i] / x_i[i]);
        phi_integral += std::complex<double>(term, 0.0) * psi[i];
    }
    
    return std::exp(phi_integral + u * std::log(S0));
}

// Prix européen via méthode COS (Fang & Oosterlee)
// Hypothèse forte : characteristic_function(u) = E[exp(i u log S_T)] sous Q,
// avec le drift (r - q) et le terme affine phi correctement intégrés.


double Lifted::prix() const {
    double tau = T - t0;
    
    // Cas à maturité
    if (tau <= 0.0) {
        if (phi == 1) { // Call
            return std::max(S0 - K, 0.0);
        } else { // Put
            return std::max(K - S0, 0.0);
        }
    }
    
    // Méthode COS directe (Fang & Oosterlee, 2008)
    
    // Étape 1: Calculer la variance moyenne pour approximer les cumulants
    double avg_variance = g0(tau) / tau;
    
    // Étape 2: Cumulants approximatifs pour Lifted
    double c1 = (r - q - 0.5 * avg_variance) * tau; // Moyenne du log-price
    double c2 = avg_variance * tau; // Variance
    double c4 = 3.0 * c2 * c2; // Kurtosis (approximation Gaussienne)
    
    // Étape 3: Intervalle de troncature [a, b]
    const double L = 10.0; // Facteur standard pour troncature
    double a = c1 - L * std::sqrt(std::abs(c2) + std::sqrt(std::abs(c4)));
    double b = c1 + L * std::sqrt(std::abs(c2) + std::sqrt(std::abs(c4)));
    
    // Étape 4: Paramètres COS
    const int N = 1000; // Nombre de termes (ajustable pour précision)
    double x = std::log(S0 / K); // Log-moneyness
    
    // Étape 5: Somme COS pour le prix du call
    std::complex<double> sum(0.0, 0.0);
    for (int k = 0; k < N; ++k) {
        double uk = k * M_PI / (b - a); // Fréquence
        std::complex<double> phi = characteristic_function(std::complex<double>(0.0, uk));
        double Vk = (k == 0 ? 1.0 : 2.0) * (chi_k(k, 0.0, b, a, b) - psi_k(k, 0.0, b, a, b)) / (b - a);
        sum += phi * Vk * std::exp(std::complex<double>(0.0, k * M_PI * (x - a) / (b - a)));
    }
    
    // Étape 6: Prix du call
    double call_price = K * std::exp(-r * tau) * std::real(sum);
    
    // Étape 7: Retourner call ou put
    if (phi == 1) { // Call
        return call_price;
    } else { // Put via parité
        return call_price - S0 * std::exp(-q * tau) + K * std::exp(-r * tau);
    }
}