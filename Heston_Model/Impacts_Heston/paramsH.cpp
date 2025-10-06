#include "Heston.hpp"
#include <complex>
#include <cmath>
#include <fstream>
#include <algorithm>
#include <vector>
#include <iostream>

double implied_volatility(double market_price, double S, double K, double T, double r, double q = 0.0) {
    
    // Prix black-scholes
    auto black_scholes_price = [&](double vol) {
        if (vol <= 1e-10) return std:: max(0.0, S * exp(- q * T) - K * exp(-r * T)); // Valeur intrinsèque pour vol proche de 0
        double d1 = (log(S / K) + (r - q + 0.5 * vol * vol) * T) / (vol * sqrt(T));
        double d2 = d1 - vol * sqrt(T);
        double call_price = S * exp(-q * T) * 0.5 * (1.0 + erf(d1 / sqrt(2.0))) - K * exp(-r * T) * 0.5 * (1.0 + erf(d2 / sqrt(2.0)));
        return call_price;
    
    };

    // Méthode de la bissection
    double tolerance = 1e-6;
    int max_iterations = 10000;
    double lower_bound = 1e-6;
    double upper_bound = 2.0; 
    double mid_vol = 0.0;

    for (int i = 0; i < max_iterations; ++i) {
        mid_vol = 0.5 * (lower_bound + upper_bound);
        double price = black_scholes_price(mid_vol);
        if (std::abs(price - market_price) < tolerance) {
            return mid_vol;
        }
        if (price < market_price) {
            lower_bound = mid_vol;
        } else {
            upper_bound = mid_vol;
        }
    }

    return mid_vol;
}

void visualize_parameter_impact() {
    std::cout << "Début de la visualisation des paramètres Heston..." << std::endl;
    
    // Paramètres de base
    double S0 = 100.0;
    double T = 1.0;
    double t0 = 0.0;
    double r = 0.05;
    double q = 0.0;
    
    // Paramètres Heston par défaut
    double v0_default = 0.04;      // Variance initiale (vol = 20%)
    double kappa_default = 1.0;    // Vitesse de retour
    double theta_default = 0.04;   // Variance long terme (vol = 20%)
    double sigma_default = 0.3;    // Vol of vol
    double rho_default = -0.7;     // Corrélation
    
    // Plage de strikes
    std::vector<double> strikes;
    for (double K = 50; K <= 175; K += 1.0) {
        strikes.push_back(K);
    }
    
    std::string kappa_path = std::string("Impacts_Heston") + "/kappa_impact.csv";
    std::string theta_path = std::string("Impacts_Heston") + "/theta_impact.csv";
    std::string sigma_path = std::string("Impacts_Heston") + "/sigma_impact.csv";
    std::string rho_path = std::string("Impacts_Heston") + "/rho_impact.csv";
    std::string v0_path = std::string("Impacts_Heston") + "/v0_impact.csv";
    
    // Fichiers de sortie
    std::ofstream file_kappa(kappa_path);
    std::ofstream file_theta(theta_path);
    std::ofstream file_sigma(sigma_path);
    std::ofstream file_rho(rho_path);
    std::ofstream file_v0(v0_path);
    
    // En-têtes
    file_kappa << "Strike,Kappa_0.5,Kappa_1.0,Kappa_2.0\n";
    file_theta << "Strike,Theta_0.01,Theta_0.04,Theta_0.09\n";
    file_sigma << "Strike,Sigma_0.1,Sigma_0.3,Sigma_0.5\n";
    file_rho << "Strike,Rho_-0.9,Rho_-0.5,Rho_0.0\n";
    file_v0 << "Strike,V0_0.01,V0_0.04,V0_0.09\n";
    
    std::cout << "Calcul des smiles de volatilité..." << std::endl;
    
    for (double K : strikes) {
        file_kappa << K;
        file_theta << K;
        file_sigma << K;
        file_rho << K;
        file_v0 << K;
        
        // === EFFET DE KAPPA ===
        std::vector<double> kappa_values = {0.5, 1.0, 2.0};
        for (double kappa : kappa_values) {
            Heston option(1, S0, K, T, t0, r, q, 0.0, 
                         kappa, theta_default, sigma_default, rho_default, v0_default, 0.0);
            double price = option.prix();
            double vol = implied_volatility(price, S0, K, T, r, q);
            file_kappa << "," << vol;
        }
        file_kappa << "\n";
        
        // === EFFET DE THETA ===
        std::vector<double> theta_values = {0.01, 0.04, 0.09}; // vol = 10%, 20%, 30%
        for (double theta : theta_values) {
            Heston option(1, S0, K, T, t0, r, q, 0.0, 
                         kappa_default, theta, sigma_default, rho_default, v0_default, 0.0);
            double price = option.prix();
            double vol = implied_volatility(price, S0, K, T, r, q);
            file_theta << "," << vol;
        }
        file_theta << "\n";
        
        // === EFFET DE SIGMA ===
        std::vector<double> sigma_values = {0.1, 0.3, 0.5};
        for (double sigma : sigma_values) {
            Heston option(1, S0, K, T, t0, r, q, 0.0, 
                         kappa_default, theta_default, sigma, rho_default, v0_default, 0.0);
            double price = option.prix();
            double vol = implied_volatility(price, S0, K, T, r, q);
            file_sigma << "," << vol;
        }
        file_sigma << "\n";
        
        // === EFFET DE RHO ===
        std::vector<double> rho_values = {0.5, -0.5, 0.0};
        for (double rho : rho_values) {
            Heston option(1, S0, K, T, t0, r, q, 0.0, 
                         kappa_default, theta_default, sigma_default, rho, v0_default, 0.0);
            double price = option.prix();
            double vol = implied_volatility(price, S0, K, T, r, q);
            file_rho << "," << vol;
        }
        file_rho << "\n";
        
        // === EFFET DE V0 ===
        std::vector<double> v0_values = {0.01, 0.04, 0.09}; // var = 10%, 20%, 30%
        for (double v0 : v0_values) {
            Heston option(1, S0, K, T, t0, r, q, 0.0, 
                         kappa_default, theta_default, sigma_default, rho_default, v0, 0.0);
            double price = option.prix();
            double vol = implied_volatility(price, S0, K, T, r, q);
            file_v0 << "," << vol;
        }
        file_v0 << "\n";
        
        if (int(K) % 10 == 0) {
            std::cout << "Progression: K=" << K << std::endl;
        }
    }
    
    file_kappa.close();
    file_theta.close();
    file_sigma.close();
    file_rho.close();
    file_v0.close();
    
    std::cout << "Fichiers générés avec succès!" << std::endl;
}

int main() {
    std::cout << "=== VISUALISATION DES PARAMÈTRES HESTON ===" << std::endl;
    
    try {
        visualize_parameter_impact();
        std::cout << "Succès! Exécutez 'python3 plot_heston_impact.py' pour voir les graphiques." << std::endl;
    } catch (const std::exception& e) {
        std::cerr << "Erreur: " << e.what() << std::endl;
        return 1;
    }
    
    return 0;
}