#ifndef LIFTED_HPP
#define LIFTED_HPP

#include "Option.hpp"
#include <vector>
#include <complex>

class Lifted : public Option {
private:
    double kappa;  // taux de retour à la moyenne
    double theta;  // variance à long terme
    double sigma;  // volatilité de la volatilité
    double rho;    // corrélation entre le spot et la variance
    double v0;     // variance initiale
    int n;         // nombre de facteurs mean-reverting
    double H;      // indice de Hurst
    double rn;     // paramètre r_N pour la paramétrisation
    
    // Paramètres calculés
    std::vector<double> c_i;  // poids
    std::vector<double> x_i;  // taux de mean-reversion
    
    // Initialisation des paramètres
    void initialize_parameters();
    
    // Fonction pour le système de Riccati
    void solve_riccati_system(double u_real, double u_imag, 
                             std::vector<std::complex<double>>& psi, 
                             double tau) const;
    
    // Fonction caractéristique
    std::complex<double> characteristic_function(std::complex<double> u) const;

    double prix(int N, double L) const;

    // Fonctions auxiliaires
    double g0(double t) const;
    std::complex<double> F(std::complex<double> u, std::complex<double> z) const;

public:
    double chi_k(int k, double c, double d, double a, double b) const;
    double psi_k(int k, double c, double d, double a, double b) const;
    Lifted(int phi_, double t0_, double T_, double S0_, double K_,
           double r_, double q_, double market_price_, double vega_,
           double kappa_, double theta_, double sigma_, double rho_, double v0_,
           int n_, double H_, double rn_);

    double prix() const override;
    
    // Implémentation de la fonction caractéristique
    std::complex<double> cf_logST(std::complex<double> u) const override {
        return characteristic_function(u);
    }
};

#endif