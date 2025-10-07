#ifndef HESTON_HPP
#define HESTON_HPP

#include "Option.hpp"
#include <complex>

class Heston : public Option {
private:
    double kappa;  // taux de retour à la moyenne
    double theta;  // variance à long terme
    double sigma;  // volatilité de la volatilité
    double rho;    // corrélation entre le spot et la variance
    double v0;     // variance initiale

    // Fonction caractéristique de Heston
    std::complex<double> characteristic_function(std::complex<double> u) const;

public:
    Heston(int phi_, double t0_, double T_, double S0_, double K_, 
           double r_, double q_, double market_price_, double vega_,
           double kappa_, double theta_, double sigma_, double rho_, double v0_);
    
    double prix() const override;
    
    // Implémentation de la fonction caractéristique
    std::complex<double> cf_logST(std::complex<double> u) const override {
        return characteristic_function(u);
    }
};

#endif