#ifndef HESTON_HPP
#define HESTON_HPP

#include "Option.hpp"
#include <complex>

using cdouble = std::complex<double>;
class Heston : public Option {
private:
    double kappa;   // Vitesse de réversion
    double theta;   // Variance à long terme
    double sigma;   // Volatilité de la variance
    double rho;     // Corrélation entre les deux processus de Wiener
    double v0;      // Variance initiale       

public:
    // constructeur par défaut
    Heston();  

    // constructeur paramétré
    Heston(int phi_, double S0_, double K_, double T_, double t0_,
           double r_, double d_, double market_price_,
           double kappa_, double theta_, double sigma_, double rho_, double v0_,
           double vega_ = 0.0); 

    virtual ~Heston();
    
    virtual double prix() const override; // implémentation de la fonction pure virtuelle

    double prix_call() const; // fonction pour calculer le prix du call européen selon le modèle de Heston

    double heston_probability(int j, double tau) const; // fonction pour calculer la probabilité P_j selon le modèle de Heston
    
    cdouble heston_characteristic(double phi, int j, double tau) const;
    // accesseurs
    double getKappa() const;
    double getTheta() const;
    double getSigma() const;
    double getRho() const;
    double getV0() const;
};

#endif // HESTON_HPP