#ifndef OPTION_HPP
#define OPTION_HPP

#include <complex>

class Option {
protected:
    int phi;           // +1 pour Call, -1 pour Put
    double t0;         // date d'évaluation
    double T;          // maturité
    double S0;         // spot en l'instant t0
    double K;          // strike
    double r;          // taux sans risque
    double q;          // taux de dividende
    double market_price; // prix observé sur le marché en t0
    double vega;       // vega observé sur le marché en t0

public:
    // Constructeurs
    Option();
    Option(int phi_, double S0_, double K_, double T_, double t0_,
           double r_, double q_, double market_price_, double vega_);
    virtual ~Option();

    // Getters
    double getS0() const;
    double getK() const;
    double getT() const;
    double getR() const;
    double getQ() const;
    int getPhi() const;
    double getMarketPrice() const;
    double getVega() const;

    // Méthode virtuelle pure pour le prix
    virtual double prix() const = 0;
    
    // Méthode virtuelle pour la fonction caractéristique (optionnelle)
    virtual std::complex<double> cf_logST(std::complex<double> u) const {
        return std::complex<double>(0,0); // Implémentation par défaut
    }
};

#endif