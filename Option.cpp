#include "Option.hpp"

// ===== Constructeur par défaut =====
Option::Option() 
    : phi(1), t0(0.0), T(1.0), S0(100.0), K(100.0), r(0.05), q(0.0), 
      market_price(0.0), vega(0.0) {}

// ===== Constructeur paramétré =====
// Correction de l'ordre d'initialisation pour éviter les warnings
Option::Option(int phi_, double S0_, double K_, double T_, double t0_,
               double r_, double q_, double market_price_, double vega_)
    : phi(phi_), t0(t0_), T(T_), S0(S0_), K(K_), r(r_), q(q_), 
      market_price(market_price_), vega(vega_) {}

Option::~Option() {}

double Option::getS0() const { return S0; }
double Option::getK() const { return K; }
double Option::getT() const { return T; }
double Option::getR() const { return r; }
double Option::getQ() const { return q; }
int Option::getPhi() const { return phi; }
double Option::getMarketPrice() const { return market_price; }
double Option::getVega() const { return vega; }