#include "Option.hpp"

    /*
    Nous implementons la classe abstraite Option de laquelle heriteront les classes Modeles qui seront devellopees par la suite.
    params : 
    phi -> +1 pour Call, -1 pour Put
    t0 -> date d'evaluation
    T -> maturite
    S0 -> spot en l'intant t0
    K -> strike
    r -> taux sans risque
    q -> taux de dividende
    market_price -> prix observe sur le marche en t0    
    vega -> vega observe sur le marche en t0
    */


// ===== Constructeur par défaut =====
Option::Option() 

    /*
    Nous initialisons les attributs avec des valeurs par defaut.
    Par defaut ce constructeur cree un Call avec les parametres suivants :
    t0 = 0.0
    T = 1.0 
    S0 = 100.0              
    K = 100.0
    r = 0.05 (5%)
    q = 0.0
    market_price = 0.0
    vega = 0.0
    */
    : phi(1), t0(0.0), T(1.0), S0(100.0), K(100.0), r(0.05), q(0.0), market_price(0.0), vega(0.0) {}


// ===== Constructeur paramétré =====
Option::Option(int phi_, double S0_, double K_, double T_, double t0_,
               double r_, double q_, double market_price_, double vega_)
    /*
    Nous initialisons les attributs avec les valeurs passees en parametres.
    */
    : phi(phi_), t0(t0_), T(T_), S0(S0_), K(K_), r(r_), q(q_), market_price(market_price_), vega(vega_) {}

Option::~Option() {
    // Destructeur
    // Pas de ressources dynamiques à libérer dans cette classe de base
}

double Option::getS0() const {
    return S0;
}
double Option::getK() const {
    return K;
}
double Option::getT() const {
    return T;
}
double Option::getR() const {
    return r;
}
double Option::getQ() const {
    return q;
}
int Option::getPhi() const {
    return phi;
}
double Option::getMarketPrice() const {
    return market_price;
}
double Option::getVega() const {
    return vega;
}
