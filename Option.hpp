#pragma once   // évite les inclusions multiples

class Option {
protected:
    int phi;            // +1 = Call, -1 = Put
    double t0;          // Date d’évaluation
    double T;           // Maturité
    double S0;          // Spot
    double K;           // Strike
    double r;           // Taux sans risque
    double q;           // Taux de dividende
    double market_price;// Prix observé sur le marché
    double vega;        // Vega observé

public:
    Option();  // constructeur par défaut
    Option(int phi_, double S0_, double K_, double T_, double t0_,
           double r_, double q_, double market_price_, double vega_ = 0.0);

    virtual ~Option();

    virtual double prix() const = 0; // fonction pure virtuelle

    // accesseurs
    double getS0() const;
    double getK() const;
    double getT() const;
    double getR() const;
    double getQ() const;
    int    getPhi() const;
    double getMarketPrice() const;
    double getVega() const;
};
