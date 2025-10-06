#include "../Option.hpp"
#include <vector>
#include <random>
#include <iostream>
#include <cmath>
#include <numeric>
#include <algorithm>

class LiftedHeston : public Option {
private:
    double kappa;   // Vitesse de réversion
    double theta;   // Variance à long terme
    double sigma;   // Volatilité de la variance
    double rho;     // Corrélation entre les deux processus de Wiener
    double v0;      // Variance initiale        
    double H;       // Paramètre de Hurst
    int N;          // Nombre de facteurs
    double alpha;   // alpha = H + 0.5
    double rn;      // Facteur d'espacement

public:
    // constructeur par défaut
    LiftedHeston(); 

    // constructeur paramétré
    LiftedHeston(int phi_, double S0_, double K_, double T_, double t0_,
                 double r_, double q_, double market_price_,
                 double kappa_, double theta_, double sigma_, double rho_, double v0_, double H_, int N_,
                 double vega_ , double rn_) ;

    virtual ~LiftedHeston();  

    double prix() const override {
        return 0.0; // On attend Felix
    }
    void simulate_paths(int N_steps, double T_hori, std::vector<double>& S_paths, std::vector<double>& V_paths, 
                        int N_paths, double kappa_, double theta_, double sigma_, double rho_, double v0_, 
                        double S0_, double H_, int N_, double rn_, double r_);
    
    void simulate_compare_N(
        const std::vector<int> &N_levels,
        const std::vector<double> &rn_levels,
        int N_steps,
        double T_hori,
        double kappa_, double theta_, double sigma_, double rho_,
        double v0_, double S0_,
        double r_, double q_,
        unsigned int seed = 42
    );
    struct Param_sup {
        std::vector<double> c;
        std::vector<double> x;
    };
   
    Param_sup fix_params(int N, double H, double alpha, double rn) const;

    double g0(double t, Param_sup params) const;  
    
    // accesseurs
    double getKappa() const;
    double getTheta() const;
    double getSigma() const;
    double getRho() const;
    double getV0() const;
    double getH() const;
    int    getN() const;
};