// Calibration.hpp
#pragma once
#include <vector>
#include <array>
#include <functional>
#include <random>
#include <limits>
#include "Option.hpp"
#include "Heston.hpp"

struct HestonParams {
    double kappa;
    double theta;
    double sigma;
    double rho;
    double v0;
};

struct CalibConfig {
    bool   use_implied_vol = false; // false = calib sur prix; true = sur vols impl.
    bool   vega_weights    = true;
    double feller_penalty  = 100.0; // lambda pour pénaliser max(0, sigma^2 - 2 kappa theta)
    double price_eps       = 1e-10; // plancher prix pour erreurs relatives
    double vol_bp_target   = 1e-4;  // tolérence cible ~ 1bp
    int    max_iters_local = 500;
    int    multistart      = 8;     // nb redémarrages
    unsigned seed          = 42;
};

// -------- utilitaires Black-Scholes pour l’inversion de vol --------
double norm_cdf(double x);
double norm_pdf(double x);
double bs_price(bool is_call, double F, double K, double T, double vol, double df=1.0);
double bs_vega (double F, double K, double T, double vol, double df=1.0);

// robuste : Brent + newton safe
double implied_vol_bs(bool is_call, double F, double K, double T, double price, double df=1.0);



// ------------------ Calibrateur principal --------------------------
class HestonCalibrator {
public:
    HestonCalibrator(const std::vector<Option*>& market, CalibConfig cfg = {});
    // lance la calibration et renvoie les paramètres + erreur RMS
    std::pair<HestonParams,double> calibrate();
    
    // utilitaire: prix/vol modèle pour une option avec params p
    double model_price_for(const Option& opt, const HestonParams& p) const;
    double model_impl_vol_for(const Option& opt, const HestonParams& p) const;
private:
    const std::vector<Option*>& market_;
    CalibConfig cfg_;

    // mapping x (R^5 non borné) -> params bornés
    static HestonParams transform(const std::array<double,5>& x);
    static std::array<double,5> inv_transform(const HestonParams& p);

    // fonction objectif (RMS pondéré)
    double objective_from_x(const std::array<double,5>& x) const;



    // simple Nelder-Mead (ou branche vers ton LBFGS si tu en as un)
    std::array<double,5> nelder_mead(const std::array<double,5>& x0, double scale, int max_iter) const;
};
