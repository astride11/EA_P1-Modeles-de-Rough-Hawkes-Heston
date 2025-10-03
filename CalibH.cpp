// Calibration.cpp
#include "CalibH.hpp"
#include <cmath>
#include <algorithm>
#include <cassert>
#include <iostream>

// ------------------ BS utils -------------------
static inline double sgn(bool is_call){ return is_call ? +1.0 : -1.0; }

double norm_pdf(double x) { static const double inv_sqrt_2pi=0.3989422804014327; return inv_sqrt_2pi*std::exp(-0.5*x*x); }
double norm_cdf(double x) { return 0.5*std::erfc(-x*M_SQRT1_2); }

double bs_price(bool is_call, double F, double K, double T, double vol, double df){
    if (vol<=0.0 || T<=0.0){ // limite: prix intrinsèque sous-forward
        double intrinsic = std::max(sgn(is_call)*(F-K), 0.0);
        return df*intrinsic;
    }
    double s = vol*std::sqrt(T);
    double d1 = (std::log(F/K) + 0.5*s*s) / s;
    double d2 = d1 - s;
    if (is_call) return df*( F*norm_cdf(d1) - K*norm_cdf(d2) );
    else         return df*( K*norm_cdf(-d2) - F*norm_cdf(-d1) );
}

double bs_vega(double F, double K, double T, double vol, double df){
    if (vol<=0.0 || T<=0.0) return 0.0;
    double s = vol*std::sqrt(T);
    double d1 = (std::log(F/K)+0.5*s*s)/s;
    return df * F * std::sqrt(T) * norm_pdf(d1);
}

// Brent + safeguarded Newton
double implied_vol_bs(bool is_call, double F, double K, double T, double price, double df){
    // bornes raisonnables
    double p = price/df;
    double lo = 1e-6, hi = 5.0; // [~0%, 500%]
    // ajuste hi pour contenir le prix
    auto pr = [&](double v){ return bs_price(is_call,F,K,T,v,1.0); };
    if (pr(hi) < p) hi = 10.0;
    if (pr(hi) < p) hi = 20.0;

    // bracketing
    for(int i=0;i<50 && (pr(lo)>p);++i) lo*=0.5;
    for(int i=0;i<50 && (pr(hi)<p);++i) hi*=2.0;

    double v = 0.5*(lo+hi);
    for(int it=0; it<100; ++it){
        double pv = pr(v);
        double diff = pv - p;
        if (std::abs(diff) < 1e-12) break;
        // newton step
        double vega = bs_vega(F,K,T,v,1.0);
        double v_new = v - diff/(vega>1e-12?vega:1e-12);
        // garde dans [lo,hi]
        if (v_new<=lo || v_new>=hi) v_new = 0.5*(lo+hi);
        // update bracket
        if (diff>0) hi=v; else lo=v;
        v = v_new;
    }
    return std::max(1e-8, v);
}

// ------------------ Calibrator -------------------
HestonCalibrator::HestonCalibrator(const std::vector<Option*>& market, CalibConfig cfg)
: market_(market), cfg_(cfg) {}

// non borné -> borné
HestonParams HestonCalibrator::transform(const std::array<double,5>& x){
    HestonParams p;
    p.kappa = std::exp(x[0]);          // >0
    p.theta = std::exp(x[1]);          // >0
    p.sigma = std::exp(x[2]);          // >0
    p.rho   = std::tanh(x[3]);         // (-1,1)
    p.v0    = std::exp(x[4]);          // >=0
    return p;
}

std::array<double,5> HestonCalibrator::inv_transform(const HestonParams& p){
    return { std::log(std::max(1e-12,p.kappa)),
             std::log(std::max(1e-12,p.theta)),
             std::log(std::max(1e-12,p.sigma)),
             std::atanh(std::clamp(p.rho,-0.999,0.999)),
             std::log(std::max(1e-12,p.v0)) };
}

double HestonCalibrator::model_price_for(const Option& opt, const HestonParams& p) const{
    // On crée un Heston temporaire avec les mêmes données marché que opt
    const bool is_call = (opt.getPhi()==+1);
    Heston h(is_call?+1:-1,
             opt.getS0(), opt.getK(), opt.getT(), /*t0*/0.0,
             opt.getR(), opt.getQ(), /*market_price*/0.0,
             p.kappa, p.theta, p.sigma, p.rho, p.v0);
    // On appelle le priceur
    double price_call = h.prix_call();
    if (is_call) return price_call;
    // put via parité : P = C - F + K*df
    double T = opt.getT();
    double r = opt.getR(), q=opt.getQ();
    double df = std::exp(-r*T);
    double fwd= opt.getS0()*std::exp((r-q)*T);
    return price_call - fwd*df + opt.getK()*df;
}

double HestonCalibrator::model_impl_vol_for(const Option& opt, const HestonParams& p) const{
    double T = opt.getT();
    double r = opt.getR(), q = opt.getQ();
    double df = std::exp(-r*T);
    double F  = opt.getS0()*std::exp((r-q)*T);
    double price = model_price_for(opt,p);
    bool is_call = (opt.getPhi()==+1);
    return implied_vol_bs(is_call, F, opt.getK(), T, price, df);
}

double HestonCalibrator::objective_from_x(const std::array<double,5>& x) const{
    HestonParams p = transform(x);

    // pénalisation Feller douce si 2 kappa theta <= sigma^2
    double feller_violation = std::max(0.0, p.sigma*p.sigma - 2.0*p.kappa*p.theta);
    double penalty = cfg_.feller_penalty * feller_violation*feller_violation;

    double num=0.0, den=0.0;

    for (const Option* opt : market_){
        double w = 1.0;
        if (cfg_.vega_weights){
            // vega BS au vol de marché si disponible via opt->getVega(), sinon approx ATM
            double vega_mkt = opt->getVega();
            if (vega_mkt<=0.0){
                double T = opt->getT();
                double r = opt->getR(), q=opt->getQ();
                double df = std::exp(-r*T);
                double F  = opt->getS0()*std::exp((r-q)*T);
                // heuristique vol: à défaut, 20%
                double vol0 = 0.2;
                vega_mkt = bs_vega(F, opt->getK(), T, vol0, df);
            }
            w = 1.0 / std::max(1e-8, vega_mkt);
        }

        if (!cfg_.use_implied_vol){
            // erreur sur PRIX, relative (stabilisée)
            double model_p = model_price_for(*opt, p);
            double mkt_p   = opt->getMarketPrice();
            double scale   = std::max(cfg_.price_eps, std::abs(mkt_p));
            double e = (model_p - mkt_p)/scale;
            num += w*e*e;
            den += w;
        }else{
            // erreur sur VOL implicite (en points)
            double model_iv = model_impl_vol_for(*opt, p);
            double T = opt->getT();
            double r = opt->getR(), q=opt->getQ();
            double df = std::exp(-r*T);
            double F  = opt->getS0()*std::exp((r-q)*T);
            bool is_call = (opt->getPhi()==+1);
            double mkt_iv = implied_vol_bs(is_call, F, opt->getK(), T, opt->getMarketPrice(), df);

            double e = (model_iv - mkt_iv);
            num += w*e*e;
            den += w;
        }
    }

    double mse = (den>0 ? num/den : num);
    return std::sqrt(mse) + penalty;
}

// --------- Simplex (Nelder–Mead) minimaliste, suffisant pour un premier jet -----
std::array<double,5> HestonCalibrator::nelder_mead(const std::array<double,5>& x0, double scale, int max_iter) const{
    const int n=5;
    const double alpha=1.0, gamma=2.0, rho=0.5, sigma=0.5;

    std::array<std::array<double,5>,6> X;
    X[0]=x0;
    for(int i=1;i<=n;++i){
        X[i]=x0; X[i][i-1]+=scale; // simplexe initial
    }
    auto f=[&](const std::array<double,5>& x){ return objective_from_x(x); };
    std::array<double,6> FX;
    for(int i=0;i<=n;++i) FX[i]=f(X[i]);

    for(int it=0; it<max_iter; ++it){
        // tri par valeur
        std::array<int,6> idx{0,1,2,3,4,5};
        std::sort(idx.begin(), idx.end(), [&](int a,int b){return FX[a]<FX[b];});
        // centre des meilleurs (sauf le pire)
        std::array<double,5> xc{}; 
        for(int k=0;k<n;++k) for(int j=0;j<n;++j) xc[j]+=X[idx[k]][j];
        for(int j=0;j<n;++j) xc[j]/=n;

        // réflexion
    std::array<double,5> xr;
    for(int j=0;j<n;++j)
        xr[j] = xc[j] + alpha * (xc[j] - X[idx[n]][j]);
    double fr=f(xr);

    if (fr<FX[idx[0]]){
        // expansion
        std::array<double,5> xe;
        for(int j=0;j<n;++j)
            xe[j] = xc[j] + gamma * (xr[j] - xc[j]);
        double fe=f(xe);
        if (fe<fr){ X[idx[n]]=xe; FX[idx[n]]=fe; }
        else       { X[idx[n]]=xr; FX[idx[n]]=fr; }
    } else if (fr<FX[idx[n-1]]){
        X[idx[n]]=xr; FX[idx[n]]=fr;
    } else {
        // contraction
        std::array<double,5> xk;
        if (fr<FX[idx[n]]){ // outside
            for(int j=0;j<n;++j)
                xk[j] = xc[j] + rho * (xr[j] - xc[j]);
        } else {            // inside
            for(int j=0;j<n;++j)
                xk[j] = xc[j] + rho * (X[idx[n]][j] - xc[j]);
        }
        double fk=f(xk);
        if (fk<FX[idx[n]]){ X[idx[n]]=xk; FX[idx[n]]=fk; }
        else {
            // shrink
            for(int i=1;i<=n;++i){
                for(int j=0;j<n;++j)
                    X[idx[i]][j] = X[idx[0]][j] + sigma * (X[idx[i]][j] - X[idx[0]][j]);
                FX[idx[i]]=f(X[idx[i]]);
            }
        }
    }
        // critère d’arrêt simple
        double fbest=FX[idx[0]], fworst=FX[idx[n]];
        if (std::abs(fbest-fworst) < (cfg_.use_implied_vol ? cfg_.vol_bp_target : 1e-6))
            break;
    }
    // renvoie le meilleur
    int best=0; for(int i=1;i<=n;++i) if (FX[i]<FX[best]) best=i;
    return X[best];
}

std::pair<HestonParams,double> HestonCalibrator::calibrate(){
    // heuristiques d’initialisation
    // v0 : var ATM courte; theta : var LT ~ moyenne; rho négatif; sigma 0.5; kappa 2.0
    double v0   = 0.04, theta=0.04, rho=-0.5, sigma=0.5, kappa=2.0;
    if (!market_.empty()){
        // si tu as une maturité mini proche ATM, on peut raffiner ici
    }
    HestonParams p0 {kappa, theta, sigma, rho, v0};
    std::array<double,5> x0 = inv_transform(p0);

    // multi-start autour d’un nuage raisonnable
    std::mt19937 rng(cfg_.seed);
    std::normal_distribution<double> N(0.0, 0.5); // écart-type sur l’espace non borné

    double best_val = std::numeric_limits<double>::infinity();
    std::array<double,5> best_x = x0;

    for(int s=0; s<std::max(1, cfg_.multistart); ++s){
        std::array<double,5> xs = x0;
        for(double& xi: xs) xi += N(rng);

        // Nelder-Mead (si tu as LBFGS dispo, remplace ici)
        auto xl = nelder_mead(xs, /*scale*/0.5, cfg_.max_iters_local);
        double val = objective_from_x(xl);
        if (val < best_val){ best_val=val; best_x=xl; }
    }
    return { transform(best_x), best_val };
}
