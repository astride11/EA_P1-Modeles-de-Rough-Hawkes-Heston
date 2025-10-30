#include <iostream>
#include <vector>
#include <memory>
#include <cmath>
#include <random>
#include <iomanip>
#include <fstream>
#include <sstream>
#include <string>


#include "../../Option.hpp"
#include "../Heston.hpp"
#include "../CalibH.hpp" 

struct TrueParams {
    double kappa, theta, sigma, rho, v0;
};
 
static void build_synthetic_calls(
    std::vector<std::unique_ptr<Option>>& storage,
    std::vector<Option*>& market,
    double S0, double r, double q,
    const std::vector<double>& maturities,
    const std::vector<double>& moneyness,   // ex: {-0.2,-0.1,0,0.1,0.2} pour K = F*exp(m)
    const TrueParams& tp,
    double price_noise_bps = 0.0,           // bruit additif en bps *du prix* (0 = noiseless)
    unsigned seed = 123
){
    std::mt19937 rng(seed);
    std::normal_distribution<double> N01(0.0, 1.0);

    for (double T : maturities){
        double F = S0 * std::exp((r - q) * T);
        for (double m : moneyness){
            double K = F * std::exp(m);

            // On crée un Heston "temporaire" pour pricer le prix de marché synthétique
            Heston h(+1, S0, K, T, /*t0*/0.0, r, q, /*market_price*/0.0,
                     tp.kappa, tp.theta, tp.sigma, tp.rho, tp.v0);

            double model_price = h.prix_call(); // prix call modèle

            // Ajoute un (petit) bruit additif si demandé (en bps du prix)
            if (price_noise_bps > 0.0 && model_price > 0.0){
                double sd = (price_noise_bps * 1e-4) * model_price; // bps du prix
                model_price = std::max(0.0, model_price + sd * N01(rng));
            }

            // Stocke une Option/Heston avec market_price = model_price
            storage.emplace_back(std::make_unique<Heston>(+1, S0, K, T, 0.0, r, q, model_price,
                                                          /*kappa*/1.0, /*theta*/0.04, /*sigma*/0.5, /*rho*/-0.5, /*v0*/0.04));
            market.push_back(storage.back().get());
        }
    }
}

static void print_params(const char* title, double kappa, double theta, double sigma, double rho, double v0){
    std::cout << std::fixed << std::setprecision(6);
    std::cout << title
              << "kappa=" << kappa << "  "
              << "theta=" << theta << "  "
              << "sigma=" << sigma << "  "
              << "rho="   << rho   << "  "
              << "v0="    << v0    << "\n";
}

int main(){
    // renomer en main() pour exécuter ce test
    // Test de calibration sur une surface synthétique

    // ---------- 1) “Vrais” paramètres ----------
    TrueParams tp { /*kappa*/1.80, /*theta*/0.035, /*sigma*/0.30, /*rho*/-0.55, /*v0*/0.030 };

    // ---------- 2) Marché (synthétique) ----------
    double S0 = 100.0;
    double r  = 0.02;
    double q  = 0.00;

    // Grille de maturités et moneyness
    std::vector<double> maturities {0.50, 1.00, 1.5, 2.0};
    std::vector<double> moneyness  { -0.4, -0.30, -0.15, -0.01, -0.02, -0.03, -0.05, -0.07, -0.1, 0.0, 0.01, 0.02, 0.03, 0.05, 0.07, 0.15, 0.30, 0.4}; // m = log(K/F)

    // Construis les options synthétiques
    std::vector<std::unique_ptr<Option>> storage;
    std::vector<Option*> market;

    // -> commence “sans bruit” pour un test de récupération exact
    double noise_bps = 0.0; // mets 2.0 ou 5.0 pour tester la robustesse
    build_synthetic_calls(storage, market, S0, r, q, maturities, moneyness, tp, noise_bps);

    // ---------- 3) Calibration ----------
    CalibConfig cfg;
    cfg.use_implied_vol = false;   // on calibre sur PRIX (plus robuste; ici noiseless, c’est parfait)
    cfg.vega_weights    = true;    // bon équilibrage strikes/maturités
    cfg.multistart      = 10;
    cfg.max_iters_local = 600;
    cfg.feller_penalty  = 50.0;    // pénalité douce

    HestonCalibrator calib(market, cfg);
    auto [p, err] = calib.calibrate(0.04, 0.04, -0.7, 0.5, 2.0);
    /*
        // -> ajout d'un bruit 
    double noise_bps2 = 2.0; 
    build_synthetic_calls(storage, market, S0, r, q, maturities, moneyness, tp, noise_bps2);

    // ---------- 3) Calibration ----------
    CalibConfig cfg2;
    cfg2.use_implied_vol = false;   // on calibre sur PRIX (plus robuste; ici noiseless, c’est parfait)
    cfg2.vega_weights    = true;    // bon équilibrage strikes/maturités
    cfg2.multistart      = 10;
    cfg2.max_iters_local = 600;
    cfg2.feller_penalty  = 50.0;    // pénalité douce

    HestonCalibrator calib2(market, cfg);
    auto [p2, err2] = calib2.calibrate(0.04, 0.04, -0.7, 0.5, 2.0);
    */

    // ---------- 4) Résultats ----------

    std::cout << "\n=== Test de récupération sur surface synthétique ===\n";
    print_params("VRAI     : ", tp.kappa, tp.theta, tp.sigma, tp.rho, tp.v0);
    print_params("ESTIME   : ", p.kappa,  p.theta,  p.sigma,  p.rho,  p.v0);

    /*
    std::cout << "\n=== Test de récupération sur surface synthétique avec " << noise_bps2 << " ===\n";
    print_params("VRAI     : ", tp.kappa, tp.theta, tp.sigma, tp.rho, tp.v0);
    print_params("ESTIME   : ", p2.kappa,  p2.theta,  p2.sigma,  p2.rho,  p2.v0);
    */

    // Export des smiles “vrai” vs “modèle calibré” dans un CSV


    // ---------- 4) Export du smile calibré (marché vs modèle) ----------
    std::ofstream out("Heston_Model/test_Heston/smile_data_calib.csv");
    out << "T,K,moneyness,iv_market,iv_model\n";

    for (auto* opt : market) {
        double T = opt->getT();
        double K = opt->getK();
        double S0 = opt->getS0();
        double r = opt->getR();
        double q = opt->getQ();

        // --- Prix marché et prix modèle ---
        double price_market = opt->getMarketPrice(); // ton "mid" lu dans le CSV
        double price_model = Heston(+1, S0, K, T, 0.0, r, q, 0.0,
                                    p.kappa, p.theta, p.sigma, p.rho, p.v0).prix_call();

        // --- Calcul des vols implicites ---
        double df = std::exp(-r * T);
        double Fwd = S0 * std::exp((r - q) * T);
        double iv_market = implied_vol_bs(true, Fwd, K, T, price_market, df);
        double iv_model  = implied_vol_bs(true, Fwd, K, T, price_model, df);

        // --- Moneyness (log(K/F)) ---
        double m = std::log(K / Fwd);

        // --- Sauvegarde dans le CSV ---
        out << T << "," << K << "," << m << ","
            << iv_market << "," << iv_model << "\n";
    }

    out.close();
    std::cout << "\n✅ Fichier 'Heston_Model/test_Heston/smile_data_calib.csv' exporté (marché vs modèle).\n";
        // std::cout << "\n✅ Fichier 'smile_data_pert.csv' exporté.\n";

    return 0;
}





int func() {
    // ---------- 1) Lecture des données réelles ----------
    std::vector<std::unique_ptr<Option>> storage;
    std::vector<Option*> market;

    std::ifstream file("Heston_Model/test_Heston/AAPL_calls_filtered.csv");
    std::string line;
    std::getline(file, line); // ignore header

    while (std::getline(file, line)) {
        std::stringstream ss(line);
        std::string contractSymbol, expiry_str;
        double S0, strike, mid, T;

        std::getline(ss, contractSymbol, ',');
        ss >> S0; ss.ignore(1);
        std::getline(ss, expiry_str, ',');
        T = std::stod(expiry_str);
        ss >> strike; ss.ignore(1);
        ss >> mid;

        storage.emplace_back(std::make_unique<Heston>(
            +1, S0, strike, T, 0.0,
            /*r*/0.02, /*q*/0.0, mid,
            /*kappa*/1.0, /*theta*/0.04, /*sigma*/0.5, /*rho*/-0.5, /*v0*/0.04
        ));
        market.push_back(storage.back().get());
    }

    std::cout << "Loaded " << market.size() << " options.\n";

    // ---------- 2) Calibration ----------
    CalibConfig cfg;
    cfg.use_implied_vol = true;
    cfg.vega_weights    = true;
    cfg.multistart      = 10;
    cfg.max_iters_local = 1000;
    cfg.feller_penalty  = 50.0; 

    HestonCalibrator calib(market, cfg);
    auto [p, err] = calib.calibrate(0.04, 0.04, -0.6, 0.5, 1.8);

    // ---------- 3) Résultats ----------
    std::cout << "\n=== Calibration sur données réelles ===\n";
    print_params("Paramètres estimés : ", p.kappa, p.theta, p.sigma, p.rho, p.v0);
    std::cout << "Erreur de calibration : " << err << "\n";

    // ---------- 4) Export du smile calibré (marché vs modèle) ----------
    std::ofstream out("Heston_Model/test_Heston/smile_data_calib.csv");
    out << "T,K,iv_market,iv_model\n";

    for (auto* opt : market) {
        double T = opt->getT();
        double K = opt->getK();
        double S0 = opt->getS0();
        double r = opt->getR();
        double q = opt->getQ();

        // --- Prix marché et prix modèle ---
        double price_market = opt->getMarketPrice(); //  "mid" lu dans le CSV
        double price_model = Heston(+1, S0, K, T, 0.0, r, q, 0.0,
                                    p.kappa, p.theta, p.sigma, p.rho, p.v0).prix_call();

        // --- Calcul des vols implicites ---
        double df = std::exp(-r * T);
        double Fwd = S0 * std::exp((r - q) * T);
        double iv_market = implied_vol_bs(true, Fwd, K, T, price_market, df);
        double iv_model  = implied_vol_bs(true, Fwd, K, T, price_model, df);

        // --- Sauvegarde dans le CSV ---
        out << T << "," << K << ","
            << iv_market << "," << iv_model << "\n";
    }

    out.close();
    std::cout << "\n✅ Fichier 'Heston_Model/test_Heston/smile_data_calib.csv' exporté (marché vs modèle).\n";

    return 0;
}