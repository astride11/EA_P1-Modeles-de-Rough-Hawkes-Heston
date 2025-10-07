#include <iostream>
#include <iomanip>
#include <chrono>
#include <cmath>
#include "Option.hpp"
#include "Heston.hpp"

static void test_heston_pricing() {
    std::cout << "==========================================\n";
    std::cout << "       TEST DU MODÈLE DE HESTON\n";
    std::cout << "==========================================\n";

    // Paramètres du marché
    const double S0 = 100.0;   // Prix spot
    const double K  = 100.0;   // Strike
    const double T  = 1.0;     // Maturité (1 an)
    const double t0 = 0.0;     // Date initiale
    const double r  = 0.05;    // Taux sans risque (5%)
    const double q  = 0.0;     // Dividende

    // Paramètres Heston
    const double kappa = 1.0;  // Vitesse de retour
    const double theta = 0.04; // Variance LT (vol 20%)
    const double sigma = 0.2;  // Vol de la vol
    const double rho   = -0.7; // Corrélation
    const double v0    = 0.04; // Variance initiale (vol 20%)

    // Test 1: Call At-the-Money
    std::cout << "\n--- TEST 1: Call ATM ---\n";
    Heston call_atm(
        +1, t0, T, S0, K, r, q, /*market*/0.0, /*vega*/0.0,
        kappa, theta, sigma, rho, v0
    );
    const double prix_call = call_atm.prix();
    std::cout << "S0=" << S0 << ", K=" << K << ", T=" << T << " an\n";
    std::cout << "Paramètres Heston: kappa=" << kappa << ", theta=" << theta
              << ", sigma=" << sigma << ", rho=" << rho << ", v0=" << v0 << "\n";
    std::cout << "Prix du call: " << std::fixed << std::setprecision(6) << prix_call << "\n";

    // Test 2: Put ATM (via parité)
    std::cout << "\n--- TEST 2: Put ATM ---\n";
    Heston put_atm(
        -1, t0, T, S0, K, r, q, /*market*/0.0, /*vega*/0.0,
        kappa, theta, sigma, rho, v0
    );
    const double prix_put = put_atm.prix();
    std::cout << "Prix du put: " << std::fixed << std::setprecision(6) << prix_put << "\n";

    // Vérification de la parité call-put (avec dividende q)
    const double parity_check = prix_call - prix_put - S0 * std::exp(-q * T) + K * std::exp(-r * T);
    std::cout << "Parité call-put (doit ~ 0): " << std::scientific << parity_check << "\n";

    // Test 3: Call In-the-Money (K=90)
    std::cout << "\n--- TEST 3: Call ITM (K=90) ---\n";
    Heston call_itm(
        +1, t0, T, S0, 90.0, r, q, 0.0, 0.0,
        kappa, theta, sigma, rho, v0
    );
    std::cout << "Prix call ITM (K=90): " << std::fixed << std::setprecision(6) << call_itm.prix() << "\n";

    // Test 4: Call Out-of-the-Money (K=110)
    std::cout << "\n--- TEST 4: Call OTM (K=110) ---\n";
    Heston call_otm(
        +1, t0, T, S0, 110.0, r, q, 0.0, 0.0,
        kappa, theta, sigma, rho, v0
    );
    std::cout << "Prix call OTM (K=110): " << std::fixed << std::setprecision(6) << call_otm.prix() << "\n";

    // Test 5: Impact de la volatilité (theta,v0)
    std::cout << "\n--- TEST 5: Impact de la volatilité ---\n";
    // Volatilité plus élevée (≈30%)
    Heston high_vol(+1, t0, T, S0, K, r, q, 0.0, 0.0,
                    kappa, 0.09, sigma, rho, 0.09);
    std::cout << "Prix call (vol~30%): " << std::fixed << std::setprecision(6) << high_vol.prix() << "\n";
    // Volatilité plus basse (≈10%)
    Heston low_vol(+1, t0, T, S0, K, r, q, 0.0, 0.0,
                   kappa, 0.01, sigma, rho, 0.01);
    std::cout << "Prix call (vol~10%): " << std::fixed << std::setprecision(6) << low_vol.prix() << "\n";

    // Test 6: Impact de la maturité
    std::cout << "\n--- TEST 6: Impact de la maturité ---\n";
    Heston short_term(+1, t0, 0.25, S0, K, r, q, 0.0, 0.0,
                      kappa, theta, sigma, rho, v0);
    Heston long_term (+1, t0, 2.00, S0, K, r, q, 0.0, 0.0,
                      kappa, theta, sigma, rho, v0);
    std::cout << "Prix call (T=0.25): " << short_term.prix() << "\n";
    std::cout << "Prix call (T=2.00): " << long_term.prix()  << "\n";
}

static void test_cas_limites() {
    std::cout << "\n==========================================\n";
    std::cout << "          TESTS CAS LIMITES\n";
    std::cout << "==========================================\n";

    // Maturité nulle (T=0)
    std::cout << "\n--- Maturité nulle (T=0) ---\n";
    Heston zero_maturity(+1, 0.0, 0.0, 100.0, 100.0, 0.05, 0.0, 0.0, 0.0,
                         1.0, 0.04, 0.2, -0.7, 0.04);
    std::cout << "Prix call (T=0): " << zero_maturity.prix() << "\n";

    // Spot élevé
    std::cout << "\n--- Spot élevé (S0=200, K=100) ---\n";
    Heston high_spot(+1, 0.0, 1.0, 200.0, 100.0, 0.05, 0.0, 0.0, 0.0,
                     1.0, 0.04, 0.2, -0.7, 0.04);
    std::cout << "Prix call (S0=200): " << high_spot.prix() << "\n";

    // Spot bas
    std::cout << "\n--- Spot bas (S0=50, K=100) ---\n";
    Heston low_spot(+1, 0.0, 1.0, 50.0, 100.0, 0.05, 0.0, 0.0, 0.0,
                    1.0, 0.04, 0.2, -0.7, 0.04);
    std::cout << "Prix call (S0=50): " << low_spot.prix() << "\n";
}

static void test_performance() {
    std::cout << "\n==========================================\n";
    std::cout << "          TEST PERFORMANCE\n";
    std::cout << "==========================================\n";

    Heston test_option(+1, 0.0, 1.0, 100.0, 100.0, 0.05, 0.0, 0.0, 0.0,
                       1.0, 0.04, 0.2, -0.7, 0.04);

    // Mesure du temps de calcul
    auto start = std::chrono::high_resolution_clock::now();

    const int nb_calculs = 10;
    double acc = 0.0;
    for (int i = 0; i < nb_calculs; ++i) {
        acc += test_option.prix();
    }

    auto end = std::chrono::high_resolution_clock::now();
    auto duration_ms = std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count();

    std::cout << "Temps pour " << nb_calculs << " calculs: " << duration_ms << " ms\n";
    std::cout << "Temps moyen par calcul: " << (double)duration_ms / nb_calculs << " ms\n";
    std::cout << "(checksum) = " << acc << "\n";
}

int main() {
    try {
        std::cout << std::fixed << std::setprecision(6);
        test_heston_pricing();
        test_cas_limites();
        test_performance();
        std::cout << "\n==========================================\n";
        std::cout << "         TESTS TERMINÉS AVEC SUCCÈS\n";
        std::cout << "==========================================\n";
    } catch (const std::exception& e) {
        std::cerr << "ERREUR: " << e.what() << "\n";
        return 1;
    }
    return 0;
}
