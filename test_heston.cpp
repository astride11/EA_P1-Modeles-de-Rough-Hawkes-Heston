#include <iostream>
#include <iomanip>
#include <chrono>
#include "Option.hpp"
#include "Heston.hpp"

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

void test_heston_pricing() {
    std::cout << "==========================================" << std::endl;
    std::cout << "       TEST DU MODÈLE DE HESTON" << std::endl;
    std::cout << "==========================================" << std::endl;
    
    // Paramètres du marché
    double S0 = 100.0;      // Prix spot
    double K = 100.0;       // Strike
    double T = 1.0;         // Maturité (1 an)
    double t0 = 0.0;        // Date initiale
    double r = 0.05;        // Taux sans risque (5%)
    double q = 0.0;         // Dividende
    
    // Paramètres Heston 
    double kappa = 1.0;     // Vitesse de retour
    double theta = 0.04;    // Variance long terme (vol = 20%)
    double sigma = 0.2;     // Vol de la vol
    double rho = -0.7;      // Corrélation
    double v0 = 0.04;       // Variance initiale (vol = 20%)
    
    // Test 1: Call At-the-Money
    std::cout << "\n--- TEST 1: Call ATM ---" << std::endl;
    Heston call_atm(
        1, S0, K, T, t0, r, q, 0.0,  // phi=1 pour call
        kappa, theta, sigma, rho, v0, 0.0
    );
    
    double prix_call = call_atm.prix();
    std::cout << "S0: " << S0 << ", K: " << K << ", T: " << T << " an" << std::endl;
    std::cout << "Paramètres Heston: kappa=" << kappa << ", theta=" << theta 
              << ", sigma=" << sigma << ", rho=" << rho << ", v0=" << v0 << std::endl;
    std::cout << "Prix du call: " << std::fixed << std::setprecision(4) << prix_call << std::endl;
    
    // Test 2: Put ATM (via parité)
    std::cout << "\n--- TEST 2: Put ATM ---" << std::endl;
    Heston put_atm(
        -1, S0, K, T, t0, r, q, 0.0,  // phi=-1 pour put
        kappa, theta, sigma, rho, v0, 0.0
    );
    
    double prix_put = put_atm.prix();
    std::cout << "Prix du put: " << prix_put << std::endl;
    
    // Vérification de la parité call-put
    double parity_check = prix_call - prix_put - S0 + K * exp(-r * T);
    std::cout << "Vérification parité call-put: " << std::scientific << parity_check << std::endl;
    
    // Test 3: Call In-the-Money
    std::cout << "\n--- TEST 3: Call ITM ---" << std::endl;
    Heston call_itm(
        1, S0, 90.0, T, t0, r, q, 0.0,  // K=90 (ITM)
        kappa, theta, sigma, rho, v0, 0.0
    );
    std::cout << "Prix call ITM (K=90): " << call_itm.prix() << std::endl;
    
    // Test 4: Call Out-of-the-Money
    std::cout << "\n--- TEST 4: Call OTM ---" << std::endl;
    Heston call_otm(
        1, S0, 110.0, T, t0, r, q, 0.0,  // K=110 (OTM)
        kappa, theta, sigma, rho, v0, 0.0
    );
    std::cout << "Prix call OTM (K=110): " << call_otm.prix() << std::endl;
    
    // Test 5: Différents paramètres de volatilité
    std::cout << "\n--- TEST 5: Impact de la volatilité ---" << std::endl;
    
    // Volatilité plus élevée
    Heston high_vol(1, S0, K, T, t0, r, q, 0.0,
                   kappa, 0.09, sigma, rho, 0.09, 0.0); // vol = 30%
    std::cout << "Prix call (vol=30%): " << high_vol.prix() << std::endl;
    
    // Volatilité plus basse
    Heston low_vol(1, S0, K, T, t0, r, q, 0.0,
                  kappa, 0.01, sigma, rho, 0.01, 0.0); // vol = 10%
    std::cout << "Prix call (vol=10%): " << low_vol.prix() << std::endl;
    
    // Test 6: Maturité différente
    std::cout << "\n--- TEST 6: Impact de la maturité ---" << std::endl;
    Heston short_term(1, S0, K, 0.25, t0, r, q, 0.0,  // T=3 mois
                     kappa, theta, sigma, rho, v0, 0.0);
    std::cout << "Prix call (T=0.25 an): " << short_term.prix() << std::endl;
    
    Heston long_term(1, S0, K, 2.0, t0, r, q, 0.0,  // T=2 ans
                    kappa, theta, sigma, rho, v0, 0.0);
    std::cout << "Prix call (T=2.0 ans): " << long_term.prix() << std::endl;
}

void test_cas_limites() {
    std::cout << "\n==========================================" << std::endl;
    std::cout << "          TESTS CAS LIMITES" << std::endl;
    std::cout << "==========================================" << std::endl;
    
    // Test: Maturité nulle
    std::cout << "\n--- Maturité nulle ---" << std::endl;
    Heston zero_maturity(1, 100.0, 100.0, 0.0, 0.0, 0.05, 0.0, 0.0,
                        1.0, 0.04, 0.2, -0.7, 0.04, 0.0);
    std::cout << "Prix call (T=0): " << zero_maturity.prix() << std::endl;
    
    // Test: Prix spot très élevé
    std::cout << "\n--- Prix spot élevé ---" << std::endl;
    Heston high_spot(1, 200.0, 100.0, 1.0, 0.0, 0.05, 0.0, 0.0,
                    1.0, 0.04, 0.2, -0.7, 0.04, 0.0);
    std::cout << "Prix call (S0=200, K=100): " << high_spot.prix() << std::endl;
    
    // Test: Prix spot très bas
    std::cout << "\n--- Prix spot bas ---" << std::endl;
    Heston low_spot(1, 50.0, 100.0, 1.0, 0.0, 0.05, 0.0, 0.0,
                   1.0, 0.04, 0.2, -0.7, 0.04, 0.0);
    std::cout << "Prix call (S0=50, K=100): " << low_spot.prix() << std::endl;
}

void test_performance() {
    std::cout << "\n==========================================" << std::endl;
    std::cout << "          TEST PERFORMANCE" << std::endl;
    std::cout << "==========================================" << std::endl;
    
    Heston test_option(1, 100.0, 100.0, 1.0, 0.0, 0.05, 0.0, 0.0,
                      1.0, 0.04, 0.2, -0.7, 0.04, 0.0);
    
    // Mesure du temps de calcul
    auto start = std::chrono::high_resolution_clock::now();
    
    const int nb_calculs = 10;
    for (int i = 0; i < nb_calculs; ++i) {
        double prix = test_option.prix();
    }
    
    auto end = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(end - start);
    
    std::cout << "Temps pour " << nb_calculs << " calculs: " 
              << duration.count() << " ms" << std::endl;
    std::cout << "Temps moyen par calcul: " 
              << duration.count() / double(nb_calculs) << " ms" << std::endl;
}

int main() {
    try {
        std::cout << std::fixed << std::setprecision(6);
        
        test_heston_pricing();
        test_cas_limites();
        test_performance();
        
        std::cout << "\n==========================================" << std::endl;
        std::cout << "         TESTS TERMINÉS AVEC SUCCÈS" << std::endl;
        std::cout << "==========================================" << std::endl;
        
    } catch (const std::exception& e) {
        std::cerr << "ERREUR: " << e.what() << std::endl;
        return 1;
    }
    
    return 0;
}