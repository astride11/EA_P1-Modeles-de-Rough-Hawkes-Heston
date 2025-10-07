#include <iostream>
#include <iomanip>
#include <chrono>
#include "Lifted.hpp"
#include "Heston.hpp"

static void test_lifted_pricing() {
    std::cout << "==========================================\n";
    std::cout << "       TEST DU MODELE LIFTED HESTON\n";
    std::cout << "==========================================\n";

    const double S0 = 100.0, K = 100.0, T = 1.0, t0 = 0.0, r = 0.05, q = 0.00;
    const double kappa = 1.0, theta = 0.04, sigma = 0.2, rho = -0.7, v0 = 0.04;
    const int n = 20; const double H = 0.1; const double rn = 10.0;

    // Compare with Heston
    Heston heston_ref(+1, t0, T, S0, K, r, q, 0.0, 0.0, kappa, theta, sigma, rho, v0);
    std::cout << "Reference Heston price: " << heston_ref.prix() << "\n";

    std::cout << "\n--- TEST 1: Call ATM ---\n";
    Lifted call_atm(+1, t0, T, S0, K, r, q, 0.0, 0.0, kappa, theta, sigma, rho, v0, n, H, rn);
    const double price_call = call_atm.prix();
    std::cout << "S0=" << S0 << ", K=" << K << ", T=" << T << "\n";
    std::cout << "Params: kappa=" << kappa << ", theta=" << theta
              << ", sigma=" << sigma << ", rho=" << rho << ", v0=" << v0
              << ", n=" << n << ", H=" << H << ", rn=" << rn << "\n";
    std::cout << "Prix call: " << std::setprecision(6) << price_call << "\n";

    std::cout << "\n--- TEST 2: Put ATM ---\n";
    Lifted put_atm(-1, t0, T, S0, K, r, q, 0.0, 0.0, kappa, theta, sigma, rho, v0, n, H, rn);
    const double price_put = put_atm.prix();
    std::cout << "Prix put:  " << std::setprecision(6) << price_put << "\n";

    const double parity = price_call - price_put - S0 * std::exp(-q * T) + K * std::exp(-r * T);
    std::cout << "Parité call-put (doit ~ 0): " << std::scientific << parity << "\n";

    std::cout << "\n--- TEST 3: Call ITM (K=90) ---\n";
    Lifted call_itm(+1, t0, T, S0, 90.0, r, q, 0.0, 0.0, kappa, theta, sigma, rho, v0, n, H, rn);
    std::cout << "Prix call ITM: " << std::fixed << std::setprecision(6) << call_itm.prix() << "\n";

    std::cout << "\n--- TEST 4: Call OTM (K=110) ---\n";
    Lifted call_otm(+1, t0, T, S0, 110.0, r, q, 0.0, 0.0, kappa, theta, sigma, rho, v0, n, H, rn);
    std::cout << "Prix call OTM: " << std::fixed << std::setprecision(6) << call_otm.prix() << "\n";

    std::cout << "\n--- TEST 5: Impact (theta,v0) ---\n";
    Lifted high_vol(+1, t0, T, S0, K, r, q, 0.0, 0.0, kappa, 0.09, sigma, rho, 0.09, n, H, rn);
    Lifted low_vol (+1, t0, T, S0, K, r, q, 0.0, 0.0, kappa, 0.01, sigma, rho, 0.01, n, H, rn);
    std::cout << "Prix call (vol~30%): " << high_vol.prix() << "\n";
    std::cout << "Prix call (vol~10%): " << low_vol.prix()  << "\n";

    std::cout << "\n--- TEST 6: Impact maturité ---\n";
    Lifted short_term(+1, t0, 0.25, S0, K, r, q, 0.0, 0.0, kappa, theta, sigma, rho, v0, n, H, rn);
    Lifted long_term (+1, t0, 2.00, S0, K, r, q, 0.0, 0.0, kappa, theta, sigma, rho, v0, n, H, rn);
    std::cout << "Prix call (T=0.25): " << short_term.prix() << "\n";
    std::cout << "Prix call (T=2.00): "  << long_term.prix()  << "\n";
}

static void test_lifted_cas_limites() {
    std::cout << "\n==========================================\n";
    std::cout << "          TESTS CAS LIMITES LIFTED\n";
    std::cout << "==========================================\n";

    const int n = 10; const double H = 0.1; const double rn = 10.0;

    std::cout << "\n--- Maturite nulle (T=0) ---\n";
    Lifted zero_T(+1, 0.0, 0.0, 100.0, 100.0, 0.05, 0.0, 0.0, 0.0, 1.0, 0.04, 0.2, -0.7, 0.04, n, H, rn);
    std::cout << "Prix call (T=0): " << zero_T.prix() << "\n";

    std::cout << "\n--- Spot eleve (S0=200, K=100) ---\n";
    Lifted high_S(+1, 0.0, 1.0, 200.0, 100.0, 0.05, 0.0, 0.0, 0.0, 1.0, 0.04, 0.2, -0.7, 0.04, n, H, rn);
    std::cout << "Prix call (S0=200): " << high_S.prix() << "\n";

    std::cout << "\n--- Spot bas (S0=50, K=100) ---\n";
    Lifted low_S(+1, 0.0, 1.0, 50.0, 100.0, 0.05, 0.0, 0.0, 0.0, 1.0, 0.04, 0.2, -0.7, 0.04, n, H, rn);
    std::cout << "Prix call (S0=50): " << low_S.prix() << "\n";
}

static void test_lifted_performance() {
    std::cout << "\n==========================================\n";
    std::cout << "          TEST PERFORMANCE LIFTED\n";
    std::cout << "==========================================\n";

    Lifted opt(+1, 0.0, 1.0, 100.0, 100.0, 0.05, 0.0, 0.0, 0.0, 1.0, 0.04, 0.2, -0.7, 0.04, 10, 0.1, 10.0);

    const int N = 10;
    auto t0 = std::chrono::high_resolution_clock::now();
    double acc = 0.0;
    for (int i = 0; i < N; ++i) {
        double price = opt.prix();
        acc += price;
        std::cout << "Iteration " << i << ": price = " << price << "\n";
    }
    auto t1 = std::chrono::high_resolution_clock::now();

    auto ms = std::chrono::duration_cast<std::chrono::milliseconds>(t1 - t0).count();
    std::cout << "Temps pour " << N << " prix: " << ms << " ms\n";
    std::cout << "Temps moyen: " << (double)ms / N << " ms\n";
    std::cout << "(checksum interne) = " << acc << "\n";
}

int main() {
    try {
        std::cout << std::fixed << std::setprecision(6);
        test_lifted_pricing();
        test_lifted_cas_limites();
        test_lifted_performance();
        std::cout << "\n==========================================\n";
        std::cout << "         TESTS TERMINES AVEC SUCCES\n";
        std::cout << "==========================================\n";
    } catch (const std::exception& e) {
        std::cerr << "ERREUR: " << e.what() << "\n";
        return 1;
    }
    return 0;
}