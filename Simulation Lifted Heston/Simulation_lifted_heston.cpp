#include <random>
#include <iostream>
#include <fstream>
#include "Eigen/Dense"

using namespace std;

class Simulation_lifted_heston {
public:
    double kappa;   // Vitesse de réversion
    double theta;   // Variance à long terme
    double sigma;   // Volatilité de la variance
    double rho;     // Corrélation entre les deux processus de Wiener
    double V0;      // Variance initiale 
    double H;       // Paramètre de Hurst
    int N;          // Nombre de facteurs
    double S0;      // Prix initial de l'actif sous-jacent
    double r;      // Taux d'intérêt sans risque
    double T;   // Horizon de temps
    double alpha;
    double rn= 2.5;
    double N_steps=1000;
    double eps = 1e-8; // Facteur d'espacement

    Simulation_lifted_heston(double kappa_, double theta_, double sigma_, double rho_, double V0_, double H_, int N_, double S0_, double r_, double T_)
        : kappa(kappa_), theta(theta_), sigma(sigma_), rho(rho_), V0(V0_), H(H_), N(N_), S0(S0_), r(r_), T(T_), alpha(H_ + 0.5) {}

struct Result {
    double g0;
    std::vector<double> c;
    std::vector<double> x;
};

// Génération de deux variables aléatoires normales corrélées
Eigen::Vector2d generateZ() {
    Eigen::Matrix2d Sigma;
    Sigma << 1, rho,
             rho, 1;
    static std::mt19937 gen(42); // Graine fixe pour la reproductibilité
    std::normal_distribution<> dist(0.0, 1.0);
    Eigen::Vector2d Z;
    Z << dist(gen), dist(gen);
    Eigen::Matrix2d L = Sigma.llt().matrixL();
    Eigen::Vector2d Z_final = L * Z;
    return Z_final;
}

// calcule de la fonction g et de c et x
Result g0(double t) {
    Result result;
    result.c.resize(N);
    result.x.resize(N);
    for (int i = 0; i < N; ++i) {
        result.c[i] = ((pow(rn, 1 - alpha) - 1) * pow(rn, (alpha - 1) * (1 + N/2.0)))* pow(rn, (1 - alpha) * (i+1)) / 
               (tgamma(alpha) * tgamma(2 - alpha)) ;

        result.x[i] = (((1 - alpha) * (pow(rn, 2 - alpha) - 1))* pow(rn, i - N/2.0) / ((2 - alpha) * (pow(rn, 1 - alpha) - 1))); 
    }
    double sum = 0.0;
    for (int i = 0; i < N; ++i) {
        sum += result.c[i] * (1 - exp(-result.x[i] * t)) / result.x[i];
    }
    result.g0 = V0 + kappa * theta * sum;
    return result;
}

//fonction qui simule le modèle 
void Simulate(std::vector<double>& S, std::vector<double>& V) {
        double dt = T/N_steps;
        S.resize(N_steps+1);
        V.resize(N_steps+1);
        std::vector<std::vector<double>> U(N, std::vector<double>(N_steps+1, 0.0));
        S[0] = S0;
        V[0] = V0;

        for (int j = 1; j <= N_steps; ++j) {
            Result res = g0(j*dt);
            double Vplus = std::max(V[j-1], 0.0);
            Eigen::Vector2d Z = generateZ();

            // Mise à jour du prix S
            S[j] = S[j-1] * std::exp((r - 0.5 * Vplus) * dt + std::sqrt(Vplus * dt) * Z(0));

            // Mise à jour des facteurs U
            for (int i = 0; i < N; ++i) {
            double x_i = std::max(abs(res.x[i]), eps);
            U[i][j] = (U[i][j-1] - kappa * Vplus * dt + sigma * std::sqrt(Vplus * dt) * Z(1)) / 
                      (1 + x_i * dt);
            }
            // Mise à jour de la variance V
            double sumU = 0.0;
            for (int i = 0; i < N; ++i) {
                sumU += res.c[i] * U[i][j];
            };
            V[j] = std::max(res.g0 + sumU, eps);

        }
    }
};
int main() {
     std::vector<int> Nvals = {10, 20, 50, 100, 200, 500};
    for (int idx = 0; idx < Nvals.size(); ++idx) {
        int N = Nvals[idx];
        Simulation_lifted_heston sim(0.3, 0.05, 0.1, -0.5, 0.05, 0.1, N, 100.0, 0.01, 1.0);
        std::vector<double> S, V;
        sim.Simulate(S, V);

        std::ofstream file("simulation_N" + std::to_string(N) + ".csv");
        file << "S,V\n";
        for (size_t j = 0; j < S.size(); ++j) {
            file << S[j] << "," << V[j] << "\n";
        }
        file.close();
    }
    std::cout << "Fichiers CSV créés pour chaque N." << std::endl;
    return 0;
}