#include "LiftedH.hpp"
#include <random>
#include <vector>
#include <fstream>
#include <sstream>
#include <iostream>
#include <iomanip>
#include <cmath>
#include <algorithm>
#include <string>

using namespace std;

// ===== Constructeur par défaut =====
LiftedHeston::LiftedHeston()
    /*
    Nous initialisons les attributs avec des valeurs par defaut.
    Par defaut ce constructeur cree un Call avec les parametres suivants :
    */
    : Option(), kappa(0.0), theta(0.0), sigma(0.0), rho(0.0), v0(0.0), H(0.1), N(1), alpha(H + 0.5), rn(1 + 10 * std::pow(N, -0.9)) {}

// ===== Constructeur paramétré =====
LiftedHeston::LiftedHeston(int phi_, double S0_, double K_, double T_, double t0_,
                           double r_, double q_, double market_price_,
                           double kappa_, double theta_, double sigma_, double rho_, double v0_, double H_, int N_,
                           double vega_, double rn_)
    /*
    Nous initialisons les attributs avec les valeurs passees en parametres.
    */
    : Option(phi_, S0_, K_, T_, t0_, r_, q_, market_price_, vega_),
      kappa(kappa_), theta(theta_), sigma(sigma_), rho(rho_), v0(v0_), H(H_), N(N_), alpha(H_ + 0.5), rn(rn_) {}    

LiftedHeston::~LiftedHeston() {
    // Destructeur
    // Pas de ressources dynamiques à libérer dans cette classe de base
}

LiftedHeston::Param_sup LiftedHeston::fix_params(int N, double H, double alpha, double rn) const {
    Param_sup params;
    params.c.resize(N);
    params.x.resize(N);
    for (int i = 0; i < N; ++i) {
        params.c[i] = ((pow(rn, 1 - alpha) - 1) * pow(rn, (alpha - 1) * (1 + N/2.0)))* pow(rn, (1 - alpha) * (i+1)) / 
               (tgamma(alpha) * tgamma(2 - alpha)) ;

        params.x[i] = (((1 - alpha) * (pow(rn, 2 - alpha) - 1))* pow(rn, i - N/2.0) / ((2 - alpha) * (pow(rn, 1 - alpha) - 1))); 
    }
    return params;
}

double LiftedHeston::g0(double t, Param_sup params) const {
    double sum = 0.0;
    for (int i = 0; i < N; ++i) {
        sum += params.c[i] * (1 - exp(-params.x[i] * t)) / params.x[i];
    }
    return v0 + kappa * theta * sum;
}

void LiftedHeston::simulate_paths(int N_steps, double T_hori, std::vector<double>& S_paths, std::vector<double>& V_paths,
                                  int N_paths, double kappa_, double theta_, double sigma_, double rho_,
                                  double v0_, double S0_, double H_, int N_, double rn_, double r_){

    double dt = T_hori / N_steps;
    S_paths[0] = S0_;
    V_paths[0] = v0_;

    Param_sup params = fix_params(N_, H_, H_ + 0.5, rn_);
    std::vector<std::vector<double>> U(N_, std::vector<double>(N_steps + 1, 0.0));

    std::mt19937 gen(42);
    std::normal_distribution<> dist(0.0, 1.0);

    for (int i = 1; i <= N_steps; ++i) {
        double t = i * dt;
        double g0_t = g0(t, params);
        double Vplus = std::max(V_paths[i - 1], 0.0);

        double Z1 = dist(gen);
        double Z2 = dist(gen);

        double X1 = Z1;
        double X2 = rho_ * Z1 + std::sqrt(1 - rho_ * rho_) * Z2;

        // Mise à jour des processus auxiliaires U_j
        for (int j = 0; j < N_; ++j) {
            U[j][i] = (U[j][i - 1] - kappa_ * U[j][i - 1] * dt + sigma_ * std::sqrt(Vplus * dt) * X2)
                       / (1 + params.x[j] * dt);
        }

        // Mise à jour de V
        double sumU = 0.0;
        for (int j = 0; j < N_; ++j)
            sumU += params.c[j] * U[j][i];
        double V_new = std::max(g0_t + sumU, 1e-8);
        V_paths[i] = V_new;

        // Mise à jour de S (Euler en log)
        S_paths[i] = S_paths[i - 1] * std::exp((r - 0.5 * Vplus) * dt + std::sqrt(Vplus * dt) * X1);
    }
}

// Simuler une trajectoire par niveau N (comparaison visuelle de la convergence)
void LiftedHeston::simulate_compare_N(
    const std::vector<int> &N_levels,
    const std::vector<double> &rn_levels,
    int N_steps,
    double T_hori,
    double kappa_, double theta_, double sigma_, double rho_,
    double v0_, double S0_,
    double r_, double q_,
    unsigned int seed)
{
    if (N_levels.size() != rn_levels.size()) {
        std::cerr << "N_levels and rn_levels must have same length\n";
        return;
    }
    const size_t L = N_levels.size();
    double dt = T_hori / N_steps;
    std::vector<double> times(N_steps + 1);
    for (int i = 0; i <= N_steps; ++i) times[i] = i * dt;

    // RNG et tirages Z1,Z2 (une seule suite de bruits pour tous les modèles)
    std::mt19937 gen(seed);
    std::normal_distribution<> dist(0.0, 1.0);
    std::vector<double> Z1(N_steps), Z2(N_steps);
    for (int k = 0; k < N_steps; ++k) {
        Z1[k] = dist(gen);
        Z2[k] = dist(gen);
    }
    const double sqrt_rho_term = 0.0; // not used, compute per step if needed

    // Préparer données par niveau
    std::vector<Param_sup> params_list(L);
    std::vector<std::vector<double>> U_list(L); // U_list[l] size = Nf for level l
    std::vector<std::vector<double>> S_cols(L, std::vector<double>(N_steps + 1));
    std::vector<std::vector<double>> V_cols(L, std::vector<double>(N_steps + 1));
    std::vector<std::string> col_names(L);

    for (size_t l = 0; l < L; ++l) {
        int Nf = N_levels[l];
        double rn_ = rn_levels[l];
        // alpha local = H + 0.5 (si tu veux utiliser member H, passe-le en paramètre)
        double alpha_local = this->alpha; // ou calcule H + 0.5 si besoin
        params_list[l] = fix_params(Nf, this->H, alpha_local, rn_);
        U_list[l].assign(Nf, 0.0);
        S_cols[l][0] = S0_;
        V_cols[l][0] = v0_;
        std::ostringstream ss; ss << "N=" << Nf;
        col_names[l] = ss.str();
    }

    const double eps = 1e-12;

    // boucle temps (une seule trajectoire par niveau)
    for (int k = 1; k <= N_steps; ++k) {
        double z1 = Z1[k-1];
        double z2 = Z2[k-1];
        double X1 = z1;
        double X2 = rho_ * z1 + std::sqrt(std::max(0.0, 1.0 - rho_ * rho_)) * z2;
        double t = times[k];

        for (size_t l = 0; l < L; ++l) {
            int Nf = N_levels[l];
            LiftedHeston::Param_sup &params = params_list[l];
            std::vector<double> &U = U_list[l];

            double Vprev = std::max(V_cols[l][k-1], 0.0);
            double sqrtVdt = std::sqrt(Vprev * dt);

            // mise à jour des U_j (semi-implicite comme dans ton schéma)
            for (int j = 0; j < Nf; ++j) {
                double x_j = params.x[j];
                double Uprev = U[j];
                double numerator = Uprev - kappa_ * Uprev * dt + sigma_ * sqrtVdt * X2;
                double denom = 1.0 + x_j * dt;
                U[j] = numerator / denom;
            }

            // reconstruction V: g0(t) + sum c_j * U_j
            double sumU = 0.0;
            for (int j = 0; j < Nf; ++j) sumU += params.c[j] * U[j];

            // calcule g0 localement (sécurité x_j ~ 0)
            double sum_g = 0.0;
            for (int j = 0; j < Nf; ++j) {
                double x = params.x[j];
                if (std::abs(x) < eps) sum_g += params.c[j] * t;
                else sum_g += params.c[j] * (1.0 - std::exp(-x * t)) / x;
            }
            double g0_t = v0_ + kappa_ * theta_ * sum_g;

            double Vnew = std::max(g0_t + sumU, 1e-12);
            V_cols[l][k] = Vnew;

            // mise à jour S (Euler log)
            double Sprev = S_cols[l][k-1];
            double dW_s = std::sqrt(std::max(0.0, Vprev) * dt) * X1;
            S_cols[l][k] = Sprev * std::exp((r_ - q_ - 0.5 * Vprev) * dt + dW_s);
        } // end loop levels
    } // end loop time

    // Construire matrices "row = time" pour écrire CSV
    std::vector<std::vector<double>> S_mat(N_steps + 1, std::vector<double>(L));
    std::vector<std::vector<double>> V_mat(N_steps + 1, std::vector<double>(L));
    for (int k = 0; k <= N_steps; ++k) {
        for (size_t l = 0; l < L; ++l) {
            S_mat[k][l] = S_cols[l][k];
            V_mat[k][l] = V_cols[l][k];
        }
    }

    // Ecrivons deux CSV: time, N=..., N=..., ...
    auto write_named_csv = [&](const std::string &filename,
                               const std::vector<double> &times,
                               const std::vector<std::string> &col_names,
                               const std::vector<std::vector<double>> &mat) {
        std::ofstream ofs(filename);
        ofs << "time";
        for (const auto &cn : col_names) ofs << "," << cn;
        ofs << "\n";
        for (size_t i = 0; i < times.size(); ++i) {
            ofs << std::setprecision(12) << times[i];
            for (size_t j = 0; j < col_names.size(); ++j) ofs << "," << std::setprecision(12) << mat[i][j];
            ofs << "\n";
        }
        ofs.close();
    };

    write_named_csv("S_compare_N.csv", times, col_names, S_mat);
    write_named_csv("V_compare_N.csv", times, col_names, V_mat);

    std::cout << "Ecrit fichiers : S_compare_N.csv et V_compare_N.csv\n";
}

int main() {
    // paramètres
    double kappa = 1.5, theta = 0.04, sigma = 0.3, rho = -0.7;
    double v0 = 0.04, S0 = 100.0, r = 0.0, q = 0.0;
    double H = 0.1; // utilisé par fix_params via this->H
    int N_steps = 500;

    std::vector<int> N_levels = {2, 20, 100, 500, 1000};
    std::vector<double> rn_levels;
    for (int N : N_levels) rn_levels.push_back(1 + 10 * std::pow(N, -0.9));

    LiftedHeston model; 
    model.simulate_compare_N(N_levels, rn_levels, N_steps, 1.0,
                             kappa, theta, sigma, rho,
                             v0, S0, r, q, 42);
    return 0;
}
