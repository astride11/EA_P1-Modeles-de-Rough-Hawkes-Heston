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

static inline cd F_fun(cd u, cd v, double kappa, double sigma, double rho) {
    return 0.5 * (u * u - u)
         + (rho * sigma * u - kappa) * v
         + 0.5 * sigma * sigma * v * v;
}
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
    
    double gamma_alpha = tgamma(alpha);
    double gamma_2_minus_alpha = tgamma(2.0 - alpha);
    
    for (int i = 0; i < N; ++i) {
        // Correction du calcul de c[i]
        double term1 = pow(rn, 1.0 - alpha) - 1.0;
        double term2 = pow(rn, (alpha - 1.0) * (1.0 + N/2.0));
        double term3 = pow(rn, (1.0 - alpha) * (i + 1));
        
        params.c[i] = (term1 * term2 * term3) / (gamma_alpha * gamma_2_minus_alpha);
        
        // Correction du calcul de x[i]
        double num = (1.0 - alpha) * (pow(rn, 2.0 - alpha) - 1.0);
        double denom = (2.0 - alpha) * (pow(rn, 1.0 - alpha) - 1.0);
        params.x[i] = (num / denom) * pow(rn, i - N/2.0);
        
        // Protection contre les valeurs aberrantes
        if (std::isnan(params.c[i]) || std::isinf(params.c[i])) {
            params.c[i] = 0.0;
        }
        if (std::isnan(params.x[i]) || std::isinf(params.x[i])) {
            params.x[i] = 1.0; // valeur par défaut
        }
    }
    return params;
}

double LiftedHeston::g0(double t, Param_sup params) const {
    double sum = 0.0;
    const double eps = 1e-12;
    
    for (int i = 0; i < N; ++i) {
        double x_i = params.x[i];
        if (std::abs(x_i) < eps) {
            // Développement de Taylor pour x petit
            sum += params.c[i] * (t - 0.5 * x_i * t * t);
        } else {
            sum += params.c[i] * (1.0 - exp(-x_i * t)) / x_i;
        }
    }
    
    double result = v0 + kappa * theta * sum;
    return std::max(result, 1e-8); // Éviter les valeurs négatives
}
void LiftedHeston::simulate_paths(int N_steps, double T_hori, std::vector<double>& S_paths, std::vector<double>& V_paths){


    double S0_ = this->S0;
    double v0_ = this->v0;
    double kappa_ = this->kappa;
    double theta_ = this->theta;
    double sigma_ = this->sigma;
    double rho_ = this->rho;
    double H_ = this->H;
    int N_ = this->N;
    double rn_ = this->rn;
    double r = this->r;
    double q = this->q;
    
    double dt = T_hori / N_steps;
    S_paths[0] = S0_;
    V_paths[0] = v0_;

    Param_sup params = fix_params(N_, H_, H_ + 0.5, rn_);
    std::vector<std::vector<double>> U(N_, std::vector<double>(N_steps + 1, 0.0));

    std::random_device rd;
    std::mt19937 gen(rd());
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

double LiftedHeston::prix_mc(int N_steps, int N_paths, double T_hori) {

    double K = this->K;
    double phi = this->phi;
    double r = this->r;
    double q = this->q;
    
    std::vector<double> S_paths(N_steps + 1);
    std::vector<double> V_paths(N_steps + 1);
    double payoff_sum = 0.0;

    for (int i = 0; i < N_paths; ++i) {
        simulate_paths(N_steps, T_hori, S_paths, V_paths);
        double payoff = std::max(phi * (S_paths.back() - K), 0.0);
        payoff_sum += payoff;
    }

    double discount_factor = std::exp(-r * T_hori);
    return (payoff_sum / N_paths) * discount_factor;
}

// ======== Pricing avec Cosinus ========

double LiftedHeston::Ksi(double c, double d, int k, double a, double b) const{

    if (k == 0) {
        return std::exp(d) - std::exp(c);
    }

    double delta1 = b - a; 
    double const1 = (k * M_1_PI) / delta1;
    double denum = 1 + const1 * const1;
    cout << "denum: " << denum << endl;

    cout << "le pb c'est exp b"<< std::exp(b) << endl;
    cout << "je pense y a pas de suci avec exp c"<< std::exp(c) << endl;
    double num = std::cos(const1 * (d - a)) * std::exp(d) - std::cos(const1 * (c - a)) * std::exp(c);
    num += const1 * (std::sin(const1 * (d - a)) * std::exp(d) - std::sin(const1 * (c - a)) * std::exp(c));

    return  num / denum;
}

double LiftedHeston::Chi(double c, double d, int k, double a, double b) const{

    if (k == 0) {
        return d-c;
    }

    double delta1 = b - a; 
    double const1 = (k * M_1_PI) / delta1;
    //cout << "const1: " << const1 << endl;
    double num = std::sin(const1 * (d - a)) - std::sin(const1 * (c - a));

    return num / const1;
}

double LiftedHeston::Vk(int k, double a, double b) const{

    double K = this->K;
    //cout <<"strike: " << K << endl;

    //cout << "Calcul de Vk pour k=" << k << ", a=" << a << ", b=" << b << endl;
    double term1 = Ksi(0, b, k, a, b);
    //cout << "Term1 pour k=" << k << ": " << term1 << endl;
    double term2 = Chi(0, b, k, a, b);
    //cout << "Term2 pour k=" << k << ": " << term2 << endl;

    return (2.0 / (b - a)) * K * (term1 - term2);
}

double LiftedHeston::prix_cos(int N_cos, double L) const{
    
    double K = this->K;
    double S0 = this->S0;
    double r = this->r;
    double q = this->q;
    double T = this->T;

    // Étape 1: Calculer la variance moyenne pour approximer les cumulants
    double avg_variance = g0(T, fix_params(this->N, this->H, this->alpha, this->rn)) / T;
    
    // Étape 2: Cumulants approximatifs pour Lifted
    double c1 = (r - q - 0.5 * avg_variance) * T; // Moyenne du log-price
    double c2 = avg_variance * T; // Variance
    double c4 = 3.0 * c2 * c2; // Kurtosis (approximation Gaussienne)
    

    double a = c1 - L * std::sqrt(std::abs(c2) + std::sqrt(std::abs(c4)));
    double b = c1 + L * std::sqrt(std::abs(c2) + std::sqrt(std::abs(c4)));

    // Calcul des coefficients Vk
    std::vector<double> Vk_vec(N_cos);
    for (int k = 0; k < N_cos; ++k) {
        Vk_vec[k] = Vk(k, a, b);
        cout << "Calcul de Vk pour k=" << k << ", Vk=" << Vk_vec[k] << endl;
    }

    // Calcul des coefficients Re(φ(k))
    std::vector<double> Re_phi_vec(N_cos);
    for (int k = 0; k < N_cos; ++k) {
        double u = (k * M_2_PI) / (b - a); 

        std::complex<double> phi_u = char_func(u, T);
        std::complex<double> exp_term = std::exp(std::complex<double>(0, -u  * a));
        std::complex<double> temp = phi_u * exp_term;
        Re_phi_vec[k] = temp.real();
        if (k == 0) Re_phi_vec[k] *= 0.5; // Premier terme divisé par 2
    }

    // Calcul du prix
    double price = 0.0;
    for (int k = 0; k < N_cos; ++k) {
        price += Re_phi_vec[k] * Vk_vec[k];
    }
    price *= std::exp(-r * T);
    return price;
}

vector<vector<cd>> LiftedHeston::riccati_solver_all( const LiftedHeston::Param_sup& params, cd u, double T, int N_steps) const {
    
    //params 
    double rho = this->rho; double kappa = this->kappa; double sigma = this->sigma;

    double dt = T / N_steps;

    vector<vector<cd>> PSI(N_steps + 1, vector<cd>(this->N + 1, 0.0));

    vector<cd> v_sum_grid(N_steps + 1, 0.0);

    for (int i = 1; i < this->N + 1; ++i) {
        PSI[0][i] = 0.0; // condition initiale
    }

    for (int n = 0; n < N_steps; ++n) {
        cd v_sum = 0.0;
        for (int j = 1; j < this->N + 1; ++j) {
            v_sum += params.c[j-1] * PSI[n][j];
        }
        v_sum_grid[n + 1] = v_sum;

        for (int i = 1; i < this->N + 1; ++i) {
            cd F_val = F_fun(u, v_sum, kappa, sigma, rho);

            PSI[n + 1][i] = PSI[n][i] + dt * F_val;

            PSI[n + 1][i] /= (1.0 + dt * params.x[i - 1]); 
        }
    }
    return PSI;
}

cd LiftedHeston::char_func(double u, double T, int N_steps) const{
    
    //params 
    double v0 = this->v0; double r = this->r; double q = this->q; double S0 = this->S0;

    //fix params
    LiftedHeston::Param_sup params = fix_params(this->N, this->H, this->alpha, this->rn);
    
    //solve riccati
    vector<vector<cd>> PSI = riccati_solver_all(params, u, T, N_steps);
    //compute phi
    double dt = T / N_steps;
    cd intergral = 0.0;

    cd v_sum_T = 0.0;
    for (int k = 0; k<N_steps; ++k){
        for (int j = 1; j < this->N + 1; ++j) {
            v_sum_T += params.c[j-1] * PSI[k][j];
        }
        intergral += F_fun(u, v_sum_T, this->kappa, this->sigma, this->rho) * dt;
        intergral *= g0(T - k * dt, params);

        v_sum_T = 0.0;
    }

    cd phi = std::exp(cd(0, u * std::log(S0)) + intergral);

    return phi;
}

int main() {

    LiftedHeston model(1, 100.0, 100.0, 1.0, 0.0, 0.0, 0.0, 10.0,
                       1.5, 0.04, 0.3, -0.7, 0.04, 0.1, 2, 
                       0.0, 2.5);

    int N_steps = 500;  
    int N_paths = 10000;
    double T_hori = 1.0;
    //double prix_mc = model.prix_mc(N_steps, N_paths, T_hori);
    double prix_cos = model.prix_cos(200, 10.0);
    //std::cout << "Prix MC: " << prix_mc << std::endl;
    std::cout << "Prix COS: " << prix_cos << std::endl;
    return 0;
}
