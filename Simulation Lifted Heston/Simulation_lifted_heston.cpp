#include <random>
#include <iostream>
#include <fstream>
#include <cmath> 
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
    double delta =30.0/365.0; // 30 jours en années
    int M=1000; // nombre de simulations pour le prix de l'option VIX
    int n=4000; // nombre de points pour l'intégration numérique
    double A;
    double B;
    double b;
   
    Simulation_lifted_heston(double kappa_, double theta_, double sigma_, double rho_, double V0_, double H_, int N_, double S0_, double r_, double T_)
        : kappa(kappa_), theta(theta_), sigma(sigma_), rho(rho_), V0(V0_), H(H_), N(N_), S0(S0_), r(r_), T(T_), alpha(H_ + 0.5) {
     delta = 30.0 / 365.0;
    B = (1.0 - std::exp(-kappa * delta)) / (kappa * delta);
    A = theta * (1.0 - B);
    b = kappa * kappa / (2.0 * sigma * sigma * B) * 0.999;
   
        }
    
   

 // fonctions spécifiques au pricing de l'option VIX
    double d(double u) {
        double discriminant = kappa * kappa -2 * sigma * sigma * u;
        return std::sqrt(std::max(discriminant, 1e-12));

    }
    

    double Bphi(double u) {
        double term=(1-std::exp(-d(u)*T));
        return 2*u*term/((kappa-d(u))*term+2*d(u));
    }
    double Aphi(double u) {
        double term=2*kappa*theta/(sigma*sigma);
        double denom=(kappa+d(u))*(1-std::exp(-d(u)*T))+2*d(u);
        return term*(std::log(2*d(u)*std::exp((kappa+d(u))*T/2)/denom));
    }

    double erfc(double x) {
        return 1.0 - std::erf(x);
    }

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
    if (N==1)
    {
        result.c[0] = 1.0;
        result.x[0] = 0.0;
        result.g0 = V0 + kappa * theta * t;
        return result;
    }
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
void Simulate(std::vector<double>& S, std::vector<double>& V, std::vector<double>& VIX) {
        double dt = T/N_steps;
        S.resize(N_steps+1);
        V.resize(N_steps+1);
        VIX.resize(N_steps+1);
        std::vector<std::vector<double>> U(N, std::vector<double>(N_steps+1, 0.0));
        S[0] = S0;
        V[0] = V0;
        VIX[0]=std::sqrt(theta + (V[0] - theta) * (1 - std::exp(-kappa * delta)) / (kappa * delta));

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
            VIX[j]= std::sqrt(theta + (V[j] - theta) * (1 - std::exp(-kappa * delta)) / (kappa * delta));

        }
    }
    
double call_price() {
    std::ofstream paths("vix_paths.csv");
    paths << "path,step,t,VIX\n";
    double sum_payoff = 0.0; 
    double dt = T / N_steps;
    double K=0.1;
    std::cout << "Valeur du Strike " <<K<<std::endl;
    for (int m = 0; m < M; ++m) {
        std::vector<double> S, V, VIX;
        Simulate(S, V, VIX);
        for (int j = 0; j <= N_steps; ++j) {
            double t = j * dt;
            paths << m << "," << j << "," << t << "," << VIX[j] << "\n";
        }
        sum_payoff += std::max(VIX.back() - K, 0.0);
    }
    paths.close();
    double C0 = std::exp(-r * T) * (sum_payoff / M);
    std::cout << std::fixed << "C0 (prix du call VIX via Monte Carlo) = " << C0 << std::endl;
    return C0;
}
double integrate_simpson(std::function<double(double)> f, double a, double b, int n) {
    if (n % 2 != 0) n++;  // n doit être pair
    double h = (b - a) / n;
    double s = f(a) + f(b);

    for (int i = 1; i < n; ++i) {
        double x = a + i * h;
        s += f(x) * ( (i % 2 == 0) ? 2.0 : 4.0 );
    }

    return s * h / 3.0;
}

std::complex<double> erfc_complex(const std::complex<double>& z) {
    // Approximation numérique stable de l'erreur complémentaire complexe
    // basée sur la fonction de Dawson et sur la série de Abramowitz-Stegun 7.1.23
    const std::complex<double> i(0.0, 1.0);
    std::complex<double> t = 1.0 / (1.0 + 0.5 * std::abs(z));
    std::complex<double> ans = t * std::exp(-z * z - 1.26551223 +
                                           t * (1.00002368 +
                                           t * (0.37409196 +
                                           t * (0.09678418 +
                                           t * (-0.18628806 +
                                           t * (0.27886807 +
                                           t * (-1.13520398 +
                                           t * (1.48851587 +
                                           t * (-0.82215223 +
                                           t * 0.17087277)))))))));
    if (std::real(z) < 0)
        ans = 2.0 - ans;
    return ans;
}
std::complex<double> laplace_CIR(std::complex<double> u) {
    std::complex<double> exp_kappa_t = std::exp(-kappa * T);
    std::complex<double> a = (sigma * sigma / (2.0 * kappa)) * (1.0 - exp_kappa_t);
    std::complex<double> denom = 1.0 + a * u;
    std::complex<double> pow_term = -2.0 * kappa * theta / (sigma * sigma);
    std::complex<double> part1 = std::pow(denom, pow_term);
    std::complex<double> part2 = std::exp(- (u * V0 * exp_kappa_t) / denom);
    return part1 * part2;
}
double call_price_laplace() {
    double K=0.1;
    double phi_R = 0.1;
    double Imax  = 100.0;
    auto integrand = [&](double phi_I) {
    std::complex<double> phi(phi_R, phi_I);  // φ = φR + iφI
    std::complex<double> erfc_term = erfc_complex(K * std::sqrt(phi));
    std::complex<double> expo = std::exp(phi * A )* laplace_CIR(-phi * B);
    return real( (erfc_term * expo / std::pow(phi, 1.5)));
    
};

    double integral = integrate_simpson(integrand, 1e-6, 500, n);
    double C0 = std::exp(-r * T) * (1.0/(2.0*std::sqrt(M_PI)) * integral);
    std::cout << std::fixed << "C0 (prix du call VIX via formule analytique) = " << C0 << std::endl;
    return C0;
}

};
int main() {
    Simulation_lifted_heston sim(0.3, 0.05, 0.1, -0.5, 0.05, 0.1, 1, 100.0, 0.01, 1.0);
    sim.call_price(); 
    sim.call_price_laplace();     
    return 0;
}