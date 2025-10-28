// LiftedHeston_CarrMadan.cpp
// Traduction C++ du script "Call_Price_Lifted_Heston.py"
// Compile: g++ -std=c++17 -O3 -march=native -o lifted_heston LiftedHeston_CarrMadan.cpp

#include <bits/stdc++.h>
using namespace std;

static inline double norm_cdf(double x) {
    // CDF N(0,1)
    return 0.5 * erfc(-x / M_SQRT2);
}

static inline double bs_call(double K, double S0, double T, double r, double sigma) {
    if (sigma <= 0.0) return max(0.0, S0 - K * exp(-r * T));
    double vol = sigma * sqrt(T);
    double d1 = (log(S0 / K) + (r + 0.5 * sigma * sigma) * T) / vol;
    double d2 = d1 - vol;
    return S0 * norm_cdf(d1) - K * exp(-r * T) * norm_cdf(d2);
}

// Brent pour volatilité implicite (prix de call)
static inline double implied_vol_brent(double K, double S0, double r, double T, double target_price,
                                       double lo=1e-6, double hi=5.0, double tol=1e-8, int maxit=200) {
    auto f = [&](double s){ return bs_call(K,S0,T,r,s) - target_price; };
    double a=lo, b=hi, fa=f(a), fb=f(b);
    // élargit la borne haute si nécessaire
    int grow = 0;
    while (fa * fb > 0.0 && grow < 20) { b*=2.0; fb=f(b); grow++; }
    if (fa * fb > 0.0) return numeric_limits<double>::quiet_NaN();

    double c=a, fc=fa, d=b-a, e=d;
    for (int it=0; it<maxit; ++it) {
        if (fabs(fc) < fabs(fb)) { a=b; b=c; c=a; fa=fb; fb=fc; fc=fa; }
        double tol1 = 2.0 * numeric_limits<double>::epsilon() * fabs(b) + 0.5*tol;
        double xm = 0.5*(c-b);
        if (fabs(xm) <= tol1 || fb == 0.0) return b;

        if (fabs(e) >= tol1 && fabs(fa) > fabs(fb)) {
            // interpolation inverse quadratique / sécante
            double s, p, q;
            if (a == c) {
                // sécante
                s = fb/fa;
                p = 2.0 * xm * s;
                q = 1.0 - s;
            } else {
                // IQI
                double s1 = fb / fa;
                double s2 = fb / fc;
                double s3 = fa / fc;
                p = s2*(2.0*xm*s1*(s1 - s3) - (b - a)*(s1 - 1.0));
                q = (s1 - 1.0)*(s2 - 1.0)*(s3 - 1.0);
            }
            if (p > 0) q = -q;
            p = fabs(p);
            double min1 = 3.0*xm*q - fabs(tol1*q);
            double min2 = fabs(e*q);
            if (2.0*p < (min1 < min2 ? min1 : min2)) { e=d; d=p/q; }
            else { d=xm; e=d; }
        } else { d=xm; e=d; }

        a = b; fa = fb;
        if (fabs(d) > tol1) b += d; else b += (xm > 0 ? tol1 : -tol1);
        fb = f(b);
        if ((fb>0 && fc>0) || (fb<0 && fc<0)) { c=a; fc=fa; e=d=b-a; }
    }
    return b; // tolérance atteinte
}

// Simpson 1D sur [a,b] avec n intervalles pairs
template<class F>
double simpson(F&& f, double a, double b, int n){
    if (n % 2) n++;
    double h = (b - a) / n;
    double s = f(a) + f(b);
    for (int i = 1; i < n; ++i) {
        double x = a + i * h;
        s += f(x) * (i % 2 ? 4.0 : 2.0);
    }
    return s * (h / 3.0);
}

// Fonction caractéristique Lifted Heston (approche "lifted" par N facteurs)
std::complex<double> Ch_Lifted_Heston(
    std::complex<double> omega, // i*ω déjà appliqué à l'entrée
    double S0, double T,
    double rho, double lamb, double theta, double nu, double V0,
    int N, double rN, double alpha, int M)
{
    using cd = complex<double>;
    const cd i(0.0, 1.0);

    // Grille pour c_j et gamma_j
    vector<double> h(N);
    for (int j=0;j<N;++j) h[j]=j;
    vector<double> rpowerN(N);
    for (int j=0;j<N;++j) rpowerN[j] = pow(rN, h[j] - 0.5*N);

    // poids c_j et mean-reversions gamma_j
    vector<double> c(N), gammas(N);
    double ga = tgamma(alpha);
    double g2a = tgamma(2.0 - alpha);
    double coeff_c = (pow(rN, 1.0 - alpha) - 1.0) / (ga * g2a);
    for (int j=0;j<N;++j) c[j] = coeff_c * pow(rpowerN[j], 1.0 - alpha);

    double ratio = ((1.0 - alpha) / (2.0 - alpha))
                 * ((pow(rN, 2.0 - alpha) - 1.0) / (pow(rN, 1.0 - alpha) - 1.0));
    for (int j=0;j<N;++j) gammas[j] = ratio * rpowerN[j];

    // Courbe initiale g(t)
    auto g = [&](double t){
        // V0 + λθ * sum_j (c_j / gamma_j) * (1 - e^{-gamma_j t})
        double sum = 0.0;
        for (int j=0;j<N;++j) sum += (c[j] / gammas[j]) * (1.0 - exp(-gammas[j] * t));
        return V0 + lamb * theta * sum;
    };

    // Discrétisation en temps
    double delta = T / M;
    vector<double> t(M+1);
    for (int k=0;k<=M;++k) t[k] = k * delta;

    auto F = [&](cd u, cd v) {
        return 0.5*(u*u - u) + (rho*nu*u - lamb)*v + 0.5*nu*nu*v*v;
    };

    // itération pour psi (M+1 x N)
    vector<vector<cd>> psi(M+1, vector<cd>(N, cd(0.0,0.0)));
    vector<double> denom(N);
    for (int j=0;j<N;++j) denom[j] = 1.0 + delta * gammas[j];

    for (int k=1; k<=M; ++k) {
        // v = dot(c, psi[k-1,:])
        cd v = 0.0;
        for (int j=0;j<N;++j) v += c[j] * psi[k-1][j];
        cd incr = delta * F(omega, v);
        for (int j=0;j<N;++j) {
            psi[k][j] = (psi[k-1][j] + incr) / denom[j];
        }
    }

    // g_0[k] = g(T - t[k])
    vector<double> g0(M+1);
    for (int k=0;k<=M;++k) g0[k] = g(T - t[k]);

    // Y[k] = F(omega, dot(c, psi[k,:])) * g0[k]
    vector<cd> Y(M+1);
    for (int k=0;k<=M;++k) {
        cd v=0.0;
        for (int j=0;j<N;++j) v += c[j] * psi[k][j];
        Y[k] = F(omega, v) * g0[k];
    }

    // Trapèzes pour phi = ∫ Y dt
    cd phi_int = 0.0;
    for (int k=0;k<=M;++k) {
        double w = (k==0 || k==M) ? 0.5*delta : delta;
        phi_int += w * Y[k];
    }

    cd phi = exp(omega*log(S0) + phi_int);
    return phi;
}

// Prix du call via Carr–Madan
double Call_Price_Lifted_Heston(
    double S0, double K, double T, double r,
    double rho, double lamb, double theta, double nu, double V0,
    int N, double rN, double alpha, int M,
    double alpha2, double L, int n_int = 20000)
{
    using cd = complex<double>;
    const cd i(0.0, 1.0);

    auto phi = [&](double w_real){
        // En Python: i*omega est passé à Ch_Lifted_Heston
        cd w = i * w_real;
        return Ch_Lifted_Heston(w, S0, T, rho, lamb, theta, nu, V0, N, rN, alpha, M);
    };

    auto integrand = [&](double w){
        // integrand = Re{ [phi(w - i*(a+1)) / (a^2 + a - w^2 + i(2a+1)w)] * e^{-i w logK} }
        double a = alpha2;
        cd wshift = w - (-1.0)*(a+1.0); // minus i*(a+1)  -> add i*(a+1) in imaginary part
        // Plus soigneusement:
        cd s = cd(w, -(a+1.0)); // ω - i(α+1)
        cd denom = (a*a + a - w*w) + cd(0.0, (2.0*a + 1.0) * w);
        cd val = Ch_Lifted_Heston(cd(0.0,1.0)*s, S0, T, rho, lamb, theta, nu, V0, N, rN, alpha, M);
        cd e = exp(cd(0.0, -w * log(K)));
        cd frac = val / denom * e;
        return frac.real();
    };

    // Intégration Simpson sur [0, L]
    double I = simpson(integrand, 0.0, L, n_int);

    double P = exp(-r*T - alpha2 * log(K)) * I / M_PI;
    return P;
}

int main(){
    // Paramètres (comme dans le script Python)
    double S0=1.0, r=0.0, rho=-0.7, lamb=0.3, theta=0.02, nu=0.3, V0=0.02;
    int N=20;
    double rN = 1.0 + 10.0 * pow(N, -0.9);
    double alpha = 0.6;
    int M = 100;
    double L = 1000.0;
    double alpha2 = 1.0;
    double T = 1.0/26.0;

    // Grille log-strike
    int Nk = 20;
    vector<double> logK(Nk);
    double a = -0.15, b = 0.05;
    for (int i=0;i<Nk;++i) logK[i] = a + (b-a) * i / (Nk-1);

    cout.setf(std::ios::fixed); cout<<setprecision(6);
    cout << "idx,logK,K,price,impvol\n";

    for (int i=0;i<Nk;++i){
        double K = exp(logK[i]);
        double price = Call_Price_Lifted_Heston(
            S0, K, T, r, rho, lamb, theta, nu, V0,
            N, rN, alpha, M, alpha2, L, 20000 // n_int
        );

        double iv = implied_vol_brent(K, S0, r, T, price);
        cout << i << "," << logK[i] << "," << K << "," << price << "," << iv << "\n";
    }

    return 0;
}
