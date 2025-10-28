#include <random>
#include <iostream>
#include <fstream>
#include "../LiftedH.hpp"

using namespace std;


int main() {
    std::vector<int> Nvals = {10, 20, 50, 100, 200, 500};

    for (int idx = 0; idx < Nvals.size(); ++idx) {
        int N = Nvals[idx];
        LiftedHeston lifted_heston_model = LiftedHeston(1, 100.0, 100.0, 1.0, 0.0, 0.01, 0.0, 0.0, 0.3, 0.05, 0.1, -0.5, 0.05, 0.1, N, 0.0, 2.5);
        int N_steps = 1000;
        double T_hori = 1.0;   
        std::vector<double> S(N_steps + 1);
        std::vector<double> V(N_steps + 1);
        lifted_heston_model.simulate_paths(N_steps, T_hori, S, V);

        // Write to CSV
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