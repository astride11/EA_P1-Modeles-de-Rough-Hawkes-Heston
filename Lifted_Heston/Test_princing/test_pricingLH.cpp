#include <iostream> 
#include <iomanip>
#include <vector>
#include <random>
#include <algorithm>
#include "../LiftedH.hpp"
#include <cmath>


using namespace std;

int main() {
    LiftedHeston model(1, 100.0, 100.0, 1.0, 0.0, 0.01, 0.0, 0.0, 1.5, 0.04, 0.3, -0.7, 0.04, 0.1, 50, 0.0, 2.5);
    
    int N_steps = 1000;
    vector<int> N_paths_list = {1000, 5000, 10000, 50000, 100000};
    double T_hori = 1.0;

    for (int N_paths : N_paths_list) {
        double price = model.prix_mc(N_steps, N_paths, T_hori);
        std::cout << "Prix MC Lifted Heston (N_paths=" << N_paths << "): " << std::setprecision(6) << price << std::endl;
    }

    return 0;
}