#include <iostream>
#include <vector>
#include <random>
#include <cmath>
#include <algorithm>

float GARCHEuropeanCallMCPricer(float S0, float r, const std::vector<float>& sigma, float T, int n, int m, float K) {
    float dt = T / n;
    std::vector<std::vector<float>> gbm_paths(m, std::vector<float>(n + 1, 0.0));

    for (int i = 0; i < m; i++) {
        gbm_paths[i][0] = S0;
    }
    
    std::default_random_engine generator;
    std::normal_distribution<float> standard_normal(0.0, 1.0);

    for (int i = 1; i <= n; i++) {
        std::vector<float> Z(m);
        for (int j = 0; j < m; j++) {
            Z[j] = standard_normal(generator);
        }

        for (int j = 0; j < m; j++) {
            gbm_paths[j][i] = gbm_paths[j][i - 1] * (1 + r * dt + sigma[i - 1] * std::sqrt(dt) * Z[j]);
        }
    }

    float sum_payoffs = 0.0;
    for (int j = 0; j < m; j++) {
        float payoff = std::max(gbm_paths[j][n] - K, 0.0f);
        sum_payoffs += payoff;
    }

    float discount_factor = std::exp(-r * T);
    float call_price = discount_factor * (sum_payoffs / m);

    return call_price;
}

int main() {
    float S0 = 100.0f;
    float K = 110.0f;
    float r = 0.05f;

    std::vector<float> sigma_arr(30, 0.2f);  

    float T = 30.0f / 252.0f;
    int n = 30;
    int m = 100000;

    float call_price = GARCHEuropeanCallMCPricer(S0, r, sigma_arr, T, n, m, K);

    std::cout << "Estimated Call Option Price: " << call_price << std::endl;

    return 0;
}
