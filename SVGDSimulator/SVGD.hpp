#pragma once
#include <vector>
#include <functional>
#include <cmath>
#include <algorithm>

template <typename VecType>
class SVGD {
public:
    using ScoreFunction = std::function<VecType(const VecType&)>;

    static void step(std::vector<VecType>& particles, ScoreFunction score_func, double step_size) {
        size_t n = particles.size();
        if (n == 0) return;
        double h = compute_median_heuristic(particles);
        std::vector<VecType> phi(n);
        std::vector<VecType> scores(n);
        for(size_t i=0; i<n; ++i) scores[i] = score_func(particles[i]);

        for (size_t i = 0; i < n; ++i) {
            VecType sum_term;
            for (size_t j = 0; j < n; ++j) {
                VecType diff = particles[i] - particles[j];
                double sq_dist = diff.normSq();
                double k_val = std::exp(-sq_dist / h);
                VecType grad_k = diff * (k_val * 2.0 / h);
                sum_term = sum_term + (scores[j] * k_val) + grad_k;
            }
            phi[i] = sum_term * (1.0 / n);
        }
        for (size_t i = 0; i < n; ++i) particles[i] = particles[i] + phi[i] * step_size;
    }

private:
    static double compute_median_heuristic(const std::vector<VecType>& particles) {
        size_t n = particles.size();
        if(n < 2) return 1.0;
        std::vector<double> sq_dists;
        sq_dists.reserve(n * n);
        for (size_t i = 0; i < n; ++i) {
            for (size_t j = i + 1; j < n; ++j) {
                sq_dists.push_back((particles[i] - particles[j]).normSq());
            }
        }
        if (sq_dists.empty()) return 1.0;
        size_t mid_idx = sq_dists.size() / 2;
        std::nth_element(sq_dists.begin(), sq_dists.begin() + mid_idx, sq_dists.end());
        double median = sq_dists[mid_idx];
        return (median < 1e-6) ? 1.0 : median / std::log(n + 1.0);
    }
};