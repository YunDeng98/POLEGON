//
//  Scaler.cpp
//  arg_branch_length
//
//  Created by Yun Deng on 10/31/23.
//  Modified by Wonseop Lim on 08/15/26.
//

#include <omp.h>
#include "Scaler.hpp"

Scaler::Scaler() {}

void Scaler::compute_deltas(DAG &dag) {
    int n = (int)dag.nodes.size();
    node_deltas.assign(n, 0.0);
    for (int j = 0; j < n; j++) {
        double d = 0;
        for (int k = dag.child_start[j]; k < dag.child_start[j+1]; k++)
            d -= dag.child_span[k];
        for (int k = dag.parent_start[j]; k < dag.parent_start[j+1]; k++)
            d += dag.parent_span[k];
        node_deltas[j] = d;
    }
}

void Scaler::compute_accumulated_arg_length(const vector<vector<double>> &samples) {
    int n = (int)node_deltas.size();
    int K = (int)samples.size();
    size_t m = (size_t)K*n;
    vector<double> pooled(m);
    for (int s = 0; s < K; s++)
        for (int j = 0; j < n; j++)
            pooled[(size_t)s*n + j] = samples[s][j];

    if (sorted_order.size() != m) {
        sorted_order.resize(m);
        for (size_t i = 0; i < m; i++) sorted_order[i] = i;
        sort(sorted_order.begin(), sorted_order.end(),
             [&pooled](size_t a, size_t b) { return pooled[a] < pooled[b]; });
    }

    sorted_times.resize(m);
    rates.resize(m);
    accumulated_arg_length.assign(m, 0.0);
    double r = 0;
    for (size_t i = 0; i < m; i++) {
        size_t idx = sorted_order[i];
        sorted_times[i] = pooled[idx];
        r += node_deltas[idx % n]/K;
        rates[i] = r;
    }
    for (size_t i = 1; i < m; i++)
        accumulated_arg_length[i] = accumulated_arg_length[i-1]
            + rates[i-1]*(sorted_times[i] - sorted_times[i-1]);
}

// Partitions the time axis into num_bins equal-ARG-length windows
void Scaler::compute_old_grid() {
    size_t m = sorted_times.size();
    old_grid.assign(num_bins + 1, 0.0);
    for (int b = 1; b <= num_bins; b++) {
        double partial_arg_length = accumulated_arg_length.back()*b/num_bins;
        auto it = upper_bound(accumulated_arg_length.begin(),
                              accumulated_arg_length.end(), partial_arg_length);
        size_t i = min(m - 1, (size_t)distance(accumulated_arg_length.begin(), it));
        double residue = max(0.0, accumulated_arg_length[i] - partial_arg_length);
        old_grid[b] = sorted_times[i] - residue/rates[i-1];
    }
    old_grid[num_bins] = nextafter(sorted_times[m-1], (double)INT_MAX);
    for (int b = 1; b <= num_bins; b++)
        old_grid[b] = max(old_grid[b], nextafter(old_grid[b-1], (double)INT_MAX));

    observed_arg_length.assign(num_bins, 0.0);
    expected_arg_length.assign(num_bins, accumulated_arg_length.back()/num_bins);
}

void Scaler::map_mutations(DAG &dag, const vector<double> &times) {
    int nb = (int)dag.branches.size();
    const vector<double> &og = old_grid;
    #pragma omp parallel num_threads(num_cores)
    {
        vector<double> local(num_bins, 0.0);
        vector<double> slope(num_bins + 1, 0.0);
        #pragma omp for schedule(static)
        for (int bi = 0; bi < nb; bi++) {
            Branch *b = dag.branches[bi];
            double lb = times[b->lower_node->index];
            double ub = times[b->upper_node->index];
            double w = b->mutation_count;
            double l = ub - lb;
            int x = (int)(upper_bound(og.begin(), og.end(), lb) - og.begin()) - 1;
            int y = (int)(upper_bound(og.begin(), og.end(), ub) - og.begin()) - 1;
            x = min(max(x, 0), num_bins - 1);
            y = min(max(y, 0), num_bins - 1);
            if (x == y) {
                local[x] += w;
                continue;
            }
            double p = w/l;
            local[x] += p*(og[x+1] - lb);
            local[y] += p*(ub - og[y]);
            slope[x+1] += p;
            slope[y] -= p;
        }
        double r = 0;
        for (int k = 0; k < num_bins; k++) {
            r += slope[k];
            local[k] += r*(og[k+1] - og[k]);
        }
        #pragma omp critical
        for (int k = 0; k < num_bins; k++)
            observed_arg_length[k] += local[k];
    }
}

void Scaler::compute_new_grid(double theta) {
    scaling_factors.assign(num_bins, 1.0);
    new_grid.assign(num_bins + 1, 0.0);
    for (int k = 0; k < num_bins; k++) {
        scaling_factors[k] = observed_arg_length[k]/(theta*expected_arg_length[k]);
        new_grid[k+1] = new_grid[k] + (old_grid[k+1] - old_grid[k])*scaling_factors[k];
    }
}

void Scaler::rescale(DAG &dag, vector<vector<double>> &samples, int subsample,
                     double theta) {
    int n_samples = (int)samples.size();
    int K = min(subsample, n_samples);
    vector<vector<double>> sub(K);
    for (int i = 0; i < K; i++)
        sub[i] = samples[(long long)i*n_samples/K];
    compute_accumulated_arg_length(sub);
    compute_old_grid();
    for (int i = 0; i < K; i++)
        map_mutations(dag, sub[i]);
    for (int k = 0; k < num_bins; k++)
        observed_arg_length[k] /= K;
    compute_new_grid(theta);
    #pragma omp parallel for schedule(static) num_threads(num_cores)
    for (int s = 0; s < n_samples; s++)
        apply_scaling_factors(dag, samples[s]);
}

void Scaler::apply_scaling_factors(DAG &dag, vector<double> &times) const {
    int n = (int)dag.nodes.size();
    const vector<double> &og = old_grid;
    for (int i = 0; i < n; i++) {
        if (dag.nodes[i]->is_sample) continue;
        int k = (int)(upper_bound(og.begin(), og.end(), times[i]) - og.begin()) - 1;
        k = min(max(k, 0), num_bins - 1);
        times[i] = new_grid[k] + scaling_factors[k]*(times[i] - og[k]);
    }
}
