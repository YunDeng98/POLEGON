//
//  Scaler.cpp
//  arg_branch_length
//
//  Created by Yun Deng on 10/31/23.
//  Modified by Wonseop Lim on 08/13/26.
//

#include <omp.h>
#include "Scaler.hpp"

Scaler::Scaler() {}

void Scaler::compute_deltas(DAG &dag) {
    if (sorted_nodes.empty()) {
        sorted_nodes.resize(dag.nodes.size());
        copy(dag.nodes.begin(), dag.nodes.end(), sorted_nodes.begin());
        sort(sorted_nodes.begin(), sorted_nodes.end(), [this](const Node* a, const Node* b) {
            if (local_times[a->index] != local_times[b->index])
                return local_times[a->index] < local_times[b->index];
            return a->index < b->index;
        });
    }
    int n = (int)sorted_nodes.size();
    node_deltas.assign(n, 0.0);
    #pragma omp parallel for schedule(static) num_threads(num_cores)
    for (int i = 0; i < n; i++) {
        int j = sorted_nodes[i]->index;
        double d = 0;
        for (int k = dag.child_start[j]; k < dag.child_start[j+1]; k++)
            d -= dag.child_span[k];
        for (int k = dag.parent_start[j]; k < dag.parent_start[j+1]; k++)
            d += dag.parent_span[k];
        node_deltas[i] = d;
    }
}

void Scaler::compute_accumulated_arg_length() {
    rates.resize(sorted_nodes.size());
    accumulated_arg_length.resize(sorted_nodes.size());

    partial_sum(node_deltas.begin(), node_deltas.end(), rates.begin());

    for (int i = 1; i < (int)sorted_nodes.size(); i++) {
        accumulated_arg_length[i] = accumulated_arg_length[i-1]
            + rates[i-1] * (local_times[sorted_nodes[i]->index] - local_times[sorted_nodes[i-1]->index]);
    }
}

// Partitions the time axis into num_bins equal-ARG-length windows
void Scaler::compute_old_grid() {
    expected_arg_length.resize(num_bins);
    compute_accumulated_arg_length();

    double unit_arg_length = accumulated_arg_length.back() / num_bins;
    double partial_arg_length = 0;
    int new_index = 0;
    double rate = 0;
    double residue = 0;

    for (int i = 1; i <= num_bins; i++) {
        partial_arg_length = accumulated_arg_length.back() * i / num_bins;
        auto it = upper_bound(accumulated_arg_length.begin(), accumulated_arg_length.end(), partial_arg_length);
        new_index = (int) distance(accumulated_arg_length.begin(), it);
        new_index = min((int) sorted_nodes.size() - 1, new_index);
        rate = rates[new_index - 1];
        residue = accumulated_arg_length[new_index] - partial_arg_length;
        residue = max(0.0, residue);
        expected_arg_length[i-1] = unit_arg_length;
        old_grid.push_back(local_times[sorted_nodes[new_index]->index] - residue / rate);
    }
    old_grid.back() = nextafter(local_times[sorted_nodes.back()->index], INT_MAX);
    assert((int)old_grid.size() == num_bins + 1);
}

void Scaler::compute_new_grid(double theta) {
    for (auto &x : observed_arg_length) {
        x /= theta;
    }
    double base_time = 0;
    double old_window_width = 0, scaling_factor = 0;
    new_grid.reserve(old_grid.size());
    for (int i = 1; i < (int)old_grid.size(); i++) {
        old_window_width = old_grid[i] - old_grid[i-1];
        scaling_factor = observed_arg_length[i-1] / expected_arg_length[i-1];
        assert(!isnan(scaling_factor));
        scaling_factors.push_back(scaling_factor);
        base_time += old_window_width * scaling_factor;
        new_grid.push_back(base_time);
    }
}

void Scaler::map_mutations(DAG &dag) {
    observed_arg_length.assign(num_bins, 0.0);
    int nb = (int)dag.branches.size();
    const vector<double> &og = old_grid;
    #pragma omp parallel num_threads(num_cores)
    {
        vector<double> local(num_bins, 0.0);
        #pragma omp for schedule(static)
        for (int bi = 0; bi < nb; bi++) {
            Branch *b = dag.branches[bi];
            double lb = local_times[b->lower_node->index];
            double ub = local_times[b->upper_node->index];
            double w  = b->mutation_count;
            auto it = upper_bound(og.begin(), og.end(), lb);
            --it;
            int idx = (int)(it - og.begin());
            while (og[idx] < ub) {
                double x = og[idx], y = og[idx+1];
                double l = min(ub, y) - max(lb, x);
                double p = (ub == lb) ? 1.0 : min(l / (ub - lb), 1.0);
                local[idx] += w * p;
                ++idx;
            }
        }
        #pragma omp critical
        for (int k = 0; k < num_bins; k++)
            observed_arg_length[k] += local[k];
    }
}

void Scaler::rescale(DAG &dag, double theta) {
    old_grid = {0};
    new_grid = {0};
    scaling_factors.clear();

    compute_deltas(dag);
    compute_old_grid();
    map_mutations(dag);
    compute_new_grid(theta);

    int k = 0;
    int node_index = 0;
    double t;

    for (int i = 0; i < (int)sorted_nodes.size(); i++) {
        while (local_times[sorted_nodes[i]->index] > old_grid[k+1]) {
            k++;
        }
        node_index = sorted_nodes[i]->index;
        if (!sorted_nodes[i]->is_sample) {
            t = scaling_factors[k] * (local_times[node_index] - old_grid[k]) + new_grid[k];
            local_times[node_index] = t;
        }
    }

    for (int i = 0; i < (int)dag.nodes.size(); i++) {
        if (!dag.nodes[i]->is_sample) {
            double lb = dag.lower_bound(i, local_times);
            if (local_times[i] <= lb) {
                local_times[i] = lb + 1e-6;
            }
        }
    }

    sort(sorted_nodes.begin(), sorted_nodes.end(), [this](const Node* a, const Node* b) {
        if (local_times[a->index] != local_times[b->index])
            return local_times[a->index] < local_times[b->index];
        return a->index < b->index;
    });
}
