//
//  Scaler.cpp
//  arg_branch_length
//
//  Created by Yun Deng on 10/31/23.
//  Modified by Wonseop Lim on 03/16/26.
//

#include "Scaler.hpp"

Scaler::Scaler() {}

// Computes node_deltas: net change in total ARG span at each node's time
void Scaler::compute_deltas(DAG &dag) {
    sorted_nodes.resize(dag.nodes.size());
    copy(dag.nodes.begin(), dag.nodes.end(), sorted_nodes.begin());
    sort(sorted_nodes.begin(), sorted_nodes.end(), compare_node);
    node_deltas.resize(sorted_nodes.size());

    int j = 0;

    for (int i = 0; i < (int)dag.nodes.size(); i++) {
        j = sorted_nodes[i]->index;
        for (int k = dag.child_start[j]; k < dag.child_start[j+1]; k++) {
            node_deltas[i] -= dag.child_data[k]->span;
        }
    }

    for (int i = 0; i < (int)dag.nodes.size(); i++) {
        j = sorted_nodes[i]->index;
        for (int k = dag.parent_start[j]; k < dag.parent_start[j+1]; k++) {
            node_deltas[i] += dag.parent_data[k]->span;
        }
    }

    for (int i = 0; i < (int)sorted_nodes.size() - 1; i++) {
        assert(sorted_nodes[i]->time <= sorted_nodes[i+1]->time);
    }
    /*
    double delta_sum = 0;
    for (int i = 0; i < sorted_nodes.size(); i++) {
        delta_sum += node_deltas[i];
    }
    double delta_sum_b = accumulate(node_deltas.begin(), node_deltas.end(), 0.0);
    assert(delta_sum == delta_sum_b);
     */
}

// Integrates the piecewise-constant ARG rate over time
// rates[i] = cumulative ARG span at sorted_nodes[i]->time (partial_sum of node_deltas)
// accumulated_arg_length[i] = integral of rates from t=0 to sorted_nodes[i]->time
void Scaler::compute_accumulated_arg_length() {
    rates.resize(sorted_nodes.size());
    accumulated_arg_length.resize(sorted_nodes.size());

    partial_sum(node_deltas.begin(), node_deltas.end(), rates.begin());

    for (int i = 1; i < (int)sorted_nodes.size(); i++) {
        accumulated_arg_length[i] = accumulated_arg_length[i-1]
            + rates[i-1] * (sorted_nodes[i]->time - sorted_nodes[i-1]->time);
    }
}

/*
void Scaler::compute_old_grid() {
    double arg_length = 0;
    double rate = 0;
    double base_time = 0;
    for (int i = 0; i < node_deltas.size(); i++) {
        rate += node_deltas[i];
        arg_length += (sorted_nodes[i]->time - base_time)*rate;
        base_time = sorted_nodes[i]->time;
    }
    double unit_arg_length = arg_length/num_windows;
    double partial_length = 0;
    double window_width = 0;
    base_time = 0;
    rate = 0;
    for (int i = 0; i < node_deltas.size(); i++) {
        rate += node_deltas[i];
        partial_length += (sorted_nodes[i]->time - base_time)*rate;
        while (partial_length > unit_arg_length) {
            window_width = (partial_length - unit_arg_length)/rate;
            old_grid.push_back(sorted_nodes[i]->time - window_width);
            base_time = old_grid.back();
            partial_length -= unit_arg_length;
        }
        base_time = sorted_nodes[i]->time;
    }
    if (old_grid.size() < num_windows + 1) {
        old_grid.push_back(sorted_nodes.back()->time);
    } else {
        old_grid.back() = sorted_nodes.back()->time;
    }
    assert(old_grid.size() == num_windows + 1);
    expected_arg_length.resize(num_windows);
    fill(expected_arg_length.begin(), expected_arg_length.end(), unit_arg_length);
}
*/

// Partitions the time axis into num_windows equal-ARG-length windows
void Scaler::compute_old_grid() {
    expected_arg_length.resize(num_windows);
    compute_accumulated_arg_length();

    double unit_arg_length = accumulated_arg_length.back() / num_windows;
    double partial_arg_length = 0;
    int new_index = 0;
    double rate = 0;
    double residue = 0;

    for (int i = 1; i <= num_windows; i++) {
        partial_arg_length = accumulated_arg_length.back() * i / num_windows;
        auto it = upper_bound(accumulated_arg_length.begin(), accumulated_arg_length.end(), partial_arg_length);
        new_index = (int) distance(accumulated_arg_length.begin(), it);
        new_index = min((int) sorted_nodes.size() - 1, new_index);
        rate = rates[new_index - 1];
        residue = accumulated_arg_length[new_index] - partial_arg_length;
        residue = max(0.0, residue);
        expected_arg_length[i-1] = unit_arg_length;
        old_grid.push_back(sorted_nodes[new_index]->time - residue / rate);
    }
    old_grid.back() = nextafter(sorted_nodes.back()->time, INT_MAX);
    assert((int)old_grid.size() == num_windows + 1);
}

// Derives per-window scaling factors and the new time grid
// scaling_factor[k] = observed_arg_length[k] / expected_arg_length[k]
// new_grid is accumulated as: new_width[k] = old_width[k] * scaling_factor[k]
void Scaler::compute_new_grid(double theta) {
    for (auto &x : observed_arg_length) {
        x /= theta; // convert mutation counts to ARG length units
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

// Distributes observed mutations from all branches into windows
// proportionally to each branch's overlap with each window
void Scaler::map_mutations(DAG &dag) {
    observed_arg_length.resize(num_windows);
    for (Branch *b : dag.branches) {
        add_mutation(b->mutation_count, b->lower_node->time, b->upper_node->time);
    }
}

// Adds w mutations from a branch spanning [lb, ub] to windows in proportion
// to their overlap with [lb, ub]. Zero-length branches are treated as point masses
void Scaler::add_mutation(double w, double lb, double ub) {
    double x, y, l, p;
    int index;
    auto it = upper_bound(old_grid.begin(), old_grid.end(), lb);
    it--;
    index = (int) distance(old_grid.begin(), it);
    while (old_grid[index] < ub) {
        x = old_grid[index];
        y = old_grid[index + 1];
        l = min(ub, y) - max(lb, x);
        if (ub - lb == 0) {
            p = 1.0;
        } else {
            p = l / (ub - lb);
            p = min(p, 1.0);
        }
        observed_arg_length[index] += w * p;
        index++;
    }
}

void Scaler::rescale(DAG &dag, double theta) {
    compute_deltas(dag);
    compute_old_grid();
    map_mutations(dag);
    compute_new_grid(theta);

    int k = 0;
    int node_index = 0;
    double t;

    for (int i = 0; i < (int)sorted_nodes.size(); i++) {
        while (sorted_nodes[i]->time > old_grid[k+1]) {
            k++;
        }
        node_index = sorted_nodes[i]->index;
        if (sorted_nodes[i]->is_sample) {
            dag.scaling_factors[node_index] = 1;
        } else {
            t = scaling_factors[k] * (sorted_nodes[i]->time - old_grid[k]) + new_grid[k];
            dag.scaling_factors[node_index] = t / dag.nodes[node_index]->time;
            sorted_nodes[i]->time = t;
        }
    }

    for (int i = 0; i < (int)dag.node_times.size(); i++)
        dag.node_times[i] = dag.nodes[i]->time;

    // Forward topological pass: correct any internal node whose rescaled time
    // now falls at or below its oldest child
    for (int i = 0; i < (int)dag.nodes.size(); i++) {
        if (!dag.nodes[i]->is_sample) {
            double lb = dag.lower_bound(i);
            if (dag.nodes[i]->time <= lb) {
                dag.nodes[i]->time = lb + 1e-6;
                dag.node_times[i] = dag.nodes[i]->time;
            }
        }
    }

    sort(sorted_nodes.begin(), sorted_nodes.end(), compare_node);
    for (int i = 0; i < (int)sorted_nodes.size() - 1; i++) {
        assert(sorted_nodes[i]->time <= sorted_nodes[i+1]->time);
    }
}

// private:

/*
double Scaler::arg_length(double t) {
    if (t >= sorted_nodes.back()->time) {
        return accumulated_arg_length.back();
    } else if (t == 0) {
        return 0.0;
    }
    auto it = lower_bound(sorted_nodes.begin(), sorted_nodes.end(), t, [](const Node* node, double t) {return node->time < t;});
    int index = (int) (it - sorted_nodes.begin());
    double p = (t - sorted_nodes[index]->time)/(sorted_nodes[index+1]->time - sorted_nodes[index]->time);
    double l = p*accumulated_arg_length[index+1] + (1 - p)*accumulated_arg_length[index];
    return l;
}

void Scaler::add_branch_length(Scaler &scaler) {
    expected_arg_length.resize(num_windows);
    for (int i = 0; i < num_windows; i++) {
        expected_arg_length[i] += scaler.arg_length(old_grid[i+1]) - scaler.arg_length(old_grid[i]);
    }
}

void Scaler::all_sample_rescale(DAG &dag, double theta) {
    compute_new_grid(theta);
    int k = 0;
    int node_index = 0;
    double t;
    for (int i = 0; i < sorted_nodes.size(); i++) {
        while (sorted_nodes[i]->time > old_grid[k+1]) {
            k++;
        }
        t = scaling_factors[k]*(sorted_nodes[i]->time - old_grid[k]) + new_grid[k];
        node_index = sorted_nodes[i]->index;
        if (sorted_nodes[i]->is_sample) {
            dag.scaling_factors[node_index] = 1;
        } else {
            dag.scaling_factors[node_index] = t/dag.nodes[node_index]->time;
        }
        sorted_nodes[i]->time = t;
    }
    for (int i = 0; i < sorted_nodes.size() - 1; i++) {
        assert(sorted_nodes[i]->time <= sorted_nodes[i+1]->time);
    }
}

double Scaler::rescale_time(double t) {
    double new_time = 0;
    if (t == 0) {
        return 0;
    } else if (t >= old_grid.back()) {
        new_time = scaling_factors.back()*(t - old_grid[num_windows - 1]) + new_grid[num_windows - 1];
    }
    auto it = lower_bound(old_grid.begin(), old_grid.end(), t);
    int index = (int) (it - old_grid.begin());
    new_time = scaling_factors[index]*(t - old_grid[index]) + new_grid[index];
    return new_time;
}

void Scaler::rescale_samples(DAG &dag) {
    double new_time = 0;
    int num_samples = (int) dag.node_age_samples.front().size();
    for (int i = 0; i < dag.nodes.size(); i++) {
        for (int k = 0; k < num_samples; k++) {
            new_time = rescale_time(dag.node_age_samples[i][k]);
            dag.scaled_node_age_samples[i].push_back(new_time);
        }
    }
}
*/
