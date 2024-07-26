//
//  Scaler.cpp
//  arg_branch_length
//
//  Created by Yun Deng on 10/31/23.
//

#include "Scaler.hpp"

Scaler::Scaler() {}

void Scaler::compute_deltas(DAG &dag) {
    sorted_nodes.resize(dag.nodes.size());
    copy(dag.nodes.begin(), dag.nodes.end(), sorted_nodes.begin());
    sort(sorted_nodes.begin(), sorted_nodes.end(), compare_node);
    node_deltas.resize(sorted_nodes.size());
    int j = 0;
    for (int i = 0; i < dag.children.size(); i++) {
        j = sorted_nodes[i]->index; // the position in the dag vectors
        for (Branch *b : dag.children[j]) {
            node_deltas[i] -= b->span;
        }
    }
    for (int i = 0; i < dag.parents.size(); i++) {
        j = sorted_nodes[i]->index;
        for (Branch *b : dag.parents[j]) {
            node_deltas[i] += b->span;
        }
    }
    for (int i = 0; i < sorted_nodes.size() - 1; i++) {
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

void Scaler::compute_old_grid() {
    expected_arg_length.resize(num_windows);
    rates.resize(sorted_nodes.size());
    accumulated_arg_length.resize(sorted_nodes.size());
    partial_sum(node_deltas.begin(), node_deltas.end(), rates.begin());
    for (int i = 1; i < sorted_nodes.size(); i++) {
        accumulated_arg_length[i] = accumulated_arg_length[i-1] + rates[i-1]*(sorted_nodes[i]->time - sorted_nodes[i-1]->time);
    }
    double unit_arg_length = accumulated_arg_length.back()/num_windows;
    double partial_arg_length = 0;
    int new_index = 0;
    double rate = 0;
    double residue = 0;
    for (int i = 1; i <= num_windows; i++) {
        partial_arg_length = accumulated_arg_length.back()*i/num_windows;
        auto it = upper_bound(accumulated_arg_length.begin(), accumulated_arg_length.end(), partial_arg_length);
        new_index = (int) distance(accumulated_arg_length.begin(), it);
        new_index = min((int) sorted_nodes.size() - 1, new_index);
        rate = rates[new_index - 1];
        residue = accumulated_arg_length[new_index] - partial_arg_length;
        residue = max(0.0, residue);
        expected_arg_length[i-1] = unit_arg_length;
        old_grid.push_back(sorted_nodes[new_index]->time - residue/rate);
    }
    old_grid.back() = nextafter(sorted_nodes.back()->time, INT_MAX);
    assert(old_grid.size() == num_windows + 1);
}


void Scaler::compute_new_grid(double theta) {
    for (auto &x : observed_arg_length) {
        x /= theta; // convert mutation counts back to arg length
    }
    double base_time = 0;
    double old_window_width = 0, scaling_factor = 0;
    new_grid.reserve(old_grid.size());
    for (int i = 1; i < old_grid.size(); i++) {
        old_window_width = old_grid[i] - old_grid[i-1];
        scaling_factor = observed_arg_length[i-1]/expected_arg_length[i-1];
        assert(!isnan(scaling_factor));
        scaling_factors.push_back(scaling_factor);
        base_time += old_window_width*scaling_factor;
        new_grid.push_back(base_time);
    }
}

void Scaler::map_mutations(DAG &dag) {
    observed_arg_length.resize(num_windows);
    for (Branch *b : dag.branches) {
        add_mutation(b->mutation_count, b->lower_node->time, b->upper_node->time);
    }
}

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
            p = l/(ub - lb);
            p = min(p, 1.0);
        }
        observed_arg_length[index] += w*p;
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
