//
//  Mutation_map.cpp
//  POLEGON
//
//  Created by Yun Deng on 1/19/24.
//

#include "Mutation_map.hpp"

Mutation_map::Mutation_map() {}

void Mutation_map::load_map(string mut_map_file) {
    ifstream fin(mut_map_file);
    if (!fin.good()) {
        cerr << "input mutation map file not found" << endl;
        exit(1);
    }
    double curr_pos = 0;
    double prev_pos = 0;
    double prev_rate = 0;
    double curr_rate = 0;
    while (fin >> curr_pos >> curr_rate) {
        coordinates.push_back(curr_pos);
        mutation_distance.push_back((curr_pos - prev_pos)*prev_rate);
        prev_pos = curr_pos;
        prev_rate = curr_rate;
    }
}

int Mutation_map::find_index(double x) {
    auto it = upper_bound(coordinates.begin(), coordinates.end(), x);
    it--;
    int index = (int) distance(coordinates.begin(), it);
    assert(index >= 0 and index < coordinates.size() - 1);
    return index;
}

double Mutation_map::mut_distance(double x) {
    int index = find_index(x);
    double p = (x - coordinates[index])/(coordinates[index+1] - coordinates[index]);
    double d = mutation_distance[index]*(1 - p) + mutation_distance[index + 1]*p;
    return d;
}

double Mutation_map::mut_rate_sum(double x, double y) {
    assert(x < y);
    double sum = mut_distance(y) - mut_distance(x);
    assert(sum > 0);
    return sum;
}
