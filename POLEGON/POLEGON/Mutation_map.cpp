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
    mutation_distances.push_back(0);
    double left;
    double right = 0;
    double rate;
    double mut_dist;
    while (fin >> left >> right >> rate) {
        coordinates.push_back(left);
        mut_dist = mutation_distances.back() + rate*(right - left);
        mutation_distances.push_back(mut_dist);
    }
    if (coordinates.empty()) {
        cerr << "mutation map file has no intervals: " << mut_map_file << endl;
        exit(1);
    }
    sequence_length = right;
    coordinates.push_back(sequence_length);
}

int Mutation_map::find_index(double x) {
    auto it = upper_bound(coordinates.begin(), coordinates.end(), x);
    int index = (int) distance(coordinates.begin(), it) - 1;
    return min(max(index, 0), (int) coordinates.size() - 2);
}

double Mutation_map::mutation_distance(double x) {
    int index = find_index(x);
    double prev_dist = mutation_distances[index];
    double next_dist = mutation_distances[index+1];
    double p = (x - coordinates[index])/(coordinates[index+1] - coordinates[index]);
    double dist = (1-p)*prev_dist + p*next_dist;
    return dist;
}

double Mutation_map::mutation_rate(double x, double y) {
    return mutation_distance(y) - mutation_distance(x);
}

double Mutation_map::mean_rate() {
    double mr = mutation_distances.back()/sequence_length;
    return mr;
}
