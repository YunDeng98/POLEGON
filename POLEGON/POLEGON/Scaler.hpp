//
//  Scaler.hpp
//  POLEGON
//
//  Created by Yun Deng on 10/31/23.
//  Modified by Wonseop Lim on 08/15/26.
//

#ifndef Scaler_hpp
#define Scaler_hpp

#include <stdio.h>
#include "DAG.hpp"

class Scaler {

public:

    int num_bins = 100;
    int num_cores = 1;

    vector<double> node_deltas = {};
    vector<size_t> sorted_order = {};
    vector<double> sorted_times = {};
    vector<double> rates = {};
    vector<double> accumulated_arg_length = {};
    vector<double> old_grid = {};
    vector<double> new_grid = {};
    vector<double> expected_arg_length = {};
    vector<double> observed_arg_length = {};
    vector<double> scaling_factors = {};

    Scaler();

    void compute_deltas(DAG &dag);
    void compute_accumulated_arg_length(const vector<vector<double>> &samples);
    void compute_old_grid();
    void compute_new_grid(double theta);
    void map_mutations(DAG &dag, const vector<double> &times);
    void apply_scaling_factors(DAG &dag, vector<double> &times) const;
    void rescale(DAG &dag, vector<vector<double>> &samples, int subsample, double theta);

};

#endif /* Scaler_hpp */
