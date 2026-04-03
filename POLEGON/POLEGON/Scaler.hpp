//
//  Scaler.hpp
//  POLEGON
//
//  Created by Yun Deng on 10/31/23.
//  Modified by Wonseop Lim on 04/03/26.
//

#ifndef Scaler_hpp
#define Scaler_hpp

#include <stdio.h>
#include "DAG.hpp"

class Scaler {

public:

    int num_bins = 100;

    vector<Node *> sorted_nodes = {};
    vector<double> node_deltas = {};
    vector<double> rates = {};
    vector<double> accumulated_arg_length = {};
    vector<double> old_grid = {0};
    vector<double> new_grid = {0};
    vector<double> expected_arg_length = {};
    vector<double> observed_arg_length = {};
    vector<double> scaling_factors = {};

    Scaler();

    void compute_deltas(DAG &dag);
    void compute_accumulated_arg_length();
    void compute_old_grid();
    void compute_new_grid(double theta);
    void map_mutations(DAG &dag);
    void add_mutation(double w, double lb, double ub);
    void rescale(DAG &dag, double theta);

};

#endif /* Scaler_hpp */
