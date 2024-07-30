//
//  Scaler.hpp
//  arg_branch_length
//
//  Created by Yun Deng on 10/31/23.
//

#ifndef Scaler_hpp
#define Scaler_hpp

#include <stdio.h>
#include "DAG.hpp"

class Scaler {

public:
    
    int num_windows = 100;
    
    vector<Node *> sorted_nodes = {}; // nodes sorted in the order of time
    vector<double> node_deltas = {}; // nodes rate changes in the same sorted order
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
    
    // void all_sample_rescale(DAG &dag, double theta);
    
};

#endif /* Scaler_hpp */
