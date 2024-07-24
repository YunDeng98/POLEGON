//
//  DAG.hpp
//  arg_branch_length
//
//  Created by Yun Deng on 10/31/23.
//

#ifndef DAG_hpp
#define DAG_hpp

#include <stdio.h>
#include <vector>
#include <map>
#include <unordered_map>
#include <set>
#include <unordered_set>
#include "Node.hpp"
#include "Branch.hpp"
#include "random_utils.hpp"
#include "Mutation_map.hpp"
#include "Distribution.hpp"

class DAG {
    
public:
    
    int num_leaf_nodes = 0;
    double Ne = 1;
    double lambda = 10;
    int num_posterior_samples = 0;
    Node *root = new Node(numeric_limits<double>::infinity(), INT_MAX);
    vector<Node *> nodes = {};
    vector<double> node_ages = {};
    vector<vector<double>> node_age_samples = {};
    vector<vector<double>> scaled_node_age_samples = {};
    vector<Branch *> branches = {};
    vector<set<Branch *>> parents = {};
    vector<set<Branch *>> children = {};
    int updates = 0;
    
    DAG(double n);
    
    void load_dag(string node_file, string branch_file);
    
    void load_dag(string node_file, string branch_file, Mutation_map &mm);
    
    void map_mutations(string mut_file);
    
    void compute_mutation_rates(double theta);
    
    void burn_in();
    
    void record_node_ages();
    
    void record_scaled_node_ages();
    
    void no_prior_MCMC();
    
    void sample(int i);
    
    void posterior_average();
    
    void scaled_sample_average();
    
    void write_node_ages(string filename);
    
// private:
    
    double lower_bound(int i);
    
    double upper_bound(int i);
    
    double log_acceptance_weight(int i, double t);
    
    double acceptance_ratio(int i, double t);
    
    double fast_acceptance_ratio(int i, double t0, double t);
    
    double no_prior_acceptance_ratio(int i, double t, double lb, double ub);
    
    void no_prior_propose(int i);
    
    void propose(int i, Distribution *d);
    
    int random_index();
    
    vector<int> get_permutation();
    
    void load_nodes(string node_file);
    
    void load_branches(string branch_file);
    
    void load_branches(string branch_file, Mutation_map &mm);
    
    // void load_mutations(string mut_file);
    
    Branch *search_branch(Node *n1, Node *n2);
    
    double random_non_root_time(double t0, double lb, double ub);
    
    double random_root_time(int i, double lb);
};

#endif /* DAG_hpp */
