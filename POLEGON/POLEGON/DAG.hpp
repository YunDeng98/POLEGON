//
//  DAG.hpp
//  POLEGON
//
//  Created by Yun Deng on 10/31/23.
//  Updated by Wonseop Lim on 03/16/26.
//

#ifndef DAG_hpp
#define DAG_hpp

#include <stdio.h>
#include <vector>
#include <map>
#include <unordered_map>
#include "Node.hpp"
#include "Branch.hpp"
#include "random_utils.hpp"
#include "Mutation_map.hpp"
#include "Distribution.hpp"

class DAG {

public:

    int num_leaf_nodes = 0;
    double Ne = 1;
    double lambda = 5;      // scale of the exponential proposal for root node times
    double time_origin = 0; // age of youngest tip in coalescent units

    // Virtual root at t=infinity; all local-tree roots connect to this node
    // so that branches leading to the root are handled uniformly without special cases
    Node *root = new Node(numeric_limits<double>::infinity(), INT_MAX);

    vector<Node *> nodes = {};
    vector<double> scaling_factors = {};
    vector<Branch *> branches = {};

    // CSR adjacency: for node i,
    //   parent_data[parent_start[i]..parent_start[i+1]) = branches where i is the lower (younger) endpoint
    //   child_data [child_start[i] ..child_start[i+1])  = branches where i is the upper (older) endpoint
    vector<Branch *> parent_data = {};
    vector<int> parent_start = {};
    vector<Branch *> child_data = {};
    vector<int> child_start = {};

    vector<int> perm_cache = {};

    int updates = 0;
    double max_step = 10.0; // maximum exponential draw for root node proposals

    DAG(double n);

    void load_dag(string node_file, string branch_file);
    void load_dag(string node_file, string branch_file, Mutation_map &mm);
    void map_mutations(string mut_file);
    void compute_mutation_rates(double theta);
    void apply_tip_ages(string tip_ages_file, double gen_time);
    void no_prior_MCMC();
    void posterior_average(string samples_file, string output_file);
    void write_node_ages(string filename, double gen_time = 1);

// private:

    double lower_bound(int i);  // max time among children of i
    double upper_bound(int i);  // min time among parents of i
    double log_acceptance_weight(int i, double t);
    double acceptance_ratio(int i, double t);
    double fast_acceptance_ratio(int i, double t0, double t1);
    double no_prior_acceptance_ratio(int i, double t, double lb, double ub);
    void no_prior_propose(int i);
    void propose(int i, Distribution *d);
    int random_index();
    void load_nodes(string node_file);
    void load_branches(string branch_file);
    void load_branches(string branch_file, Mutation_map &mm);
    Branch *search_branch(Node *n1, Node *n2);
    double random_non_root_time(double t0, double lb, double ub);
    double random_root_time(int i, double lb);
    double median(vector<double> &values);
};

#endif /* DAG_hpp */
