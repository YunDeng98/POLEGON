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
#include "Distribution.hpp"

class DAG {
    
public:
    
    int num_leaf_nodes = 0;
    float Ne = 1;
    float lambda = 10;
    int num_posterior_samples = 0;
    Node *root = new Node(numeric_limits<float>::infinity(), INT_MAX);
    vector<Node *> nodes = {};
    vector<float> node_ages = {};
    vector<Branch *> branches = {};
    vector<set<Branch *>> parents = {};
    vector<set<Branch *>> children = {};
    int updates = 0;
    
    DAG(float n);
    
    void load_dag(string node_file, string branch_file, string mut_file);
    
    void compute_mutation_rates(float theta);
    
    void burn_in();
    
    void MCMC(int n, Distribution *d);
    
    void no_prior_MCMC(int n);
    
    void posterior_average();
    
    void write_node_ages(string filename);
    
// private:
    
    float lower_bound(int i);
    
    float upper_bound(int i);
    
    float log_acceptance_weight(int i, float t);
    
    float acceptance_ratio(int i, float t);
    
    float no_prior_acceptance_ratio(int i, float t, float lb, float ub);
    
    void no_prior_propose(int i);
    
    void propose(int i, Distribution *d);
    
    int random_index();
    
    void load_nodes(string node_file);
    
    void load_branches(string branch_file);
    
    void load_mutations(string mut_file);
    
    Branch *search_branch(Node *n1, Node *n2);
    
    float random_non_root_time(float t0, float lb, float ub);
    
    float random_root_time(int i, float lb);
};

#endif /* DAG_hpp */
