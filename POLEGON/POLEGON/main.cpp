//
//  main.cpp
//  POLEGON
//
//  Created by Yun Deng on 12/13/23.
//

#include <iostream>
#include "DAG.hpp"
#include "Distribution.hpp"
#include "Scaler.hpp"
#include "Mutation_map.hpp"

int main(int argc, const char * argv[]) {
    double m = -1;
    int num_samples = -1;
    int burn_in = -1;
    int spacing = -1;
    int scaling_rep = 0;
    double max_step = 10.0;
    string input_prefix = "", output_prefix = "";
    int seed = 42;
    double Ne = 0;
    string map_file = "";
    bool write_samples = false;
    for (int i = 1; i < argc; ++i) {
        string arg = argv[i];
        if (arg == "-m") {
            if (i + 1 >= argc || argv[i+1][0] == '-') {
                cerr << "Error: -m flag cannot be empty. " << endl;
                exit(1);
            }
            try {
                m = stod(argv[++i]);
            } catch (const invalid_argument&) {
                cerr << "Error: -m flag expects a number. " << endl;
                exit(1);
            }
        }
        else if (arg == "-mutation_map") {
            if (i + 1 >= argc || argv[i+1][0] == '-') {
                cerr << "Error: -mutation_map flag cannot be empty. " << endl;
                exit(1);
            }
            try {
                map_file = argv[++i];
            } catch (const invalid_argument&) {
                cerr << "Error: -mutation_map flag expects a string. " << endl;
                exit(1);
            }
        }
        else if (arg == "-Ne") {
            if (i + 1 >= argc || argv[i+1][0] == '-') {
                cerr << "Error: -Ne flag cannot be empty. " << endl;
                exit(1);
            }
            try {
                Ne = stod(argv[++i]);
            } catch (const invalid_argument&) {
                cerr << "Error: -Ne flag expects a number. " << endl;
                exit(1);
            }
        }
        else if (arg == "-burn_in") {
            if (i + 1 >= argc || argv[i+1][0] == '-') {
                cerr << "Error: -burn_in flag cannot be empty. " << endl;
                exit(1);
            }
            try {
                burn_in = stoi(argv[++i]);
            } catch (const invalid_argument&) {
                cerr << "Error: -burn_in flag expects a number. " << endl;
                exit(1);
            }
        }
        else if (arg == "-num_samples") {
            if (i + 1 >= argc || argv[i+1][0] == '-') {
                cerr << "Error: -num_samples flag cannot be empty. " << endl;
                exit(1);
            }
            try {
                num_samples = stoi(argv[++i]);
            } catch (const invalid_argument&) {
                cerr << "Error: -num_samples flag expects a number. " << endl;
                exit(1);
            }
        }
        else if (arg == "-scaling_rep") {
            if (i + 1 >= argc || argv[i+1][0] == '-') {
                cerr << "Error: -scaling_rep flag cannot be empty. " << endl;
                exit(1);
            }
            try {
                scaling_rep = stoi(argv[++i]);
            } catch (const invalid_argument&) {
                cerr << "Error: -scaling_rep flag expects a number. " << endl;
                exit(1);
            }
        }
        else if (arg == "-input") {
            if (i + 1 > argc || argv[i+1][0] == '-') {
                cerr << "Error: -input flag cannot be empty. " << endl;
                exit(1);
            }
            input_prefix = argv[++i];
        }
        else if (arg == "-output") {
            if (i + 1 > argc || argv[i+1][0] == '-') {
                cerr << "Error: -output flag cannot be empty. " << endl;
                exit(1);
            }
            output_prefix = argv[++i];
        }
        else if (arg == "-thin") {
            if (i + 1 > argc || argv[i+1][0] == '-') {
                cerr << "Error: -thinning flag cannot be empty. " << endl;
                exit(1);
            }
            try {
                spacing = stoi(argv[++i]);
            } catch (const invalid_argument&) {
                cerr << "Error: -thinning flag expects a number. " << endl;
                exit(1);
            }
        }
        else if (arg == "-seed") {
            if (i + 1 >= argc || argv[i+1][0] == '-') {
                cerr << "Error: -seed flag cannot be empty. " << endl;
                exit(1);
            }
            try {
                seed = stoi(argv[++i]);
            } catch (const invalid_argument&) {
                cerr << "Error: -seed flag expects a number. " << endl;
                exit(1);
            }
        }
        else if (arg == "-write_samples") {
            write_samples = true;
        }
        else if (arg == "-max_step") {
            if (i + 1 >= argc || argv[i+1][0] == '-') {
                cerr << "Error: -max_step flag cannot be empty. " << endl;
                exit(1);
            }
            try {
                 max_step = stod(argv[++i]);
            } catch (const invalid_argument&) {
                cerr << "Error: -max_step flag expects a number. " << endl;
                exit(1);
            }
        }
        else {
            cerr << "Error: Unknown flag. " << arg << endl;
            exit(1);
        }
    }
    DAG dag = DAG(Ne);
    dag.max_step = max_step;
    string node_file = input_prefix + "_nodes.txt";
    string branch_file = input_prefix + "_branches.txt";
    string mut_file = input_prefix + "_muts.txt";
    if (map_file.empty()) {
        dag.load_dag(node_file, branch_file);
        dag.compute_mutation_rates(m);
    } else {
        Mutation_map mm = Mutation_map();
        mm.load_map(map_file);
        dag.load_dag(node_file, branch_file, mm);
        m = mm.mean_rate();
    }
    dag.map_mutations(mut_file);
    for (int i = 0; i < burn_in; i++) {
        cout << "Burn-in iterations: " << i << endl;
        dag.no_prior_MCMC();
    }
    dag.burn_in();
    for (int i = 0; i < num_samples; i++) {
        cout << "MCMC iteration: " << i << endl;
        for (int j = 0; j < spacing; j++) {
            dag.no_prior_MCMC();
        }
        dag.record_node_ages();
    }
    dag.posterior_average();
    for (int i = 0; i < scaling_rep; i++) {
        Scaler scaler = Scaler();
        scaler.rescale(dag, Ne*m);
    }
    string new_node_file = input_prefix + "_new_nodes.txt";
    dag.write_node_ages(new_node_file);
    for (int i = 0; i < num_samples; i++) {
        cout << "Sample index: " << i << endl;
        dag.sample(i);
        for (int k = 0; k < scaling_rep; k++) {
            Scaler scaler = Scaler();
            scaler.rescale(dag, Ne*m);
        }
        dag.record_scaled_node_ages();
    }
    if (write_samples) {
        string node_samples_file = input_prefix + "_node_samples.txt";
        dag.write_node_age_samples(node_samples_file);
    }
    return 0;
}
