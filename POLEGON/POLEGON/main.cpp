//
//  main.cpp
//  POLEGON
//
//  Created by Yun Deng on 12/13/23.
//  Modified by Wonseop Lim on 04/03/26.
//

#include <iostream>
#include <omp.h>
#include "DAG.hpp"
#include "Distribution.hpp"
#include "Scaler.hpp"
#include "Mutation_map.hpp"
#include "random_utils.hpp"

int main(int argc, const char * argv[]) {
    double m = -1;              // mutation rate per generation per bp
    double g = -1;              // generation time in years; -1 = not set (defaults to 1 if unused)
    int num_samples = -1;       // number of posterior MCMC samples
    int burn_in = -1;           // number of burn-in samples
    int spacing = -1;           // thinning interval
    int scaling_rep = 5;        // number of ARG rescaling rounds
    int scaling_bin = 100;      // number of time bins used by the Scaler
    double max_step = 10.0;     // maximum exponential draw for root node proposals
    int num_cores = 1;
    string input_prefix = "", output_prefix = "";
    int seed = 42;              // random seed
    double Ne = 0;              // effective population size
    string map_file = "";       // path to mutation map
    string tip_ages_file = "";  // path to sample ages
    bool posterior_mean = true;   // whether to compute and write the posterior mean

    for (int i = 1; i < argc; ++i) {
        string arg = argv[i];
        if (arg == "-m") {
            if (i + 1 >= argc || argv[i+1][0] == '-') {
                cerr << "Error: -m flag cannot be empty. " << endl; exit(1);
            }
            try { m = stod(argv[++i]); }
            catch (const invalid_argument&) {
                cerr << "Error: -m flag expects a number. " << endl; exit(1);
            }
        }
        else if (arg == "-m_map") {
            if (i + 1 >= argc || argv[i+1][0] == '-') {
                cerr << "Error: -m_map flag cannot be empty. " << endl; exit(1);
            }
            try { map_file = argv[++i]; }
            catch (const invalid_argument&) {
                cerr << "Error: -m_map flag expects a string. " << endl; exit(1);
            }
        }
        else if (arg == "-Ne") {
            if (i + 1 >= argc || argv[i+1][0] == '-') {
                cerr << "Error: -Ne flag cannot be empty. " << endl; exit(1);
            }
            try { Ne = stod(argv[++i]); }
            catch (const invalid_argument&) {
                cerr << "Error: -Ne flag expects a number. " << endl; exit(1);
            }
        }
        else if (arg == "-burn_in") {
            if (i + 1 >= argc || argv[i+1][0] == '-') {
                cerr << "Error: -burn_in flag cannot be empty. " << endl; exit(1);
            }
            try { burn_in = stoi(argv[++i]); }
            catch (const invalid_argument&) {
                cerr << "Error: -burn_in flag expects a number. " << endl; exit(1);
            }
        }
        else if (arg == "-n_samples") {
            if (i + 1 >= argc || argv[i+1][0] == '-') {
                cerr << "Error: -n_samples flag cannot be empty. " << endl; exit(1);
            }
            try { num_samples = stoi(argv[++i]); }
            catch (const invalid_argument&) {
                cerr << "Error: -n_samples flag expects a number. " << endl; exit(1);
            }
        }
        else if (arg == "-scaling_rep") {
            if (i + 1 >= argc || argv[i+1][0] == '-') {
                cerr << "Error: -scaling_rep flag cannot be empty. " << endl; exit(1);
            }
            try { scaling_rep = stoi(argv[++i]); }
            catch (const invalid_argument&) {
                cerr << "Error: -scaling_rep flag expects a number. " << endl; exit(1);
            }
        }
        else if (arg == "-scaling_bin") {
            if (i + 1 >= argc || argv[i+1][0] == '-') {
                cerr << "Error: -scaling_bin flag cannot be empty. " << endl; exit(1);
            }
            try { scaling_bin = stoi(argv[++i]); }
            catch (const invalid_argument&) {
                cerr << "Error: -scaling_bin flag expects a number. " << endl; exit(1);
            }
        }
        else if (arg == "-input") {
            if (i + 1 > argc || argv[i+1][0] == '-') {
                cerr << "Error: -input flag cannot be empty. " << endl; exit(1);
            }
            input_prefix = argv[++i];
        }
        else if (arg == "-output") {
            if (i + 1 > argc || argv[i+1][0] == '-') {
                cerr << "Error: -output flag cannot be empty. " << endl; exit(1);
            }
            output_prefix = argv[++i];
        }
        else if (arg == "-thin") {
            if (i + 1 > argc || argv[i+1][0] == '-') {
                cerr << "Error: -thinning flag cannot be empty. " << endl; exit(1);
            }
            try { spacing = stoi(argv[++i]); }
            catch (const invalid_argument&) {
                cerr << "Error: -thinning flag expects a number. " << endl; exit(1);
            }
        }
        else if (arg == "-seed") {
            if (i + 1 >= argc || argv[i+1][0] == '-') {
                cerr << "Error: -seed flag cannot be empty. " << endl; exit(1);
            }
            try { seed = stoi(argv[++i]); }
            catch (const invalid_argument&) {
                cerr << "Error: -seed flag expects a number. " << endl; exit(1);
            }
        }
        else if (arg == "-no_posterior_mean") {
            posterior_mean = false;
        }
        else if (arg == "-tip_ages") {
            if (i + 1 >= argc || argv[i+1][0] == '-') {
                cerr << "Error: -tip_ages flag cannot be empty. " << endl; exit(1);
            }
            tip_ages_file = argv[++i];
        }
        else if (arg == "-g") {
            if (i + 1 >= argc || argv[i+1][0] == '-') {
                cerr << "Error: -g flag cannot be empty. " << endl; exit(1);
            }
            try { g = stod(argv[++i]); }
            catch (const invalid_argument&) {
                cerr << "Error: -g flag expects a number. " << endl; exit(1);
            }
        }
        else if (arg == "-cores") {
            if (i + 1 >= argc || argv[i+1][0] == '-') {
                cerr << "Error: -cores flag cannot be empty." << endl; exit(1);
            }
            try { num_cores = stoi(argv[++i]); }
            catch (const invalid_argument&) {
                cerr << "Error: -cores flag expects a number." << endl; exit(1);
            }
        }
        else if (arg == "-max_step") {
            if (i + 1 >= argc || argv[i+1][0] == '-') {
                cerr << "Error: -max_step flag cannot be empty. " << endl; exit(1);
            }
            try { max_step = stod(argv[++i]); }
            catch (const invalid_argument&) {
                cerr << "Error: -max_step flag expects a number. " << endl; exit(1);
            }
        }
        else {
            cerr << "Error: Unknown flag. " << arg << endl; exit(1);
        }
    }

    if (!tip_ages_file.empty() && g == -1) {
        cerr << "Error: -g (generation time in years) must be provided when using -tip_ages." << endl;
        exit(1);
    }
    if (g != -1 && g <= 0) {
        cerr << "Error: -g must be positive." << endl;
        exit(1);
    }
    if (g == -1) g = 1;

    DAG dag = DAG(Ne);
    dag.max_step = max_step;
    dag.num_cores = num_cores;
    string node_file   = input_prefix + "_nodes.txt";
    string branch_file = input_prefix + "_branches.txt";
    string mut_file    = input_prefix + "_muts.txt";

    if (map_file.empty()) {
        // Constant mutation rate: load topology and set uniform branch mutation rates
        dag.load_dag(node_file, branch_file);
        dag.compute_mutation_rates(m);
    } else {
        // Spatially varying mutation rate: load mutation map, then load ARG with per-branch rates
        Mutation_map mm = Mutation_map();
        mm.load_map(map_file);
        dag.load_dag(node_file, branch_file, mm);
        m = mm.mean_rate(); // Use the genome-wide mean rate for Scaler calibration
    }

    // For heterochronous data, shift the time origin to the youngest tip, set each tip's
    // age from the supplied calendar ages, and repair any internal node whose
    // initialisation falls below its oldest child via a single DAG up-pass
    if (!tip_ages_file.empty()) {
        dag.apply_tip_ages(tip_ages_file, g);
    }

    // Map observed mutations from the mutations file
    dag.map_mutations(mut_file);

    dag.compute_coloring();
    #pragma omp parallel num_threads(num_cores)
    {
        seed_random_engine(seed, omp_get_thread_num());
    }

    // Burn-in
    cout << "Burn-in Phase..." << endl;
    for (int i = 0; i < burn_in; i++) {
        dag.no_prior_MCMC();
    }
    cout << "Done" << endl;

    // MCMC phase: collect all raw samples
    int n_nodes = (int)dag.nodes.size();
    vector<vector<double>> all_raw(num_samples, vector<double>(n_nodes));
    int total_mcmc_iters = num_samples * spacing;
    for (int i = 0; i < num_samples; i++) {
        for (int j = 0; j < spacing; j++)
            dag.no_prior_MCMC();
        int done = (i + 1) * spacing;
        if (done % 100 == 0 || i + 1 == num_samples)
            cout << "MCMC Iterations: " << done << "/" << total_mcmc_iters << endl;
        for (int j = 0; j < n_nodes; j++)
            all_raw[i][j] = dag.nodes[j]->time;
    }

    // Rescaling phase: rescale all samples in parallel
    if (scaling_rep > 0) {
        int done_count = 0;
        #pragma omp parallel for schedule(dynamic,1) num_threads(dag.num_cores)
        for (int s = 0; s < num_samples; s++) {
            Scaler scaler;
            scaler.num_bins = scaling_bin;
            scaler.local_times = all_raw[s];
            for (int k = 0; k < scaling_rep; k++)
                scaler.rescale(dag, Ne * m);
            all_raw[s] = scaler.local_times;
            int cnt;
            #pragma omp atomic capture
            cnt = ++done_count;
            if (cnt % 10 == 0 || cnt == num_samples) {
                #pragma omp critical
                cout << "ARG Rescaling: " << cnt << "/" << num_samples << endl;
            }
        }
    }

    // Output phase
    ofstream samples_file;
    string node_samples_file = output_prefix + "_node_samples.txt";
    samples_file.open(node_samples_file);
    vector<double> sums;
    if (posterior_mean) sums.assign(n_nodes, 0.0);
    for (int i = 0; i < num_samples; i++) {
        for (int j = 0; j < n_nodes; j++) {
            double t = (all_raw[i][j] + dag.time_origin) * Ne * g;
            samples_file << std::setprecision(std::numeric_limits<double>::max_digits10)
                         << t << " ";
            if (posterior_mean) sums[j] += t;
        }
        samples_file << "\n";
    }
    samples_file.close();
    
    if (posterior_mean) {
        string new_node_file = input_prefix + "_new_nodes.txt";
        ofstream fout(new_node_file);
        for (int j = 0; j < (int)dag.nodes.size(); j++)
            fout << std::setprecision(std::numeric_limits<double>::max_digits10)
                 << sums[j] / num_samples << "\n";
        fout.close();
    }

    return 0;
}
