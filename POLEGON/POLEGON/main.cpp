//
//  main.cpp
//  POLEGON
//
//  Created by Yun Deng on 12/13/23.
//  Modified by Wonseop Lim on 08/15/26.
//

#include <iostream>
#include <charconv>
#include <cstdio>
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
    int scaling_subsample = 10; // number of posterior samples used to estimate the shared rescaling grid
    int num_cores = 1;
    string input_prefix = "", output_prefix = "";
    int seed = 42;              // random seed
    double Ne = 0;              // effective population size
    string map_file = "";       // path to mutation map
    string tip_ages_file = "";  // path to sample ages
    bool posterior_mean = true;   // whether to compute and write the posterior mean
    bool memory_safe = false;     // stream unrescaled samples to disk, rescale in batches

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
        else if (arg == "-scaling_subsample") {
            if (i + 1 >= argc || argv[i+1][0] == '-') {
                cerr << "Error: -scaling_subsample flag cannot be empty. " << endl; exit(1);
            }
            try { scaling_subsample = stoi(argv[++i]); }
            catch (const invalid_argument&) {
                cerr << "Error: -scaling_subsample flag expects a number. " << endl; exit(1);
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
        else if (arg == "-no_mean") {
            posterior_mean = false;
        }
        else if (arg == "-memory_safe") {
            memory_safe = true;
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

    dag.compute_root_lambda();
    dag.compute_non_root_lambda();
    dag.compute_internal_by_time();
    dag.compute_coloring();
    init_random_streams(seed, dag.num_streams);

    // Burn-in
    cout << "Burn-in Phase..." << endl;
    for (int i = 0; i < burn_in; i++) {
        dag.no_prior_MCMC();
    }
    cout << "Done" << endl;

    int n_nodes = (int)dag.nodes.size();
    int total_mcmc_iters = num_samples * spacing;

    // ── memory-safe path: stream unrescaled samples to disk, rescale in batches ──
    if (memory_safe) {
        string unrescaled_file = output_prefix + "_unrescaled_node_samples.txt";
        string raw_file = scaling_rep == 0 ? unrescaled_file
                                           : output_prefix + "_unrescaled_node_samples.tmp";
        {
            ofstream raw_out(raw_file);
            for (int i = 0; i < num_samples; i++) {
                for (int j = 0; j < spacing; j++)
                    dag.no_prior_MCMC();
                int done = (i + 1) * spacing;
                if (done % 100 == 0 || i + 1 == num_samples)
                    cout << "MCMC Iterations: " << done << "/" << total_mcmc_iters << endl;
                for (int j = 0; j < n_nodes; j++) {
                    double t = scaling_rep == 0
                             ? dag.output_time(j, (dag.nodes[j]->time + dag.time_origin) * Ne * g)
                             : dag.nodes[j]->time;
                    raw_out << std::setprecision(std::numeric_limits<double>::max_digits10)
                            << t << " ";
                }
                raw_out << "\n";
            }
        }
        if (scaling_rep == 0) return 0;

        ofstream samples_file(output_prefix + "_node_samples.txt");
        vector<double> sums;
        if (posterior_mean) sums.assign(n_nodes, 0.0);
        vector<vector<double>> subsample;
        {
            ifstream raw_in(raw_file);
            ofstream unrescaled_out(unrescaled_file);
            vector<double> row(n_nodes);
            int K = min(scaling_subsample, num_samples);
            int subsample_spacing = max(1, num_samples/K);
            for (int s = 0; s < num_samples; s++) {
                for (int j = 0; j < n_nodes; j++)
                    raw_in >> row[j];
                for (int j = 0; j < n_nodes; j++)
                    unrescaled_out << std::setprecision(std::numeric_limits<double>::max_digits10)
                                   << dag.output_time(j, (row[j] + dag.time_origin) * Ne * g) << " ";
                unrescaled_out << "\n";
                if (s % subsample_spacing == 0 && (int)subsample.size() < K) subsample.push_back(row);
            }
        }

        vector<Scaler> scalers;
        {
            Scaler scaler;
            scaler.num_bins = scaling_bin;
            scaler.num_cores = dag.num_cores;
            scaler.compute_deltas(dag);
            for (int k = 0; k < scaling_rep; k++) {
                scaler.rescale(dag, subsample, (int)subsample.size(), Ne * m);
                Scaler round;
                round.num_bins = scaling_bin;
                round.old_grid = scaler.old_grid;
                round.new_grid = scaler.new_grid;
                round.scaling_factors = scaler.scaling_factors;
                scalers.push_back(round);
                cout << "ARG Rescaling: " << k + 1 << "/" << scaling_rep << endl;
            }
        }

        {
            ifstream raw_in(raw_file);
            vector<double> row(n_nodes);
            for (int s = 0; s < num_samples; s++) {
                for (int j = 0; j < n_nodes; j++)
                    raw_in >> row[j];
                for (Scaler &st : scalers)
                    st.apply_scaling_factors(dag, row);
                for (int j = 0; j < n_nodes; j++) {
                    double t = dag.output_time(j, (row[j] + dag.time_origin) * Ne * g);
                    samples_file << std::setprecision(std::numeric_limits<double>::max_digits10)
                                 << t << " ";
                    if (posterior_mean) sums[j] += t;
                }
                samples_file << "\n";
            }
        }
        samples_file.close();
        std::remove(raw_file.c_str());
        if (posterior_mean) {
            ofstream fout(input_prefix + "_posterior_mean_nodes.txt");
            for (int j = 0; j < n_nodes; j++)
                fout << std::setprecision(std::numeric_limits<double>::max_digits10)
                     << dag.output_time(j, sums[j] / num_samples) << "\n";
        }
        return 0;
    }

    // ── default path: hold all samples in memory, rescale fully in parallel ──
    vector<vector<double>> all_raw(num_samples, vector<double>(n_nodes));
    for (int i = 0; i < num_samples; i++) {
        for (int j = 0; j < spacing; j++)
            dag.no_prior_MCMC();
        int done = (i + 1) * spacing;
        if (done % 100 == 0 || i + 1 == num_samples)
            cout << "MCMC Iterations: " << done << "/" << total_mcmc_iters << endl;
        for (int j = 0; j < n_nodes; j++)
            all_raw[i][j] = dag.nodes[j]->time;
    }

    if (scaling_rep == 0) {
        ofstream raw_out(output_prefix + "_unrescaled_node_samples.txt");
        for (int i = 0; i < num_samples; i++) {
            for (int j = 0; j < n_nodes; j++)
                raw_out << std::setprecision(std::numeric_limits<double>::max_digits10)
                        << dag.output_time(j, (all_raw[i][j] + dag.time_origin) * Ne * g) << " ";
            raw_out << "\n";
        }
        return 0;
    }

    Scaler scaler;
    scaler.num_bins = scaling_bin;
    scaler.num_cores = dag.num_cores;
    scaler.compute_deltas(dag);
    for (int k = 0; k < scaling_rep; k++) {
        scaler.rescale(dag, all_raw, scaling_subsample, Ne * m);
        cout << "ARG Rescaling: " << k + 1 << "/" << scaling_rep << endl;
    }

    ofstream samples_file(output_prefix + "_node_samples.txt");
    vector<double> sums;
    if (posterior_mean) sums.assign(n_nodes, 0.0);
    {
        string out;
        out.reserve(1 << 22);
        char buf[64];
        for (int i = 0; i < num_samples; i++) {
            for (int j = 0; j < n_nodes; j++) {
                double t = dag.output_time(j, (all_raw[i][j] + dag.time_origin) * Ne * g);
                auto r = std::to_chars(buf, buf + sizeof buf, t, std::chars_format::general,
                                       std::numeric_limits<double>::max_digits10);
                out.append(buf, r.ptr - buf);
                out.push_back(' ');
                if (posterior_mean) sums[j] += t;
            }
            out.push_back('\n');
            if (out.size() >= (1 << 21)) {
                samples_file.write(out.data(), out.size());
                out.clear();
            }
        }
        samples_file.write(out.data(), out.size());
    }
    samples_file.close();
    
    if (posterior_mean) {
        string new_node_file = input_prefix + "_posterior_mean_nodes.txt";
        ofstream fout(new_node_file);
        for (int j = 0; j < (int)dag.nodes.size(); j++)
            fout << std::setprecision(std::numeric_limits<double>::max_digits10)
                 << dag.output_time(j, sums[j] / num_samples) << "\n";
        fout.close();
    }

    return 0;
}
