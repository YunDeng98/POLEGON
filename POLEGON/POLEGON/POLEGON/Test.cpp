//
//  Test.cpp
//  POLEGON
//
//  Created by Yun Deng on 10/31/23.
//

#include "Test.hpp"

/*
void test_load_dag() {
    DAG dag = DAG(1e4);
    dag.load_dag("/Users/yun_deng/Desktop/POLEGON/arg_files/constant_nodes.txt", "/Users/yun_deng/Desktop/POLEGON/arg_files/constant_branches.txt", "/Users/yun_deng/Desktop/POLEGON/arg_files/constant_muts.txt");
    dag.compute_mutation_rates(2e-8);
}


void test_coalescent_prior() {
    Distribution *d = new Distribution(10);
    d->load_distribution("/Users/yun_deng/Desktop/POLEGON/arg_files/distribution_50.txt");
    vector<float> random_values = vector<float>(10000);
    for (int i = 0; i< 10000; i++) {
        random_values[i] = d->propose(0, numeric_limits<float>::infinity());
    }
    ofstream output_file("/Users/yun_deng/Desktop/POLEGON/arg_files/coalescence_time_samples.txt");
    for (float x : random_values) {
        output_file << x << endl;
    }
}



void test_sampling() {
    DAG dag = DAG(2e4);
    dag.load_dag("/Users/yun_deng/Desktop/POLEGON/arg_files/constant_nodes.txt", "/Users/yun_deng/Desktop/POLEGON/arg_files/constant_branches.txt", "/Users/yun_deng/Desktop/POLEGON/arg_files/constant_muts.txt");
    dag.compute_mutation_rates(1.2e-8);
    Distribution *d = new Distribution(10);
    d->load_distribution("/Users/yun_deng/Desktop/POLEGON/arg_files/distribution_50.txt");
    for (int i = 0; i < 100; i++) {
        string node_file = "/Users/yun_deng/Desktop/POLEGON/arg_files/constant_new_nodes_" + to_string(i) + ".txt";
        dag.MCMC(10000, d);
        dag.write_node_ages(node_file);
    }
    dag.posterior_average();
    dag.write_node_ages("/Users/yun_deng/Desktop/POLEGON/arg_files/constant_new_nodes.txt");
}

void test_tsinfer_topology() {
    DAG dag = DAG(2e4);
    dag.load_dag("/Users/yun_deng/Desktop/POLEGON/arg_files/tsinfer_nodes.txt", "/Users/yun_deng/Desktop/POLEGON/arg_files/tsinfer_branches.txt", "/Users/yun_deng/Desktop/POLEGON/arg_files/tsinfer_muts.txt");
    dag.compute_mutation_rates(1.2e-8);
    Distribution *d = new Distribution(10);
    d->load_distribution("/Users/yun_deng/Desktop/POLEGON/arg_files/distribution_50.txt");
    for (int i = 0; i < 100; i++) {
        string node_file = "/Users/yun_deng/Desktop/POLEGON/arg_files/tsinfer_new_nodes_" + to_string(i) + ".txt";
        dag.MCMC(1000, d);
        dag.write_node_ages(node_file);
    }
}

void test_singer_topology() {
    DAG dag = DAG(2e4);
    dag.load_dag("/Users/yun_deng/Desktop/POLEGON/arg_files/50_0_nodes_0.txt", "/Users/yun_deng/Desktop/POLEGON/arg_files/50_0_branches_0.txt", "/Users/yun_deng/Desktop/POLEGON/arg_files/50_0_muts_0.txt");
    dag.compute_mutation_rates(2e-8);
    for (int i = 0; i < 100; i++) {
        cout << "MCMC iteration: " << i << endl;
        string node_file = "/Users/yun_deng/Desktop/POLEGON/arg_files/50_0_new_nodes_" + to_string(i) + ".txt";
        dag.no_prior_MCMC(30000);
        dag.write_node_ages(node_file);
    }
    dag.posterior_average();
    Scaler scaler = Scaler();
    scaler.rescale(dag, 4e-4);
    dag.write_node_ages("/Users/yun_deng/Desktop/POLEGON/arg_files/50_0_new_nodes.txt");
}

void test_singer_demo_topology() {
    DAG dag = DAG(2e4);
    dag.load_dag("/Users/yun_deng/Desktop/POLEGON/arg_files/50_0_demo_nodes_599.txt", "/Users/yun_deng/Desktop/POLEGON/arg_files/50_0_demo_branches_599.txt", "/Users/yun_deng/Desktop/POLEGON/arg_files/50_0_demo_muts_599.txt");
    dag.compute_mutation_rates(2e-8);
    for (int i = 0; i < 100; i++) {
        cout << "MCMC iteration: " << i << endl;
        string node_file = "/Users/yun_deng/Desktop/POLEGON/arg_files/50_0_demo_new_nodes_" + to_string(i) + ".txt";
        dag.no_prior_MCMC(30000);
        dag.write_node_ages(node_file);
    }
    dag.posterior_average();
    Scaler scaler = Scaler();
    scaler.rescale(dag, 4e-4);
    dag.write_node_ages("/Users/yun_deng/Desktop/POLEGON/arg_files/50_0_demo_new_nodes.txt");
}

void test_no_prior_sampling() {
    DAG dag = DAG(2e4);
    dag.load_dag("/Users/yun_deng/Desktop/POLEGON/arg_files/constant_nodes.txt", "/Users/yun_deng/Desktop/POLEGON/arg_files/constant_branches.txt", "/Users/yun_deng/Desktop/POLEGON/arg_files/constant_muts.txt");
    dag.compute_mutation_rates(1.2e-8);
    for (int i = 0; i < 100; i++) {
        string node_file = "/Users/yun_deng/Desktop/POLEGON/arg_files/no_prior_nodes_" + to_string(i) + ".txt";
        dag.no_prior_MCMC(10000);
        dag.write_node_ages(node_file);
    }
    dag.posterior_average();
    Scaler scaler = Scaler();
    scaler.rescale(dag, 2.4e-4);
    dag.write_node_ages("/Users/yun_deng/Desktop/POLEGON/arg_files/no_prior_nodes.txt");
}

void test_scaling() {
    DAG dag = DAG(2e4);
    dag.load_dag("/Users/yun_deng/Desktop/POLEGON/arg_files/no_prior_nodes.txt", "/Users/yun_deng/Desktop/POLEGON/arg_files/tree_branches.txt", "/Users/yun_deng/Desktop/POLEGON/arg_files/tree_muts.txt");
    dag.compute_mutation_rates(1.2e-8);
    Scaler scaler = Scaler();
    scaler.rescale(dag, 2.4e-4);
    dag.write_node_ages("/Users/yun_deng/Desktop/POLEGON/arg_files/rescaled_nodes.txt");
}

void test_demography() {
    DAG dag = DAG(2e4);
    dag.load_dag("/Users/yun_deng/Desktop/POLEGON/arg_files/demo_nodes.txt", "/Users/yun_deng/Desktop/POLEGON/arg_files/demo_branches.txt", "/Users/yun_deng/Desktop/POLEGON/arg_files/demo_muts.txt");
    dag.compute_mutation_rates(1.2e-8);
    for (int i = 0; i < 500; i++) {
        dag.no_prior_MCMC(50000);
    }
    dag.burn_in();
    for (int i = 0; i < 100; i++) {
        cout << "MCMC iteration: " << endl;
        string node_file = "/Users/yun_deng/Desktop/POLEGON/arg_files/demo_new_nodes_" + to_string(i) + ".txt";
        dag.no_prior_MCMC(100000);
        dag.write_node_ages(node_file);
    }
    dag.posterior_average();
    Scaler scaler = Scaler();
    scaler.rescale(dag, 2.4e-4);
    dag.write_node_ages("/Users/yun_deng/Desktop/POLEGON/arg_files/demo_new_nodes.txt");
}

void test_demo_scaling() {
    DAG dag = DAG(2e4);
    dag.load_dag("/Users/yun_deng/Desktop/POLEGON/arg_files/demo_new_nodes.txt", "/Users/yun_deng/Desktop/POLEGON/arg_files/demo_tree_branches.txt", "/Users/yun_deng/Desktop/POLEGON/arg_files/demo_tree_muts.txt");
    dag.compute_mutation_rates(1.2e-8);
    Scaler scaler = Scaler();
    scaler.rescale(dag, 2.4e-4);
    dag.write_node_ages("/Users/yun_deng/Desktop/POLEGON/arg_files/rescaled_demo_nodes.txt");
}

void test_bottleneck() {
    DAG dag = DAG(2e4);
    dag.load_dag("/Users/yun_deng/Desktop/POLEGON/arg_files/bottleneck_nodes.txt", "/Users/yun_deng/Desktop/POLEGON/arg_files/bottleneck_branches.txt", "/Users/yun_deng/Desktop/POLEGON/arg_files/bottleneck_muts.txt");
    dag.compute_mutation_rates(1.2e-8);
    for (int i = 0; i < 300; i++) {
        dag.no_prior_MCMC(10000);
    }
    dag.burn_in();
    for (int i = 0; i < 100; i++) {
        cout << "MCMC iteration: " << i << endl;
        string node_file = "/Users/yun_deng/Desktop/POLEGON/arg_files/bottleneck_new_nodes_" + to_string(i) + ".txt";
        dag.no_prior_MCMC(20000);
        dag.write_node_ages(node_file);
    }
    dag.posterior_average();
    Scaler scaler = Scaler();
    scaler.rescale(dag, 2.4e-4);
    dag.write_node_ages("/Users/yun_deng/Desktop/POLEGON/arg_files/bottleneck_new_nodes.txt");
}

void test_bgs() {
    DAG dag = DAG(2e3);
    dag.load_dag("/Users/yun_deng/Desktop/POLEGON/arg_files/bgs_nodes.txt", "/Users/yun_deng/Desktop/POLEGON/arg_files/bgs_branches.txt", "/Users/yun_deng/Desktop/POLEGON/arg_files/bgs_muts.txt");
    dag.compute_mutation_rates(1.5e-8);
    for (int i = 0; i < 300; i++) {
        dag.no_prior_MCMC(30000);
    }
    dag.burn_in();
    for (int i = 0; i < 100; i++) {
        cout << "MCMC iteration: " << i << endl;
        string node_file = "/Users/yun_deng/Desktop/POLEGON/arg_files/bgs_new_nodes_" + to_string(i) + ".txt";
        dag.no_prior_MCMC(100000);
        dag.write_node_ages(node_file);
    }
    dag.posterior_average();
    Scaler scaler = Scaler();
    scaler.rescale(dag, 3e-5);
    dag.write_node_ages("/Users/yun_deng/Desktop/POLEGON/arg_files/bgs_new_nodes.txt");
}

void test_pairwise_demo() {
    DAG dag = DAG(2e4);
    dag.load_dag("/Users/yun_deng/Desktop/POLEGON/arg_files/pair_nodes.txt", "/Users/yun_deng/Desktop/POLEGON/arg_files/pair_branches.txt", "/Users/yun_deng/Desktop/POLEGON/arg_files/pair_muts.txt");
    dag.compute_mutation_rates(1.2e-8);
    for (int i = 0; i < 500; i++) {
        dag.no_prior_MCMC(10000);
    }
    dag.burn_in();
    for (int i = 0; i < 100; i++) {
        cout << "MCMC iteration: " << i << endl;
        string node_file = "/Users/yun_deng/Desktop/POLEGON/arg_files/pair_new_nodes_" + to_string(i) + ".txt";
        dag.no_prior_MCMC(100000);
        dag.write_node_ages(node_file);
    }
    dag.posterior_average();
    Scaler scaler = Scaler();
    scaler.rescale(dag, 2.4e-4);
    dag.write_node_ages("/Users/yun_deng/Desktop/POLEGON/arg_files/pair_new_nodes.txt");
}

void test_migration() {
    DAG dag = DAG(2e4);
    dag.load_dag("/Users/yun_deng/Desktop/POLEGON/arg_files/migration_start_nodes_0.txt", "/Users/yun_deng/Desktop/POLEGON/arg_files/migration_start_branches_0.txt", "/Users/yun_deng/Desktop/POLEGON/arg_files/migration_start_muts_0.txt");
    // dag.load_dag("/Users/yun_deng/Desktop/POLEGON/arg_files/migration_nodes_99.txt", "/Users/yun_deng/Desktop/POLEGON/arg_files/migration_branches_99.txt", "/Users/yun_deng/Desktop/POLEGON/arg_files/migration_muts_99.txt");
    dag.compute_mutation_rates(1e-8);
    for (int i = 0; i < 500; i++) {
        dag.no_prior_MCMC(10000);
    }
    dag.burn_in();
    for (int i = 0; i < 100; i++) {
        cout << "MCMC iteration: " << i << endl;
        // string node_file = "/Users/yun_deng/Desktop/POLEGON/arg_files/migration_new_nodes_" + to_string(i) + ".txt";
        dag.no_prior_MCMC(100000);
        // dag.write_node_ages(node_file);
    }
    dag.posterior_average();
    // dag.write_node_ages("/Users/yun_deng/Desktop/POLEGON/arg_files/migration_raw_nodes_99.txt");
    Scaler scaler = Scaler();
    scaler.rescale(dag, 2e-4);
    dag.write_node_ages("/Users/yun_deng/Desktop/POLEGON/arg_files/migration_new_nodes_0.txt");
    // dag.write_node_ages("/Users/yun_deng/Desktop/POLEGON/arg_files/migration_new_nodes_99.txt");
}
*/
