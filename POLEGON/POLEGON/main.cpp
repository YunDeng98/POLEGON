//
//  main.cpp
//  POLEGON
//
//  Created by Yun Deng on 12/13/23.
//

#include <iostream>
#include "Test.hpp"

int main(int argc, const char * argv[]) {
    float m = -1;
    int num_samples = -1;
    int spacing = 1000;
    string input_filename = "", output_prefix = "";
    int seed = 42;
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
        else if (arg == "-num_samples") {
            if (i + 1 >= argc || argv[i+1][0] == '-') {
                cerr << "Error: -penalty flag cannot be empty. " << endl;
                exit(1);
            }
            try {
                num_samples = stoi(argv[++i]);
            } catch (const invalid_argument&) {
                cerr << "Error: -num_samples flag expects a number. " << endl;
                exit(1);
            }
        }
        else if (arg == "-input") {
            if (i + 1 > argc || argv[i+1][0] == '-') {
                cerr << "Error: -input flag cannot be empty. " << endl;
                exit(1);
            }
            input_filename = argv[++i];
        }
        else if (arg == "-output") {
            if (i + 1 > argc || argv[i+1][0] == '-') {
                cerr << "Error: -output flag cannot be empty. " << endl;
                exit(1);
            }
            output_prefix = argv[++i];
        }
        else if (arg == "-thinning") {
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
        else {
            cerr << "Error: Unknown flag. " << arg << endl;
            exit(1);
        }
    }
    if (m < 0) {
        cerr << "-m flag missing or invalid value. " << endl;
        exit(1);
    }
    if (input_filename.size() == 0) {
        cerr << "-input flag missing or invalid value. " << endl;
        exit(1);
    }
    if (output_prefix.size() == 0) {
        cerr << "-output flag missing or invalid value. " << endl;
        exit(1);
    }
    if (num_samples < 0) {
        cerr << "-num_samples flag is invalid. " << endl;
        exit(1);
    }
    if (spacing < 1) {
        cerr << "-thinning flag is invalid. " << endl;
        exit(1);
    }
    return 0;
}

/*
int main(int argc, const char * argv[]) {
    // insert code here...
    // test_load_dag();
    // test_coalescent_prior();
    // test_sampling();
    // test_tsinfer_topology();
    // test_singer_topology();
    // test_singer_demo_topology();
    // test_no_prior_sampling();
    // test_scaling();
    // test_demography();
    // test_demo_scaling();
    // test_bottleneck();
    // test_bgs();
    test_pairwise_demo();
    std::cout << "Hello, World!\n";
    return 0;
}
*/
