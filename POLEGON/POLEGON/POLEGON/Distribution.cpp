//
//  Distribution.cpp
//  arg_branch_length
//
//  Created by Yun Deng on 10/31/23.
//

#include "Distribution.hpp"

Distribution::Distribution(int n) {
    num_samples = n;
}

void Distribution::load_distribution(string filename) {
    ifstream fin(filename);
    if (!fin.good()) {
        cerr << "input file not found" << endl;
        exit(1);
    }
    float x;
    float y;
    float r;
    while (fin >> x >> y >> r) {
        times.push_back(x);
        probs.push_back(y);
        rates.push_back(r);
    }
}

float Distribution::propose(float lb, float ub) {
    float lq = survival(lb);
    float uq = survival(ub);
    float r = uniform_random();
    float q = lq*r + uq*(1 - r);
    float x = inverse_survival(q);
    if (x <= lb or x >= ub) {
        return 0.5*(lb + ub);
    }
    return x;
}

float Distribution::survival(float x) {
    if (isinf(x)) {
        return 0;
    }
    auto it = upper_bound(times.begin(), times.end(), x);
    int index = (int) (it - times.begin());
    float rate = rates[index - 1];
    float delta = x - times[index - 1];
    float prop = exp(-rate * delta);
    float q = probs[index - 1]*prop;
    return q;
}

float Distribution::inverse_survival(float q) {
    if (q == 0) {
        return numeric_limits<float>::infinity();
    }
    auto it = upper_bound(probs.begin(), probs.end(), q, std::greater<float>());
    int index = (int) (it - probs.begin());
    float rate = rates[index - 1];
    float prop = probs[index - 1]/q;
    float delta = log(prop)/rate;
    float x = times[index - 1] + delta;
    return x;
}


