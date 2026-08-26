//
//  Distribution.cpp
//  POLEGON
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
    double x;
    double y;
    double r;
    while (fin >> x >> y >> r) {
        times.push_back(x);
        probs.push_back(y);
        rates.push_back(r);
    }
}

double Distribution::propose(double lb, double ub) {
    double lq = survival(lb);
    double uq = survival(ub);
    double r = uniform_random();
    double q = lq*r + uq*(1 - r);
    double x = inverse_survival(q);
    if (x <= lb or x >= ub) {
        return 0.5*(lb + ub);
    }
    return x;
}

double Distribution::survival(double x) {
    if (isinf(x)) {
        return 0;
    }
    auto it = upper_bound(times.begin(), times.end(), x);
    int index = max(1, (int) (it - times.begin()));
    double rate = rates[index - 1];
    double delta = x - times[index - 1];
    double prop = exp(-rate * delta);
    double q = probs[index - 1]*prop;
    return q;
}

double Distribution::inverse_survival(double q) {
    if (q == 0) {
        return numeric_limits<double>::infinity();
    }
    auto it = upper_bound(probs.begin(), probs.end(), q, std::greater<double>());
    int index = (int) (it - probs.begin());
    double rate = rates[index - 1];
    double prop = probs[index - 1]/q;
    double delta = log(prop)/rate;
    double x = times[index - 1] + delta;
    return x;
}
