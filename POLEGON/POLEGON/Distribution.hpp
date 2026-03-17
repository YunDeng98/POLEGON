//
//  Distribution.hpp
//  POLEGON
//
//  Created by Yun Deng on 10/31/23.
//

#ifndef Distribution_hpp
#define Distribution_hpp

#include <stdio.h>
#include <map>
#include <iostream>
#include <algorithm>
#include "random_utils.hpp"

using namespace std;

class Distribution {

public:

    int num_samples = 0;
    vector<double> times = {};
    vector<double> probs = {};
    vector<double> rates = {};

    Distribution(int n);

    void load_distribution(string filename);
    double propose(double lb, double ub);

// private:

    double survival(double x);
    double inverse_survival(double q);

};

#endif /* Distribution_hpp */
