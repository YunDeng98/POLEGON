//
//  Distribution.hpp
//  arg_branch_length
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
    vector<float> times = {};
    vector<float> probs = {};
    vector<float> rates = {};
    
    Distribution(int n);
    
    void load_distribution(string filename);
    
    float propose(float lb, float ub);
    
// private:
    
    float survival(float x);
    
    float inverse_survival(float q);
    
};

#endif /* Distribution_hpp */
