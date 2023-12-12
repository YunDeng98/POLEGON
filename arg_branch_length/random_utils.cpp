//
//  random_utils.cpp
//  arg_branch_length
//
//  Created by Yun Deng on 10/31/23.
//

#include "random_utils.hpp"

std::mt19937 random_engine;
std::uniform_real_distribution<> uniform_distribution(0.0, 1.0);

float uniform_random() {
    float q = uniform_distribution(random_engine);
    if (q < 1e-5 or q > 1 - 1e-5) {
        q = uniform_distribution(random_engine);
    }
    return q;
}
