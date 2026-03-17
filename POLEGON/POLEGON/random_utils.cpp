//
//  random_utils.cpp
//  POLEGON
//
//  Created by Yun Deng on 10/31/23.
//

#include "random_utils.hpp"

std::mt19937 random_engine;
std::uniform_real_distribution<> uniform_distribution(0.0, 1.0);

void seed_random_engine(int seed) {
    random_engine.seed(seed);
}

double uniform_random() {
    double q = uniform_distribution(random_engine);
    if (q < 1e-5 or q > 1 - 1e-5) {
        q = uniform_distribution(random_engine);
    }
    return q;
}
