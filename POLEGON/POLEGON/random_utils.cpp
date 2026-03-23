//
//  random_utils.cpp
//  POLEGON
//
//  Created by Yun Deng on 10/31/23.
//  Updated by Wonseop Lim on 03/21/26.
//

#include "random_utils.hpp"

thread_local std::mt19937 random_engine;
thread_local std::uniform_real_distribution<> uniform_distribution(0.0, 1.0);

void seed_random_engine(int seed, int tid) {
    std::seed_seq seq{seed, tid};
    random_engine.seed(seq);
}

double uniform_random() {
    double q = uniform_distribution(random_engine);
    if (q < 1e-5 or q > 1 - 1e-5) {
        q = uniform_distribution(random_engine);
    }
    return q;
}
