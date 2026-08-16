//
//  random_utils.cpp
//  POLEGON
//
//  Created by Yun Deng on 10/31/23.
//  Updated by Wonseop Lim on 08/15/26.
//

#include <vector>
#include "random_utils.hpp"

static std::vector<std::mt19937> random_streams;
static thread_local int active_stream = 0;
static thread_local std::uniform_real_distribution<> uniform_distribution(0.0, 1.0);

void init_random_streams(int seed, int num_streams) {
    random_streams.resize(num_streams);
    for (int i = 0; i < num_streams; i++) {
        std::seed_seq seq{seed, i};
        random_streams[i].seed(seq);
    }
}

void bind_random_stream(int index) {
    active_stream = index;
}

double uniform_random() {
    return uniform_distribution(random_streams[active_stream]);
}
