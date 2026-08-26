//
//  random_utils.hpp
//  POLEGON
//
//  Created by Yun Deng on 10/31/23.
//  Updated by Wonseop Lim on 03/21/26.
//

#ifndef random_utils_hpp
#define random_utils_hpp

#include <stdio.h>
#include <random>
#include <chrono>
#include <iomanip>
#include <fstream>
#include <sstream>

void init_random_streams(int seed, int num_streams);
void bind_random_stream(int index);
double uniform_random();

#endif /* random_utils_hpp */
