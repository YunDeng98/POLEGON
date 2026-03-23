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

void seed_random_engine(int seed, int tid = 0);
double uniform_random(); // redraws once if result is within 1e-5 of 0 or 1

#endif /* random_utils_hpp */
