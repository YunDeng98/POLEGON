//
//  random_utils.hpp
//  arg_branch_length
//
//  Created by Yun Deng on 10/31/23.
//

#ifndef random_utils_hpp
#define random_utils_hpp

#include <stdio.h>
#include <random>
#include <chrono>
#include <iomanip>
#include <fstream>
#include <sstream>

extern std::mt19937 random_engine;
extern std::uniform_real_distribution<> uniform_distribution;

float uniform_random();

#endif /* random_utils_hpp */
