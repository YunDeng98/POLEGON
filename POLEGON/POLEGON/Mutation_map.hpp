//
//  Mutation_map.hpp
//  POLEGON
//
//  Created by Yun Deng on 1/19/24.
//

#ifndef Mutation_map_hpp
#define Mutation_map_hpp

#include <stdio.h>
#include <cassert>
#include <algorithm>
#include <memory>
#include <iostream>
#include <fstream>
#include <climits>
#include <vector>
#include <map>
#include <unordered_map>
#include <set>
#include <unordered_set>

using namespace std;

class Mutation_map {
    
public:
    
    double sequence_length = INT_MAX;
    vector<double> coordinates = {};
    vector<double> mutation_distance = {};
    
    Mutation_map();
    
    void load_map(string mut_map_file);
    
    int find_index(double x);
    
    double mut_distance(double x);
    
    double mut_rate_sum(double x, double y);

};

#endif /* Mutation_map_hpp */
