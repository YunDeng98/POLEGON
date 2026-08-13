//
//  Node.hpp
//  POLEGON
//
//  Created by Yun Deng on 10/31/23.
//

#ifndef Node_hpp
#define Node_hpp

#include <stdio.h>
#include <cassert>
#include <algorithm>
#include <memory>
#include <iostream>
#include <climits>

using namespace std;

class Node {

public:
    double time = 0;
    int index = 0;
    bool is_sample = false;
    Node(double t, int i);
};

bool compare_node(const Node *n1, const Node *n2);

#endif /* Node_hpp */
