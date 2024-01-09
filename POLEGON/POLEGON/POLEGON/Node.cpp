//
//  Node.cpp
//  arg_branch_length
//
//  Created by Yun Deng on 10/31/23.
//

#include "Node.hpp"

Node::Node(float t, int i) {
    time = t;
    index = i;
}

bool compare_node(const Node *n1, const Node *n2) {
    if (n1->time != n2->time) {
        return n1->time < n2->time;
    } else {
        return n1->index < n2->index;
    }
}

