//
//  Branch.cpp
//  POLEGON
//
//  Created by Yun Deng on 10/31/23.
//

#include "Branch.hpp"

Branch::Branch() {
    lower_node = nullptr;
    upper_node = nullptr;
}

Branch::Branch(Node *n1, Node *n2) {
    assert(n1->time <= n2->time);
    lower_node = n1;
    upper_node = n2;
}
