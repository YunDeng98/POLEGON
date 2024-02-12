//
//  Branch.hpp
//  arg_branch_length
//
//  Created by Yun Deng on 10/31/23.
//

#ifndef Branch_hpp
#define Branch_hpp

#include <stdio.h>
#include "Node.hpp"

class Branch {
    
public:
    
    Node *upper_node = nullptr;
    Node *lower_node = nullptr;
    double span = 0;
    double mutation_count = 0;
    double mutation_rate = 0;
    
    Branch();
    
    Branch(Node *n1, Node *n2);
    
};

struct compare_branch {
    
    bool operator() (const Branch *b1, const Branch *b2) const {
        if (b1->upper_node->index != b2->upper_node->index) {
            return b1->upper_node->index < b2->upper_node->index;
        } else if (b1->lower_node->index != b2->lower_node->index) {
            return b1->lower_node->index < b2->lower_node->index;
        } else {
            assert(b1->upper_node == b2->upper_node and b1->lower_node == b2->lower_node);
            return false;
        }
    }
};

#endif /* Branch_hpp */
