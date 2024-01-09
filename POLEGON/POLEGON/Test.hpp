//
//  Test.hpp
//  arg_branch_length
//
//  Created by Yun Deng on 10/31/23.
//

#ifndef Test_hpp
#define Test_hpp

#include <stdio.h>
#include "DAG.hpp"
#include "Distribution.hpp"
#include "Scaler.hpp"

void test_load_dag();

void test_coalescent_prior();

void test_sampling();

void test_tsinfer_topology();

void test_singer_topology();

void test_singer_demo_topology();

void test_no_prior_sampling();

void test_scaling();

void test_demography();

void test_demo_scaling();

void test_bottleneck();

void test_bgs();

void test_pairwise_demo();

void test_migration();

#endif /* Test_hpp */
