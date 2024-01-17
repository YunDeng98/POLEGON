//
//  DAG.cpp
//  arg_branch_length
//
//  Created by Yun Deng on 10/31/23.
//

#include "DAG.hpp"

DAG::DAG(float n) {
    Ne = n;
}

void DAG::load_dag(string node_file, string branch_file, string mut_file) {
    load_nodes(node_file);
    load_branches(branch_file);
    load_mutations(mut_file);
}

void DAG::compute_mutation_rates(float theta) {
    for (Branch *b : branches) {
        b->mutation_rate = b->span*theta*Ne;
    }
}

void DAG::burn_in() {
    fill(node_ages.begin(), node_ages.end(), 0.0);
    num_posterior_samples = 0;
}

/*
void DAG::MCMC(int n, Distribution *d) {
    for (int i = 0; i < n; i++) {
        int index = random_index();
        propose(index, d);
    }
    num_posterior_samples += 1;
    for (int i = 0; i < nodes.size(); i++) {
        node_ages[i] += nodes[i]->time;
    }
}
 */

void DAG::no_prior_MCMC() {
    /*
    int num_floating_nodes = (int) nodes.size() - num_leaf_nodes;
    for (int i = 0; i < n; i++) {
        int index = (i % num_floating_nodes) + num_leaf_nodes;
        no_prior_propose(index);
    }
     */
    vector<int> permutation = get_permutation();
    for (int index : permutation) {
        no_prior_propose(index);
    }
    num_posterior_samples += 1;
    for (int i = 0; i < nodes.size(); i++) {
        node_ages[i] += nodes[i]->time;
    }
}

void DAG::posterior_average() {
    for (int i = 0; i < node_ages.size(); i++) {
        nodes[i]->time = node_ages[i]/num_posterior_samples;
    }
}

void DAG::write_node_ages(string filename) {
    ofstream file;
    file.open(filename);
    for (Node *n : nodes) {
        file << std::setprecision(std::numeric_limits<float>::max_digits10) << n->time*Ne << "\n";
    }
    file.close();
}

// private methods:

float DAG::lower_bound(int i) {
    auto &x = children[i];
    if (x.size() == 0) {
        return nodes[i]->time;
    }
    float lb = 0;
    for (auto &b : x) {
        lb = max(b->lower_node->time, lb);
    }
    return lb;
}

float DAG::upper_bound(int i) {
    auto &x = parents[i];
    if (x.size() == 0) {
        return INT_MAX;
    }
    float ub = INT_MAX;
    for (auto &b : x) {
        ub = min(b->upper_node->time, ub);
    }
    return ub;
}

float DAG::log_acceptance_weight(int i, float t) {
    float w = 0;
    auto &xp = parents[i];
    auto &xc = children[i];
    float length = 0, count = 0, rate = 0;
    for (auto &b : xp) { // note that the mutation rate here is a product of Ne, theta, and span
        length = b->upper_node->time - t;
        count = b->mutation_count;
        rate = length*b->mutation_rate;
        if (rate > 0) {
            w += count*log(rate);
            w -= rate;
        } else {
            w = 0;
        }
        assert(!isnan(w));
    }
    for (auto &b : xc) {
        length = t - b->lower_node->time;
        count = b->mutation_count;
        rate = length*b->mutation_rate;
        if (rate > 0) {
            w += count*log(rate);
            w -= rate;
        } else {
            w = 0;
        }
        assert(!isnan(w));
    }
    return w;
}

float DAG::acceptance_ratio(int i, float t) {
    float t0 = nodes[i]->time;
    float w0 = log_acceptance_weight(i, t0);
    float w1 = log_acceptance_weight(i, t);
    float q = exp(w1 - w0);
    return q;
}

float DAG::no_prior_acceptance_ratio(int i, float t, float lb, float ub) {
    float q = acceptance_ratio(i, t);
    float t0 = nodes[i]->time;
    if (ub == INT_MAX) {
        q *= exp((t - t0)/lambda);
    }
    /*
    if (ub == INT_MAX) {
        cout << lb << " " << t0 << " " << t << " " << q << endl;
    }
     */
    return q;
}

void DAG::propose(int i, Distribution *d) {
    float lb = lower_bound(i);
    float ub = upper_bound(i);
    float t = d->propose(lb, ub);
    float ar = acceptance_ratio(i, t);
    float q = uniform_random();
    if (q < ar) {
        nodes[i]->time = t;
        updates += 1;
        cout << "Number of updates: " << updates << endl;
    }
}

void DAG::no_prior_propose(int i) {
    float lb = lower_bound(i);
    float ub = upper_bound(i);
    float t0 = nodes[i]->time;
    float t = 0;
    if (ub != INT_MAX) {
        t = random_non_root_time(t0, lb, ub);
    } else {
        t = random_root_time(i, lb);
    }
    float ar = no_prior_acceptance_ratio(i, t, lb, ub);
    float q = uniform_random();
    if (q < ar) {
        nodes[i]->time = t;
        updates += 1;
        // cout << "Number of updates: " << updates << endl;
    }
}

int DAG::random_index() {
    int index = 0;
    while (nodes[index]->is_sample) {
        index = floor(nodes.size()*uniform_random());
    }
    return index;
}

void DAG::load_nodes(string node_file) {
    ifstream fin(node_file);
    if (!fin.good()) {
        cerr << "input file not found" << endl;
        exit(1);
    }
    int count = 0;
    float x;
    while (fin >> x) {
        Node *n = new Node(x/Ne, count);
        if (x == 0) {
            n->is_sample = true;
            num_leaf_nodes += 1;
        }
        nodes.push_back(n);
        n->index = count;
        count += 1;
    }
    node_ages.resize(nodes.size());
}

void DAG::load_branches(string branch_file) {
    ifstream fin(branch_file);
    if (!fin.good()) {
        cerr << "input file not found" << endl;
        exit(1);
    }
    float x;
    float y;
    float p;
    float c;
    Node *un;
    Node *ln;
    Branch *b;
    map<pair<Node *, Node *>, float> branch_span = {};
    while (fin >> x >> y >> p >> c) {
        if (p < 0) {
            un = root;
        } else {
            un = nodes[int(p)];
        }
        ln = nodes[int(c)];
        assert(ln->index < un->index or un == root);
        branch_span[{ln, un}] += y - x;
    }
    for (auto &x : branch_span) {
        if (x.first.second != root) {
            b = new Branch(x.first.first, x.first.second);
            b->span = x.second;
            branches.push_back(b);
        }
    }
    sort(branches.begin(), branches.end(), compare_branch());
    parents.resize(nodes.size());
    children.resize(nodes.size());
    for (Branch *b : branches) {
        if (b->upper_node != root) {
            parents[b->lower_node->index].insert(b);
            children[b->upper_node->index].insert(b);
        }
    }
}

void DAG::load_mutations(string mut_file) {
    ifstream fin(mut_file);
    if (!fin.good()) {
        cerr << "input file not found" << endl;
        exit(1);
    }
    float pos;
    float n1;
    float n2;
    float s;
    Node *ln;
    Node *un;
    Branch *b;
    while (fin >> pos >> n1 >> n2 >> s) {
        if (n2 >= 0) {
            ln = nodes[n1];
            un = nodes[(int) n2];
            b = search_branch(ln, un);
            b->mutation_count += 1;
        }
    }
}

/*
Branch *DAG::search_branch(Node *n1, Node *n2) {
    int l = 0;
    int u = (int) branches.size() - 1;
    int m = 0.5*(l + u);
    Branch *b = branches[m];
    while (b->lower_node != n1 or b->upper_node != n2) {
        if (b->upper_node->index < n2->index) {
            l = m;
            m = 0.5*(l + u);
            b = branches[m];
        } else if (b->upper_node->index == n2->index and b->lower_node->index < n1->index) {
            l = m;
            m = 0.5*(l + u);
            b = branches[m];
        } else {
            u = m;
            m = 0.5*(l + u);
            b = branches[m];
        }
    }
    return branches[m];
}
*/

Branch *DAG::search_branch(Node *n1, Node *n2) {
    int l = 0;
    int u = (int) branches.size() - 1;
    while (l <= u) {
        int m = l + (u - l) / 2;
        Branch *b = branches[m];
        if (b->lower_node == n1 && b->upper_node == n2) {
            return b;
        }
        if (b->upper_node->index < n2->index ||
            (b->upper_node->index == n2->index && b->lower_node->index < n1->index)) {
            l = m + 1;
        } else {
            u = m - 1;
        }
    }
    cout << "branch search failed!" << endl;
    return nullptr;
}

float DAG::random_non_root_time(float t0, float lb, float ub) {
    assert(ub != INT_MAX);
    float t = lb + uniform_random()*(ub - lb);
    if (t <= lb or t >= ub) {
        t = 0.5*(lb + ub);
    }
    return t;
}

float DAG::random_root_time(int i, float lb) {
    float q = uniform_random();
    float delta = -lambda*log(q);
    while (delta <= 0.001 or delta > 10) {
        q = uniform_random();
        delta = -lambda*log(q);
    }
    return lb + delta;
}

vector<int> DAG::get_permutation() {
    vector<int> permutation = {};
    permutation.reserve(nodes.size() - num_leaf_nodes);
    for (int i = 0; i < nodes.size() - num_leaf_nodes; i++) {
        permutation.push_back(i + num_leaf_nodes);
    }
    shuffle(permutation.begin(), permutation.end(), random_engine);
    return permutation;
}
