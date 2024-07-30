//
//  DAG.cpp
//  arg_branch_length
//
//  Created by Yun Deng on 10/31/23.
//

#include "DAG.hpp"

DAG::DAG(double n) {
    Ne = n;
}

void DAG::load_dag(string node_file, string branch_file) {
    load_nodes(node_file);
    load_branches(branch_file);
    num_posterior_samples = 1;
}

void DAG::load_dag(string node_file, string branch_file, Mutation_map &mm) {
    load_nodes(node_file);
    load_branches(branch_file, mm);
    num_posterior_samples = 1;
}

void DAG::map_mutations(string mut_file) {
    ifstream fin(mut_file);
    if (!fin.good()) {
        cerr << "input file not found" << endl;
        exit(1);
    }
    double pos;
    double n1;
    double n2;
    double s;
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

void DAG::compute_mutation_rates(double theta) {
    for (Branch *b : branches) {
        b->mutation_rate = b->span*theta*Ne;
    }
}

void DAG::burn_in() {
    fill(node_ages.begin(), node_ages.end(), 0.0);
    num_posterior_samples = 0;
}

void DAG::record_node_ages() {
    for (int i = 0; i < node_ages.size(); i++) {
        node_age_samples[i].push_back(nodes[i]->time);
    }
}

void DAG::record_scaled_node_ages() {
    for (int i = 0; i < node_ages.size(); i++) {
        scaled_node_age_samples[i].push_back(nodes[i]->time);
    }
}

void DAG::no_prior_MCMC() {
    vector<int> permutation = get_permutation();
    for (int index : permutation) {
        no_prior_propose(index);
    }
    num_posterior_samples += 1;
    for (int i = 0; i < nodes.size(); i++) {
        node_ages[i] += nodes[i]->time;
    }
}

void DAG::sample(int i) {
    for (int j = 0; j < node_ages.size(); j++) {
        nodes[j]->time = node_age_samples[j][i];
    }
}

void DAG::posterior_average() {
    for (int i = 0; i < node_ages.size(); i++) {
        nodes[i]->time = node_ages[i]/num_posterior_samples;
    }
}

void DAG::scaled_sample_average() {
    num_posterior_samples = (int) scaled_node_age_samples[0].size();
    double sum = 0;
    for (int i = 0; i < node_ages.size(); i++) {
        sum = accumulate(scaled_node_age_samples[i].begin(), scaled_node_age_samples[i].end(), 0.0);
        nodes[i]->time = sum/num_posterior_samples;
    }
}

void DAG::write_node_ages(string filename) {
    ofstream file;
    file.open(filename);
    for (Node *n : nodes) {
        file << std::setprecision(std::numeric_limits<double>::max_digits10) << n->time*Ne << "\n";
    }
    file.close();
}

void DAG::write_node_age_samples(string filename) {
    ofstream file;
    file.open(filename);
    num_posterior_samples = (int) node_age_samples[0].size();
    for (int i = 0; i < scaled_node_age_samples.size(); i++) {
        for (int j = 0; j < num_posterior_samples; j++) {
            file << std::setprecision(std::numeric_limits<double>::max_digits10) << scaled_node_age_samples[i][j]*Ne << " ";
        }
        file << "" << endl;
    }
    file.close();
}

// private methods:

double DAG::lower_bound(int i) {
    auto &x = children[i];
    if (x.size() == 0) {
        return nodes[i]->time;
    }
    double lb = 0;
    for (auto &b : x) {
        lb = max(b->lower_node->time, lb);
    }
    return lb;
}

double DAG::upper_bound(int i) {
    auto &x = parents[i];
    if (x.size() == 0) {
        return INT_MAX;
    }
    double ub = INT_MAX;
    for (auto &b : x) {
        ub = min(b->upper_node->time, ub);
    }
    return ub;
}

double DAG::log_acceptance_weight(int i, double t) {
    double w = 0;
    auto &xp = parents[i];
    auto &xc = children[i];
    double length = 0, count = 0, rate = 0;
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

double DAG::fast_acceptance_ratio(int i, double t0, double t1) {
    double w0 = 0;
    double w1 = 0;
    auto &xp = parents[i];
    auto &xc = children[i];
    double length_0 = 0, length_1 = 0, count = 0, rate_0 = 0, rate_1 = 0;
    for (auto &b : xp) { // note that the mutation rate here is a product of Ne, theta, and span
        length_0 = b->upper_node->time - t0;
        length_1 = b->upper_node->time - t1;
        count = b->mutation_count;
        rate_0 = length_0*b->mutation_rate;
        rate_1 = length_1*b->mutation_rate;
        if (rate_0 > 0) {
            w0 += count*log(rate_0);
            w0 -= rate_0;
        } else {
            w0 = 0;
        }
        if (rate_1 > 0) {
            w1 += count*log(rate_1);
            w1 -= rate_1;
        } else {
            w1 = 0;
        }
        assert(!isnan(w0));
        assert(!isnan(w1));
    }
    for (auto &b : xc) {
        length_0 = t0 - b->lower_node->time;
        length_1 = t1 - b->lower_node->time;
        count = b->mutation_count;
        rate_0 = length_0*b->mutation_rate;
        rate_1 = length_1*b->mutation_rate;
        if (rate_0 > 0) {
            w0 += count*log(rate_0);
            w0 -= rate_0;
        } else {
            w0 = 0;
        }
        if (rate_1 > 0) {
            w1 += count*log(rate_1);
            w1 -= rate_1;
        } else {
            w1 = 0;
        }
        assert(!isnan(w0));
        assert(!isnan(w1));
    }
    return exp(w1 - w0);
}

double DAG::acceptance_ratio(int i, double t) {
    double t0 = nodes[i]->time;
    double w0 = log_acceptance_weight(i, t0);
    double w1 = log_acceptance_weight(i, t);
    double q = exp(w1 - w0);
    // double r = fast_acceptance_ratio(i, t0, t);
    // assert(r == q);
    return q;
}

double DAG::no_prior_acceptance_ratio(int i, double t, double lb, double ub) {
    // double q = acceptance_ratio(i, t);
    double t0 = nodes[i]->time;
    double q = fast_acceptance_ratio(i, t0, t);
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
    double lb = lower_bound(i);
    double ub = upper_bound(i);
    double t = d->propose(lb, ub);
    double ar = acceptance_ratio(i, t);
    double q = uniform_random();
    if (q < ar) {
        nodes[i]->time = t;
        updates += 1;
        cout << "Number of updates: " << updates << endl;
    }
}

void DAG::no_prior_propose(int i) {
    double lb = lower_bound(i);
    double ub = upper_bound(i);
    double t0 = nodes[i]->time;
    double t = 0;
    if (ub != INT_MAX) {
        t = random_non_root_time(t0, lb, ub);
    } else {
        t = random_root_time(i, lb);
    }
    double ar = no_prior_acceptance_ratio(i, t, lb, ub);
    double q = uniform_random();
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
    double x;
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
    node_age_samples.resize(nodes.size());
    scaled_node_age_samples.resize(nodes.size());
    scaling_factors.resize(nodes.size());
    for (int i = 0; i < nodes.size(); i++) {
        node_ages[i] = nodes[i]->time;
    }
}

void DAG::load_branches(string branch_file) {
    ifstream fin(branch_file);
    if (!fin.good()) {
        cerr << "input file not found" << endl;
        exit(1);
    }
    double x;
    double y;
    double p;
    double c;
    Node *un;
    Node *ln;
    Branch *b;
    map<pair<Node *, Node *>, double> branch_span = {};
    while (fin >> x >> y >> p >> c) {
        if (p < 0) {
            un = root;
        } else {
            un = nodes[int(p)];
        }
        ln = nodes[int(c)];
        // assert(ln->index < un->index or un == root);
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

void DAG::load_branches(string branch_file, Mutation_map &mm) {
    ifstream fin(branch_file);
    if (!fin.good()) {
        cerr << "input file not found" << endl;
        exit(1);
    }
    double x;
    double y;
    double p;
    double c;
    double m;
    Node *un;
    Node *ln;
    Branch *b;
    map<pair<Node *, Node *>, double> branch_span = {};
    map<pair<Node *, Node *>, double> branch_rates = {};
    while (fin >> x >> y >> p >> c) {
        if (p < 0) {
            un = root;
        } else {
            un = nodes[int(p)];
        }
        ln = nodes[int(c)];
        assert(ln->index < un->index or un == root);
        branch_span[{ln, un}] += y - x;
        m = mm.mutation_rate(x, y)*Ne;
        branch_rates[{ln, un}] += m;
    }
    for (auto &x : branch_span) {
        if (x.first.second != root) {
            b = new Branch(x.first.first, x.first.second);
            b->span = x.second;
            b->mutation_rate = branch_rates[{b->lower_node, b->upper_node}];
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

double DAG::random_non_root_time(double t0, double lb, double ub) {
    assert(ub != INT_MAX);
    double t = lb + uniform_random()*(ub - lb);
    if (t <= lb or t >= ub) {
        // t = 0.5*(lb + ub);
        t = lb + uniform_random()*(ub - lb);
    }
    return t;
}

double DAG::random_root_time(int i, double lb) {
    double q = uniform_random();
    double delta = -lambda*log(q);
    while (delta > max_step) { // max size of the exploration
        q = uniform_random();
        delta = -lambda*log(q);
    }
    return lb + delta;
}

vector<int> DAG::get_permutation() {
    vector<int> permutation = {};
    permutation.reserve(nodes.size() - num_leaf_nodes);
    for (int i = 0; i < nodes.size() - num_leaf_nodes; i++) {
        permutation.push_back((int) nodes.size() - 1 - i);
    }
    return permutation;
}

double DAG::median(std::vector<double>& vec) {
    int size = (int) vec.size();
    assert(size > 0);
    vector<double> temp(vec);
    nth_element(temp.begin(), temp.begin() + size/2, temp.end());
    return vec[size/2];
}
