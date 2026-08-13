//
//  DAG.cpp
//  POLEGON
//
//  Created by Yun Deng on 10/31/23.
//  Updated by Wonseop Lim on 05/16/26.
//

#include <cassert>
#include <sstream>
#include <queue>
#include <set>
#include <unordered_set>
#include <omp.h>
#include "DAG.hpp"

DAG::DAG(double n) {
    Ne = n;
}

void DAG::load_dag(string node_file, string branch_file) {
    load_nodes(node_file);
    load_branches(branch_file);
}

void DAG::load_dag(string node_file, string branch_file, Mutation_map &mm) {
    load_nodes(node_file);
    load_branches(branch_file, mm);
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

// Adjusts ARG for heterochronous samples
//   Convert times in coalescent units and shift all node times so the youngest tip is at t=0
//   store the shift in time_origin and correct any internal node whose time falls at or below
//   its oldest child to lb + 1e-6
void DAG::apply_tip_ages(string tip_ages_file, double gen_time) {
    ifstream fin(tip_ages_file);
    if (!fin.good()) {
        cerr << "tip_ages file not found: " << tip_ages_file << endl;
        exit(1);
    }
    vector<double> tip_ages_coal;
    double age_years;
    double min_age = numeric_limits<double>::infinity();
    for (int i = 0; i < (int)nodes.size(); i++) {
        if (nodes[i]->is_sample) {
            if (!(fin >> age_years)) {
                cerr << "tip_ages file has fewer entries than the number of tips in the ARG" << endl;
                exit(1);
            }
            double age_coal = age_years / (gen_time * Ne);
            tip_ages_coal.push_back(age_coal);
            min_age = min(min_age, age_coal);
        }
    }
    time_origin = min_age;
    for (int i = 0; i < (int)nodes.size(); i++) {
        nodes[i]->time -= min_age;
        node_times[i] = nodes[i]->time;
    }
    int s = 0;
    for (int i = 0; i < (int)nodes.size(); i++) {
        if (nodes[i]->is_sample) {
            nodes[i]->time = tip_ages_coal[s++] - min_age;
            node_times[i] = nodes[i]->time;
        }
    }
    for (int i = 0; i < (int)nodes.size(); i++) {
        if (!nodes[i]->is_sample) {
            double lb = lower_bound(i);
            if (nodes[i]->time <= lb) {
                nodes[i]->time = lb + 1e-6;
                node_times[i] = nodes[i]->time;
            }
        }
    }
}

// Chromatic decomposition via two steps:
//   1. Smallest-last ordering (Matula & Beck 1983): repeatedly remove the min-degree
//      node; reversed order guarantees greedy uses ≤ (degeneracy+1) colors,
//      where degeneracy = max min-degree over all induced sub-ARGs.
//   2. Greedy coloring: assign each node the smallest color unused by its neighbors.
void DAG::compute_coloring() {
    int n = (int)nodes.size();

    // Build deduplicated adjacency (multiple branches between same pair count once)
    vector<unordered_set<int>> adj(n);
    for (int i : perm_cache) {
        for (int k = parent_start[i]; k < parent_start[i+1]; k++) {
            int j = parent_data[k]->upper_node->index;
            if (j < n && !nodes[j]->is_sample)
                adj[i].insert(j);
        }
        for (int k = child_start[i]; k < child_start[i+1]; k++) {
            int j = child_data[k]->lower_node->index;
            if (!nodes[j]->is_sample)
                adj[i].insert(j);
        }
    }

    vector<int> deg(n, 0);
    for (int i : perm_cache) deg[i] = (int)adj[i].size();

    vector<bool> removed(n, false);
    set<pair<int,int>> pq;
    for (int i : perm_cache) pq.insert({deg[i], i});

    vector<int> order;
    order.reserve(perm_cache.size());
    while (!pq.empty()) {
        int i = pq.begin()->second;
        pq.erase(pq.begin());
        removed[i] = true;
        order.push_back(i);
        for (int j : adj[i]) {
            if (!removed[j]) {
                pq.erase({deg[j], j});
                deg[j]--;
                pq.insert({deg[j], j});
            }
        }
    }
    reverse(order.begin(), order.end());
    
    vector<int> node_color(n, -1);
    for (int i : order) {
        unordered_set<int> forbidden;
        for (int j : adj[i])
            if (node_color[j] >= 0) forbidden.insert(node_color[j]);
        int color = 0;
        while (forbidden.count(color)) color++;
        node_color[i] = color;
    }

    int num_colors = 0;
    for (int i : perm_cache) num_colors = max(num_colors, node_color[i] + 1);
    color_classes.assign(num_colors, {});
    for (int i : perm_cache) color_classes[node_color[i]].push_back(i);

    cout << "Chromatic decomposition: " << num_colors << " color classes" << endl;
    for (int c = 0; c < num_colors; c++)
        cout << "  Class " << c + 1 << ": " << color_classes[c].size() << " nodes" << endl;

    int pe = parent_start.back();
    parent_upper_idx.resize(pe);
    parent_mut_count.resize(pe);
    parent_mut_rate.resize(pe);
    parent_span.resize(pe);
    for (int i = 0; i < n; i++) {
        for (int k = parent_start[i]; k < parent_start[i+1]; k++) {
            Branch *b = parent_data[k];
            parent_upper_idx[k] = b->upper_node->index;
            parent_mut_count[k] = b->mutation_count;
            parent_mut_rate[k]  = b->mutation_rate;
            parent_span[k]      = b->span;
        }
    }
    int ce = child_start.back();
    child_lower_idx.resize(ce);
    child_mut_count.resize(ce);
    child_mut_rate.resize(ce);
    child_span.resize(ce);
    for (int i = 0; i < n; i++) {
        for (int k = child_start[i]; k < child_start[i+1]; k++) {
            Branch *b = child_data[k];
            child_lower_idx[k] = b->lower_node->index;
            child_mut_count[k] = b->mutation_count;
            child_mut_rate[k]  = b->mutation_rate;
            child_span[k]      = b->span;
        }
    }
}

void DAG::no_prior_MCMC() {
    #pragma omp parallel num_threads(num_cores)
    {
        for (const auto& class_nodes : color_classes) {
            int n = (int)class_nodes.size();
            int n_chunks = min(num_streams, n);
            #pragma omp for schedule(static)
            for (int s = 0; s < n_chunks; s++) {
                bind_random_stream(s);
                int lo = (int)((long)n*s/n_chunks);
                int hi = (int)((long)n*(s + 1)/n_chunks);
                for (int k = lo; k < hi; k++) {
                    no_prior_propose(class_nodes[k]);
                }
            }
        }
    }
}

void DAG::sync_node_times() {
    for (int i = 0; i < (int)nodes.size(); i++)
        nodes[i]->time = node_times[i];
}

void DAG::posterior_average(string samples_file, string output_file) {
    ifstream fin(samples_file);
    if (!fin.good()) {
        cerr << "Error: samples file not found: " << samples_file << endl;
        exit(1);
    }
    int n = nodes.size();
    vector<double> sums(n, 0.0);
    int count = 0;
    string line;
    while (getline(fin, line)) {
        if (line.empty()) continue;
        istringstream iss(line);
        for (int i = 0; i < n; i++) {
            double val;
            iss >> val;
            sums[i] += val;
        }
        count++;
    }
    fin.close();
    if (count == 0) {
        cerr << "Error: samples file is empty." << endl;
        exit(1);
    }
    ofstream fout(output_file);
    for (int i = 0; i < n; i++) {
        fout << setprecision(numeric_limits<double>::max_digits10) << sums[i] / count << "\n";
    }
    fout.close();
}

void DAG::write_node_ages(string filename, double gen_time) {
    ofstream file;
    file.open(filename);
    for (Node *n : nodes) {
        file << std::setprecision(std::numeric_limits<double>::max_digits10) << (n->time + time_origin)*Ne*gen_time << "\n";
    }
    file.close();
}

double DAG::lower_bound(int i) {
    if (child_start[i] == child_start[i+1]) return nodes[i]->time;
    double lb = 0;
    if (!child_lower_idx.empty()) {
        for (int k = child_start[i]; k < child_start[i+1]; k++)
            lb = max(nodes[child_lower_idx[k]]->time, lb);
    } else {
        for (int k = child_start[i]; k < child_start[i+1]; k++)
            lb = max(child_data[k]->lower_node->time, lb);
    }
    return lb;
}

double DAG::lower_bound(int i, const vector<double>& times) const {
    if (child_start[i] == child_start[i+1]) return times[i];
    double lb = 0;
    for (int k = child_start[i]; k < child_start[i+1]; k++)
        lb = max(times[child_lower_idx[k]], lb);
    return lb;
}

double DAG::upper_bound(int i) {
    if (parent_start[i] == parent_start[i+1]) return INT_MAX;
    double ub = INT_MAX;
    if (!parent_upper_idx.empty()) {
        for (int k = parent_start[i]; k < parent_start[i+1]; k++)
            ub = min(nodes[parent_upper_idx[k]]->time, ub);
    } else {
        for (int k = parent_start[i]; k < parent_start[i+1]; k++)
            ub = min(parent_data[k]->upper_node->time, ub);
    }
    return ub;
}

// Poisson log-likelihood contribution of node i
double DAG::log_acceptance_weight(int i, double t) {
    double w = 0;
    double length = 0, count = 0, rate = 0;
    for (int k = parent_start[i]; k < parent_start[i+1]; k++) {
        Branch *b = parent_data[k];
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
    for (int k = child_start[i]; k < child_start[i+1]; k++) {
        Branch *b = child_data[k];
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

// Computes exp(w(t1) - w(t0)) over incident branches
double DAG::fast_acceptance_ratio(int i, double t0, double t1) {
    double w0 = 0, w1 = 0;
    for (int k = parent_start[i]; k < parent_start[i+1]; k++) {
        double upper_t  = nodes[parent_upper_idx[k]]->time;
        double count    = parent_mut_count[k];
        double mut_rate = parent_mut_rate[k];
        double rate_0   = (upper_t - t0) * mut_rate;
        double rate_1   = (upper_t - t1) * mut_rate;
        if (rate_0 > 0) { w0 += count*log(rate_0); w0 -= rate_0; } else { w0 = 0; }
        if (rate_1 > 0) { w1 += count*log(rate_1); w1 -= rate_1; } else { w1 = 0; }
    }
    for (int k = child_start[i]; k < child_start[i+1]; k++) {
        double lower_t  = nodes[child_lower_idx[k]]->time;
        double count    = child_mut_count[k];
        double mut_rate = child_mut_rate[k];
        double rate_0   = (t0 - lower_t) * mut_rate;
        double rate_1   = (t1 - lower_t) * mut_rate;
        if (rate_0 > 0) { w0 += count*log(rate_0); w0 -= rate_0; } else { w0 = 0; }
        if (rate_1 > 0) { w1 += count*log(rate_1); w1 -= rate_1; } else { w1 = 0; }
    }
    return exp(w1 - w0);
}

double DAG::acceptance_ratio(int i, double t) {
    double t0 = nodes[i]->time;
    double w0 = log_acceptance_weight(i, t0);
    double w1 = log_acceptance_weight(i, t);
    double q = exp(w1 - w0);
    return q;
}

// For root nodes multiply by exp((t - t0)/lambda) (Hastings correction)
double DAG::no_prior_acceptance_ratio(int i, double t, double lb, double ub) {
    double t0 = nodes[i]->time;
    double q = fast_acceptance_ratio(i, t0, t);
    if (ub == INT_MAX) {
        q *= exp((t - t0)/lambda);
    }
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
    double x, sf;
    while (fin >> x >> sf) {
        Node *n = new Node(x/Ne, count);
        if (sf > 0.5) {
            n->is_sample = true;
            num_leaf_nodes += 1;
        }
        nodes.push_back(n);
        n->index = count;
        count += 1;
    }
    scaling_factors.resize(nodes.size());
    node_times.resize(nodes.size());
    for (int i = 0; i < (int)nodes.size(); i++) node_times[i] = nodes[i]->time;
    int n_internal = (int)nodes.size() - num_leaf_nodes;
    perm_cache.resize(n_internal);
    for (int i = 0; i < n_internal; i++) {
        perm_cache[i] = (int)nodes.size() - 1 - i;
    }
}

// build CSR adjacency arrays in two passes
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
    int n = (int)nodes.size();
    parent_start.assign(n + 1, 0);
    child_start.assign(n + 1, 0);
    for (Branch *b : branches) {
        if (b->upper_node != root) {
            parent_start[b->lower_node->index + 1]++;
            child_start[b->upper_node->index + 1]++;
        }
    }
    for (int i = 1; i <= n; i++) {
        parent_start[i] += parent_start[i-1];
        child_start[i] += child_start[i-1];
    }
    parent_data.resize(parent_start.back());
    child_data.resize(child_start.back());
    vector<int> ppos = parent_start, cpos = child_start;
    for (Branch *b : branches) {
        if (b->upper_node != root) {
            parent_data[ppos[b->lower_node->index]++] = b;
            child_data[cpos[b->upper_node->index]++] = b;
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
    int n = (int)nodes.size();
    parent_start.assign(n + 1, 0);
    child_start.assign(n + 1, 0);
    for (Branch *b : branches) {
        if (b->upper_node != root) {
            parent_start[b->lower_node->index + 1]++;
            child_start[b->upper_node->index + 1]++;
        }
    }
    for (int i = 1; i <= n; i++) {
        parent_start[i] += parent_start[i-1];
        child_start[i] += child_start[i-1];
    }
    parent_data.resize(parent_start.back());
    child_data.resize(child_start.back());
    vector<int> ppos = parent_start, cpos = child_start;
    for (Branch *b : branches) {
        if (b->upper_node != root) {
            parent_data[ppos[b->lower_node->index]++] = b;
            child_data[cpos[b->upper_node->index]++] = b;
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
        t = lb + uniform_random()*(ub - lb);
    }
    return t;
}

// Draws Exp(1/lambda) above lb; rejects and redraws if result exceeds lb + max_step
double DAG::random_root_time(int i, double lb) {
    double q = uniform_random();
    double delta = -lambda*log(q);
    while (delta > max_step) {
        q = uniform_random();
        delta = -lambda*log(q);
    }
    return lb + delta;
}

double DAG::median(std::vector<double>& vec) {
    int size = (int) vec.size();
    assert(size > 0);
    vector<double> temp(vec);
    nth_element(temp.begin(), temp.begin() + size/2, temp.end());
    return vec[size/2];
}
