#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdint>
#include <cstdlib>
#include <iostream>
#include <limits>
#include <numeric>
#include <random>
#include <string>
#include <utility>
#include <vector>

using namespace std;

using U64 = unsigned long long;
using I128 = __int128_t;

static constexpr U64 BAD_ID = std::numeric_limits<U64>::max();

static inline U64 comb3(U64 x) {
    return (x >= 3 ? x * (x - 1) * (x - 2) / 6 : 0);
}

struct Target {
    U64 K3 = 0;
    U64 K13 = 0;
    U64 P4 = 0;
    U64 P3K2 = 0;
    U64 K222 = 0;
};

struct Profile {
    U64 K3 = 0;
    U64 K13 = 0;
    U64 P4 = 0;
    U64 P3K2 = 0;
    U64 K222 = 0;
};

struct Graph {
    int n = 0;
    int W = 0;
    U64 m = 0;
    U64 k = 0;

    vector<int> U;
    vector<int> V;

    vector<vector<uint64_t>> adj;
    vector<int> deg;

    vector<unsigned char> inE;
    vector<U64> posE;
    vector<U64> posNE;

    vector<U64> edges;
    vector<U64> nonedges;
};

static inline bool has_edge(const Graph& G, int u, int v) {
    return ((G.adj[u][v >> 6] >> (v & 63)) & 1ULL) != 0;
}

static inline void set_edge_bit(Graph& G, int u, int v) {
    G.adj[u][v >> 6] |= (1ULL << (v & 63));
    G.adj[v][u >> 6] |= (1ULL << (u & 63));
}

static inline void clear_edge_bit(Graph& G, int u, int v) {
    G.adj[u][v >> 6] &= ~(1ULL << (v & 63));
    G.adj[v][u >> 6] &= ~(1ULL << (u & 63));
}

static inline int common_neighbors(const Graph& G, int u, int v) {
    int s = 0;
    for (int i = 0; i < G.W; ++i) {
        s += __builtin_popcountll(G.adj[u][i] & G.adj[v][i]);
    }
    return s;
}

static void build_edge_catalog(Graph& G) {
    G.m = (U64)G.n * (G.n - 1) / 2;
    G.U.resize((size_t)G.m);
    G.V.resize((size_t)G.m);

    U64 id = 0;
    for (int u = 0; u < G.n; ++u) {
        for (int v = u + 1; v < G.n; ++v) {
            G.U[(size_t)id] = u;
            G.V[(size_t)id] = v;
            ++id;
        }
    }
}

static void init_empty_graph(Graph& G, int n, U64 k) {
    G = Graph();
    G.n = n;
    G.k = k;
    G.W = (n + 63) >> 6;

    build_edge_catalog(G);

    G.adj.assign(n, vector<uint64_t>(G.W, 0));
    G.deg.assign(n, 0);

    G.inE.assign((size_t)G.m, 0);
    G.posE.assign((size_t)G.m, BAD_ID);
    G.posNE.assign((size_t)G.m, BAD_ID);

    G.edges.clear();
    G.nonedges.clear();
    G.edges.reserve((size_t)k);
    G.nonedges.reserve((size_t)(G.m - k));
}

static inline void add_edge_by_id(Graph& G, U64 eid) {
    if (G.inE[(size_t)eid]) return;
    int u = G.U[(size_t)eid];
    int v = G.V[(size_t)eid];
    set_edge_bit(G, u, v);
    ++G.deg[u];
    ++G.deg[v];
    G.inE[(size_t)eid] = 1;
}

static inline void remove_edge_by_id(Graph& G, U64 eid) {
    if (!G.inE[(size_t)eid]) return;
    int u = G.U[(size_t)eid];
    int v = G.V[(size_t)eid];
    clear_edge_bit(G, u, v);
    --G.deg[u];
    --G.deg[v];
    G.inE[(size_t)eid] = 0;
}

static void rebuild_lists(Graph& G) {
    G.edges.clear();
    G.nonedges.clear();

    for (U64 eid = 0; eid < G.m; ++eid) {
        if (G.inE[(size_t)eid]) {
            G.posE[(size_t)eid] = (U64)G.edges.size();
            G.posNE[(size_t)eid] = BAD_ID;
            G.edges.push_back(eid);
        } else {
            G.posNE[(size_t)eid] = (U64)G.nonedges.size();
            G.posE[(size_t)eid] = BAD_ID;
            G.nonedges.push_back(eid);
        }
    }
}

static inline void move_edge_to_nonedge(Graph& G, U64 eid) {
    U64 pe = G.posE[(size_t)eid];
    U64 lastE = G.edges.back();
    G.edges[(size_t)pe] = lastE;
    G.posE[(size_t)lastE] = pe;
    G.edges.pop_back();
    G.posE[(size_t)eid] = BAD_ID;

    G.posNE[(size_t)eid] = (U64)G.nonedges.size();
    G.nonedges.push_back(eid);
}

static inline void move_nonedge_to_edge(Graph& G, U64 eid) {
    U64 pn = G.posNE[(size_t)eid];
    U64 lastN = G.nonedges.back();
    G.nonedges[(size_t)pn] = lastN;
    G.posNE[(size_t)lastN] = pn;
    G.nonedges.pop_back();
    G.posNE[(size_t)eid] = BAD_ID;

    G.posE[(size_t)eid] = (U64)G.edges.size();
    G.edges.push_back(eid);
}

static inline void apply_swap(Graph& G, U64 out_eid, U64 in_eid) {
    remove_edge_by_id(G, out_eid);
    add_edge_by_id(G, in_eid);
    move_edge_to_nonedge(G, out_eid);
    move_nonedge_to_edge(G, in_eid);
}

static inline void undo_swap(Graph& G, U64 out_eid, U64 in_eid) {
    remove_edge_by_id(G, in_eid);
    add_edge_by_id(G, out_eid);
    move_edge_to_nonedge(G, in_eid);
    move_nonedge_to_edge(G, out_eid);
}

static bool graph_state_consistent(const Graph& G) {
    U64 edge_count = 0;
    vector<int> deg_check(G.n, 0);

    for (U64 eid = 0; eid < G.m; ++eid) {
        int u = G.U[(size_t)eid];
        int v = G.V[(size_t)eid];
        bool a = has_edge(G, u, v);
        bool b = (bool)G.inE[(size_t)eid];
        if (a != b) return false;
        if (a) {
            ++edge_count;
            ++deg_check[u];
            ++deg_check[v];
        }
    }

    if (edge_count != G.k) return false;
    if (G.edges.size() != (size_t)G.k) return false;
    if (G.nonedges.size() != (size_t)(G.m - G.k)) return false;

    for (int v = 0; v < G.n; ++v) {
        if (deg_check[v] != G.deg[v]) return false;
    }

    for (size_t i = 0; i < G.edges.size(); ++i) {
        U64 eid = G.edges[i];
        if (!G.inE[(size_t)eid]) return false;
        if (G.posE[(size_t)eid] != (U64)i) return false;
        if (G.posNE[(size_t)eid] != BAD_ID) return false;
    }

    for (size_t i = 0; i < G.nonedges.size(); ++i) {
        U64 eid = G.nonedges[i];
        if (G.inE[(size_t)eid]) return false;
        if (G.posNE[(size_t)eid] != (U64)i) return false;
        if (G.posE[(size_t)eid] != BAD_ID) return false;
    }

    return true;
}

static bool basic_parameter_checks(int n, U64 k) {
    if (n <= 6) return false;
    if (k < 4) return false;

    I128 lhs = (I128)4 * k;
    I128 rhs = (I128)n * (n - 1);
    if (lhs > rhs) return false;

    I128 c1 = (I128)2 * k * (k - 1);
    I128 d1 = n + 1;
    if (c1 % d1 != 0) return false;

    I128 c2 = (I128)4 * k * (k - 1) * (k - 2);
    I128 d2 = (I128)3 * (n + 1) * ((I128)n * n - n - 4);
    if (d2 == 0 || c2 % d2 != 0) return false;

    return true;
}

static bool compute_target_profile(int n, U64 k, Target& T) {
    I128 num = (I128)4 * (I128)k * (k - 1) * (k - 2);
    I128 den = (I128)3 * (n + 1) * ((I128)n * n - n - 4);
    if (den == 0 || num % den != 0) return false;

    I128 c = num / den;
    I128 s = c * (n - 3);
    I128 p = (I128)3 * c * (n - 3);

    I128 q_num = (I128)3 * c * (n - 3) * (n - 4);
    I128 r_num = c * (n - 3) * (n - 4) * (n - 5);

    if (q_num % 2 != 0) return false;
    if (r_num % 8 != 0) return false;

    I128 q = q_num / 2;
    I128 r = r_num / 8;

    if (c < 0 || s < 0 || p < 0 || q < 0 || r < 0) return false;

    I128 sum = c + s + p + q + r;
    if (sum != (I128)comb3(k)) return false;

    T.K3 = (U64)c;
    T.K13 = (U64)s;
    T.P4 = (U64)p;
    T.P3K2 = (U64)q;
    T.K222 = (U64)r;
    return true;
}

static Profile full_profile(const Graph& G) {
    Profile P{0, 0, 0, 0, 0};

    for (int v = 0; v < G.n; ++v) {
        U64 d = (U64)G.deg[v];
        if (d >= 3) P.K13 += d * (d - 1) * (d - 2) / 6;
    }

    for (U64 eid : G.edges) {
        int u = G.U[(size_t)eid];
        int v = G.V[(size_t)eid];
        int c = common_neighbors(G, u, v);
        P.K3 += (U64)c;
        P.P4 += (U64)((long long)(G.deg[u] - 1) * (long long)(G.deg[v] - 1) - c);
    }
    P.K3 /= 3;

    for (int v = 0; v < G.n; ++v) {
        vector<int> nb;
        nb.reserve((size_t)G.deg[v]);
        for (int u = 0; u < G.n; ++u) {
            if (u != v && has_edge(G, v, u)) nb.push_back(u);
        }

        int dv = G.deg[v];
        int m = (int)nb.size();
        for (int i = 0; i < m; ++i) {
            int a = nb[i];
            for (int j = i + 1; j < m; ++j) {
                int b = nb[j];
                int ab = has_edge(G, a, b) ? 1 : 0;
                P.P3K2 += (U64)((long long)G.k - dv - G.deg[a] - G.deg[b] + 2 + ab);
            }
        }
    }

    P.K222 = comb3(G.k) - P.K3 - P.K13 - P.P4 - P.P3K2;
    return P;
}

static inline bool exact_match(const Profile& P, const Target& T) {
    return P.K3 == T.K3 &&
           P.K13 == T.K13 &&
           P.P4 == T.P4 &&
           P.P3K2 == T.P3K2 &&
           P.K222 == T.K222;
}

static inline U64 l1_score4(const Profile& P, const Target& T) {
    U64 s = 0;
    s += (P.K3   >= T.K3   ? P.K3   - T.K3   : T.K3   - P.K3);
    s += (P.K13  >= T.K13  ? P.K13  - T.K13  : T.K13  - P.K13);
    s += (P.P4   >= T.P4   ? P.P4   - T.P4   : T.P4   - P.P4);
    s += (P.P3K2 >= T.P3K2 ? P.P3K2 - T.P3K2 : T.P3K2 - P.P3K2);
    return s;
}

static Graph near_regular_initial_graph(int n, U64 k, mt19937_64& rng) {
    Graph G;
    init_empty_graph(G, n, k);

    vector<U64> all_ids((size_t)G.m);
    iota(all_ids.begin(), all_ids.end(), 0ULL);
    shuffle(all_ids.begin(), all_ids.end(), rng);

    U64 sum_deg = 2 * k;
    int base = (int)(sum_deg / n);
    int extra = (int)(sum_deg % n);

    vector<int> target_deg(n, base);
    vector<int> perm(n);
    iota(perm.begin(), perm.end(), 0);
    shuffle(perm.begin(), perm.end(), rng);
    for (int i = 0; i < extra; ++i) target_deg[perm[i]]++;

    vector<int> rem = target_deg;

    U64 chosen = 0;
    for (U64 eid : all_ids) {
        int u = G.U[(size_t)eid];
        int v = G.V[(size_t)eid];
        if (rem[u] > 0 && rem[v] > 0) {
            add_edge_by_id(G, eid);
            --rem[u];
            --rem[v];
            ++chosen;
            if (chosen == k) break;
        }
    }

    if (chosen < k) {
        for (U64 eid : all_ids) {
            if (!G.inE[(size_t)eid]) {
                add_edge_by_id(G, eid);
                ++chosen;
                if (chosen == k) break;
            }
        }
    }

    rebuild_lists(G);
    return G;
}

static void print_edges(const Graph& G) {
    vector<pair<int,int>> out;
    out.reserve((size_t)G.k);
    for (U64 eid : G.edges) {
        out.push_back({G.U[(size_t)eid] + 1, G.V[(size_t)eid] + 1});
    }
    sort(out.begin(), out.end());
    for (const auto& e : out) {
        cout << e.first << " " << e.second << "\n";
    }
}

static void verified_print_and_exit(const Graph& G, const Target& T) {
    if (!graph_state_consistent(G)) {
        std::abort();
    }
    Profile R = full_profile(G);
    if (!exact_match(R, T)) {
        std::abort();
    }
    print_edges(G);
    std::exit(0);
}

int main(int argc, char** argv) {
    if (argc < 3) {
        cerr << "Usage: " << argv[0]
             << " n k [restarts] [steps] [seed] [log_every]\n";
        return 1;
    }

    int n = atoi(argv[1]);
    U64 k = strtoull(argv[2], nullptr, 10);

    int restarts = (argc > 3 ? atoi(argv[3]) : 100);
    int steps    = (argc > 4 ? atoi(argv[4]) : 200000);
    U64 seed     = (argc > 5
        ? strtoull(argv[5], nullptr, 10)
        : (U64)chrono::high_resolution_clock::now().time_since_epoch().count());
    int log_every = (argc > 6 ? atoi(argv[6]) : 0);

    if (!basic_parameter_checks(n, k)) {
        return 0;
    }

    Target T;
    if (!compute_target_profile(n, k, T)) {
        return 0;
    }

    mt19937_64 rng(seed);
    uniform_real_distribution<double> U01(0.0, 1.0);

    U64 global_best = BAD_ID;

    for (int r = 1; r <= restarts; ++r) {
        Graph G = near_regular_initial_graph(n, k, rng);
        Profile P = full_profile(G);

        if (exact_match(P, T)) {
            verified_print_and_exit(G, T);
        }

        U64 cur = l1_score4(P, T);
        U64 best_restart = cur;
        if (cur < global_best) global_best = cur;

        if (log_every > 0) {
            cerr << "[r=" << r << "] start=" << cur << "\n";
        }

        uniform_int_distribution<size_t> distE(0, G.edges.size() - 1);
        uniform_int_distribution<size_t> distN(0, G.nonedges.size() - 1);

        for (int step = 1; step <= steps; ++step) {
            U64 out_eid = G.edges[distE(rng)];
            U64 in_eid  = G.nonedges[distN(rng)];

            apply_swap(G, out_eid, in_eid);
            Profile Q = full_profile(G);

            if (exact_match(Q, T)) {
                verified_print_and_exit(G, T);
            }

            U64 nxt = l1_score4(Q, T);

            double temp = 10.0 * pow(0.99995, step);

            bool accept = false;
            if (nxt <= cur) {
                accept = true;
            } else if (temp > 1e-15) {
                double delta = (double)nxt - (double)cur;
                double prob = exp(-delta / temp);
                if (U01(rng) < prob) accept = true;
            }

            if (accept) {
                P = Q;
                cur = nxt;
                if (cur < best_restart) best_restart = cur;
                if (cur < global_best) global_best = cur;
            } else {
                undo_swap(G, out_eid, in_eid);
            }

            if (log_every > 0 && step % log_every == 0) {
                cerr << "  step " << step
                     << "/" << steps
                     << " current=" << cur
                     << " best_restart=" << best_restart
                     << " best_global=" << global_best
                     << "\n";
            }
        }
    }

    return 0;
}
