// Simulated annealing over double-edge swaps for a 3-edge-balanced graph in
// G(n,k) with a prescribed degree sequence (Corollary 3.4 of the paper).
// The degree sequence fixes a = sum C(d,2) and s = sum C(d,3); the search
// steers the triangle count c and R = sum_{uv in E} (d_u - 1)(d_v - 1)
// to their forced values c* and R* = p* + 3c*.
//
// Usage: swap3 n cstar Rstar degfile seed steps restarts weight
// Output on success: the graph in graph6 format on stdout.
#include <bits/stdc++.h>
using namespace std;
typedef long long ll; typedef unsigned long long u64;

int n, W;
vector<u64> adjb;
inline bool adj(int u, int v) { return (adjb[(size_t)u * W + (v >> 6)] >> (v & 63)) & 1; }
inline void setb(int u, int v) { adjb[(size_t)u * W + (v >> 6)] |= 1ULL << (v & 63); }
inline void clrb(int u, int v) { adjb[(size_t)u * W + (v >> 6)] &= ~(1ULL << (v & 63)); }
inline int codeg(int u, int v) {
    int c = 0; const u64 *a = &adjb[(size_t)u * W], *b = &adjb[(size_t)v * W];
    for (int i = 0; i < W; ++i) c += __builtin_popcountll(a[i] & b[i]);
    return c;
}


// graph6 encoding of the current graph (vertices 0..n-1)
static void write_graph6(FILE *out) {
    string g;
    if (n < 63) g += char(n + 63);
    else { g += char(126); for (int sh = 12; sh >= 0; sh -= 6) g += char(((n >> sh) & 63) + 63); }
    int acc = 0, nb = 0;
    for (int j = 1; j < n; ++j)
        for (int i = 0; i < j; ++i) {
            acc = (acc << 1) | (adj(i, j) ? 1 : 0);
            if (++nb == 6) { g += char(acc + 63); acc = 0; nb = 0; }
        }
    if (nb) { acc <<= (6 - nb); g += char(acc + 63); }
    fprintf(out, "%s\n", g.c_str());
}

int main(int argc, char **argv) {
    if (argc < 9) { fprintf(stderr, "usage\n"); return 1; }
    n = atoi(argv[1]); ll cstar = atoll(argv[2]); ll Qstar = atoll(argv[3]);  // R*
    vector<int> deg;
    { ifstream f(argv[4]); int x; while (f >> x) deg.push_back(x); }
    u64 seed = strtoull(argv[5], 0, 10); ll steps = atoll(argv[6]);
    int restarts = atoi(argv[7]); double wt = atof(argv[8]);
    if ((int)deg.size() != n) { fprintf(stderr, "bad degfile\n"); return 1; }
    W = (n + 63) / 64;
    mt19937_64 rng(seed);
    for (int r = 0; r < restarts; ++r) {
        // Havel-Hakimi construction
        adjb.assign((size_t)n * W, 0);
        vector<pair<int,int>> E;
        {
            vector<pair<int,int>> rem; for (int v = 0; v < n; ++v) rem.push_back({deg[v], v});
            while (true) {
                sort(rem.rbegin(), rem.rend());
                if (rem.empty() || rem[0].first == 0) break;
                int dv = rem[0].first, v = rem[0].second; rem[0].first = 0;
                for (int i = 1; i <= dv; ++i) {
                    if (i >= (int)rem.size() || rem[i].first == 0) { fprintf(stderr, "not graphical\n"); return 1; }
                    rem[i].first--; E.push_back({v, rem[i].second}); setb(v, rem[i].second); setb(rem[i].second, v);
                }
            }
        }
        size_t k = E.size();
        uniform_int_distribution<size_t> dE(0, k - 1);
        auto tryswap = [&](size_t i, size_t j, bool flip, ll &dc, ll &dQ, int &a, int &b, int &c, int &d) -> bool {
            if (i == j) return false;
            a = E[i].first; b = E[i].second; c = E[j].first; d = E[j].second;
            if (flip) swap(c, d);
            if (a == c || a == d || b == c || b == d) return false;
            if (adj(a, c) || adj(b, d)) return false;
            ll lost = codeg(a, b) + codeg(c, d);
            ll g1 = codeg(a, c) - (adj(b, c) ? 1 : 0) - (adj(d, a) ? 1 : 0);
            ll g2 = codeg(b, d) - (adj(a, d) ? 1 : 0) - (adj(c, b) ? 1 : 0);
            dc = g1 + g2 - lost;
            dQ = (ll)(deg[a] - deg[d]) * (deg[c] - deg[b]);
            return true;
        };
        auto apply = [&](size_t i, size_t j, int a, int b, int c, int d) {
            clrb(a, b); clrb(b, a); clrb(c, d); clrb(d, c);
            setb(a, c); setb(c, a); setb(b, d); setb(d, b);
            E[i] = {a, c}; E[j] = {b, d};
        };
        // randomize
        for (size_t t = 0; t < min<size_t>(20 * k, 5000000); ++t) {
            ll dc, dQ; int a, b, c, d; size_t i = dE(rng), j = dE(rng);
            if (tryswap(i, j, rng() & 1, dc, dQ, a, b, c, d)) apply(i, j, a, b, c, d);
        }
        ll C = 0, Q = 0;
        for (auto &e : E) { C += codeg(e.first, e.second); Q += (ll)(deg[e.first] - 1) * (deg[e.second] - 1); }
        C /= 3;
        auto score = [&](ll C, ll Q) { return wt * llabs(C - cstar) + (double)llabs(Q - Qstar); };
        double cur = score(C, Q), best = cur;
        double T0 = 2.0 * wt, T = T0;
        uniform_real_distribution<double> U(0, 1);
        for (ll s = 1; s <= steps; ++s) {
            ll dc, dQ; int a, b, c, d; size_t i = dE(rng), j = dE(rng);
            if (!tryswap(i, j, rng() & 1, dc, dQ, a, b, c, d)) continue;
            double nx = score(C + dc, Q + dQ);
            if (nx <= cur || U(rng) < exp(-(nx - cur) / T)) {
                apply(i, j, a, b, c, d); C += dc; Q += dQ; cur = nx;
                if (cur < best) best = cur;
                if (C == cstar && Q == Qstar) {
                    // recompute from scratch
                    ll C2 = 0, Q2 = 0;
                    for (auto &e : E) { C2 += codeg(e.first, e.second); Q2 += (ll)(deg[e.first] - 1) * (deg[e.second] - 1); }
                    C2 /= 3;
                    if (C2 != cstar || Q2 != Qstar) { fprintf(stderr, "mismatch!\n"); return 2; }
                    write_graph6(stdout);
                    fprintf(stderr, "found restart %d step %lld\n", r, s);
                    return 0;
                }
            }
            T = max(0.05, T0 * pow(0.9999995, (double)s));
            if (s % 20000000 == 0) fprintf(stderr, "r%d s%lld C-c*=%lld Q-Q*=%lld best=%.0f T=%.3f\n", r, s, C - cstar, Q - Qstar, best, T);
        }
    }
    fprintf(stderr, "not found\n");
    return 1;
}
