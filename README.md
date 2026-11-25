# t-edge-balanced

Programs and data for the paper

> Y. M. Chee, *On t-Edge-Balanced Graphs*, 2026.

A graph $G$ with $n$ vertices and $k$ edges is *t-edge-balanced* if every graph
with $n$ vertices and $t$ edges is contained in the same number of subgraphs of
$K_n$ isomorphic to $G$. Equivalently, the orbit of $G$ under the symmetric group
is the block set of a graphical $t$-design on the edges of $K_n$.

For $n \ge 6$, a graph $G$ with $n$ vertices and $k$ edges is 3-edge-balanced if
and only if
- $\sum_v \binom{\deg v}{2} = a^* = 2k(k-1)/(n+1)$,
- $\sum_v \binom{\deg v}{3} = s^* = (n-3)c^*$,
- the number of triangles is $c^* = 4k(k-1)(k-2)/(3(n+1)(n^2-n-4))$, and
- $R = \sum_{uv \in E} (\deg u - 1)(\deg v - 1) = p^* + 3c^*$, where $p^* = 3(n-3)c^*$

(Corollary 3.4 of the paper). A pair $(n,k)$ with $n \ge 6$ and
$4 \le k \le \binom{n}{2}/2$ is *admissible* if $a^*$ and $c^*$ are integers.

## Contents

- `src/admissible.py`: lists the admissible pairs $(n,k)$ with `LO <= n <= HI`.
- `src/targets.py`: the forced values $a^*$, $c^*$, $p^*$, $s^*$, $R^*$.
- `src/degseq.py`: Phase 1 of the search. Finds a graphical degree sequence with the forced values of $\sum \binom{d_i}{2}$ and $\sum \binom{d_i}{3}$.
- `src/swap3.cpp`: Phase 2 of the search. Simulated annealing over degree-preserving double-edge swaps, steering the triangle count and $R$ to their forced values. Writes the graph in graph6 format.
- `src/run.sh`: runs both phases.
- `src/verify3.cpp`: an independent check. Reads a graph in graph6 format and computes its numbers of subgraphs isomorphic to $K_3$, $P_4$, $K_{1,3}$, $P_3 \cup K_2$, and $3K_2$ from the definitions, and compares them with the forced values.
- `results/`: the 3-edge-balanced graphs in graph6 format, one for each of the 36 admissible pairs with $n \le 3000$.

## Requirements

A C++17 compiler and Python 3. `src/admissible.py` also requires `sympy`.

## Usage

```bash
g++ -O3 -std=c++17 -o src/swap3 src/swap3.cpp
g++ -O3 -std=c++17 -o src/verify3 src/verify3.cpp

python3 src/admissible.py 6 3000         # list admissible pairs
./src/run.sh 13 21                       # writes graph_13_21.g6 (seed 1 by default)
./src/verify3 13 21 graph_13_21.g6       # prints "3-edge-balanced" and the counts
```

`./src/run.sh n k seed` uses the given random seed. With seed 1, it reproduces
the graphs in `results/`. The graph6 files can be read by nauty, SageMath, and
NetworkX (`networkx.read_graph6`).

## Results

Each graph below has trivial automorphism group, so it gives a graphical
3-$(N,k,\lambda)$ design with $N = \binom{n}{2}$, consisting of $n!$ blocks, with
$\lambda = 6(n-3)!\,c^*$.

| $n$ | $k$ | $N$ | $\lambda$ | File |
|---:|---:|---:|---:|---|
| 13 | 21 | 78 | 30·10! | `results/n13_k21.g6` |
| 21 | 66 | 210 | 240·18! | `results/n21_k66.g6` |
| 34 | 260 | 561 | 3552·31! | `results/n34_k260.g6` |
| 77 | 1378 | 2926 | 45792·74! | `results/n77_k1378.g6` |
| 84 | 1275 | 3486 | 27930·81! | `results/n84_k1275.g6` |
| 114 | 3151 | 6441 | 168840·111! | `results/n114_k3151.g6` |
| 139 | 3570 | 9591 | 135456·136! | `results/n139_k3570.g6` |
| 203 | 5644 | 20503 | 171864·200! | `results/n203_k5644.g6` |
| 229 | 11546 | 26106 | 1025196·226! | `results/n229_k11546.g6` |
| 255 | 8321 | 32385 | 277890·252! | `results/n255_k8321.g6` |
| 292 | 12599 | 42486 | 642498·289! | `results/n292_k12599.g6` |
| 309 | 2976 | 47586 | 7140·306! | `results/n309_k2976.g6` |
| 311 | 10621 | 48205 | 318570·308! | `results/n311_k10621.g6` |
| 414 | 37101 | 85491 | 5757360·411! | `results/n414_k37101.g6` |
| 444 | 23941 | 98346 | 1254078·441! | `results/n444_k23941.g6` |
| 451 | 15369 | 101475 | 316536·448! | `results/n451_k15369.g6` |
| 459 | 30706 | 105111 | 2394912·456! | `results/n459_k30706.g6` |
| 531 | 29659 | 140715 | 1393926·528! | `results/n531_k29659.g6` |
| 539 | 49896 | 144991 | 6346032·536! | `results/n539_k49896.g6` |
| 571 | 61776 | 162735 | 10130400·568! | `results/n571_k61776.g6` |
| 575 | 75520 | 165025 | 18124320·572! | `results/n575_k75520.g6` |
| 620 | 72496 | 191890 | 12789240·617! | `results/n620_k72496.g6` |
| 748 | 67089 | 279378 | 5771964·745! | `results/n748_k67089.g6` |
| 948 | 170820 | 448878 | 46803240·945! | `results/n948_k170820.g6` |
| 949 | 98401 | 449826 | 8918238·946! | `results/n949_k98401.g6` |
| 1124 | 314875 | 631126 | 175874622·1121! | `results/n1124_k314875.g6` |
| 1187 | 132517 | 703891 | 11131260·1184! | `results/n1187_k132517.g6` |
| 1729 | 352575 | 1493856 | 67835430·1726! | `results/n1729_k352575.g6` |
| 2283 | 757718 | 2604903 | 292478376·2280! | `results/n2283_k757718.g6` |
| 2295 | 846077 | 2632365 | 400839450·2292! | `results/n2295_k846077.g6` |
| 2331 | 112519 | 2715615 | 899766·2328! | `results/n2331_k112519.g6` |
| 2539 | 268351 | 3221991 | 9445110·2536! | `results/n2539_k268351.g6` |
| 2595 | 619619 | 3365715 | 108905622·2592! | `results/n2595_k619619.g6` |
| 2899 | 304675 | 4200651 | 9286494·2896! | `results/n2899_k304675.g6` |
| 2929 | 1056851 | 4288056 | 375813330·2926! | `results/n2929_k1056851.g6` |
| 2984 | 221886 | 4450636 | 3289104·2981! | `results/n2984_k221886.g6` |

## Larger orders

There are 27 admissible pairs with $3000 < n \le 10^4$ (`python3 src/admissible.py 3001 10000`).
For each of them, `./src/run.sh n k` finds a 3-edge-balanced graph in at most two
minutes on a computer with two cores. These graphs have up to $1.6 \cdot 10^7$ edges
and are not included here.

## License

MIT. See `LICENSE`.
