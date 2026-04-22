# t-edge-balanced
Code and data for the first known examples of 3-edge-balanced graphs

This repository accompanies the paper:

> Y. M. Chee, "On t-Edge-Balanced Graphs", 2026.

## Overview

A graph $G$ on $n$ vertices with $k$ edges is *t-edge-balanced* if every
graph on $n$ vertices with $t$ edges is contained in exactly the same number
of subgraphs of $K_n$ isomorphic to $G$. This paper gives the first known
examples of 3-edge-balanced graphs, and proves that no nontrivial
t-edge-balanced graph exists for $t \geq 4$.

## Repository Contents

- `src/sa3eb.cpp` — C++ source code for the simulated annealing search.
- `results/` — Edge lists and graph6 encodings of the 3-edge-balanced graphs found.

## Compiling the Search Code

The code requires a C++17-compatible compiler. To compile:

```bash
g++ -O3 -std=c++17 -o sa3eb src/sa3eb.cpp
```

## Running the Search

```bash
./sa3eb n k [restarts] [steps] [seed] [log_every]
```

**Arguments:**
- `n` — number of vertices
- `k` — number of edges
- `restarts` — number of restarts (default: 100)
- `steps` — number of steps per restart (default: 200000)
- `seed` — random seed (default: time-based)
- `log_every` — log progress every this many steps (default: 0, no logging)

**Example:**

```bash
./sa3eb 13 21
./sa3eb 21 66 200 500000
./sa3eb 34 260 500 1000000 12345 10000
```

If a 3-edge-balanced graph is found, its edge list is printed to standard
output, one edge per line as `u v` with vertices numbered from 1.

## Results

The following 3-edge-balanced graphs were found. Each parameter set has its
own subfolder under `results/`, containing:
- `graph.txt` — edge list, one edge per line as `u v`
- `graph.g6` — graph6 encoding, readable by SageMath, nauty, and NetworkX

| n   | k    | Subfolder          |
|-----|------|--------------------|
| 13  | 21   | results/n13_k21/   |
| 21  | 66   | results/n21_k66/   |
| 34  | 260  | results/n34_k260/  |
| 77  | 1378 | results/n77_k1378/ |
| 84  | 1275 | results/n84_k1275/ |
| 114 | 3151 | results/n114_k3151/|
| 139 | 3570 | results/n139_k3570/|
| 203 | 5644 | results/n203_k5644/|
| 229 | 11546 | results/n229_k11546/|
| 255 | 8321 | results/n255_k8321/|
| 309 | 2976 | results/n309_k2976/|

The first 10 parameter sets are
the 10 smallest for which a 3-edge-balanced graph can exist, as proved in
the paper. All 11 parameter sets satisfy the necessary arithmetic conditions
(C1) and (C2) derived in the paper.

## Requirements

- C++17-compatible compiler (e.g. GCC 7 or later, Clang 5 or later)
- No external libraries are required for either the search or the verification.

## Citation

If you use this code or these graphs in your work, please cite:

```bibtex
@article{Chee:2026,
      author = "Y. M. Chee",
       title = "On $t$-Edge-Balanced Graphs",
        year = "2026"
}
```

## Funding

Research supported by SUTD Grant SKI 2021_07_04.

## License

MIT License. See `LICENSE` for details.
