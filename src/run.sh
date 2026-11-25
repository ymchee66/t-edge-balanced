#!/bin/bash
# Usage: ./run.sh n k [seed]
# Writes the degree sequence to deg_n_k.txt and the graph in graph6 format to graph_n_k.g6.
set -e
n=$1; k=$2; seed=${3:-1}
dir=$(dirname "$0")
read c R < <(python3 -c "
import sys; sys.path.insert(0,'$dir')
from targets import targets
t=targets($n,$k)
assert t['c'].denominator==1 and t['a'].denominator==1, 'not admissible'
print(int(t['c']), int(t['R']))")
python3 "$dir/degseq.py" $n $k $seed > deg_${n}_${k}.txt
"$dir/swap3" $n $c $R deg_${n}_${k}.txt $seed 2000000000 5 10 > graph_${n}_${k}.g6
