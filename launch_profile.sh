#! /bin/bash
export name=$1-$HOSTNAME-$(date '+%Y-%m-%d')
time python3 -m cProfile -o $name.prof $1.py --log $name.log $name.txt 
gprof2dot --colour-nodes-by-selftime -n 0.05 -e 0.01 --skew=0.05 -f pstats $name.prof -o $name.dot
dot -Tpdf $name.dot -o $name.pdf


