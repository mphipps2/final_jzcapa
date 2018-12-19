#!/bin/bash

 for t in $(seq 0 100); do
     ./build/zdc run.mac results/Event$t.root >> ~/zdc/results/out$t.txt 2>> ~/zdc/results/err$t.txt
 done
