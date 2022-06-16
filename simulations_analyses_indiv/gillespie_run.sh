##!/bin/bash

filename='gillespie_sims.py'
n_cells=5
n_mins=1200
queue_type='vshort'

for i in {1..2000} #num runs to parallelize, each has n_cells number of cells
do
    command="python ${filename} --n_cells ${n_cells} --n_mins ${n_mins} --batch ${i}"
    bsub -J ${i}[1-10] -q $queue_type -R 'select[hname!=cn001]' -R 'select[hname!=cn002]' -R 'select[hname!=cn003]' -R 'select[hname!=cn004]' -R 'select[hname!=cn005]' $command
done