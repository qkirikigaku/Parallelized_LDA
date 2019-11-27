#!/bin/bash
#$ -S /bin/bash

data_type=${1}

if [-e result/data${data_type}]; then
    rm -r result/data${data_type}
fi

#----------------------
# Data_type shows the mutation dictionary using in experiments.
# In detail, please refer to our paper.
#----------------------

Iter=100

#----------------------
# Iter shows the iteration of experiments.
# Please set Iter for large value to avoid local-minimum.
# In our paper, setting Iter to 100.
#----------------------

K=30

#-----------------------
# K shows the max number of topics.
#-----------------------


#-----------------------
# For queueing systems.
#-----------------------

# python runms.py ${data_type} ${SGE_TASK_ID}

# Comment out above line and use "qsub" command.
# Set ${SGE_TASK_ID} equal to ${K}.
# If you wont to conduct experiments with  queueing systems,
# you should set absolute path according to your environments.

#-----------------------
# For local
#-----------------------

for i in $(seq 1 ${Iter}); do
    python scripts/runms.py ${data_type} ${i} ${K}
done

bash scripts/make_figure.sh ${data_type} ${Iter} ${K}
