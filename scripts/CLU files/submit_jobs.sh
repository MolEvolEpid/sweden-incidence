#! /bin/bash
for i in {1..10823}
do
  qsub -V jobs/incidence_mbm_$i.pbs
done
