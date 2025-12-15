#! /bin/bash
for i in {1..10823}
do
  echo > jobs/incidence_mbm_$i.pbs
  echo "#! /bin/bash" > jobs/incidence_mbm_$i.pbs
  echo "#PBS -l nodes=1:ppn=1" >> jobs/incidence_mbm_$i.pbs
  echo "source ~/miniconda3/etc/profile.d/conda.sh" >> jobs/incidence_mbm_$i.pbs
  echo "conda activate" >> jobs/incidence_mbm_$i.pbs
  echo "cd ~/Documents/inference/incidence/clu" >> jobs/incidence_mbm_$i.pbs
  echo "Rscript infer_inftimes_clu_new.R $i" >> jobs/incidence_mbm_$i.pbs
done
