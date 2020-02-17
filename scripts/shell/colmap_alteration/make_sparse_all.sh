 #!/bin/bash

# script for testing all algorithms on the same data

# FILL MISSING DIRECTORIES

# FIX: set dataset_path example:
#             dataset_path=/home/emironovich/hdd_scratch/courtyard
# FIX: set output_path, example:
#             output_path=/home/emironovich/data/courtyard/ransac_16_02_snd
# FIX: set trials_num, example:
#             trials_num=1

rm -r $output_path/*
mkdir $output_path/p35p
mkdir $output_path/p4p
mkdir $output_path/real

for ((i=1;i<=trials_num;i++));
do

  mkdir $output_path/p35p/$i
  mkdir $output_path/p4p/$i
  mkdir $output_path/real/$i

  bash ./make_sparse_var.sh $dataset_path $output_path/p35p/$i 35
  bash ./make_sparse_var.sh $dataset_path $output_path/p4p/$i 4
  bash ./make_sparse_var.sh $dataset_path $output_path/real/$i 3

done
