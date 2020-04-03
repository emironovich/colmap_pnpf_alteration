 #!/bin/bash

# script for testing all algorithms on the same data

# FILL MISSING DIRECTORIES

# FIX: set dataset_path example:
#             dataset_path=/home/emironovich/hdd_scratch/courtyard
# FIX: set output_path, example:
#             output_path=/home/emironovich/data/courtyard/ransac_16_02_snd

mkdir $output_path/p35p
mkdir $output_path/p4p
mkdir $output_path/real

bash ./make_sparse_var.sh $dataset_path $output_path/p35p 35
bash ./make_sparse_var.sh $dataset_path $output_path/p4p 4
bash ./make_sparse_var.sh $dataset_path $output_path/real 3

