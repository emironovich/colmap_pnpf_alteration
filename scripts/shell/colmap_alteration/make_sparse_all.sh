 #!/bin/bash

# script for testing all algorithms on the same data

# FILL MISSING DIRECTORIES

# FIX: set dataset_path example:
#             dataset_path=/home/data/courtyard
# FIX: set output_path, example:
#             output_path=/home/data/courtyard/output

mkdir -p $output_path/p35p
mkdir -p $output_path/p4p
mkdir -p $output_path/real

bash ./make_sparse_var.sh $dataset_path $output_path/p35p 35
bash ./make_sparse_var.sh $dataset_path $output_path/p4p 4
bash ./make_sparse_var.sh $dataset_path $output_path/real 3

