#!/bin/bash
# FILL MISSING DIRECTORIES

# usage:
# bash make_sparse_var.sh dataset_path output_path algo_num

# FIX: cd colmap_install_path/bin, example:
#         cd /usr/local/bin

 DATASET_PATH=$1
 OUTPUT_PATH=$2

 mkdir $DATASET_PATH/sparse

 ./colmap mapper \
    --database_path $DATASET_PATH/database.db \
    --image_path $DATASET_PATH/images \
    --output_path $DATASET_PATH/sparse \
    --Mapper.abs_pose_algo $3

 ./colmap model_converter --input_path $DATASET_PATH/sparse/0 --output_path $OUTPUT_PATH --output_type TXT

 mv *.txt $OUTPUT_PATH

 rm -r $DATASET_PATH/sparse

