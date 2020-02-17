# FILL MISSING DIRECTORIES

# FIX: cd colmap_install_path/bin, example:
#         cd ~/colmap_install/all/bin

# FIX: set DATASET_PATH, example:
#         DATASET_PATH=/home/emironovich/hdd_scratch/courtyard

 ./colmap feature_extractor \
   --database_path $DATASET_PATH/database.db \
   --image_path $DATASET_PATH/images \
   --ImageReader.camera_model "PINHOLE"

