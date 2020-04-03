# FILL MISSING DIRECTORIES

# FIX: cd colmap_install_path/bin, example:
#         cd ~/colmap_install/all/bin

# FIX: set dataset_path, example:
#         dataset_path=/home/emironovich/hdd_scratch/courtyard

 ./colmap feature_extractor \
   --database_path $dataset_path/database.db \
   --image_path $dataset_path/images \
   --ImageReader.camera_model "PINHOLE"

