# COLMAP

This repository contains implementations of the algorithm described in "P3,5P: Estimation with Unknown Focal Length"
by Changchang Wu and "Efficient Intersection of Three Quadrics and Applications in Computer Vision" by Zuzana Kukelova
and others. It was done as a part of bachelor's thesis "Fast localisation of uncalibrated cameras".

### Dependancies
Use COLMAP's [original instructions](https://colmap.github.io/install.html) for installation, but make sure to use
ceres-solver version 1.14.0.

### Usage
This colmap alteration contains pnpf library as a submodule. Before building run this:
```
bash
git submodule update --init --recursive
```

To use specific algorithms in `colmap mapper` choose the version in the parameters:
```
colmap mapper \
    --Mapper.abs_pose_algo 3 # 3 for the original algo
                             # 35 for P3.5Pf by Wu
                             # 4 for P4Pf by Kukelova
```

To test creating the reconstruction using every type of absolute pose estimation (colmap, p3.5pf, p4pf) you can use
template of the scripts in `scripts/shell/colmap_alteration`.

To use these scripts you will need to edit them with the COLMAP installation directory, your project dataset directory and
the output directory.
These scripts assume that the dataset directory contains an `images` directory with images inside it.

After *manually* changing the paths in all the scripts (look for `#FIX: set ...` comments) run:
```
bash extract_features.sh
bash make_sparse_all.sh
```

In the output directory there will be created a corresponding directory for every used algorithm.
