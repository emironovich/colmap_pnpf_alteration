# COLMAP

This repository containts implementations of the algorithm described in "P3,5P: Estimation with Unknown Focal Length"
by Changchang Wu and "Efficient Intersection of Three Quadrics and Applications in Computer Vision" by Zuzana Kukelova
and others. It was done as a part of bachelor's thesis "Fast localisation of uncalibrated cameras".

### Dependancies
Use COLMAP's [original instructions](https://colmap.github.io/install.html) for installation, but make sure to use
ceres-solver version 1.14.0.

### Usage
This colmap alteration containts pnpf library as a submodule. Before building run this:
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

