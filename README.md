# proto-Nucleic Acid Builder (pNAB)
## Try the Package: [![Binder](https://mybinder.org/badge_logo.svg)](https://mybinder.org/v2/gh/alenaizan/pnab/v1.0.0?filepath=binder%2Fdriver.ipynb)
## Read the documentation: [![docs latest](https://img.shields.io/badge/docs-latest-5077AB.svg?logo=read%20the%20docs)](https://alenaizan.github.io/pnab/html/index.html)

[![Build Status](https://travis-ci.com/alenaizan/pnab.svg?branch=master)](https://travis-ci.com/alenaizan/pnab)
[![Build status](https://ci.appveyor.com/api/projects/status/9uew76hr6ftlwso4/branch/master?svg=true)](https://ci.appveyor.com/project/alenaizan/pnab/branch/master)
[![Language grade: C/C++](https://img.shields.io/lgtm/grade/cpp/g/alenaizan/pnab.svg?logo=lgtm&logoWidth=18)](https://lgtm.com/projects/g/alenaizan/pnab/context:cpp)
[![Language grade: Python](https://img.shields.io/lgtm/grade/python/g/alenaizan/pnab.svg?logo=lgtm&logoWidth=18)](https://lgtm.com/projects/g/alenaizan/pnab/context:python)

## Overview
The proto-Nucleic Acid Builder is a program for modeling the 3D strucutres of DNA, RNA, and nucleic acid analogs. Nucleic acids with alternative backbones or nucleobases can be constructed by the program by supplying the 3D structure of isolated backbones or nucleobases. The program can perform a helical parameter search and backbone conformation search and find reasonable nucleic acid structures.
![image](Doxygen/images/output.png)

## Installing the Package Using the conda package manager
To install the conda package, first install miniconda/anaconda. Then,
```
conda config --add channels conda-forge
conda install -c alenaizan pnab
```

## Compiling the Package in Linux
If you prefer to compile the package, install miniconda or anaconda. Then,

```
git clone https://github.com/alenaizan/pnab.git
cd pnab
sh install.sh 
```
