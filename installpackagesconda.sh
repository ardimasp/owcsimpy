#!/bin/bash

conda install python=3.7.3
conda install gcc==4.8.5
conda install -c conda-forge openmp==8.0.0
conda install numpy=1.16.3
conda install scipy=1.2.1
conda install matplotlib=3.0.3
conda install cython=0.29.7
conda install jupyter=1.0.0
conda install sphinx=2.0.1
conda install numpydoc=0.9.1
conda install -c conda-forge nbsphinx=0.4.2
conda install sphinx_rtd_theme=0.4.3

pip install JSAnimation
pip install tqdm