# MARS via LASSO

This repository contains the code for the numerical experiments in the 
paper [MARS via LASSO](https://arxiv.org/abs/2111.11694). You can find 
all R scripts in the folder `code`. To run those R scripts, you need to 
install the R package `regmdc`. Please refer to this 
[page](https://github.com/DohyeongKi/regmdc) for how to install it. We 
used the version 0.5.0 for our experiments.

We put the datasets we used for the experiments into the folder `data`. 
We obtained the [Airfoil Self-Noise](data/airfoil_self_noise.dat) 
dataset from 
[here](https://archive.ics.uci.edu/dataset/291/airfoil+self+noise),
the [Concrete](data/Concrete_Data.xls) dataset from 
[here](https://archive.ics.uci.edu/dataset/165/concrete+compressive+strength),
and the [wine](data/winequality-red.csv) dataset from 
[here](https://archive.ics.uci.edu/dataset/186/wine+quality).