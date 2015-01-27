#!/bin/bash

cd R_scripts
Rscript hiHMM_Pre1.r
cd ..

cd MATLAB_scripts
matlab < driv_hihmm.m
cd ..

