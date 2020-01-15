#!/bin/bash

gfortran functions.f95 -c
gfortran cmcalc.f95 -c
gfortran reducekgrid.f95 -c
gfortran read_xanes.f95 -c
gfortran count_xanes.f95 -c
gfortran predict_xanes.f95 -c
gfortran read_setup.f95 -c
gfortran read_xyz.f95 -c
gfortran x2_gradient_descent.f95 -c

gfortran -g -fcheck=all -Wall -o mlxanes mlxanes.f95 functions.o cmcalc.o reducekgrid.o read_xyz.o predict_xanes.o read_xanes.o read_setup.o count_xanes.o x2_gradient_descent.o 
