#!/bin/bash

cd src/
mpiexec --quiet -np 4 ../build/cluster
cd ..
