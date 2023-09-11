#!/usr/bin/env bash

#@ job_name = n12
#@ job_type = mpich
#@ output = mujscan.out
#@ error  = mujscan.err
#@ notification = complete
#@ notify_user = hbraun@fkf.mpg.de
#@ node = 1
#@ class = 128core_new
#@ tasks_per_node=256
#@ queue

module load julia/1.7.3
which julia
julia --startup-file=no ./start.jl
