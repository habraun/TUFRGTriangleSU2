#!/usr/bin/env bash 

#@ job_name = julia_testrun 
#@ job_type = mpich                         
#@ output = julia_testrun.out  
#@ error  = julia_testrun.err                                 
#@ notification = complete                  
#@ notify_user = hbraun@fkf.mpg.de 
#@ node = 1
#@ class = 32core
#@ tasks_per_node=40
#@ queue 

module load julia
which julia
julia --startup-file=no ./start.jl
