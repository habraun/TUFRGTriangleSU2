#!/usr/bin/env bash 

#@ job_name = julia_testrun 
#@ job_type = mpich                         
#@ output = julia_testrun.out  
#@ error  = julia_testrun.err              
#@ environment = COPY_ALL                   
#@ notification = complete                  
#@ notify_user = hbraun@fkf.mpg.de 
#@ node = 1
#@ tasks per node=40
#@ queue 

./start.jl
