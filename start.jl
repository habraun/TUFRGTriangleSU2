
using Pkg

#==
If used for the first time, choose package_install=1 to install
necessary packages
==#

package_install=1

if package_install==1
    Pkg.add("HDF5")
    Pkg.add("LoopVectorization")
    Pkg.add("StaticArrays")
    Pkg.add("TensorOperations")
    Pkg.add("JLD")
    Pkg.add("Roots")
end

#Load Tufrg code
include("main.jl")



#Set parameters, example parameters are for simple t-U Hubbard model at Van Hove filling

N       = 4                     #Selects momentum resolution, e.g. N=4 equals 180 points, N=5 equals 336 points, N=6 equals 768 points
U       = 4.0                   #Hubbard Interaction
V1      = 0.0                   #Nearest, n-nearest and n-n-nearest neigbhbour Interactions
V2      = 0.0
V3      = 0.0
J       = 0.0                   #Magnetic exchange coupling, PROTOTYPICAL, better keep at J=0.0
t       = 1.0                   #Tight-binding hoppings
t2      = 0.0
t3      = 0.0
mu      = 2.0                   #Chemical potential
Gamma   = false                 #Set to "true" for higher momentum resolution at Gamma,M or K point.
M       = false
K       = false
shell   = 2                     #Hexagon shell of form factors. shell=2 equaios 19 form factors, which is the minimum one should use


#start flow
init_flow!(N,U,V1,V2,V3,J,t,t2,t3,mu,Gamma,M,K,shell)
