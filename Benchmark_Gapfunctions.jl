using Pkg
using Distributed
rmprocs(workers())
addprocs(8)
nworkers()

@everywhere using HDF5
@everywhere using LinearAlgebra
@everywhere using LoopVectorization
@everywhere using StaticArrays
@everywhere using TensorOperations
@everywhere using JLD
@everywhere using Roots
@everywhere using SharedArrays
@everywhere using ProgressBars

include("BubbleStructure.jl")
include("GridStructure.jl")
include("BosonicGrid.jl")
include("VertexStructure.jl")
include("Formfactors.jl")
include("BareHamiltonian.jl")
include("BubbleIntegration.jl")
include("Fourier.jl")
include("Projection.jl")
include("Increment.jl")
include("DiffEq.jl")
include("GapEquations.jl")

N           = 8                     #Selects momentum resolution, e.g. N=4 equals 180 points, N=5 equals 336 points, N=6 equals 768 points
U           = 3.0                   #Hubbard Interaction
V1          = 0.0                   #Nearest, n-nearest and n-n-nearest neigbhbour Interactions
V2          = 0.0
V3          = 0.0
J           = 0.4                   #Magnetic exchange coupling, PROTOTYPICAL, better keep at J=0.0
t           = 1.0                   #Tight-binding hoppings
t2          = 0.0
t3          = 0.0
mu          = 2.0                   #Chemical potential
Gamma       = false                 #Set to "true" for higher momentum resolution at Gamma,M or K point.
M           = false
K           = false
shell       = 4                     #Hexagon shell of form factors. shell=2 equaios 19 form factors, which is the minimum one should use
Vh          = 0.0
Kh          = 0.0
phifaktor   = 3
NM          = [2,2]

println("initialization complete")
data=h5open("triangle_N_"*string(N)*"_U_"*string(U)*"_V1_"*string(V1)*"_V2_"*string(V2)*"_V3_"*string(V3)*"_J_"*string(J)*"_t_"*string(round(t,digits=3))*"_t2_"*string(round(t2,digits=3))*"_t3_"*string(round(t3,digits=3))*"_mu_"*string(round(mu,digits=3))*"_Gam_"*string(Gamma)*"_M_"*string(M)*"_K_"*string(K)*"_Sh_"*string(shell)*"_phifaktor_"*string(phifaktor)*".h5", "cw")
println("Opened file")
pv = data["pv"]
dv = data["dv"]
grid_bosons = kgrid_initialization(N,shell,200,Gamma,M,K)
println("")
#leadingvalsvp,leadingvecsvp,scv,vPbuffv= gapper(pv,grid_bosons)
leadingvalsvd,leadingvecsvd,dv,vDbuffv=@sync dgapper(dv,grid_bosons)

#if haskey(data,"leadingvecs1111vp")
#    delete_object(data,"leadingvecs1111vp")
#    delete_object(data,"leadingvals1111vp")
#end
if  haskey(data,"leadingvecs1111vd")
    delete_object(data,"leadingvecs1111vd")
    delete_object(data,"leadingvals1111vd")
end
if  haskey(data,"leadingvecs1111v")
    delete_object(data,"leadingvecs1111v")
    delete_object(data,"leadingvals1111v")
end
#write(data, "leadingvecs1111vp", leadingvecsvp)
#write(data, "leadingvals1111vp", leadingvalsvp)
write(data, "leadingvecs1111vd", leadingvecsvd)
write(data, "leadingvals1111vd", leadingvalsvd)

close(data)
