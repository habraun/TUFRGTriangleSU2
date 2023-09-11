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

using Plots;gr()

t=1.0
t2=1/9
t3=0.0
mu=1.0*2*(t+t2)


####prepare objects###########
shell           = 2
L               = 1+3*shell + 3*shell^2
grid_bosons     = kgrid_initialization(4,shell,200,false,false,false)
bubbles         = bubbles_initialization(L,grid_bosons.N,shell)

sumph           = zeros(Float64,length(bubbles.formfactorgrid))*im
sumpp           = zeros(Float64,length(bubbles.formfactorgrid))*im
init0           = zeros(Float64,length(bubbles.formfactorgrid))*im
ff              = zeros(Float64,length(bubbles.formfactorgrid))*im
res             = zeros(Float64,length(bubbles.formfactorgrid))*im
buff1           = zeros(Float64,length(bubbles.formfactorgrid))*im
####################################


testoX          = []
testoY1         = []
testoY2         = []

testopp0        = zeros(50,4)
testoph0        = zeros(50,4)
testoppM        = zeros(50,4)
testophM        = zeros(50,4)

for i in 1:1:50
    println("")
    Lambda=10.0*exp(-i*0.25)
    faktor=Int64(ceil(log(10,10/Lambda)))

    for j in 1:4
        println(j)
        ph1,pp1=get_bubbles_qiadaptive(j*96,3*(2^faktor),grid_bosons,1,Lambda,t,t2,t3,mu,ff,init0,sumph,sumpp,res,buff1,bubbles)
        phM,ppM=get_bubbles_qiadaptive(j*96,3*(2^faktor),grid_bosons,11,Lambda,t,t2,t3,mu,ff,init0,sumph,sumpp,res,buff1,bubbles)

        testoph0[i,j]   = real(ph1[1])
        testopp0[i,j]   = real(pp1[1])

        testophM[i,j]   = real(phM[1])
        testoppM[i,j]   = real(ppM[1])
    end

    push!(testoX,Lambda)

    println(i)
    println(Lambda)
end

plot( testoX,testoX.*testoph0[:,1] ,lab =  "ph(0), 120" ,xaxis=:log)
plot!(testoX,-testoX.*testopp0[:,1] ,lab =  "pp(0), 120" ,xaxis=:log)
plot!(testoX,testoX.*testophM[:,1] ,lab =  "ph(M), 120" ,xaxis=:log)
plot!(testoX,testoX.*testoppM[:,1] ,lab =  "pp(M), 120" ,xaxis=:log)

plot!(testoX,testoX.*testoph0[:,2] ,lab =  "ph(0), 240" ,xaxis=:log)
plot!(testoX,-testoX.*testopp0[:,2] ,lab =  "pp(0), 240" ,xaxis=:log)
plot!(testoX,testoX.*testophM[:,2] ,lab =  "ph(M), 240" ,xaxis=:log)
plot!(testoX,testoX.*testoppM[:,2] ,lab =  "pp(M), 240" ,xaxis=:log)

plot!(testoX,testoX.*testoph0[:,3] ,lab =  "ph(0), 360" ,xaxis=:log)
plot!(testoX,-testoX.*testopp0[:,3] ,lab =  "pp(0), 360" ,xaxis=:log)
plot!(testoX,testoX.*testophM[:,3] ,lab =  "ph(M), 360" ,xaxis=:log)
plot!(testoX,testoX.*testoppM[:,3] ,lab =  "pp(M), 360" ,xaxis=:log)

plot!(testoX,testoX.*testoph0[:,4] ,lab =  "ph(0), 480" ,xaxis=:log)
plot!(testoX,-testoX.*testopp0[:,4] ,lab =  "pp(0), 480" ,xaxis=:log)
plot!(testoX,testoX.*testophM[:,4] ,lab =  "ph(M), 480" ,xaxis=:log)
plot!(testoX,testoX.*testoppM[:,4] ,lab =  "pp(M), 480" ,xaxis=("T/t",:log), legend = :outertopright, yaxis=("T*Bubble"))


savefig("HOVH.svg")
