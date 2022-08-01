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



function init_flow!(
        N       ::Int64,
        U       ::Float64,
        V1      ::Float64,
        V2      ::Float64,
        V3      ::Float64,
        J       ::Float64,
        t       ::Float64,
        t2      ::Float64,
        t3      ::Float64,
        mu      ::Float64,
        Gamma   ::Bool,
        M       ::Bool,
        K       ::Bool,
        shell   ::Int64,
        )

        L               = 1+3*shell + 3*shell^2
        grid_bosons     = kgrid_initialization(N,shell,200,Gamma,M,K) # gamma,k,m
        bubbles         = bubbles_initialization(L,grid_bosons.N,shell)
        grid_r          = rgrid_initialization(Int64(ceil(sqrt(grid_bosons.N))))
        fv              = fouriervertex_initialization(L,grid_r,shell) #fw
        v               = vertex_initialization(L,grid_bosons.N) 


        println("Momenta:")
        println(grid_bosons.N)
        println("Form factors:")
        println(L)


        #Calculate Flow
        LambdaArr,pmaxv,pmaxw,cmaxv,cmaxw,dmaxv,dmaxw,BubblesGamma,BubblesM = start_flow(t,t2,t3,mu,U,V1,V2,V3,J,grid_bosons,bubbles,grid_r,v,fv)# add two arguments w and fw

        ################################################################################
        #prepare supplemental information for plotting
        ################################################################################
        idxs    = rand(L,L,grid_bosons.N)
        idxlist = []
        for i in 1:length(idxs)
            push!(idxlist, CartesianIndices(idxs)[i])
        end

        idxlist = (idxlist)
        coords  = (symhexagon(N,Gamma,M,K))

        idxarr          = Array{Int64,2}(undef,length(idxlist),3)
        coordsarr       = Array{Float64,2}(undef,grid_bosons.N,2)

        for i in 1:length(idxlist)
                for j in 1:3
                        idxarr[i,j]     = Int64(idxlist[i][j])
                end

        end

        for i in 1:length(coords)
                coordsarr[i,1]  = coords[i][1]
                coordsarr[i,2]  = coords[i][2]
        end
        ################################################################################
        ################################################################################


        #Use superconducting vertex c_sc to calculate leading gap functions

        leadingvalsv,leadingvecsv,scv,vPbuffv=gapper(v,grid_bosons)

        jldopen("data/triangle"*string(N)*string(U)*string(V1)*string(V2)*string(V3)*string(J)*string(round(t,digits=3))*string(round(t2,digits=3))*string(round(t3,digits=3))*string(round(mu,digits=3))*string(Gamma)*string(M)*string(K)*string(shell)*".jld", "w") do file
            write(file, "SCv", scv)  # alternatively, say "@write file A"
            write(file, "VPv", v.P)
            write(file, "VPbuffv", vPbuffv)
        end

        h5open("data/triangle_N_"*string(N)*"_U_"*string(U)*"_V1_"*string(V1)*"_V2_"*string(V2)*"_V3_"*string(V3)*"_J_"*string(J)*"_t_"*string(round(t,digits=3))*"_t2_"*string(round(t2,digits=3))*"_t3_"*string(round(t3,digits=3))*"_mu_"*string(round(mu,digits=3))*"_Gam_"*string(Gamma)*"_M_"*string(M)*"_K_"*string(K)*"_Sh_"*string(shell)*".h5", "w") do file
                write(file, "pv", real.(v.P)) # save also w etc.
                write(file, "cv", real.(v.C))
                write(file, "dv", real.(v.D))

                write(file, "pmaxv", abs.(pmaxv))
                write(file, "cmaxv", abs.(cmaxv))
                write(file, "dmaxv", abs.(dmaxv))

                write(file, "leadingvecs1111v", leadingvecsv)
                write(file, "leadingvals1111v", leadingvalsv)

                write(file, "Lambda", abs.(LambdaArr))

                write(file,"idxlist",idxarr)
                write(file,"coords",coordsarr)

                write(file,"BubblesGamma",BubblesGamma)
                write(file,"BubblesM",BubblesM)
                create_group(file,"Parameters")
                write(file["Parameters"],"N",N)
                write(file["Parameters"],"U",U)
                write(file["Parameters"],"V1",V1)
                write(file["Parameters"],"V2",V2)
                write(file["Parameters"],"V3",V3)
                write(file["Parameters"],"J",J)
                write(file["Parameters"],"t",t)
                write(file["Parameters"],"t2",t2)
                write(file["Parameters"],"t3",t3)
                write(file["Parameters"],"mu",mu)
                write(file["Parameters"],"Shell",shell)
        end

end
