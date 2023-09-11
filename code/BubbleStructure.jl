###Structures storing the particle-particle and particle-hole bubbles.
# the struct "bubble" has the following fields:
# .ph and .pp, a complex array with 3 indices, for two form factors and one momentum index
# .L the amount of form factors used in a calculation
# .ffidxs an array of tuples, indicating which form factor index belongs to which integer combination
# of the real space vectors, e.g. [2,3] corresponds to 2*R_1 + 3*R_2 with R_1,2 the lattice vectors of the real space lattice
# .formfactorgrid
# .formfactorgridMut
#  .mn_arr



@everywhere struct bubble
    ph::Array{Complex{Float64},3}
    pp::Array{Complex{Float64},3}
    L::Int64
    ffidxs::Array{Tuple{Int64,Int64},1}
    formfactorgrid::Array{SArray{Tuple{2},Int64,1,2},1}
    formfactorgridMut ::Array{Array{Int64,1},1}
    mn_arr::Array{Int64,2}
end

@everywhere function bubbles_initialization(
    L           ::Int64,
    N           ::Int64,
    shell       ::Int64
    )           ::bubble

    XPH=Array{Complex{Float64},3}(undef,L,L,N)
    XPP=Array{Complex{Float64},3}(undef,L,L,N)
    XPH.=0.0+0.0
    XPP.=0.0+0.0

    ffidxs=[(1,1)]

    for i in 1:L
        for j in i+1:L
            push!(ffidxs,(i,j))
        end
    end

    for i in 3:2:L
        push!(ffidxs,(i,i-1))
    end

    formfactorgrid=formfactor_grid(shell*2)
    formfactorgrid_mutable=formfactor_grid_mutable(shell*2)
    arr=generate_mn_to_site_arr(L,formfactorgrid_mutable)

    grid_bosons=bubble(XPH,XPP, L,ffidxs,formfactorgrid,formfactorgrid_mutable,arr)

    return grid_bosons
end
