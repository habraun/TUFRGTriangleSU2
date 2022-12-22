@everywhere function get_formfactor_explicit(
        kx  ::Float64,
        ky  ::Float64,
        site::SArray,
        R1  ::SArray,
        R2  ::SArray
    )       ::Complex{Float64}

    Rarr=@SArray[R1,R2]
    R=site'*Rarr
    k=SVector(kx,ky)

    return exp(im*k'*R)
end

@everywhere function gapper(
    pv::HDF5.Dataset,
    grid_bosons::kgrid
    )
    println("begins pairing gap equation")
    c_sc    = Array{Complex{Float64},2}(undef,grid_bosons.N,grid_bosons.N)
    c_sc    .= 0.0+0.0*im

    sites=grid_bosons.formfactorgrid

    R1  = SVector(3.0/2.0,sqrt(3.0)/2.0)
    R2  = SVector(3.0/2.0,-sqrt(3.0)/2.0)

    Np      = Int64(grid_bosons.N/12)

    vPbuff  = similar(pv[:,:,1])

    vPbuff  .=0.0+0.0*im

    for n in 1:size(pv)[1]
        for m in 1:size(pv)[1]
            for i in 1:12
                vPbuff[m,n] += pv[m,n,1+(i-1)*Np]./12
            end
        end
    end

    Midx=1
    maxQ = findmax(abs.(pv[:,:,:]))[2][3]

    println("maxQidx: ",maxQ," | Gammaidx: ",Midx)
    println("maxQ: " ,grid_bosons.grid[maxQ]," | Gamma: ",grid_bosons.grid[Midx])
    println("PV[maxQ]:",findmax(abs.(pv[:,:,maxQ]))[1]," | PV[Gamma]: ", findmax(abs.(pv[:,:,1]))[1])

    println("vertex at q calculated")
    for k in tqdm(1:grid_bosons.N)
        for n in 1:size(pv)[1]
            ffk = get_formfactor_explicit(grid_bosons.grid[k][1],grid_bosons.grid[k][2],sites[n],R1,R2)
            for p in 1:grid_bosons.N
                for m in 1:size(pv)[1]
                    ffp         = get_formfactor_explicit(grid_bosons.grid[p][1],grid_bosons.grid[p][2],sites[m],R1,R2)
                    c_sc[p,k]   += vPbuff[m,n]*(conj(ffp)*ffk)
                end
            end
        end
    end

    println("gap calculation")
    eigensystem     = eigen(-c_sc[:,:])
    vals            = eigensystem.values
    u               = eigensystem.vectors
    vals            = real.(vals)
    leadingvals     = real.(vals[end-3:end])
    leadingvecs     = real.(u[:,end-3:end])

    return leadingvals,leadingvecs, -c_sc,vPbuff
end

@everywhere function dgapper(
    dv::HDF5.Dataset,
    grid_bosons::kgrid
    )

    println("begins denisty gap equation")
    c_d    = Array{Complex{Float64},2}(undef,grid_bosons.N,grid_bosons.N)
    c_d    .= 0.0+0.0*im

    sites=grid_bosons.formfactorgrid

    R1  = SVector(3.0/2.0,sqrt(3.0)/2.0)
    R2  = SVector(3.0/2.0,-sqrt(3.0)/2.0)

    Np      = Int64(grid_bosons.N/6)
    vDbuff  = similar(dv[:,:,1])

    Midx        = findmin(norm.([grid_bosons.grid[i]-SVector(2*pi/3.0, 0.0) for i in 1:Int64((grid_bosons.N-0)/12)]))[2]

    vDbuff  .=0.0+0.0*im


    maxQ = findmax(abs.(dv[:,:,:]))[2][3]
    println("maxQidx: ",maxQ," | Midx: ",Midx)
    println("maxQ: " ,grid_bosons.grid[maxQ]," | M: ",grid_bosons.grid[Midx])
    println("DV[maxQ]:",findmax(abs.(dv[:,:,maxQ]))[1]," | DV[M]: ", findmax(abs.(dv[:,:,Midx]))[1])

    for n in 1:size(dv)[1]
       for m in 1:size(dv)[1]
           vDbuff[m,n] += dv[m,n,maxQ]
       end
    end

    println("vertex at q")

    for k in tqdm(1:grid_bosons.N)
        for n in 1:size(dv)[1]
            ffk = get_formfactor_explicit(grid_bosons.grid[k][1],grid_bosons.grid[k][2],sites[n],R1,R2)
            for p in 1:grid_bosons.N
                for m in 1:size(dv)[1]
                    ffp         = get_formfactor_explicit(grid_bosons.grid[p][1],grid_bosons.grid[p][2],sites[m],R1,R2)
                    c_d[p,k]   += vDbuff[m,n]*(conj(ffp)*ffk)
                end
            end
        end
    end
    println("direct particle hole channel gap calculation")
    eigensystem     = eigen(-c_d[:,:])
    vals            = eigensystem.values
    u               = eigensystem.vectors
    vals            = real.(vals)
    leadingvals     = real.(vals[end-3:end])
    leadingvecs     = real.(u[:,end-3:end])

    return leadingvals,leadingvecs, -c_d,vDbuff
end
