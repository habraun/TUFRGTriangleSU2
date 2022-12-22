@everywhere function formfactor_grid(
        s::Int64
    )    :: Array{SArray{Tuple{2},Int64,1,2},1}

    N=3*(s-1)^2 + 9*(s-1) + 7

    sites=[SVector(0,0)]

    edge1=SVector(1,0)
    edge2=SVector(0,1)
    edge3=SVector(-1,1)
    T1=SVector(-1,1)
    T2=SVector(-1,0)
    T3=SVector(0,-1)

    for i in 1:s
        for j in 0:i-1
            vec=i*edge1 + T1*j
            push!(sites,vec)
            push!(sites,-vec)
        end
        for j in 0:i-1
            vec=i*edge2 + T2*j
            push!(sites,vec)
            push!(sites,-vec)
        end
        for j in 0:i-1
            vec=i*edge3 + T3*j
            push!(sites,vec)
            push!(sites,-vec)
        end
    end
    return sites
end


@everywhere function formfactor_grid_mutable(
        s::Int64
    )    ::Array{Array{Int64,1},1}

    N=3*(s-1)^2 + 9*(s-1) + 7

    sites=[[0,0]]

    edge1=[1,0]
    edge2=[0,1]
    edge3=[-1,1]
    T1=[-1,1]
    T2=[-1,0]
    T3=[0,-1]

    for i in 1:s
        for j in 0:i-1
            vec=i*edge1 + T1*j
            push!(sites,vec)
            push!(sites,-vec)
        end
        for j in 0:i-1
            vec=i*edge2 + T2*j
            push!(sites,vec)
            push!(sites,-vec)
        end
        for j in 0:i-1
            vec=i*edge3 + T3*j
            push!(sites,vec)
            push!(sites,-vec)
        end
    end
    return sites
end

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


@everywhere function fill_formfactors!(
    kx      ::Float64,
    ky      ::Float64,
    ff      ::Array{Complex{Float64}},
    sites   ::Array
    )

    R1= SVector(3.0/2.0,sqrt(3.0)/2.0)
    R2= SVector(3.0/2.0,-sqrt(3.0)/2.0)

    for i in 1 : size(ff)[1]
        ff[i] = get_formfactor_explicit(kx ,ky ,sites[i],R1,R2)
    end

end

####

@everywhere function mn_to_site(
    m       ::Int64,
    n       ::Int64,
    sites   ::Array
    )       ::Int64

    return  findall(x->x==(sites[m]-sites[n]), sites)[1]
end

@everywhere function generate_mn_to_site_arr(
    L       ::Int64,
    sites   ::Array
    )       ::Matrix{Int64}

    arr=Matrix{Int64}(undef, L, L)
    for m in 1:L
        for n in 1:L
            arr[m,n]=mn_to_site(m,n,sites)
        end
    end
    return arr
end

@everywhere function vec_to_site(
    vec     ::SArray,
    sites   ::Array
    )       ::Int64

    ret= findall(x->x==(vec), sites)
    if size(ret)[1]>0
        site=ret[1]
    else
        site=9999
    end
    return site
end


@everywhere function Mx_symmetrizer(
    sites   ::Array
    )       ::Vector{Int64}


    R1= SVector(3.0/2.0,sqrt(3.0)/2.0)
    R2= SVector(3.0/2.0,-sqrt(3.0)/2.0)

    #Real space lattice
    sites_real=[R1*sites[i][1]+R2*sites[i][2] for i in 1:length(sites)]

    #Buffer
    sites_symmetry=Array{Int64,1}(undef,length(sites_real))

    #Apply y->  -y mirror symmetry
    sites_real_Mx=[[1,-1].*sites_real[i] for i in 1:length(sites_real)]

    for j in 1:length(sites_real)
        v=sites_real[j]
        identifier=[norm(sites_real_Mx[i]-v) for i in 1:length(sites_real_Mx)]
        minval,minpos=findmin(identifier)

        if minval<0.001
            sites_symmetry[j]=Int64(minpos)
        else
            sites_symmetry[j]=Int64(9999)
        end
    end
    return sites_symmetry
end


@everywhere function Rot6_symmetrizer(
    sites   ::Array,
    L       ::Int64
    )       ::Vector{Int64}
    #sites=bubbles.formfactorgrid
    R1= SVector(3.0/2.0,sqrt(3.0)/2.0)
    R2= SVector(3.0/2.0,-sqrt(3.0)/2.0)

    #Real space lattice
    sites_real=[R1*sites[i][1]+R2*sites[i][2] for i in 1:length(sites)]

    #Buffer
    sites_symmetry=Array{Int64,1}(undef,length(sites_real))

    #Apply Rotation according to orbital combination

    for j in 1:length(sites_real)
        v=sites_real[j]

        identifier=[norm(sites_real[i]-(rot(2*pi/6)*v)) for i in 1:length(sites_real)]
        minval,minpos=findmin(identifier)
        if minval<0.001
            sites_symmetry[j]=Int64(minpos)
        else
            sites_symmetry[j]=Int64(9999)
        end


    end
    return sites_symmetry
end


@everywhere function symmetrizer(
    sites   ::Array,
    L       ::Int64
    )

    my=Mx_symmetrizer(sites)
    rottrafo=Rot6_symmetrizer(sites,L::Int64)


    return my,rottrafo
end
