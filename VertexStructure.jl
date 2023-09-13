@everywhere struct vertices
    P::Array{Complex{Float64},3}
    C::Array{Complex{Float64},3}
    D::Array{Complex{Float64},3} #p -> VP and WP

    p0::Array{Complex{Float64},3}
    c0::Array{Complex{Float64},3}
    d0::Array{Complex{Float64},3}

    pc::Array{Complex{Float64},3}
    pd::Array{Complex{Float64},3}

    cp::Array{Complex{Float64},3}
    cd::Array{Complex{Float64},3}

    dp::Array{Complex{Float64},3}
    dc::Array{Complex{Float64},3} # c on d! from right to left!, see projection.jl

    L   ::Int64
    N   ::Int64


end

@everywhere mutable struct fouriervertices
    P::Array{Complex{Float64},3}
    C::Array{Complex{Float64},3}
    D::Array{Complex{Float64},3}

    L   ::Int64 #amount of ff
    N   ::Int64 # amount of momenta

    vec2idx::Dict
end

@everywhere function vertex_initialization(
        L   ::Int64,
        N   ::Int64
    )       ::vertices

    matrixP=Array{Complex{Float64},3}(undef,L,L,N)
    matrixP.=0.0+0.0*im
    matrixC=Array{Complex{Float64},3}(undef,L,L,N)
    matrixC.=0.0+0.0*im
    matrixD=Array{Complex{Float64},3}(undef,L,L,N)
    matrixD.=0.0+0.0*im

    matrixP0=Array{Complex{Float64},3}(undef,L,L,N)
    matrixP0.=0.0+0.0*im
    matrixC0=Array{Complex{Float64},3}(undef,L,L,N)
    matrixC0.=0.0+0.0*im
    matrixD0=Array{Complex{Float64},3}(undef,L,L,N)
    matrixD0.=0.0+0.0*im

    matrixPC=Array{Complex{Float64},3}(undef,L,L,N)
    matrixPC.=0.0+0.0*im
    matrixPD=Array{Complex{Float64},3}(undef,L,L,N)
    matrixPD.=0.0+0.0*im

    matrixCP=Array{Complex{Float64},3}(undef,L,L,N)
    matrixCP.=0.0+0.0*im
    matrixCD=Array{Complex{Float64},3}(undef,L,L,N)
    matrixCD.=0.0+0.0*im

    matrixDP=Array{Complex{Float64},3}(undef,L,L,N)
    matrixDP.=0.0+0.0*im
    matrixDC=Array{Complex{Float64},3}(undef,L,L,N)
    matrixDC.=0.0+0.0*im

    v=vertices(matrixP,matrixC,matrixD,matrixP0,matrixC0,matrixD0,matrixPC,matrixPD,matrixCP,matrixCD,matrixDP,matrixDC,L,N)

    return v
end

@everywhere function fouriervertex_initialization(
    L           ::Int64,
    grid_r      ::rgrid,
    shell       ::Int64,
    )           ::fouriervertices

    formfactorgrid=formfactor_grid(shell*2)
    D=Dict(zip(formfactorgrid, 1:length(formfactorgrid)))
    matrix1=Array{Complex{Float64},3}(undef,L,L,length(D))
    matrix1.=0.0+0.0*im

    matrix2=Array{Complex{Float64},3}(undef,L,L,length(D))
    matrix2.=0.0+0.0*im

    matrix3=Array{Complex{Float64},3}(undef,L,L,length(D))
    matrix3.=0.0+0.0*im

    v=fouriervertices(matrix1,matrix2,matrix3,L,grid_r.Nreal,D)
    return v
end


@everywhere struct incrementvertices
    P::Array{Complex{Float64},3}
    C::Array{Complex{Float64},3}
    D::Array{Complex{Float64},3}

    L   ::Int64
    N   ::Int64

end

@everywhere function incrementvertices_initialization(
    L   ::Int64,
    N   ::Int64
    )           ::incrementvertices

    matrixP=Array{Complex{Float64},3}(undef,L,L,N)
    matrixP.=0.0+0.0*im
    matrixC=Array{Complex{Float64},3}(undef,L,L,N)
    matrixC.=0.0+0.0*im
    matrixD=Array{Complex{Float64},3}(undef,L,L,N)
    matrixD.=0.0+0.0*im

    return incrementvertices(matrixP,matrixC,matrixD,L,N)
end
