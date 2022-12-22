#==
Structures for the bosonic momentum grid and the real space grid used for other functions
kgrid is a structure containing the kvectors on the positions of the 1st Bz (grid), the
amount of points (N) the reciprocal lattice vectors (b1,b2) and the weights of the points
(weights)

formfactorgrid and form2idx store the form factor plane waves and a dictionary
relating these to the indices.

rgrid is saves the information about the real space lattice, i.e. the two lattice vectors.
Note that the lattice vectors are hard-coded here and moreover are the lattice vectors which
one usually uses for graphene (which is also a triangular lattice with a two-atom basis)

==#


@everywhere struct kgrid
    grid            ::Array{SArray{Tuple{2},Float64,1,2},1}
    N               ::Int64
    b1              ::SArray{Tuple{2},Float64,1,2}
    b2              ::SArray{Tuple{2},Float64,1,2}
    weights         ::Array{Float64,1}
    formfactorgrid  ::Array{SArray{Tuple{2},Int64,1,2},1}
    form2idx        ::Dict{SArray{Tuple{2},Int64,1,2},Int64}
end

@everywhere function kgrid_initialization(
    N           ::Int64,
    shell       ::Int64,
    wpoints     ::Int64,
    Gammapoint  ::Bool,
    Mpoint      ::Bool,
    Kpoint      ::Bool
    )           ::kgrid

    gridpoints          = symhexagon(N,Gammapoint::Bool,Mpoint::Bool,Kpoint::Bool)
    gridpoints_weights  = weighter(gridpoints,wpoints)*8*pi*pi/(3*sqrt(3))

    formfactorgrid      = formfactor_grid(shell*2)
    formfactorgrid2     = formfactor_grid(shell*6)

    d                   = Dict(zip(formfactorgrid2, 1:length(formfactorgrid2)))
    grid_bosons         = kgrid(gridpoints, length(gridpoints),SVector(2*pi/3,2*pi*sqrt(3)/3),SVector(2*pi/3,-2*pi*sqrt(3)/3),gridpoints_weights,formfactorgrid,d)

    return grid_bosons
end


@everywhere struct rgrid
    Nreal   ::Int64
    r1      ::SArray{Tuple{2},Float64,1,2}
    r2      ::SArray{Tuple{2},Float64,1,2}
end

@everywhere function rgrid_initialization(
    Nreal   ::Int64
    )       ::rgrid
    r1      = @SVector [3/2,sqrt(3)/2]
    r2      = @SVector [3/2,-sqrt(3)/2]
    grid_r  = rgrid(Nreal,r1,r2)

    return grid_r
end
