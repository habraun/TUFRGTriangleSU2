#== Differential equation solver for the actual flow equations. The routine  is
adaptive such that the steps scale logarithmic as the flow continues. In case of an
onsetting divergence/instability the solver adapts and shrinkes steps even more.

Right now the solver stops if one entry of the vertices surpasses 18.0 or if 350
steps have been calcul or if the scale 10.0^-4 is reached.
==#

@everywhere function start_flow(
    t           ::Float64,
    t2          ::Float64,
    t3          ::Float64,
    mu          ::Float64,
    U           ::Float64,
    V1          ::Float64,
    V2          ::Float64,
    V3          ::Float64,
    J           ::Float64,
    grid_bosons ::kgrid,
    bubbles     ::bubble,
    grid_r      ::rgrid,
    v           ::vertices,
    fv          ::fouriervertices,
    phifaktor      ::Int64
    )

    #Set initial conditions#####################################################
    R2= SVector(3.0/2.0,sqrt(3.0)/2.0)
    R4= SVector(3.0/2.0,-sqrt(3.0)/2.0)
    R6= R2-R4


    println("Set initial conditions")
    v.p0[1,1,:] .=U
    v.c0[1,1,:] .=U
    v.d0[1,1,:] .=U

    v2_Arr  = [10,11,14,15,18,19]
    v3_Arr  = [8,9,12,13,16,17]

    for m in 2:7
        v.p0[m,m,:].+=(-J/4)
        v.c0[m,m,:].+=(-J/4)
        v.d0[m,m,:].+=(-J/2)

        v.p0[m,m,:].+=V1
        v.c0[m,m,:].+=V1

        v.p0[v2_Arr[m-1],v2_Arr[m-1],:].+=V2
        v.c0[v2_Arr[m-1],v2_Arr[m-1],:].+=V2

        v.p0[v3_Arr[m-1],v3_Arr[m-1],:].+=V3
        v.c0[v3_Arr[m-1],v3_Arr[m-1],:].+=V3
    end

    #==
    #spinless
    for qi in 1:grid_bosons.N
        kx,ky=grid_bosons.grid[qi]
        v.p0[2,3,qi]+=(-V1/1)*get_formfactor_explicit(kx,ky,grid_bosons.formfactorgrid[3],grid_r.r1,grid_r.r2)
        v.p0[3,2,qi]+=(-V1/1)*get_formfactor_explicit(kx,ky,grid_bosons.formfactorgrid[2],grid_r.r1,grid_r.r2)
        v.p0[4,5,qi]+=(-V1/1)*get_formfactor_explicit(kx,ky,grid_bosons.formfactorgrid[5],grid_r.r1,grid_r.r2)
        v.p0[5,4,qi]+=(-V1/1)*get_formfactor_explicit(kx,ky,grid_bosons.formfactorgrid[4],grid_r.r1,grid_r.r2)
        v.p0[6,7,qi]+=(-V1/1)*get_formfactor_explicit(kx,ky,grid_bosons.formfactorgrid[5],grid_r.r1,grid_r.r2)
        v.p0[7,6,qi]+=(-V1/1)*get_formfactor_explicit(kx,ky,grid_bosons.formfactorgrid[4],grid_r.r1,grid_r.r2)
    end
    ==#

    for qi in 1:grid_bosons.N
        kx,ky=grid_bosons.grid[qi]
        v.p0[2,3,qi]+=-(J/2)*get_formfactor_explicit(kx,ky,grid_bosons.formfactorgrid[3],grid_r.r1,grid_r.r2)
        v.p0[3,2,qi]+=-(J/2)*get_formfactor_explicit(kx,ky,grid_bosons.formfactorgrid[2],grid_r.r1,grid_r.r2)
        v.p0[4,5,qi]+=-(J/2)*get_formfactor_explicit(kx,ky,grid_bosons.formfactorgrid[5],grid_r.r1,grid_r.r2)
        v.p0[5,4,qi]+=-(J/2)*get_formfactor_explicit(kx,ky,grid_bosons.formfactorgrid[4],grid_r.r1,grid_r.r2)
        v.p0[6,7,qi]+=-(J/2)*get_formfactor_explicit(kx,ky,grid_bosons.formfactorgrid[7],grid_r.r1,grid_r.r2)
        v.p0[7,6,qi]+=-(J/2)*get_formfactor_explicit(kx,ky,grid_bosons.formfactorgrid[6],grid_r.r1,grid_r.r2)
    end


    for qi in 1:grid_bosons.N
        kx,ky=grid_bosons.grid[qi]
        for m in 2:7
            #v.c0[1,1,qi]+=(-J/2)*get_formfactor_explicit(kx,ky,grid_bosons.formfactorgrid[m],grid_r.r1,grid_r.r2)
            #v.d0[1,1,qi]+=(-J/4)*get_formfactor_explicit(kx,ky,grid_bosons.formfactorgrid[m],grid_r.r1,grid_r.r2)
            #v.d0[1,1,qi]+=V1*get_formfactor_explicit(kx,ky,grid_bosons.formfactorgrid[m],grid_r.r1,grid_r.r2)

            v.d0[1,1,qi]+=V1*get_formfactor_explicit(kx,ky,grid_bosons.formfactorgrid[m],grid_r.r1,grid_r.r2)
            v.d0[1,1,qi]+=V2*get_formfactor_explicit(kx,ky,grid_bosons.formfactorgrid[v2_Arr[m-1]],grid_r.r1,grid_r.r2)
            v.d0[1,1,qi]+=V3*get_formfactor_explicit(kx,ky,grid_bosons.formfactorgrid[v3_Arr[m-1]],grid_r.r1,grid_r.r2)

            #spinlful
            #v.c0[1,1,qi]-=(V1/1)*get_formfactor_explicit(kx,ky,grid_bosons.formfactorgrid[m],grid_r.r1,grid_r.r2)
        end
    end


    ################################################################################

    #Find upper bandwidth###########################################################
    energiesC=[]
    for q in grid_bosons.grid
        e = energies(q[1],q[2],t,t2,t3,0.0)
        push!(energiesC,e)
    end

    println("bandwidth: ",findmax(energiesC)[1]-findmin(energiesC)[1])
    println("Initial scale is now :", (findmax(energiesC)[1]-findmin(energiesC)[1])*1.05)
    println(" ")
    ################################################################################

    #Buffer entities################################################################
    println("Construct buffers")

    Iter        = 0
    IterMax     = 400
    LambdaMin   = 10.0^-6
    vstop       = (findmax(energiesC)[1]-findmin(energiesC)[1])*3.0

    Lambda      = (findmax(energiesC)[1]-findmin(energiesC)[1])*1.05
    dLambda     = 0.05*Lambda

    increments  = incrementvertices_initialization(bubbles.L,grid_bosons.N)

    LambdaArr   = []
    pmax        = []
    cmax        = []
    dmax        = []

    inc         = []

    Midx        = findmin(norm.([grid_bosons.grid[i]-SVector(2*pi/3.0, 0.0) for i in 1:Int64((grid_bosons.N-0)/12)]))[2]
    Gammaidx    = 1

    BubblesGamma    = zeros(2,IterMax) #1111,1122,1221,1112,2111
    BubblesM        = zeros(2,IterMax)



    ################################################################################

    #Start flow#####################################################################
    println("Start flow")
    while (Iter<IterMax)&(Lambda>LambdaMin)

        println(" ")
        println(" ")
        println("Iteration: ",Iter)

        println("Set flow parameters")
        Lambda  = Lambda-dLambda
        println("Lambda: ",Lambda)
        println("dLambda: ",dLambda)

        println(" ")
        println("Calculate Increment")


        increment!(increments,bubbles,Lambda,t,t2,t3,mu,v,fv,grid_bosons,grid_r,phifaktor)

        v.P .+= increments.P.*dLambda
        v.C .+= increments.C.*dLambda
        v.D .+= increments.D.*dLambda

        Vsymmetrizer!(v,bubbles,grid_bosons)

        vmaxp   = maximum(real.(v.P))
        vminp   = maximum(real.(-v.P))
        vmaxc   = maximum(real.(v.C))
        vminc   = maximum(real.(-v.C))
        vmaxd   = maximum(real.(v.D))
        vmind   = maximum(real.(-v.D))


        println("Check for instabilities")
        println("vmax P: ", vmaxp, " at ", findmax(real.(v.P))[2])
        println("vmin P: ", vminp, " at ", findmin(real.(-v.P))[2])
        println("vmax C: ", vmaxc, " at ", findmax(real.(v.C))[2])
        println("vmin C: ", vminc, " at ", findmin(real.(-v.C))[2])
        println("vmax D: ", vmaxd, " at ", findmax(real.(v.D))[2])
        println("vmin D: ", vmind, " at ", findmin(real.(-v.D))[2])


        vmaxp   = maximum(abs.(v.P))
        vmaxc   = maximum(abs.(v.C))
        vmaxd   = maximum(abs.(v.D))

        if vmaxp>vstop
            Iter    = IterMax
            println("P INSTABILITY!!!")
        elseif vmaxc>vstop
            Iter    = IterMax
            println("C INSTABILITY!!!")
        elseif vmaxd>vstop
            println("D INSTABILITY!!!")
            Iter    = IterMax
        elseif Lambda<LambdaMin
            Iter    = IterMax
        else
            Iter    += 1
        end

        push!(LambdaArr,Lambda)
        push!(pmax,vmaxp)
        push!(cmax,vmaxc)
        push!(dmax,vmaxd)


        pL      = dLambda*findmax(abs.((increments.P)))[1]
        cL      = dLambda*findmax(abs.((increments.C)))[1]
        dL      = dLambda*findmax(abs.((increments.D)))[1]
        dLambda = min(0.05*Lambda, dLambda/(2*pL),dLambda/(2*cL),dLambda/(2*dL) )


        BubblesGamma[1:2,Iter]  .= real(bubbles.ph[1,1,1]), real(bubbles.pp[1,1,1])
        BubblesM[1:2,Iter]      .= real(bubbles.ph[1,1,Midx]), real(bubbles.pp[1,1,Midx])

    end

    return LambdaArr,pmax,cmax,dmax, BubblesGamma, BubblesM,increments
end

@everywhere function gapper(v::vertices,grid_bosons::kgrid)

    c_sc    = Array{Complex{Float64},2}(undef,grid_bosons.N,grid_bosons.N)
    c_sc    .= 0.0+0.0*im

    sites=grid_bosons.formfactorgrid

    R1  = SVector(3.0/2.0,sqrt(3.0)/2.0)
    R2  = SVector(3.0/2.0,-sqrt(3.0)/2.0)


    Np      = Int64(grid_bosons.N/12)
    vPbuff  = similar(v.P[:,:,1])

    vPbuff  .=0.0+0.0*im
    for m in 1:v.L
        for n in 1:v.L
            for i in 1:12
                vPbuff[m,n] += v.P[m,n,1+(i-1)*Np]./12
            end
        end
    end



    println("sc channel calculation")

    for k in 1:grid_bosons.N
        for p in 1:grid_bosons.N
            for n in 1:v.L
                for m in 1:v.L
                    ffk         = get_formfactor_explicit(grid_bosons.grid[k][1],grid_bosons.grid[k][2],sites[n],R1,R2)
                    ffp         = get_formfactor_explicit(grid_bosons.grid[p][1],grid_bosons.grid[p][2],sites[m],R1,R2)
                    c_sc[p,k]   += vPbuff[m,n].*(conj(ffp)*ffk)
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
