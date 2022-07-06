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
    w           ::vertices,
    fw          ::fouriervertices,
    Vh          ::Float64,
    Kh          ::Float64
    )

    #Set initial conditions#####################################################
    R2= SVector(3.0/2.0,sqrt(3.0)/2.0)
    R4= SVector(3.0/2.0,-sqrt(3.0)/2.0)
    R6= R2-R4


    println("Set initial conditions")
    v.p0[1,1,:] .=U+Vh/2+K/2
    v.c0[1,1,:] .=U+Vh/2+K/2
    v.d0[1,1,:] .=U+Vh/2+K/2
    w.p0[1,1,:] .=Vh-K
    w.c0[1,1,:] .=Vh-K
    w.d0[1,1,:] .=Vh-K

    v2_Arr  = [10,11,14,15,18,19]
    v3_Arr  = [8,9,12,13,16,17]

    for m in 2:7
        v.p0[m,m,:].+=(-J/8)
        v.c0[m,m,:].+=(-J/8)
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
            v.c0[1,1,qi]+=(-J/2)*get_formfactor_explicit(kx,ky,grid_bosons.formfactorgrid[m],grid_r.r1,grid_r.r2)
            v.d0[1,1,qi]+=(-J/8)*get_formfactor_explicit(kx,ky,grid_bosons.formfactorgrid[m],grid_r.r1,grid_r.r2)
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

    println("bandwidth:")
    println(findmax(energiesC)[1]-findmin(energiesC)[1])
    println("Initial scale is now:")
    println((findmax(energiesC)[1]-findmin(energiesC)[1])*1.05  )
    println(" ")
    ################################################################################

    #Buffer entities################################################################
    println("Construct buffers")

    Iter        = 0
    IterMax     = 400
    LambdaMin   = 10.0^-6
    vstop       = (findmax(energiesC)[1]-findmin(energiesC)[1])*3.0
    wstop       = vstop

    Lambda      = (findmax(energiesC)[1]-findmin(energiesC)[1])*1.05
    dLambda     = 0.05*Lambda

    incrementsv  = incrementvertices_initialization(bubbles.L,grid_bosons.N)
    incrementsw  = incrementvertices_initialization(bubbles.L,grid_bosons.N)

    LambdaArr   = []
    pmaxv        = []
    cmaxv        = []
    dmaxv        = []
    pmaxw        = []
    cmaxw        = []
    dmaxw        = []

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
        println(Iter)

        println("Set flow parameters")
        Lambda  = Lambda-dLambda
        println(Lambda)
        println(dLambda)

        println(" ")
        println("Calculate Increment")


        increment!(incrementsv,incrementsw,bubbles,Lambda,t,t2,t3,mu,v,fv,w,fw,grid_bosons,grid_r)

        v.P .+= incrementsv.P.*dLambda
        v.C .+= incrementsv.C.*dLambda
        v.D .+= incrementsv.D.*dLambda

        w.P .+= incrementsw.P.*dLambda
        w.C .+= incrementsw.C.*dLambda
        w.D .+= incrementsw.D.*dLambda

        Vsymmetrizer!(v,bubbles,grid_bosons)
        Vsymmetrizer!(w,bubbles,grid_bosons)


        vmaxp   = maximum(real.(v.P))
        vminp   = maximum(real.(-v.P))
        vmaxc   = maximum(real.(v.C))
        vminc   = maximum(real.(-v.C))
        vmaxd   = maximum(real.(v.D))
        vmind   = maximum(real.(-v.D))

        wmaxp   = maximum(real.(w.P))
        wminp   = maximum(real.(-w.P))
        wmaxc   = maximum(real.(w.C))
        wminc   = maximum(real.(-w.C))
        wmaxd   = maximum(real.(w.D))
        wmind   = maximum(real.(-w.D))


        println("Check for instabilities")
        println("vmax P:")
        println(vmaxp)
        println(findmax(real.(v.P))[2])
        println("vmin P:")
        println(vminp)
        println(findmax(real.(-v.P))[2])
        println("vmax C")
        println(vmaxc)
        println(findmax(real.(v.C))[2])
        println("vmin C:")
        println(vminc)
        println(findmax(real.(-v.C))[2])
        println("vmax D:")
        println(vmaxd)
        println(findmax(real.(v.D))[2])
        println("vmin D:")
        println(vmind)
        println(findmax(real.(-v.D))[2])

        println("wmax P:")
        println(wmaxp)
        println(findmax(real.(w.P))[2])
        println("wmin P:")
        println(wminp)
        println(findmax(real.(-w.P))[2])
        println("wmax C")
        println(wmaxc)
        println(findmax(real.(w.C))[2])
        println("wmin C:")
        println(wminc)
        println(findmax(real.(-w.C))[2])
        println("wmax D:")
        println(wmaxd)
        println(findmax(real.(w.D))[2])
        println("wmin D:")
        println(wmind)
        println(findmax(real.(-w.D))[2])

        vmaxp   = maximum(abs.(v.P))
        vmaxc   = maximum(abs.(v.C))
        vmaxd   = maximum(abs.(v.D))

        wmaxp   = maximum(abs.(w.P))
        wmaxc   = maximum(abs.(w.C))
        wmaxd   = maximum(abs.(w.D))

        if vmaxp>vstop
            Iter    = IterMax
            println("P INSTABILITY IN V!!!")
        elseif vmaxc>vstop
            Iter    = IterMax
            println("C INSTABILITY IN V!!!")
        elseif vmaxd>vstop
            println("D INSTABILITY IN V!!!")
            Iter    = IterMax
        elseif wmaxp>wstop
            Iter    = IterMax
            println("P INSTABILITY IN W!!!")
        elseif wmaxc>wstop
            Iter    = IterMax
            println("C INSTABILITY IN W!!!")
        elseif wmaxd>wstop
            println("D INSTABILITY IN W!!!")
            Iter    = IterMax
        elseif Lambda<LambdaMin
            Iter    = IterMax
        else
            Iter    += 1
        end

        push!(LambdaArr,Lambda)
        push!(pmaxv,vmaxp)
        push!(cmaxv,vmaxc)
        push!(dmaxv,vmaxd)
        push!(pmaxw,wmaxp)
        push!(cmaxw,wmaxc)
        push!(dmaxw,wmaxd)


        pL      = dLambda*findmax(abs.((incrementsv.P)))[1]
        cL      = dLambda*findmax(abs.((incrementsv.C)))[1]
        dL      = dLambda*findmax(abs.((incrementsv.D)))[1]
        dLambda = min(0.05*Lambda, dLambda/(2*pL),dLambda/(2*cL),dLambda/(2*dL) )


        BubblesGamma[1:2,Iter]  .= real(bubbles.ph[1,1,1]), real(bubbles.pp[1,1,1])
        BubblesM[1:2,Iter]      .= real(bubbles.ph[1,1,Midx]), real(bubbles.pp[1,1,Midx])

    end

    return LambdaArr,pmaxv,pmaxw,cmaxv,cmaxw,dmaxv,dmaxw, BubblesGamma, BubblesM,incrementsv, incrementsw
end


#gapper f
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
