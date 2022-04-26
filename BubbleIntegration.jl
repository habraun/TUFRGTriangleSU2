# Functions for searching zeros/fermi surface in bz

@everywhere function lineenergy(
    rho ::Float64,
    phi ::Float64,
    t   ::Float64,
    t2  ::Float64,
    t3  ::Float64,
    mu  ::Float64
    )   ::Float64

    kx=rho*cos(phi)
    ky=rho*sin(phi)

    return energies(kx,ky,t,t2,t3,mu)
end

@everywhere function fermi_rho_finder(
    phi ::Float64,
    t   ::Float64,
    t2  ::Float64,
    t3  ::Float64,
    mu  ::Float64
    )   ::Float64

    if (lineenergy(0.0, phi,t,t2,t3,mu)*lineenergy(rmaxfinder(phi), phi,t,t2,t3,mu))<=0.0
        rhomu=find_zero(x -> lineenergy(x,phi,t,t2,t3,mu), (0.0,rmaxfinder(phi)), Bisection())
    else
        rhomu=0.5*rmaxfinder(phi)
    end

    return rhomu
end


@everywhere function lineenergyDirac(
    j   ::Int64,
    rho ::Float64,
    phi ::Float64,
    t   ::Float64,
    t2  ::Float64,
    t3  ::Float64,
    mu  ::Float64
    )   ::Float64

    Kdirac=rot(j*pi/3)*SVector(2*pi/3,2*pi/(3*sqrt(3)))
    kx,ky=SVector(rho*cos(phi), rho*sin(phi)) + Kdirac
    return energies(kx,ky,t,t2,t3,mu)
end

@everywhere function fermi_rho_finderDirac(
    j   ::Int64,
    phi ::Float64,
    t   ::Float64,
    t2  ::Float64,
    t3  ::Float64,
    mu  ::Float64
    )   ::Float64

    if (lineenergyDirac(j,0.0, phi,t,t2,t3,mu)*lineenergyDirac(j, 2*pi/(3*sqrt(3)) , phi,t,t2,t3,mu))<=0.0
        rhomu=find_zero(x -> lineenergyDirac(j,x, phi,t,t2,t3,mu), (0.0,2*pi/(3*sqrt(3))), Bisection())
    else
        rhomu=rmaxfinder(phi)
    end

    return rhomu
end

###############################################################################
#functions for integrand kernel

@everywhere function lambda(
    e   ::Float64,
    T   ::Float64
    )   ::Float64

    return (e/(2.0*T^2))/(cosh(e/T)+1)
end

@everywhere function dlambda(
    e   ::Float64,
    T   ::Float64
    )   ::Float64

    inc= (T-e*sinh(e/T)+T*cosh(e/T))/((2*T^3)*(1+cosh(e/T))^2)
    if isnan(inc)==false
        inc=inc
    else
        inc=0.0
    end

    return inc
end

@everywhere function lambdabandsPP(
    value1  :: Float64,
    value3  :: Float64,
    Lambda  :: Float64
    )       ::Float64


    f1=lambda(value1,Lambda)
    f3=lambda(value3,Lambda)

    if abs(value1+value3)<=1e-10
        a=dlambda(-value3,Lambda)
    else
        a=(f1+f3)/(value1+value3)
    end

    return a
end


@everywhere function lambdabandsPH(
    value1  :: Float64,
    value2  :: Float64,
    Lambda  :: Float64
    )       ::Float64


    f1=lambda(value1,Lambda)
    f2=lambda(value2,Lambda)


    if abs(value1-value2)<=1e-10
        a=dlambda(-value2,Lambda)
    else
        a=(f1-f2)/(value1-value2)
    end

    return a
end

@everywhere function rot(
    phi ::Float64
    )   ::SArray{Tuple{2,2},Float64,2,4}

    m = @SMatrix [ cos(phi)  -sin(phi) ;
                    sin(phi)  cos(phi) ]

    return m
end

################################################################################
############ordinary integration routine, from gamma point######################
################################################################################

@everywhere function integrandPH!(
    grid_bosons :: kgrid,
    rho         :: Float64,
    phi         :: Float64,
    qb          :: SArray{Tuple{2},Float64,1,2},
    Lambda      :: Float64,
    t           :: Float64,
    t2          :: Float64,
    t3          :: Float64,
    mu          :: Float64,
    ff          :: Array{Complex{Float64}},
    bubbles     :: bubble,
    sites       :: Array,
    init        :: Array{Complex{Float64},1}
    )

    # init return value
    qxph = 0.0 + 0.0 * im
    #init.=0.0+0.0*im

    qx=rho*cos(phi)
    qy=rho*sin(phi)

    #compute folded vectors
    qx1=qx + qb[1]
    qy1=qy + qb[2]

    qx2=qx
    qy2=qy

    values1 = energies(qx1, qy1, t, t2, t3, mu)
    values2 = energies(qx2, qy2, t, t2, t3, mu)


    Mlambda=lambdabandsPH(values1,values2,Lambda)
    qxph += Mlambda

    # fill result array inplace
    fill_formfactors!(qx,qy, ff,sites)
    ff.*=(qxph*rho)

    for i in 1:length(sites)
        init[i]                    += ff[i]
        #for integration-check, keep commented
        #init[i]+=rho
    end
end


@everywhere function integrandPP!(
    grid_bosons :: kgrid,
    rho         :: Float64,
    phi         :: Float64,
    qb          :: SArray{Tuple{2},Float64,1,2},
    Lambda      :: Float64,
    t           :: Float64,
    t2          :: Float64,
    t3          :: Float64,
    mu          :: Float64,
    ff          :: Array{Complex{Float64}},
    bubbles     :: bubble,
    sites       :: Array,
    init        ::Array{Complex{Float64},1}
    )

    # init return value
    qxpp = 0.0 + 0.0 * im

    qx=rho*cos(phi)
    qy=rho*sin(phi)

    #compute folded vectors
    qx1=qx + qb[1]
    qy1=qy + qb[2]

    qx2=-qx
    qy2=-qy

    # compute dispersion
    values1 = energies(qx1, qy1, t, t2, t3, mu)
    values3 = energies(qx2, qy2, t, t2, t3, mu)


    Mlambda=lambdabandsPP(values1,values3,Lambda)
    qxpp -= Mlambda

    # fill result array inplace
    fill_formfactors!(qx,qy, ff,sites)
    ff.*=(qxpp*rho)

    for i in 1:length(sites)
        init[i]                    += ff[i]
        #for integration-check, keep commented
        #init[i]+=rho
    end
end

@everywhere function integralinitPH!(
    grid_bosons :: kgrid,
    a           :: Float64,
    b           :: Float64,
    phi         :: Float64,
    qb          :: SArray{Tuple{2},Float64,1,2},
    Lambda      :: Float64,
    t           :: Float64,
    t2          :: Float64,
    t3          :: Float64,
    mu          :: Float64,
    ff          :: Array{Complex{Float64}},
    bubbles     :: bubble,
    sites       :: Array,
    init        :: Array{Complex{Float64},1},
    res         :: Array{Complex{Float64},1}
    )

    init.=0.0
    res .=0.0

    integrandPH!(grid_bosons,a,phi,qb,Lambda,t,t2,t3,mu,ff,bubbles,sites,init)
    integrandPH!(grid_bosons,b,phi,qb,Lambda,t,t2,t3,mu,ff,bubbles,sites,init)

    res.+=(0.5*((b-a))).*init
    init.=0.0
end

@everywhere function integralinitPP!(
    grid_bosons :: kgrid,
    a           :: Float64,
    b           :: Float64,
    phi         :: Float64,
    qb          :: SArray{Tuple{2},Float64,1,2},
    Lambda      :: Float64,
    t           :: Float64,
    t2          :: Float64,
    t3          :: Float64,
    mu          :: Float64,
    ff          :: Array{Complex{Float64}},
    bubbles     :: bubble,
    sites       :: Array,
    init        :: Array{Complex{Float64},1},
    res         ::  Array{Complex{Float64},1}
    )

    init.=0.0
    res .=0.0

    integrandPP!(grid_bosons,a,phi,qb,Lambda,t,t2,t3,mu,ff,bubbles,sites,init)
    integrandPP!(grid_bosons,b,phi,qb,Lambda,t,t2,t3,mu,ff,bubbles,sites,init)

    res.+=(0.5*((b-a))).*init
    init.=0.0
end


@everywhere function nextintegralPH!(
    Nint        :: Int64,
    grid_bosons :: kgrid,
    bubbles     :: bubble,
    a           :: Float64,
    b           :: Float64,
    phi         :: Float64,
    qb          :: SArray{Tuple{2},Float64,1,2},
    Lambda      :: Float64,
    t           :: Float64,
    t2          :: Float64,
    t3          :: Float64,
    mu          :: Float64,
    ff          :: Array{Complex{Float64},1},
    sites       :: Array,
    init        :: Array{Complex{Float64},1},
    res         ::  Array{Complex{Float64},1}
    )

    res.*=0.5

    w=((b-a)/Nint)

    init.=0.0

    for i in 1:2:Nint-1
        integrandPH!(grid_bosons,a+i*w,phi,qb,Lambda,t,t2,t3,mu,ff,bubbles,sites,init)
    end
    init.*=w

    res.+=init
    init.=0.0
end

@everywhere function nextintegralPP!(
    Nint        :: Int64,
    grid_bosons :: kgrid,
    bubbles     :: bubble,
    a           :: Float64,
    b           :: Float64,
    phi         :: Float64,
    qb          :: SArray{Tuple{2},Float64,1,2},
    Lambda      :: Float64,
    t           :: Float64,
    t2          :: Float64,
    t3          :: Float64,
    mu          :: Float64,
    ff          :: Array{Complex{Float64}},
    sites       :: Array,
    init        :: Array{Complex{Float64},1},
    res         :: Array{Complex{Float64},1}
    )

    res.*=0.5

    w=((b-a)/Nint)

    init.=0.0

    for i in 1:2:Nint-1
        integrandPP!(grid_bosons,a+i*w,phi,qb,Lambda,t,t2,t3,mu,ff,bubbles,sites,init)
    end
    init.*=w

    res.+=init
    init.=0.0
end


@everywhere function LineIntegrationPH!(
    grid_bosons :: kgrid,
    bubbles     :: bubble,
    a           :: Float64,
    b           :: Float64,
    phi         :: Float64,
    kmax        :: Int64,
    res         :: Array{Complex{Float64},1},
    init        :: Array{Complex{Float64},1},
    buff1       :: Array{Complex{Float64},1},
    rtol        :: Float64,
    atol        :: Float64,
    qb          :: SArray{Tuple{2},Float64,1,2},
    Lambda      :: Float64,
    t           :: Float64,
    t2          :: Float64,
    t3          :: Float64,
    mu          :: Float64,
    ff          :: Array{Complex{Float64}},
    sites       :: Array
    )

    buff1  .=0.0
    res    .=0.0
    init   .=0.0
    k       =1
    breakk  =kmax

    integralinitPH!(grid_bosons,a,b,phi,qb,Lambda,t,t2,t3,mu,ff,bubbles,sites,init,res)
    integralinitPH!(grid_bosons,a,b,phi,qb,Lambda,t,t2,t3,mu,ff,bubbles,sites,init,buff1)

    while k <kmax

        nextintegralPH!(2^k,grid_bosons,bubbles,a,b,phi,qb,Lambda,t,t2,t3,mu,ff,sites,init,res)
        check  = isapprox(buff1,res; rtol = rtol)
        check2 = isapprox(buff1,res; atol = atol)

        if (all(check == true) || all(check2==true))
            res     .=buff1
            breakk  =k
            k       =kmax+1
        else
            buff1   .=res
            k       +=1
        end
    end
end


@everywhere function LineIntegrationPP!(
    grid_bosons :: kgrid,
    bubbles     :: bubble,
    a           :: Float64,
    b           :: Float64,
    phi         :: Float64,
    kmax        :: Int64,
    res         :: Array{Complex{Float64},1},
    init        :: Array{Complex{Float64},1},
    buff1       :: Array{Complex{Float64},1},
    rtol        :: Float64,
    atol        :: Float64,
    qb          :: SArray{Tuple{2},Float64,1,2},
    Lambda      :: Float64,
    t           :: Float64,
    t2          :: Float64,
    t3          :: Float64,
    mu          :: Float64,
    ff          :: Array{Complex{Float64}},
    sites       :: Array
    )

    buff1.=0.0
    res.=0.0
    init.=0.0
    k=1
    breakk=kmax

    integralinitPP!(grid_bosons,a,b,phi,qb,Lambda,t,t2,t3,mu,ff,bubbles,sites,init,res)
    integralinitPP!(grid_bosons,a,b,phi,qb,Lambda,t,t2,t3,mu,ff,bubbles,sites,init,buff1)

    while k <kmax

        nextintegralPP!(2^k,grid_bosons,bubbles,a,b,phi,qb,Lambda,t,t2,t3,mu,ff,sites,init,res)
        check  = isapprox(buff1,res; rtol = rtol)
        check2 = isapprox(buff1,res; atol = atol)

        if (all(check == true) || all(check2==true))
            res     .=buff1
            breakk  =k
            k       =kmax+1
        else
            buff1   .=res
            k       +=1
        end
    end
end


@everywhere function phiklapper(
    phi ::Float64
    )   ::Float64
	if (phi<=2*pi/12)&(phi>=-2*pi/12)
		phi = phi
	else
		while ((phi<=2*pi/12)&(phi>=-2*pi/12))==false
			phi= phi - 2*pi/6
		end
	end
	return phi
end

@everywhere function rmaxfinder(
    phi ::Float64
    )   ::Float64
	phi=phiklapper(phi)
	return (2*pi/3)/cos(phi)
end

################################################################################
############poacket integration routine, from K point###########################
################################################################################

@everywhere function integrandPHPocket!(
    grid_bosons :: kgrid,
    j           :: Int64,
    rho         :: Float64,
    phi         :: Float64,
    qb          :: SArray{Tuple{2},Float64,1,2},
    Lambda      :: Float64,
    t           :: Float64,
    t2          :: Float64,
    t3          :: Float64,
    mu          :: Float64,
    ff          :: Array{Complex{Float64}},
    bubbles     :: bubble,
    sites       :: Array,
    init        :: Array{Complex{Float64},1}
    )


    Kdirac=rot(j*pi/3)*SVector(2*pi/3,2*pi/(3*sqrt(3)))
    qx,qy=SVector(rho*cos(phi), rho*sin(phi)) + Kdirac
    # init return value
    qxph = 0.0 + 0.0 * im
    #init.=0.0+0.0*im

    #compute folded vectors
    qx1=qx + qb[1]
    qy1=qy + qb[2]

    qx2=qx
    qy2=qy

    values1 = energies(qx1, qy1, t, t2, t3, mu)
    values2 = energies(qx2, qy2, t, t2, t3, mu)


    Mlambda=lambdabandsPH(values1,values2,Lambda)

    qxph += Mlambda
    fill_formfactors!(qx, qy, ff,sites)
    ff.*=(qxph*rho)

    # fill result array inplace
    for i in 1:length(sites)
        init[i]                    += ff[i]
        #init[i]+=rho
    end
end


@everywhere function integrandPPPocket!(
    grid_bosons :: kgrid,
    j           :: Int64,
    rho         :: Float64,
    phi         :: Float64,
    qb          :: SArray{Tuple{2},Float64,1,2},
    Lambda      :: Float64,
    t           :: Float64,
    t2          :: Float64,
    t3          :: Float64,
    mu          :: Float64,
    ff          :: Array{Complex{Float64}},
    bubbles     :: bubble,
    sites       :: Array,
    init        :: Array{Complex{Float64},1}
    )


    Kdirac=rot(j*pi/3)*SVector(2*pi/3,2*pi/(3*sqrt(3)))
    qx,qy=SVector(rho*cos(phi), rho*sin(phi)) + Kdirac
    # init return value
    qxpp = 0.0 + 0.0 * im

    #compute folded vectors
    qx1=qx + qb[1]
    qy1=qy + qb[2]

    qx2=-qx
    qy2=-qy

    # compute dispersion
    values1 = energies(qx1, qy1, t, t2, t3, mu)
    values3 = energies(qx2, qy2, t, t2, t3, mu)

    Mlambda=lambdabandsPP(values1,values3,Lambda)
    qxpp -= Mlambda
    fill_formfactors!(qx, qy, ff,sites)
    ff.*=(qxpp*rho)


    for i in 1:length(sites)
        init[i]                    += ff[i]
        #init[i]+=rho
    end
end

@everywhere function integralinitPHPocket!(
    grid_bosons :: kgrid,
    j           :: Int64,
    a           :: Float64,
    b           :: Float64,
    phi         :: Float64,
    qb          :: SArray{Tuple{2},Float64,1,2},
    Lambda      :: Float64,
    t           :: Float64,
    t2          :: Float64,
    t3          :: Float64,
    mu          :: Float64,
    ff          :: Array{Complex{Float64}},
    bubbles     :: bubble,
    sites       :: Array,
    init        :: Array{Complex{Float64},1},
    res         :: Array{Complex{Float64},1}
    )

    init.=0.0
    res .=0.0

    integrandPHPocket!(grid_bosons,j,a,phi,qb,Lambda,t,t2,t3,mu,ff,bubbles,sites,init)
    integrandPHPocket!(grid_bosons,j,b,phi,qb,Lambda,t,t2,t3,mu,ff,bubbles,sites,init)

    res.+=(0.5*((b-a))).*init
    init.=0.0
end

@everywhere function integralinitPPPocket!(
    grid_bosons :: kgrid,
    j           :: Int64,
    a           :: Float64,
    b           :: Float64,
    phi         :: Float64,
    qb          :: SArray{Tuple{2},Float64,1,2},
    Lambda      :: Float64,
    t           :: Float64,
    t2          :: Float64,
    t3          :: Float64,
    mu          :: Float64,
    ff          :: Array{Complex{Float64}},
    bubbles     :: bubble,
    sites       :: Array,
    init        :: Array{Complex{Float64},1},
    res         ::  Array{Complex{Float64},1}
    )

    init.=0.0
    res .=0.0

    integrandPPPocket!(grid_bosons,j,a,phi,qb,Lambda,t,t2,t3,mu,ff,bubbles,sites,init)
    integrandPPPocket!(grid_bosons,j,b,phi,qb,Lambda,t,t2,t3,mu,ff,bubbles,sites,init)

    res.+=(0.5*((b-a))).*init
    init.=0.0
end


@everywhere function nextintegralPHPocket!(
    Nint        :: Int64,
    grid_bosons :: kgrid,
    bubbles     :: bubble,
    j           :: Int64,
    a           :: Float64,
    b           :: Float64,
    phi         :: Float64,
    qb          :: SArray{Tuple{2},Float64,1,2},
    Lambda      :: Float64,
    t           :: Float64,
    t2          :: Float64,
    t3          :: Float64,
    mu          :: Float64,
    ff          :: Array{Complex{Float64},1},
    sites       :: Array,
    init        :: Array{Complex{Float64},1},
    res         :: Array{Complex{Float64},1}
    )

    res.*=0.5

    w=((b-a)/Nint)

    init.=0.0

    for i in 1:2:Nint-1
        integrandPHPocket!(grid_bosons,j,a+i*w,phi,qb,Lambda,t,t2,t3,mu,ff,bubbles,sites,init)
    end
    init.*=w

    res.+=init
    init.=0.0
end

@everywhere function nextintegralPPPocket!(
    Nint        :: Int64,
    grid_bosons :: kgrid,
    bubbles     :: bubble,
    j           :: Int64,
    a           :: Float64,
    b           :: Float64,
    phi         :: Float64,
    qb          :: SArray{Tuple{2},Float64,1,2},
    Lambda      :: Float64,
    t           :: Float64,
    t2          :: Float64,
    t3          :: Float64,
    mu          :: Float64,
    ff          :: Array{Complex{Float64}},
    sites       :: Array,
    init        :: Array{Complex{Float64},1},
    res         :: Array{Complex{Float64},1}
    )

    res.*=0.5

    w=((b-a)/Nint)

    init.=0.0

    for i in 1:2:Nint-1
        integrandPPPocket!(grid_bosons,j,a+i*w,phi,qb,Lambda,t,t2,t3,mu,ff,bubbles,sites,init)
    end
    init.*=w

    res.+=init
    init.=0.0
end


@everywhere function LineIntegrationPHPocket!(
    grid_bosons :: kgrid,
    bubbles     :: bubble,
    j           :: Int64,
    a           :: Float64,
    b           :: Float64,
    phi         :: Float64,
    kmax        :: Int64,
    res         :: Array{Complex{Float64},1},
    init        :: Array{Complex{Float64},1},
    buff1       :: Array{Complex{Float64},1},
    rtol        :: Float64,
    atol        :: Float64,
    qb          :: SArray{Tuple{2},Float64,1,2},
    Lambda      :: Float64,
    t           :: Float64,
    t2          :: Float64,
    t3          :: Float64,
    mu          :: Float64,
    ff          :: Array{Complex{Float64}},
    sites       :: Array
    )

    buff1   .=0.0
    res     .=0.0
    init    .=0.0
    k       =1
    breakk  =kmax

    integralinitPHPocket!(grid_bosons,j,a,b,phi,qb,Lambda,t,t2,t3,mu,ff,bubbles,sites,init,res)
    integralinitPHPocket!(grid_bosons,j,a,b,phi,qb,Lambda,t,t2,t3,mu,ff,bubbles,sites,init,buff1)

    while k <kmax

        nextintegralPHPocket!(2^k,grid_bosons,bubbles,j,a,b,phi,qb,Lambda,t,t2,t3,mu,ff,sites,init,res)
        check   =isapprox(buff1,res; rtol = rtol)
        check2  =isapprox(buff1,res; atol = atol)

        if (all(check == true) || all(check2==true))
            res     .=buff1
            breakk  =k
            k       =kmax+1
        else
            buff1   .=res
            k       +=1
        end
    end
end


@everywhere function LineIntegrationPPPocket!(
    grid_bosons :: kgrid,
    bubbles     :: bubble,
    j           :: Int64,
    a           :: Float64,
    b           :: Float64,
    phi         :: Float64,
    kmax        :: Int64,
    res         :: Array{Complex{Float64},1},
    init        :: Array{Complex{Float64},1},
    buff1       :: Array{Complex{Float64},1},
    rtol        :: Float64,
    atol        :: Float64,
    qb          :: SArray{Tuple{2},Float64,1,2},
    Lambda      :: Float64,
    t           :: Float64,
    t2          :: Float64,
    t3          :: Float64,
    mu          :: Float64,
    ff          :: Array{Complex{Float64}},
    sites       :: Array
    )

    buff1   .=0.0
    res     .=0.0
    init    .=0.0
    k       =1
    breakk  =kmax

    integralinitPPPocket!(grid_bosons,j,a,b,phi,qb,Lambda,t,t2,t3,mu,ff,bubbles,sites,init,res)
    integralinitPPPocket!(grid_bosons,j,a,b,phi,qb,Lambda,t,t2,t3,mu,ff,bubbles,sites,init,buff1)

    while k <kmax

        nextintegralPPPocket!(2^k,grid_bosons,bubbles,j,a,b,phi,qb,Lambda,t,t2,t3,mu,ff,sites,init,res)
        check   =isapprox(buff1,res; rtol = rtol)
        check2  =isapprox(buff1,res; atol = atol)

        if (all(check == true) || all(check2==true))
            res     .=buff1
            breakk  =k
            k       =kmax+1
        else
            buff1   .=res
            k       +=1
        end
    end
end

########################Total integration##########################################



@everywhere function TotalIntegrationPH!(
    grid_bosons :: kgrid,
    bubbles     :: bubble,
    kmax        :: Int64,
    res         :: Array{Complex{Float64},1},
    init        :: Array{Complex{Float64},1},
    buff1       :: Array{Complex{Float64},1},
    rtol        :: Float64,
    atol        :: Float64,
    restot      :: Array,
    Nges        :: Int64,
    phiges      :: Int64,
    qb          :: SArray{Tuple{2},Float64,1,2},
    Lambda      :: Float64,
    t           :: Float64,
    t2          :: Float64,
    t3          :: Float64,
    mu          :: Float64,
    ff          :: Array{Complex{Float64}},
    sites       :: Array
    )


    if mu<=2*(t+t2-3*t3)
        restot  .=0.0+0.0*im
        deltaphi =2*pi/phiges

        for j in 0:phiges-1
            phi     =j*deltaphi
            rhomax  =rmaxfinder(phi)
            rhoF    =fermi_rho_finder(phi,t,t2,t3,mu)
            deltamax=rhoF/Nges

            for i in 0:Nges-1
                LineIntegrationPH!(grid_bosons,bubbles,i*deltamax,(i+1)*deltamax,phi, kmax,res,init,buff1,rtol,atol,qb,Lambda,t,t2,t3,mu,ff,sites)
                restot.+=(res)
            end

            deltamax=(rhomax-rhoF)/Nges
            for i in 0:Nges-1
                LineIntegrationPH!(grid_bosons,bubbles,rhoF+i*deltamax,rhoF+(i+1)*deltamax,phi, kmax,res,init,buff1,rtol,atol,qb,Lambda,t,t2,t3,mu,ff,sites)
                restot.+=(res)
            end

        end
        restot.*=deltaphi
        restot./=(8*pi^2/(3*sqrt(3)))
    else
        restot.=0.0+0.0*im

        #interior integration
        for j in 0:phiges-1
            deltaphi=2*pi/phiges
            phi     =j*deltaphi
            rhoF    =fermi_rho_finder(phi,t,t2,t3,2*(t+t2-3*t3))
            deltamax=rhoF/Nges

            for i in 0:Nges-1
                LineIntegrationPH!(grid_bosons,bubbles,i*deltamax,(i+1)*deltamax,phi, kmax,res,init,buff1,rtol,atol,qb,Lambda,t,t2,t3,mu,ff,sites)
                restot.+=(res.*deltaphi)
            end
        end

        #exterior integration
        offsets=[(pi/6)+2*pi/3, (pi/6)+pi, (pi/6)+2*2*pi/3, (2*2*pi/3)+pi/2, (pi/6), pi/2]

        for angle in 0:5
            for j in 0:Int64(phiges/6)
                deltaphi    =2*pi*(1/3)/(Int64(phiges/6))
                phi_weight  =2*pi*(1/3)/(Int64(phiges/6)+1)
                phi         =j*deltaphi + offsets[angle+1]

                rhomax      =fermi_rho_finderDirac(angle,phi,t,t2,t3,2*(t+t2-3*t3))
                rhoF        =fermi_rho_finderDirac(angle,phi,t,t2,t3,mu)
                deltamax    =rhoF/Nges

                for i in 0:Nges-1
                    LineIntegrationPHPocket!(grid_bosons,bubbles,angle,i*deltamax,(i+1)*deltamax,phi, kmax,res,init,buff1,rtol,atol,qb,Lambda,t,t2,t3,mu,ff,sites)
                    restot.+=(res.*phi_weight)
                end


                deltamax=(rhomax-rhoF)/Nges
                for i in 0:Nges-1
                    LineIntegrationPHPocket!(grid_bosons,bubbles,angle,rhoF+i*deltamax,rhoF+(i+1)*deltamax,phi, kmax,res,init,buff1,rtol,atol,qb,Lambda,t,t2,t3,mu,ff,sites)
                    restot.+=(res.*phi_weight)
                end
            end
        end

        restot./=(8*pi^2/(3*sqrt(3)))
    end
end

@everywhere function TotalIntegrationPP!(
    grid_bosons :: kgrid,
    bubbles     :: bubble,
    kmax        :: Int64,
    res         :: Array{Complex{Float64},1},
    init        :: Array{Complex{Float64},1},
    buff1       :: Array{Complex{Float64},1},
    rtol        :: Float64,
    atol        :: Float64,
    restot      :: Array,
    Nges        :: Int64,
    phiges      :: Int64,
    qb          :: SArray{Tuple{2},Float64,1,2},
    Lambda      :: Float64,
    t           :: Float64,
    t2          :: Float64,
    t3          :: Float64,
    mu          :: Float64,
    ff          :: Array{Complex{Float64}},
    sites       :: Array
    )


    if mu<=2*(t+t2-3*t3)
        restot  .=0.0+0.0*im
        deltaphi =2*pi/phiges

        for j in 0:phiges-1
            phi     =j*deltaphi
            rhomax  =rmaxfinder(phi)
            rhoF    =fermi_rho_finder(phi,t,t2,t3,mu)
            deltamax=rhoF/Nges

            for i in 0:Nges-1
                LineIntegrationPP!(grid_bosons,bubbles,i*deltamax,(i+1)*deltamax,phi, kmax,res,init,buff1,rtol,atol,qb,Lambda,t,t2,t3,mu,ff,sites)
                restot.+=(res)
            end

            deltamax=(rhomax-rhoF)/Nges
            for i in 0:Nges-1
                LineIntegrationPP!(grid_bosons,bubbles,rhoF+i*deltamax,rhoF+(i+1)*deltamax,phi, kmax,res,init,buff1,rtol,atol,qb,Lambda,t,t2,t3,mu,ff,sites)
                restot.+=(res)
            end

        end
        restot.*=deltaphi
        restot./=(8*pi^2/(3*sqrt(3)))
    else
        restot.=0.0+0.0*im

        #interior integration
        for j in 0:phiges-1
            deltaphi=2*pi/phiges
            phi     =j*deltaphi
            rhoF    =fermi_rho_finder(phi,t,t2,t3,2*(t+t2-3*t3))
            deltamax=rhoF/Nges

            for i in 0:Nges-1
                LineIntegrationPP!(grid_bosons,bubbles,i*deltamax,(i+1)*deltamax,phi, kmax,res,init,buff1,rtol,atol,qb,Lambda,t,t2,t3,mu,ff,sites)
                restot.+=(res.*deltaphi)
            end
        end

        #exterior integration
        offsets=[(pi/6)+2*pi/3, (pi/6)+pi, (pi/6)+2*2*pi/3, (2*2*pi/3)+pi/2, (pi/6), pi/2]

        for angle in 0:5
            for j in 0:Int64(phiges/6)
                deltaphi    =2*pi*(1/3)/(Int64(phiges/6))
                phi_weight  =2*pi*(1/3)/(Int64(phiges/6)+1)
                phi         =j*deltaphi + offsets[angle+1]

                rhomax      =fermi_rho_finderDirac(angle,phi,t,t2,t3,2*(t+t2-3*t3))
                rhoF        =fermi_rho_finderDirac(angle,phi,t,t2,t3,mu)
                deltamax    =rhoF/Nges

                for i in 0:Nges-1
                    LineIntegrationPPPocket!(grid_bosons,bubbles,angle,i*deltamax,(i+1)*deltamax,phi, kmax,res,init,buff1,rtol,atol,qb,Lambda,t,t2,t3,mu,ff,sites)
                    restot.+=(res.*phi_weight)
                end


                deltamax=(rhomax-rhoF)/Nges
                for i in 0:Nges-1
                    LineIntegrationPPPocket!(grid_bosons,bubbles,angle,rhoF+i*deltamax,rhoF+(i+1)*deltamax,phi, kmax,res,init,buff1,rtol,atol,qb,Lambda,t,t2,t3,mu,ff,sites)
                    restot.+=(res.*phi_weight)
                end
            end
        end

        restot./=(8*pi^2/(3*sqrt(3)))
    end
end
@everywhere function get_bubbles_qiadaptive(
    phiges      :: Int64,
    Nges        :: Int64,
    grid_bosons :: kgrid,
    qi          :: Int64,
    Lambda      :: Float64,
    t           :: Float64,
    t2          :: Float64,
    t3          :: Float64,
    mu          :: Float64,
    ff          :: Array{Complex{Float64},1},
    init        :: Array{Complex{Float64},1},
    sumph       :: Array{Complex{Float64},1},
    sumpp       :: Array{Complex{Float64},1},
    res         :: Array,
    buff1       :: Array,
    bubbles     :: bubble
    )

    sumph   .=0.0+0.0*im
    sumpp   .=0.0+0.0*im
    res     .=0.0+0.0*im
    init    .=0.0+0.0*im
    buff1   .=0.0+0.0*im
    ff      .=0.0+0.0*im

    qb = grid_bosons.grid[qi]

    L       = bubbles.L
    sites   = bubbles.formfactorgrid
    # allocate buffers
    qxph = zeros(Complex{Float64},L,L)
    qxpp = zeros(Complex{Float64},L,L)

    kmax    = 15
    atol    = 10.0^(-10)
    rtol    = 10.0^(-3)

    TotalIntegrationPH!(grid_bosons,bubbles,kmax,res,init,buff1,rtol,atol,sumph,Nges,phiges,qb,Lambda,t,t2,t3,mu,ff,sites)
    TotalIntegrationPP!(grid_bosons,bubbles,kmax,res,init,buff1,rtol,atol,sumpp,Nges,phiges,qb,Lambda,t,t2,t3,mu,ff,sites)

    for n in 1:L
        for m in 1:L
            idx=bubbles.mn_arr[m,n]
            qxph[m, n] = sumph[idx]
            qxpp[m, n] = sumpp[idx]
        end
    end

    return qxph, qxpp
end

@everywhere function fill_bubblesadaptive!(
    bubbles     :: bubble,
    grid_bosons :: kgrid,
    Lambda      :: Float64,
    t           :: Float64,
    t2          :: Float64,
    t3          :: Float64,
    mu          :: Float64,
    phiges      :: Int64,
    Nges        :: Int64
    )

    # allocate buffers
    sumph = zeros(Float64,length(bubbles.formfactorgrid))*im
    sumpp = zeros(Float64,length(bubbles.formfactorgrid))*im
    init  = zeros(Float64,length(bubbles.formfactorgrid))*im
    ff    = zeros(Float64,length(bubbles.formfactorgrid))*im
    res   = zeros(Float64,length(bubbles.formfactorgrid))*im
    buff1 = zeros(Float64,length(bubbles.formfactorgrid))*im
    # compute integrals in parallel
    num   = Int64((grid_bosons.N)/12)
    L     = bubbles.L

    bubblesres = pmap(qi -> get_bubbles_qiadaptive(phiges,Nges,grid_bosons,qi,Lambda,t,t2,t3,mu,ff,init,sumph,sumpp,res,buff1,bubbles), 1 : num, batch_size = ceil(Int64, num / nworkers()) )

    for qi in 1 : num
        bubbles.ph[:,:,qi] .= bubblesres[qi][1]
        bubbles.pp[:,:,qi] .= bubblesres[qi][2]
    end

end
