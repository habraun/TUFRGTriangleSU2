@everywhere function ConPprojection!(
    v           ::vertices,
    fv          ::fouriervertices,
    grid_bosons ::kgrid,
    grid_r      ::rgrid
    )

    v.pc    .=0.0+0.0*im
    sites   =grid_bosons.formfactorgrid
    N_m     =Int64((grid_bosons.N)/12)

        for n in 1:v.L
            for m in 1:v.L
                for ridx in 1:v.L

                    m1=ridx
                    Rl2=sites[m]+sites[n]-sites[ridx]
                    m2=grid_bosons.form2idx[Rl2]

                    if m2<=v.L
                        Rarg=sites[n]-sites[ridx]

                        argidx=fv.vec2idx[Rarg]::Int64
                        Rexp=sites[n]-sites[ridx]

                        for qidx in 1:N_m+0


                            factor=get_formfactor_explicit(grid_bosons.grid[qidx][1],grid_bosons.grid[qidx][2],Rexp,grid_r.r1,grid_r.r2)

                            v.pc[m,n,qidx] += fv.C[m1,m2,argidx].*factor

                        end
                    end
                end
            end
        end
end



@everywhere function DonPprojection!(
    v           ::vertices,
    fv          ::fouriervertices,
    grid_bosons ::kgrid,
    grid_r      ::rgrid
    )           ::Nothing

    v.pd.=0.0+0.0*im
    sites=grid_bosons.formfactorgrid
    N_m=Int64((grid_bosons.N)/12)

        for n in 1:v.L
            for m in 1:v.L
                for ridx in 1:v.L

                    m1=ridx
                    Rl2=sites[m]-sites[n]-sites[ridx]
                    m2=grid_bosons.form2idx[Rl2]

                    if m2<=v.L
                        Rarg=-(sites[n]-sites[ridx])

                        argidx=fv.vec2idx[Rarg]::Int64
                        Rexp=-(sites[ridx])

                        for qidx in 1:N_m+0

                            factor=get_formfactor_explicit(grid_bosons.grid[qidx][1],grid_bosons.grid[qidx][2],Rexp,grid_r.r1,grid_r.r2)

                            v.pd[m,n,qidx]+= fv.D[m1,m2,argidx].*factor

                        end

                    end
                end
           end
        end

end




@everywhere function PonCprojection!(
    v           ::vertices,
    fv          ::fouriervertices,
    grid_bosons ::kgrid,
    grid_r      ::rgrid
    )           ::Nothing

    v.cp.=0.0+0.0*im
    sites=grid_bosons.formfactorgrid
    N_m=Int64((grid_bosons.N)/12)
        for n in 1:v.L
            for m in 1:v.L
                for ridx in 1:v.L


                    m1=ridx
                    Rl2=sites[m]+sites[n]-sites[ridx]
                    m2=grid_bosons.form2idx[Rl2]

                    if m2<=v.L
                        Rarg=sites[n]-sites[ridx]

                        argidx=fv.vec2idx[Rarg]::Int64
                        Rexp=sites[n]-sites[ridx]

                        for qidx in 1:N_m+0


                            factor=get_formfactor_explicit(grid_bosons.grid[qidx][1],grid_bosons.grid[qidx][2],Rexp,grid_r.r1,grid_r.r2)

                            v.cp[m,n,qidx]+= fv.P[m1,m2,argidx].*factor

                        end
                        #end
                    end
                end
           end
        end

end


@everywhere function DonCprojection!(
    v           ::vertices,
    fv          ::fouriervertices,
    grid_bosons ::kgrid,
    grid_r      ::rgrid
    )           ::Nothing

    v.cd.=0.0+0.0*im
    sites=grid_bosons.formfactorgrid
    N_m=Int64((grid_bosons.N)/12)
        for n in 1:v.L
            for m in 1:v.L
                for ridx in 1:v.L


                    m1=ridx
                    Rl2=-sites[m]+sites[n]+sites[ridx]
                    m2=grid_bosons.form2idx[Rl2]

                    if m2<=v.L
                        Rarg=-sites[m]#?
                            argidx=fv.vec2idx[Rarg]::Int64
                            Rexp=(-sites[ridx])


                            for qidx in 1:N_m+0

                                #factor=exp(im*Rexp_explicit'*grid_bosons.grid[qidx])
                                factor=get_formfactor_explicit(grid_bosons.grid[qidx][1],grid_bosons.grid[qidx][2],Rexp,grid_r.r1,grid_r.r2)

                                v.cd[m,n,qidx]+= fv.D[m1,m2,argidx].*factor

                        end
                        #end
                    end
                end
           end
        end

end


@everywhere function PonDprojection!(
    v           ::vertices,
    fv          ::fouriervertices,
    grid_bosons ::kgrid,
    grid_r      ::rgrid
    )           ::Nothing

    v.dp.=0.0+0.0*im
    sites=grid_bosons.formfactorgrid
    N_m=Int64((grid_bosons.N)/12)
        for n in 1:v.L
            for m in 1:v.L
                for ridx in 1:v.L


                    m1=ridx
                    Rl2=-(sites[m])-(sites[n])+(sites[ridx])
                    m2=grid_bosons.form2idx[Rl2]


                    if m2<=v.L
                        Rarg=-(sites[m])
                            argidx=fv.vec2idx[Rarg]::Int64
                            Rexp=(sites[n])-(sites[ridx])


                            for qidx in 1:N_m+0

                                factor=get_formfactor_explicit(grid_bosons.grid[qidx][1],grid_bosons.grid[qidx][2],Rexp,grid_r.r1,grid_r.r2)

                                v.dp[m,n,qidx]+= fv.P[m1,m2,argidx].*factor

                            end
                        #end
                    end
                end
            end
        end

end


@everywhere function ConDprojection!(
    v           ::vertices,
    fv          ::fouriervertices,
    grid_bosons ::kgrid,
    grid_r      ::rgrid
    )           ::Nothing

    v.dc.=0.0+0.0*im
    sites=grid_bosons.formfactorgrid
    N_m=Int64((grid_bosons.N)/12)
        for n in 1:v.L
            for m in 1:v.L
                for ridx in 1:v.L


                    m1=ridx
                    Rl2=-(sites[m])+(sites[n])+(sites[ridx])
                    m2=grid_bosons.form2idx[Rl2]


                    if m2<=v.L
                        Rarg=-(sites[m])

                            argidx=fv.vec2idx[Rarg]::Int64
                            Rexp=-(sites[ridx])


                            for qidx in 1:N_m+0

                                factor=get_formfactor_explicit(grid_bosons.grid[qidx][1],grid_bosons.grid[qidx][2],Rexp,grid_r.r1,grid_r.r2)


                                v.dc[m,n,qidx]+= fv.C[m1,m2,argidx].*factor

                            end
                        #end
                    end
                end
           end
        end

end


@everywhere function projection!(
    v           ::vertices,
    fv          ::fouriervertices,
    grid_bosons ::kgrid,
    grid_r      ::rgrid
    )           ::Nothing

    println("Fourier transform...")
    factor_buffer=zeros(grid_bosons.N, size(fv.P)[end])*im
    vfourier!(fv,v,grid_bosons,grid_r,factor_buffer)

    println("Project")
    ConPprojection!(v,fv,grid_bosons,grid_r)
    DonPprojection!(v,fv,grid_bosons,grid_r)
    PonCprojection!(v,fv,grid_bosons,grid_r)
    DonCprojection!(v,fv,grid_bosons,grid_r)
    ConDprojection!(v,fv,grid_bosons,grid_r)
    PonDprojection!(v,fv,grid_bosons,grid_r)

end
