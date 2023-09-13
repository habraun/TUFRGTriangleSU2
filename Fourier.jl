@everywhere function get_factor!(
    grid_bosons ::kgrid,
    grid_r      ::rgrid,
    buffer      ::Array
    )

    sitesf=grid_bosons.formfactorgrid
    buffer.=0.0+0.0*im
    for ri in 1:length(sitesf)
        for ki in 1:grid_bosons.N
            f=get_formfactor_explicit(grid_bosons.grid[ki][1] ,grid_bosons.grid[ki][2] ,sitesf[ri],grid_r.r1,grid_r.r2)
            buffer[ki,ri]=conj(f) *grid_bosons.weights[ki]/sum(grid_bosons.weights)
        end
    end
end

@everywhere function vfourier!(
    fv              ::fouriervertices,
    v               ::vertices,
    grid_bosons     ::kgrid,
    grid_r          ::rgrid,
    buffer          ::Array)

    get_factor!(grid_bosons,grid_r,buffer)

    @tensor begin
        fv.P[n,m,ri]=v.P[n,m,ki]*buffer[ki,ri]
        fv.C[n,m,ri]=v.C[n,m,ki]*buffer[ki,ri]
        fv.D[n,m,ri]=v.D[n,m,ki]*buffer[ki,ri]
    end

end
