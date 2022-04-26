@everywhere function s1(
    kx  ::Float64,
    ky  ::Float64
    )   ::Float64
    d1=@SArray [1.5,sqrt(3)/2]
    d2=@SArray [1.5,-sqrt(3)/2]
    d3=d1-d2
    k=SVector(kx,ky)

    return real(exp(im*sum(k.*d1))+exp(im*sum(k.*d2))+exp(im*sum(k.*d3))+exp(-im*sum(k.*d1))+exp(-im*sum(k.*d2))+exp(-im*sum(k.*d3)))
end

@everywhere function s2(
    kx  ::Float64,
    ky  ::Float64
    )   ::Float64
    d1=@SArray [1.5,sqrt(3)/2]
    d2=@SArray [1.5,-sqrt(3)/2]
    d3=d1-d2

    s1=d1+d2
    s2=d3+d1
    s3=d3-d2

    k=SVector(kx,ky)

    return real(exp(im*sum(k.*s1))+exp(im*sum(k.*s2))+exp(im*sum(k.*s3))+exp(-im*sum(k.*s1))+exp(-im*sum(k.*s2))+exp(-im*sum(k.*s3)))
end

@everywhere function s3(
    kx  ::Float64,
    ky  ::Float64
    )   ::Float64

    s1=2*@SArray [1.5,sqrt(3)/2]
    s2=2*@SArray [1.5,-sqrt(3)/2]
    s3=s1-s2

    k=SVector(kx,ky)

    return real(exp(im*sum(k.*s1))+exp(im*sum(k.*s2))+exp(im*sum(k.*s3))+exp(-im*sum(k.*s1))+exp(-im*sum(k.*s2))+exp(-im*sum(k.*s3)))
end

@everywhere function energies(
    kx  ::Float64,
    ky  ::Float64,
    t   ::Float64,
    t2  ::Float64,
    t3  ::Float64,
    mu  ::Float64
    )   ::Float64

    return -t*s1(kx,ky) -t2*s2(kx,ky)-t3*s3(kx,ky) -mu
end
