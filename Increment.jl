@everywhere function fill_main_patch!(
	increments	:: incrementvertices,
	bubbles		:: bubble,
	VP			:: Array{Complex{Float64},3},
	VC			:: Array{Complex{Float64},3},
	VD			:: Array{Complex{Float64},3},
	bufferP		:: Array{Complex{Float64}, 2},
	bufferC		:: Array{Complex{Float64}, 2},
	bufferD 	:: Array{Complex{Float64}, 2}
	)

	xpp	= bubbles.pp
	xph	= bubbles.ph

	N_m	= Int64((increments.N)/12)

	for qi in 1:N_m
		xpp_qi	= xpp[:,:,qi]
		xph_qi	= xph[:,:,qi]

		VP_qi	= VP[:,:,qi]
		VC_qi	= VC[:,:,qi]
		VD_qi	= VD[:,:,qi]

		bufferP.= 0.0+0.0*im
		bufferC.= 0.0+0.0*im
		bufferD.= 0.0+0.0*im

		@tensor begin
			bufferP[L1,L4]=  VP_qi[L1,L2]*xpp_qi[L2,L3]*VP_qi[L3,L4]
			bufferC[L1,L4]=  VC_qi[L1,L2]*xph_qi[L2,L3]*VC_qi[L3,L4]
			bufferD[L1,L4]=((VD_qi[L1,L2]*xph_qi[L2,L3]*VD_qi[L3,L4])*(-2.0)
			        	+  1.0*(VC_qi[L1,L2]*xph_qi[L2,L3]*VD_qi[L3,L4])
				    	+  1.0*(VD_qi[L1,L2]*xph_qi[L2,L3]*VC_qi[L3,L4]))
		end

		increments.P[:,:,qi].=bufferP.*1.0
		increments.C[:,:,qi].=bufferC.*1.0
		increments.D[:,:,qi].=bufferD.*1.0

	end
end


@everywhere function increment!(
    increments  :: incrementvertices,
    bubbles     :: bubble,
    Lambda      :: Float64,
    t           :: Float64,
    t2          :: Float64,
    t3          :: Float64,
    mu          :: Float64,
    v           :: vertices,
    fv          :: fouriervertices,
    grid_bosons :: kgrid,
    grid_r      :: rgrid
    )

	L			= bubbles.L
	sites		= bubbles.formfactorgrid
	my,rottrafo	= symmetrizer(sites,L)


	println("Calculate Bubbles")
	faktor=Int64(ceil(log(10,10/Lambda)))

	if mu<=2*(t+t2-3*t3)
		phifaktor=1
	else
		phifaktor=3
	end

	println("Resolution:")
	println("Radial: "*string(3*(2^faktor)))
	println("Angular: "*string(phifaktor*120))
	fill_bubblesadaptive!(bubbles,grid_bosons,Lambda,t,t2,t3,mu,phifaktor*96,3*(2^faktor))


	println("Calculate Projections")
    projection!(v,fv,grid_bosons,grid_r)


    xpp	= bubbles.pp
    xph	= bubbles.ph
	VP	= v.p0+v.pc+v.pd+v.P
	VC	= v.c0+v.cp+v.cd+v.C
	VD	= v.d0+v.dc+v.dp+v.D

	bufferP	= similar(VP[:,:,1])
	bufferC	= similar(VC[:,:,1])
	bufferD	= similar(VD[:,:,1])

    increments.P.=0.0+0.0*im
    increments.C.=0.0+0.0*im
    increments.D.=0.0+0.0*im

	println("calculate main patch...")
	fill_main_patch!(increments,bubbles,VP,VC,VD,bufferP,bufferC,bufferD)

end


@everywhere function Vsymmetrizer!(
	v			::vertices,
	bubbles		::bubble,
	grid_bosons	::kgrid
	)

	L			=bubbles.L
	sites		=bubbles.formfactorgrid
	my,rottrafo	=symmetrizer(sites,L)

	N_m	= Int64((grid_bosons.N)/12)
 	N_p	= Int64((grid_bosons.N)/6)

	############################################################################
	#PHS on first main patch
	for qi in 1:N_m
		for L1 in 1:bubbles.L, L4 in L1+1:bubbles.L
			for L1 in 1:bubbles.L, L4 in L1+1:bubbles.L
				v.P[L1,L4,qi]	= conj(v.P[L4,L1,qi])
				v.C[L1,L4,qi]	= conj(v.C[L4,L1,qi])
				v.D[L1,L4,qi]	= conj(v.D[L4,L1,qi])
			end
		end
	end
	############################################################################
	#Mirror
	println("apply mirror symmetry...")

	for qi in 1:N_m
		for f2 in 1:L, f1 in 1:L
			f1t=my[f1]
			f2t=my[f2]
			if (f1t==9999)|(f2t==9999)
				#nothing
			else
				v.P[f1t,f2t,qi+N_m] = v.P[f1,f2,qi]
				v.C[f1t,f2t,qi+N_m] = v.C[f1,f2,qi]
				v.D[f1t,f2t,qi+N_m] = v.D[f1,f2,qi]
			end
		end
	end

	for qi in N_m+1:N_p
		for L1 in 1:bubbles.L, L4 in L1+1:bubbles.L
			for L1 in 1:bubbles.L, L4 in L1+1:bubbles.L
				v.P[L1,L4,qi]	= conj(v.P[L4,L1,qi])
				v.C[L1,L4,qi]	= conj(v.C[L4,L1,qi])
				v.D[L1,L4,qi]	= conj(v.D[L4,L1,qi])
			end
		end
	end
	##############################################################################
	#Rotation

	println("Apply rotational symmetry...")
	for i in 1:5
		for qi in 1:N_p
			for f1 in 1:bubbles.L, f2 in 1:bubbles.L
				f1t	= rottrafo[f1]
				f2t	= rottrafo[f2]
				if (f1t==9999)|(f2t==9999)
					#blub
				else
					v.P[f1t,f2t,qi+i*N_p]	= v.P[f1,f2,qi+(i-1)*N_p]
					v.C[f1t,f2t,qi+i*N_p]	= v.C[f1,f2,qi+(i-1)*N_p]
					v.D[f1t,f2t,qi+i*N_p]	= v.D[f1,f2,qi+(i-1)*N_p]
				end
			end
		end

		for qi in i*N_p+1:N_p*(i+1)
			for L1 in 1:bubbles.L, L4 in L1+1:bubbles.L
				v.P[L1,L4,qi]	= conj(v.P[L4,L1,qi])
				v.C[L1,L4,qi]	= conj(v.C[L4,L1,qi])
				v.D[L1,L4,qi]	= conj(v.D[L4,L1,qi])
			end
		end
    end
end
