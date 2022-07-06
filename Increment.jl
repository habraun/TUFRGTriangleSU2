@everywhere function fill_main_patch!(
	incrementsv	:: incrementvertices,
	incrementsw :: incrementvertices,
	bubbles		:: bubble,
	VP			:: Array{Complex{Float64},3},
	VC			:: Array{Complex{Float64},3},
	VD			:: Array{Complex{Float64},3},
	WP			:: Array{Complex{Float64},3},
	WC			:: Array{Complex{Float64},3},
	WD			:: Array{Complex{Float64},3},
	bufferVP	:: Array{Complex{Float64}, 2},
	bufferVC	:: Array{Complex{Float64}, 2},
	bufferVD 	:: Array{Complex{Float64}, 2},
	bufferWP	:: Array{Complex{Float64}, 2},
	bufferWC	:: Array{Complex{Float64}, 2},
	bufferWD 	:: Array{Complex{Float64}, 2}
	)
	# SU(Nsymm)xSU(Msymm)
	Nsymm=2
	Msymm=2

	xpp	= bubbles.pp
	xph	= bubbles.ph

	N_m	= Int64((incrementsv.N)/12)

	for qi in 1:N_m
		xpp_qi	= xpp[:,:,qi]
		xph_qi	= xph[:,:,qi]

		VP_qi	= VP[:,:,qi]
		VC_qi	= VC[:,:,qi]
		VD_qi	= VD[:,:,qi]

		bufferVP.= 0.0+0.0*im
		bufferVC.= 0.0+0.0*im
		bufferVD.= 0.0+0.0*im

		WP_qi	= WP[:,:,qi]
		WC_qi	= WC[:,:,qi]
		WD_qi	= WD[:,:,qi]

		bufferWP.= 0.0+0.0*im
		bufferWC.= 0.0+0.0*im
		bufferWD.= 0.0+0.0*im

		@tensor begin
			bufferVP[L1,L4]=  VP_qi[L1,L2]*xpp_qi[L2,L3]*VP_qi[L3,L4]+WP_qi[L1,L2]*xpp_qi[L2,L3]*WP_qi[L3,L4]
			bufferVC[L1,L4]=  VC_qi[L1,L2]*xph_qi[L2,L3]*VC_qi[L3,L4]
			bufferVD[L1,L4]=((VD_qi[L1,L2]*xph_qi[L2,L3]*VD_qi[L3,L4])*(-Nsymm*Msymm)#4.0 is n*m 2 statt 4!!!
			        	+  (VC_qi[L1,L2]*xph_qi[L2,L3]*VD_qi[L3,L4])
				    	+  (VD_qi[L1,L2]*xph_qi[L2,L3]*VC_qi[L3,L4])
						+  (WC_qi[L1,L2]*xph_qi[L2,L3]*WD_qi[L3,L4])
				    	+  (WD_qi[L1,L2]*xph_qi[L2,L3]*WC_qi[L3,L4])
						+  (VD_qi[L1,L2]*xph_qi[L2,L3]*WD_qi[L3,L4])*(-Nsymm) #n
						+  (WD_qi[L1,L2]*xph_qi[L2,L3]*VD_qi[L3,L4])*(-Nsymm) #n
						+  (WC_qi[L1,L2]*xph_qi[L2,L3]*VD_qi[L3,L4])*(Msymm) #m
						+  (VD_qi[L1,L2]*xph_qi[L2,L3]*WC_qi[L3,L4])*(Msymm) #m
						) # add W equations only here it makes a difference

			bufferWP[L1,L4]=  VP_qi[L1,L2]*xpp_qi[L2,L3]*WP_qi[L3,L4]+WP_qi[L1,L2]*xpp_qi[L2,L3]*VP_qi[L3,L4]
			bufferWC[L1,L4]=  (VC_qi[L1,L2]*xph_qi[L2,L3]*WC_qi[L3,L4]
						+  WC_qi[L1,L2]*xph_qi[L2,L3]*VC_qi[L3,L4]
						+  WC_qi[L1,L2]*xph_qi[L2,L3]*WC_qi[L3,L4]*(-Msymm)) #m
			bufferWD[L1,L4]=(( WD_qi[L1,L2]*xph_qi[L2,L3]*WD_qi[L3,L4])*(-Nsymm) #n
						+  (WD_qi[L1,L2]*xph_qi[L2,L3]*VC_qi[L3,L4])
						+  (VC_qi[L1,L2]*xph_qi[L2,L3]*WD_qi[L3,L4])
						) # add W equations only here it makes a difference
		end

		incrementsv.P[:,:,qi].=bufferVP.*1.0
		incrementsv.C[:,:,qi].=bufferVC.*1.0
		incrementsv.D[:,:,qi].=bufferVD.*1.0
		#use here 0.0 if you want to neglect the "w-vertices"
		incrementsw.P[:,:,qi].=bufferWP.*1.0
		incrementsw.C[:,:,qi].=bufferWC.*1.0
		incrementsw.D[:,:,qi].=bufferWD.*1.0

	end
end


@everywhere function increment!(
    incrementsv  :: incrementvertices,
	incrementsw  :: incrementvertices,
    bubbles     :: bubble,
    Lambda      :: Float64,
    t           :: Float64,
    t2          :: Float64,
    t3          :: Float64,
    mu          :: Float64,
    v           :: vertices,
    fv          :: fouriervertices,
	w           :: vertices,
    fw          :: fouriervertices,
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
	projection!(w,fw,grid_bosons,grid_r)


    xpp	= bubbles.pp
    xph	= bubbles.ph

	VP	= v.p0+v.pc+v.pd+v.P
	VC	= v.c0+v.cp+v.cd+v.C
	VD	= v.d0+v.dc+v.dp+v.D

	WP	= (w.p0+w.pc+w.pd+w.P)
	WC	= (w.c0+w.cp+w.cd+w.C)
	WD	= (w.d0+w.dc+w.dp+w.D)

	bufferVP	= similar(VP[:,:,1])
	bufferVC	= similar(VC[:,:,1])
	bufferVD	= similar(VD[:,:,1])

	bufferWP	= similar(WP[:,:,1])
	bufferWC	= similar(WC[:,:,1])
	bufferWD	= similar(WD[:,:,1])

    incrementsv.P.=0.0+0.0*im
    incrementsv.C.=0.0+0.0*im
    incrementsv.D.=0.0+0.0*im
	incrementsw.P.=0.0+0.0*im
    incrementsw.C.=0.0+0.0*im
    incrementsw.D.=0.0+0.0*im

	println("calculate main patch...")
	fill_main_patch!(incrementsv, incrementsw,bubbles,VP,VC,VD,WP,WC,WD,bufferVP,bufferVC,bufferVD,bufferWP,bufferWC,bufferWD)

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
