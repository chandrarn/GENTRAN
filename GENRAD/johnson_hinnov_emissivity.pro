FUNCTION johnson_hinnov_emissivity,ene,te,n0,debug=debug,quiet=quiet
;print,' '
;print,'the call to this procedure is johnson_hinnov_h_emission,Ne (cm-3),Te (in eV) ,n0 (cm-3), Ha_emiss'
;print,' where the H_alpha emissivity is returned as the variable ha_emiss (ph/sec/ster/cm^3)'
;print,' '
;print,'if you want emission or populations other than H_alpha, then, after completion, type'
;print,' '
;print,'common h0_pred, balmer_int,lyman_recomb,balmer_recomb,tot_rec,tot_ion,pop_dens_exc,pop_dens_rec'
;print,'common einstein_coefs,a_lyman,a_balmer'
;print,' '
;print,'the variable POP_DENS_EXC(0:4) holds the population densities (in cm^-3)
;print,'resulting from excitation from the ground state'
;print,'for the levels with  n > or = 2 (in cm-3) 
;print,'(where POP_DENS_EXC(0) gives the n=2 density, 1 gives the n=3 density, etc)'
;print,' '
;print,'the population densities for the levels with  n > or = 2 (in cm-3) 
;print,'which are the result of recombination are given as the variable:'
;print,' '
;print,'POP_DENS_REC(0:13) (where 0 gives the n=2 density, etc)'
;print,' '
;print,'To get emission rates (photons/cm-3/s)'
;print,'(e.g. for n=4 to 2) construct E=(pop_dens_exc(2)+pop_dens_rec(2))*A_balmer(1)'
;print,'or to get, say, Lyman beta, take E=(pop_dens_exc(1)+pop_dens_rec(1))*A_Lyman(1)'
;print,' '
;print,'Thus (pop_dens_exc(n)+pop_dens_rec(n))*A_Lyman(n) gives the emissivity of the n=n to 1 Lyman line,'
;print,'while (pop_dens_exc(n)+pop_dens_rec(n))*A_Balmer(n-1) gives the emissivity of the n=n to 2 Balmer line.'
;print,' '
;print,'*************************************'

;C-PROGRAM TO USE THE JOHNSON & HINNOV RESULTS TO CALCULATE THE
;C-population DENSITIES AND VOL EMISSION RATE AS A FUNCTION OF T(EV),NE, AND N1
;C-THE CONSTANTS R1 AND R0 ARE STORED IN A 4-DIMENSIONAL ARRAY. THE 
;C-FIRST DIMENSION IS DENSITY, THE SECOND TEMPERATURE, THE THIRD
;;C-SPECIFIES EITHER R0 (P), OR R1(P), AND THE FOURTH SPECIFIES
;C-ENERGY LEVEL,STARTING WITH N=2.
;C
; 	COMMON C,R,NP,DENS,TEMP
; R has indices from R(6,11,0:1,2:6)
	R=fltarr(7,11,2,5)
	DENS=[1.E10,1.E11,1.E12,1.E13,1.E14,1.E15,1.e16]
	TEMP=[.345,.69,1.38,2.76,5.52,11.0,22.1,44.1,88.0,176.5,706.]
        s=fltarr(n_elements(dens),n_elements(temp))
        al=fltarr(n_elements(dens),n_elements(temp))
;R02
	R(*,*,0,0)=[7.6E-6,1.1E-5,1.9E-5,4.9E-5,2.4E-4,2.2E-3,1.8e-2,$
     	1.5E-3,1.8E-3,2.5E-3,4.5E-3,1.3E-2,7.1E-2,3.7e-1,$
     	2.6E-2,2.9E-2,3.5E-2,4.9E-2,9.6E-2,3.2E-1,7.8e-1,$
	1.3E-1,1.4E-1,1.5E-1,1.9E-1,2.8E-1,6.1E-1,9.2e-1,$
	3.6E-1,3.7E-1,3.8E-1,4.2E-1,5.2E-1,8.0E-1,9.6e-1,$
	        6.9E-1,6.9E-1,7.0E-1,7.3E-1,7.9E-1,9.2E-1,9.8e-1,$
	1.1,1.1,1.1,1.1,1.1,1.0,1.0,$
      	1.5,1.5,1.5,1.5,1.4,1.1,1.0,$
	2.0,2.0,1.9,1.9,1.7,1.3,1.0,$
	2.4,2.4,2.4,2.3,2.1,1.4,1.1,$
        3.4,3.4,3.3,3.2,2.9,2.0,1.2]
;R12
        R(*,*,1,0)=[2.5E-7,2.5E-6,2.5E-5,2.5E-4,2.5E-3,2.4E-2,2.0e-1,$
     	1.9E-7,1.9E-6,1.9E-5,1.9E-4,1.9E-3,1.8E-2,1.0e-1,$
     	1.6E-7,1.6E-6,1.6E-5,1.6E-4,1.5E-3,1.1E-2,3.2e-2,$
     	1.5E-7,1.5E-6,1.5E-5,1.5E-4,1.3E-3,7.2E-3,1.3e-2,$
        1.6E-7,1.6E-6,1.6E-5,1.5E-4,1.3E-3,5.4E-3,8.0e-3,$
	1.8E-7,1.8E-6,1.8E-5,1.7E-4,1.4E-3,5.1E-3,7.0e-3,$
     	2.1E-7,2.1E-6,2.1E-5,2.0E-4,1.6E-3,5.6E-3,7.5e-3,$
	2.3E-7,2.3E-6,2.3E-5,2.2E-4,1.7E-3,6.3E-3,8.7e-3,$
     	2.3E-7,2.3E-6,2.3E-5,2.2E-4,1.8E-3,7.0E-3,1.0e-2,$
	2.2e-7,2.2e-6,2.1e-5,2.1e-4,1.7e-3,7.4e-3,1.1e-2,$
        1.6E-7,1.6E-6,1.6E-5,1.6E-4,1.4E-3,7.2E-3,1.3e-2]
;R03
        R(*,*,0,1)=[2.2E-3,3.1E-3,6.0E-3,2.2E-2,1.3E-1,3.5E-1,4.2e-1,$
     	2.6E-2,3.3E-2,5.0E-2,1.2E-1,4.3E-1,7.2E-1,8.5e-1,$
     	1.1E-1,1.3E-1,1.6E-1,3.0E-1,6.8E-1,8.9E-1,9.7e-1,$
	2.7E-1,2.9E-1,3.4E-1,5.0E-1,8.2E-1,9.5E-1,9.9e-1,$
     	4.8E-1,5.0E-1,5.4E-1,6.8E-1,9.0E-1,9.8E-1,1.0,$
        7.3E-1,7.4E-1,7.7E-1,8.5E-1,9.5E-1,9.9E-1,1.0,$
	1.0,1.0,1.0,1.0,1.0,1.0,1.0,$
     	1.3,1.3,1.3,1.2,1.1,1.0,1.0,$
	1.6,1.6,1.5,1.4,1.1,1.0,1.0,$
        1.9,1.9,1.8,1.6,1.2,1.1,1.0,$
        2.5,2.4,2.4,2.1,1.5,1.1,1.0]
;R13
        R(*,*,1,1)=[1.0E-7,1.0E-6,1.0E-5,1.0E-4,1.1E-3,1.3E-2,1.1e-1,$
     	8.2E-8,8.1E-7,8.0E-6,7.7E-5,6.1E-4,4.5E-3,4.2e-2,$
	7.1E-8,7.0E-7,6.8E-6,5.9E-5,3.3E-4,1.6E-3,4.1e-3,$
	6.8E-8,6.7E-7,6.3E-6,4.9E-5,2.1E-4,7.3E-4,1.2e-3,$
     	7.2E-8,7.0E-7,6.5E-6,4.7E-5,1.8E-4,4.9E-4,6.8e-4,$
     	8.1E-8,7.8E-7,7.2E-6,5.1E-5,1.9E-4,4.5E-4,5.8e-4,$
     	9.1E-8,8.9E-7,8.2E-6,5.8E-5,2.1E-4,5.0E-4,6.4e-4,$
	9.7E-8,9.5E-7,8.8E-6,6.5E-5,2.5E-4,6.0E-4,7.6e-4,$
     	9.7E-8,9.4E-7,8.8E-6,6.7E-5,2.7E-4,6.9E-4,9.1e-4,$
	8.9e-8,8.7e-7,8.2e-6,6.5e-5,2.8e-4,7.6e-4,1.0e-3,$
        6.5E-8,6.4E-7,6.1E-6,5.2E-5,2.7E-4,8.0E-4,1.2e-3]
;R04
        R(*,*,0,2)=[1.8E-2,2.8E-2,7.3E-2,3.1E-1,6.0E-1,7.4E-1,7.7e-1,$
     	8.2E-2,1.1E-1,2.2E-1,5.6E-1,8.3E-1,9.3E-1,9.6e-1,$
	2.0E-1,2.4E-1,3.9E-1,7.4E-1,9.2E-1,.98,.99,$
	3.7E-1,4.1E-1,5.7E-1,8.4E-1,9.6E-1,.99,1.0,$
     	5.6E-1,6.0E-1,7.2E-1,9.0E-1,9.8E-1,1.0,1.0,$
     	7.7E-1,7.9E-1,8.5E-1,9.5E-1,9.9E-1,1.0,1.0,$
	9.9E-1,9.9E-1,9.9E-1,1.0,1.0,1.0,1.0,$
        1.2,1.2,1.1,1.1,1.0,1.0,1.0,$
	1.4,1.4,1.3,1.1,1.0,1.0,1.0,$
	1.7,1.6,1.5,1.2,1.1,1.0,1.0,$
	2.1,2.1,1.9,1.5,1.1,1.0,1.0]
;R14
        R(*,*,1,2)=[7.2E-8,7.1E-7,6.9E-6,5.7E-5,4.8E-4,5.3E-3,4.5e-2,$
     	5.9E-8,5.7E-7,5.1E-6,3.1E-5,1.7E-4,1.1E-3,5.9e-3,$
     	5.1E-8,4.9E-7,4.0E-6,1.9E-5,7.1E-5,3.0E-4,7.8e-4,$
	4.8E-8,4.5E-7,3.4E-6,1.4E-5,4.2E-5,1.3E-4,2.1e-4,$
     	5.0E-8,4.7E-7,3.4E-6,1.3E-5,3.5E-5,8.6E-5,1.2e-4,$
     	5.6E-8,5.2E-7,3.7E-6,1.4E-5,3.6E-5,8.1E-5,1.0e-4,$
     	6.3E-8,5.9E-7,4.3E-6,1.6E-5,4.3E-5,9.3E-5,1.2e-4,$
	6.7E-8,6.3E-7,4.8E-6,1.9E-5,5.2E-5,1.1E-4,1.4e-4,$
     	6.6E-8,6.3E-7,4.9E-6,2.1E-5,5.9E-5,1.3E-4,1.7e-4,$
	6.1e-8,5.8e-7,4.7e-6,2.2e-5,6.3e-5,1.4e-4,2.0e-4,$
        4.4E-8,4.2E-7,3.7E-6,2.0E-5,6.4E-5,1.6E-4,2.5e-4]
;R05
        R(*,*,0,3)=[5.5E-2,1.E-1,3.3E-1,6.8E-1,8.5E-1,.9,.92,$
     	1.5E-1,2.4E-1,5.5E-1,8.4E-1,9.5E-1,.98,.99,$
     	2.9E-1,4.0E-1,7.0E-1,9.1E-1,9.8E-1,.99,1.0,$
	4.5E-1,5.5E-1,8.E-1,9.5E-1,9.9E-1,1.0,1.0,$
     	6.2E-1,7.0E-1,8.7E-1,9.7E-1,9.9E-1,1.0,1.0,$
     	8.0E-1,8.4E-1,9.3E-1,9.8E-1,1.0,1.0,1.0,$
	9.8E-1,9.8E-1,9.9E-1,1.0,1.0,1.0,1.0,$
        1.2,1.1,1.1,1.0,1.0,1.0,1.0,$
	1.4,1.3,1.2,1.0,1.0,1.0,1.0,$
	1.5,1.5,1.3,1.1,1.0,1.0,1.0,$
	1.9,1.9,1.6,1.2,1.0,1.0,1.0]
;R15
        R(*,*,1,3)=[6.0E-8,5.7E-7,4.4E-6,2.5E-5,1.8E-4,2.0E-3,1.6e-2,$
     	4.8E-8,4.4E-7,2.7E-6,1.1E-5,5.0E-5,3.2E-4,1.7e-3,$
     	4.1E-8,3.5E-7,1.8E-6,5.9E-6,1.9E-5,8.1E-5,2.0e-4,$
	3.8E-8,3.2E-7,1.5E-6,4.2E-6,1.1E-5,3.4E-5,5.5e-5,$
     	4.0E-8,3.2E-7,1.4E-6,3.8E-6,9.4E-6,2.3E-5,3.1e-5,$
     	4.4E-8,3.6E-7,1.6E-6,4.3E-6,1.0E-5,2.2E-5,2.8e-5,$
     	5.0E-8,4.1E-7,1.9E-6,5.2E-6,1.2E-5,2.5E-5,3.2e-5,$
	5.3E-8,4.5E-7,2.2E-6,6.2E-6,1.5E-5,3.1E-5,3.9e-5,$
        5.2E-8,4.5E-7,2.4E-6,7.1E-6,1.7E-5,3.7E-5,4.8e-5,$
	4.8e-8,4.3e-7,2.4e-6,7.6e-6,1.9e-5,4.3e-5,5.7e-5,$
        3.5E-8,3.2E-7,2.1E-6,7.5E-6,2.0E-5,4.8E-5,7.2e-5]
;R06
        R(*,*,0,4)=[1.1E-1,2.7E-1,6.4E-1,8.6E-1,9.4E-1,.96,.97,$
     	2.4E-1,4.5E-1,7.9E-1,9.4E-1,9.8E-1,.99,1.0,$
     	3.8E-1,6.E-1,8.7E-1,9.7E-1,9.9E-1,1.0,1.0,$
	5.3E-1,7.2E-1,9.1E-1,9.8E-1,1.0,1.0,1.0,$
     	6.8E-1,8.1E-1,9.4E-1,9.9E-1,1.0,1.0,1.0,$
     	8.2E-1,9.0E-1,9.7E-1,9.9E-1,1.0,1.0,1.0,$
	9.7E-1,9.9E-1,1.0,1.0,1.0,1.0,1.0,$
     	1.1,1.1,1.0,1.0,1.0,1.0,1.0,$
	1.3,1.2,1.1,1.0,1.0,1.0,1.0,$
	1.5,1.3,1.1,1.0,1.0,1.0,1.0,$
	1.8,1.7,1.3,1.1,1.0,1.0,1.0]
;R16
       	R(*,*,1,4)=[5.2E-8,4.3E-7,2.3E-6,1.0E-5,7.2E-5,7.7E-4,6.5e-3,$
       	4.1E-8,3.0E-7,1.2E-6,4.0E-6,1.7E-5,1.1E-4,5.9e-4,$
       	3.4E-8,2.2E-7,7.7E-7,2.1E-6,6.6E-6,2.7E-5,6.8e-5,$
	3.1E-8,1.9E-7,6.0E-7,1.5E-6,3.8E-6,1.1E-5,1.8e-5,$
     	3.2E-8,1.9E-7,5.9E-7,1.4E-6,3.2E-6,7.7E-6,1.0e-5,$
     	3.6E-8,2.1E-7,6.7E-7,1.6E-6,3.5E-6,7.5E-6,9.5e-6,$
     	4.1E-8,2.5E-7,8.2E-7,1.9E-6,4.3E-6,8.9E-6,1.1e-5,$
	4.4E-8,2.9E-7,9.8E-7,2.4E-6,5.3E-6,1.1E-5,1.4e-5,$
     	4.4E-8,3.0E-7,1.1E-6,2.7E-6,6.3E-6,1.3E-5,1.7e-5,$
	4.0e-8,3.0e-7,1.2e-6,3.0e-6,7.0e-6,1.5e-5,2.1e-5,$
       	3.0E-8,2.4E-7,1.1E-6,3.1E-6,7.5E-6,1.8E-5,2.6e-5]

        s(*,0)=[2.1e-26,3.2e-26,6.5e-26,2.1e-25,1.3e-24,1.4e-23,1.2e-22]
        s(*,1)=[1.0e-17,1.3e-17,2.0e-17,4.3e-17,1.5e-16,9.4e-16,5.0e-15]
        s(*,2)=[3.0e-13,3.4e-13,4.4e-13,7.1e-13,1.7e-12,6.1e-12,1.5e-11]
        s(*,3)=[6.7e-11,7.3e-11,8.6e-11,1.1e-10,2.0e-10,4.9e-10,7.6e-10]
        s(*,4)=[1.3e-9,1.4e-9,1.5e-9,1.9e-9,2.7e-9,5.0e-9,6.4e-9]
        s(*,5)=[6.9e-9,7.2e-9,7.7e-9,8.9e-9,1.2e-8,1.9e-8,2.2e-8]
        s(*,6)=[1.8e-8,1.8e-8,1.9e-8,2.1e-8,2.7e-8,4.0e-8,4.5e-8]
        s(*,7)=[2.8e-8,2.9e-8,3.0e-8,3.3e-8,4.1e-8,5.8e-8,6.7e-8]
        s(*,8)=[3.4e-8,3.5e-8,3.6e-8,3.9e-8,4.8e-8,6.7e-8,7.7e-8]
        s(*,9)=[3.4e-8,3.4e-8,3.6e-8,3.9e-8,4.7e-8,6.5e-8,7.7e-8]
        s(*,10)=[2.5e-8,2.6e-8,2.6e-8,2.8e-8,3.3e-8,4.6e-8,5.8e-8]

        al(*,0)=[1.2e-12,1.7e-12,2.9e-12,7.1e-12,2.7e-11,1.6e-10,1.4e-9]
        al(*,1)=[6.1e-13,7.3e-13,1.0e-12,1.7e-12,3.9e-12,1.4e-11,7.1e-11]
        al(*,2)=[3.3e-13,3.6e-13,4.3e-13,5.7e-13,9.2e-13,2.0e-12,4.8e-12]
        al(*,3)=[1.8e-13,1.9e-13,2.1e-13,2.4e-13,3.1e-13,4.8e-13,7.0e-13]
        al(*,4)=[1.0e-13,1.0e-13,1.1e-13,1.2e-13,1.3e-13,1.6e-13,1.9e-13]
        al(*,5)=[5.6e-14,5.7e-14,5.7e-14,5.9e-14,6.1e-14,6.5e-14,7.2e-14]
        al(*,6)=[3.0e-14,3.0e-14,3.0e-14,3.0e-14,3.0e-14,3.0e-14,3.2e-14]
        al(*,7)=[1.5e-14,1.5e-14,1.5e-14,1.5e-14,1.5e-14,1.4e-14,1.5e-14]
        al(*,8)=[7.3e-15,7.3e-15,7.2e-15,7.1e-15,6.9e-15,6.6e-15,6.7e-15]
        al(*,9)=[3.4e-15,3.4e-15,3.3e-15,3.3e-15,3.2e-15,3.0e-15,3.0e-15]
        al(*,10)=[6.5e-16,6.5e-16,6.4e-16,6.4e-16,6.2e-16,5.8e-16,5.7e-16]


;
; the following are the spontaneous emission coeffs for n=2 to 1, 3 to 1, ... 16 to 1
;
	A_lyman=[4.699E8,5.575E7,1.278E7,4.125E6,1.644E6,7.568E5,3.869E5,2.143E5,1.263E5,7.834E4,$
                 5.066E4,3.393E4,2.341E4,1.657E4,1.200E4]
	A121=5.066E4
	A131=3.393E4
	A141=2.341E4
	A151=1.657E4
	A161=1.200E4

;
; the following are the spontaneous emission coeffs for n=3 to 2, 4 to 2, ... 16 to 2
	A_balmer=[4.41E7,8.42E6,2.53E6,9.732E5,4.389e5,2.215e5,1.216e5,7.122e4,4.397e4,2.83e4,$
                  18288.8,12249.1,8451.26,5981.95,4332.13]


        ;check ne/te arrays and form indices for use in interpolate
        num_te=n_elements(te)
        num_ne=n_elements(ene)
        num_n0=n_elements(n0)
        IF num_te NE num_ne OR num_te NE num_n0 OR num_ne NE num_n0 THEN RETURN,-1
        te_index=interp_vec_reform(alog10(temp),alog10(te))
        ne_index=interp_vec_reform(alog10(dens),alog10(ene))

        ;interpolate s and a tables to get rec/ion rates
        act_s_i=interpolate(alog10(s),ne_index,te_index,missing=0.0)
        act_s=10.^act_s_i
        act_a_i=interpolate(alog10(al),ne_index,te_index,missing=0.0)
        act_a=10.^act_a_i
        tot_rec=act_a*ene^2
        tot_ion=act_s*n0*ene

        n_lev_offset=2	;sets n[0]=2
	n=indgen(14)+n_lev_offset
	ryd=13.605
	POT1=(ryd/n^2)
	PDEL1=POT1-ryd

        ;initialize lyman/balmer emiss
	lyman_emiss=dblarr(num_te, 14)
	balmer_emiss=dblarr(num_te, 13)
        
        ;intialize other crap
	balmer_recomb=dblarr(num_te,13)
	lyman_recomb=dblarr(num_te,14)
	pop_dens_exc=dblarr(num_te,14)
	pop_dens_rec=dblarr(num_te,14)
        tot_exc=dblarr(num_te)
        grn1=n0

        FOR i=0,n_elements(n)-1 DO BEGIN
            	SAHA=ENE^2*N[i]^2*4.14E-16*EXP(POT1[i]/TE)/((TE*11605.)^1.5)  ;saha density for level n
                SAHAR=N[i]^2*EXP(PDEL1[i]/TE)  ;ratio of saha density for level n to ground state

		; NOTE that i is the n level index and that we are considering levels with principle
                ; quantum #'s from n=2 thru n=15. NOTE that n=1 the gnd state is not calcuated, 
		; rather it is given and it is assumed that recombination population directly to
		; the gnd state doesn't influence things.
            
		; NOW INTERPOLATE TO GET THE RO AND R1 VALUES
;
		IF(n[i] LE 6) THEN BEGIN
                	r0=interpolate(alog10(R(*,*,0,n[i]-n_lev_offset)),ne_index,te_index,missing=0.0)
			r0=10.^r0
			r1=interpolate(alog10(R(*,*,1,n[i]-N_lev_offset)),ne_index,te_index,missing=0.0)
			r1=10.^r1

			IF(i EQ 0) THEN BEGIN
				lyman_emiss[*,i]=(R0*SAHA+R1*SAHAR*GRN1)*A_lyman[i] 
				tot_exc=tot_exc+(R0*SAHA+R1*SAHAR*GRN1)
				lyman_recomb[*,i]=R0*SAHA*A_lyman[i]
				pop_dens_exc[*,i]=R1*SAHAR*GRN1
				pop_dens_rec[*,i]=R0*SAHA		
                        ENDIF ELSE BEGIN
                        	tot_exc=tot_exc+(R0*SAHA+R1*SAHAR*GRN1)
				lyman_emiss[*,i]=(R0*SAHA+R1*SAHAR*GRN1)*A_lyman[i]  ; the units of photons/s/cm^3
				lyman_recomb[*,i]=R0*SAHA*A_lyman[i]
				pop_dens_exc[*,i]=R1*SAHAR*GRN1
				pop_dens_rec[*,i]=R0*SAHA
	       			balmer_emiss[*,i-1]=(R0*SAHA+R1*SAHAR*GRN1)*A_balmer[i-1] ; the units of photons/s/cm^3
				balmer_recomb[*,i-1]=R0*SAHA*A_balmer[i-1]

                       ENDELSE
                ENDIF ELSE BEGIN
			tot_exc=tot_exc+SAHA
			lyman_emiss[*,i]=SAHA*A_lyman[i]
			balmer_emiss[*,i-1]=SAHA*A_balmer[i-1]
			lyman_recomb[*,i]=R0*SAHA*A_lyman[i]
			pop_dens_exc[*,i]=0.1   ;no idea where this comes from (MLR)
			pop_dens_rec[*,i]=R0*SAHA
			balmer_recomb[*,i-1]=R0*SAHA*A_balmer[i-1]
                ENDELSE
        	IF keyword_set(debug) THEN stop
        ENDFOR
        
        IF NOT keyword_set(quiet) THEN BEGIN
        	print, 'For THE n0 given, the total excited/ground density is'
                print, tot_exc/grn1
        ENDIF

        emiss_str={lyman:lyman_emiss,balmer:balmer_emiss}
        recomb_str={lyman:lyman_recomb, balmer:balmer_recomb}
        misc_str={tot_exc:tot_exc,pop_dens_exc:pop_dens_exc,pop_dens_rec:pop_dens_rec,a_lyman:a_lyman,a_balmer:a_balmer}


        output={ene:ene, te:te, n0:n0, emiss:emiss_str,recomb:recomb_str,misc:misc_str}
        IF keyword_set(debug) THEN stop
        RETURN,output

END


FUNCTION lyman_ratio,fo,fo_lam,fwhm,te,ene,f_n0,zeff=zeff,debug=debug
	IF NOT keyword_set(zeff) THEN zeff=1.0

	spec_emiss_brem=4.60*zeff^0.945*ene^2/te^0.318*exp(-10.2/te)  	;[W/m^3/Ang]
        out=johnson_hinnov_emissivity(ene*1.0e14,te,f_n0*ene*1.0e14,debug=debug)
        emiss_lyman=out.emiss.lyman[0]*1.0e6*ang2ev(1215.3)*1.6e-19	;[W/m^3]

        ratio=fo/fo_lam*sqrt(!pi)/(2*sqrt(alog(2)))*fwhm*spec_emiss_brem/emiss_lyman

        RETURN,ratio
END

FUNCTION neutral_density,n_e,t_e,ly_em,debug=debug,jhdebug=jhdebug
	;n_e in 10^20 m^-3
	;T_e in eV
	;ly_em in W/m^3

        IF n(n_e) EQ n(t_e) AND n(n_e) EQ n(ly_em) THEN npts=n(n_e)+1 ELSE RETURN,-1

        n_0=fltarr(npts)
        n0_test=[1.0e-5,5.0e-5,1.0e-4,5.0e-4,0.001,0.005,0.01,0.05,0.1]*1.0e14

        n_test=n(n0_test)+1
        ly_test=fltarr(n_test)
        FOR i=0,npts-1 DO BEGIN
            	FOR j=0,n_test-1 DO BEGIN
        		out=johnson_hinnov_emissivity(n_e[i]*1.0e14,t_e[i],n0_test[j],debug=jhdebug,/quiet)
                	ly_test[j]=out.emiss.lyman[0]*1.0e6*ang2ev(1215.3)*1.6e-19 ;[W/m^3]
                ENDFOR
                n_0[i]=interpol(n0_test,ly_test,ly_em[i])/1.0e14
        ENDFOR

        output=n_0
        IF keyword_set(debug) THEN stop
        RETURN,output
END

FUNCTION neutral_rate,n_e,t_e,n0=n0,jhdebug=jhdebug
	IF NOT keyword_set(n0) THEN n0=1.0e-3*1.0e14
	;n_e in 10^20 m^-3
	;T_e in eV
	;n0 in 10^20 m^-3

	n_te=n(t_e)+1
        n_ne=n(n_e)+1
        

        rate=fltarr(n_te,n_ne)
        FOR i=0,n_ne-1 DO BEGIN
            	FOR j=0,n_te-1 DO BEGIN
        		out=johnson_hinnov_emissivity(n_e[i]*1.0e14,t_e[j],n0,debug=jhdebug,/quiet)
                	rate[j,i]=out.emiss.lyman[0]*1.0e6*ang2ev(1215.3)*1.6e-19/n0 ;[W]
                ENDFOR
        ENDFOR

        output=rate
        IF keyword_set(debug) THEN stop
        RETURN,output
END
