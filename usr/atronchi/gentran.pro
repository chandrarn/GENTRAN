;+
;NAME:
;	CALC_FRACABUND
;
;PURPOSE:
;	This function calculates the temperature dependent fractional abundance of impurities 
;	given an ionization and recombination rate structure
;
;CALLING SEQUECTION:
;	result=CALC_FRACABUND(ion,rec)
;	
;INPUTS:
;	ion:	STRUC	ionization rate structure (see READ_ION_DATA for format)
;	rec:	STRUC	recombination rate structure (see READ_REC_DATA for format)
;
;OPTIONAL INPUTS:
;	te:	FLTARR	[nt] of the temperature values at which to calculate the data [eV] DEFAULT: rec.temp
;	nel	FLT	of the electron density to use for the rates if density dependent [m^-3] DEFAULT: 1.0e20
;
;OUTPUTS:
;	result:	STRUC	containing the fractional abundance data
;		*.fq	FLTARR	[z+1,nt] of the fractional abundances [i,*] for ionization stage i
;		*.temp	FLTARR	[nt] of the electron temperatures [eV]
;		*.dens	FLTARR 	[1] of the electron density [m^-3]
;		*.z	INT	it's just the z folks
;		*.rec	STR	path from rec.path
;		*.ion	STR	path from ion.path
;
;MODIFICATION HISTORY:
;	Written by:	ML Reinke - 2/09
;
;-

FUNCTION calc_fracabund,ion,rec,te=te,nel=nel,debug=debug

	z=rec.z	
	temp=rec.temp
	IF NOT keyword_set(te) THEN te=temp
	ntemp=n(te)+1

	IF NOT keyword_set(nel) THEN nel=1.0e20
	IF n(ion.dens) EQ 0 THEN idpt=0 ELSE idpt=ipt(ion.dens,nel)
	IF n(rec.dens) EQ 0 THEN rdpt=0 ELSE rdpt=ipt(rec.dens,nel)
	
	;both are at the same temp
	fq=fltarr(z+1,ntemp)
	A=dblarr(z+1,z+1)
	b=fltarr(z+1)
	b[0]=1.0

;	MATRIX HELP
;	matrix=[[1,2,3],[4,5,6],[7,8,9]]
;	print, matrix
;		1       2       3
;      		4       5       6
;		7       8       9
;	print, matrix[0,1]
;		4
;	print, matrix[2,1]
;		6	
;	print, matrix##[1,2,1]
;		6
;		20
;		24

	itemp=alog10(ion.temp)
	rtemp=alog10(rec.temp)
	logte=alog10(te)
	FOR i=0,ntemp-1 DO BEGIN
		FOR j=0,z DO BEGIN
			A[j,j]=-interpol(ion.rates[j,*,idpt],itemp,logte[i])-interpol(rec.rates[j,*,rdpt],rtemp,logte[i])
			IF j NE 0 THEN A[j-1,j]=interpol(ion.rates[j-1,*,idpt],itemp,logte[i])
			IF j NE z THEN A[j+1,j]=interpol(rec.rates[j+1,*,rdpt],rtemp,logte[i])
			;A[j,j]=-interpol(ion.rates[j,*,idpt],ion.temp,te[i])-interpol(rec.rates[j,*,rdpt],rec.temp,te[i])
			;IF j NE 0 THEN A[j-1,j]=interpol(ion.rates[j-1,*,idpt],ion.temp,te[i])
			;IF j NE z THEN A[j+1,j]=interpol(rec.rates[j+1,*,rdpt],rec.temp,te[i])
		ENDFOR
		A[*,0]=1.0
		x=la_invert(A,/double)##b
		tmp=where(x LT 0)
		IF tmp[0] NE -1 THEN x[tmp]=0.0
		fq[*,i]=x
	ENDFOR

	output={fq:fq,temp:te,dens:[nel],z:z,rec:rec.path,ion:ion.path}
	IF keyword_set(debug) THEN stop
	RETURN,output	
END

;+
;NAME:
;	GENTRAN_FQ
;
;MODIFICATION HISTORY:
;	Written by:	M.L. Reinke - 2/10/11
;
;-

FUNCTION gentran_fq,z,fq,temp=temp,load=load,adas=adas
	rec=read_rec_data(z,load=load,adas=adas)
	ion=read_ion_data(z,load=load,adas=adas)
	IF NOT keyword_set(temp) THEN temp=10^(make(0,4.3,200))
	fq=calc_fracabund(ion,rec,te=temp)
	RETURN,fq
END

;+
;NAME:
;	CALC_VOLTERM_PROFILE
;
;PURPOSE:
;	This function computes the volume weighting term from EFIT reconstruction
;
;CALLING SEQUENCE:
;	result=CALC_VOLTERM_PROFILE(shot,rhovec)
;
;INPUTS:
;	shot	LONG	shot number
;	rhovec	FLTARR	[nrho,ntau] of the radial grid points to compute volume weighting [r/a]
;
;OPTIONAL INPUTS:
;	tr	FLTARR	[t1,t2] to truncate the time range of the calculation DEFAULT: EFIT time base
;	tree	STRING	identifies EFIT tree to use DEFAULT: 'ANALYSIS'
;
;KEYWORD PARAMETERS"
;	circular	/circular will assume circular flux surfaces, returning result=1/rhovec
;
;OUTPUTS:
;	result	FLTARR [nrho,ntau] of the "volume weighting term" (see PROCEDURE)	
;
;OPTIONAL OUTPUTS:
;	etree	STRING	will be filled with 'ANALYSIS' if left unspecified
;
;PROCEDURE:
;	Changes in radial particle flux occur due to changes in the area as function of minor radius.
;	In the continuity equation (dV/drho)^-1*d/drho(dV/drho*FLUX) term reduces to dFLUX/drho+(d2V/drho2)/(dV/drho)*FLUX
;	The "volume weighting term calculated here is (d2V/drho2)/(dV/drho), which for a cylinder is simply 1/rho.
;
;MODIFICATION HISTORY:
;	Written by:	M.L. Reinke - spring 2009
;	6/6/12		M.L. Reinke - added the ability to specify the EFIT tree
;
;-

FUNCTION calc_volterm_profile,shot,rhovec,tr=tr,debug=debug,circular=circular,etree=etree

	IF NOT keyword_set(etree) THEN etree='ANALYSIS'
	IF etree EQ 'ANALYSIS' THEN BEGIN
		mdsopen,'analysis',shot
		volp=mdsvalue('\ANALYSIS::TOP:EFIT.RESULTS.FITOUT:VOLP')
		time=mdsvalue('dim_of(\ANALYSIS::TOP:EFIT.RESULTS.FITOUT:VOLP,0)')
		rmid=mdsvalue('\efit_rmid')
		mdsclose,'analysis',shot
	ENDIF ELSE BEGIN
		mdsopen,etree,shot
		volp=mdsvalue('\'+etree+'::TOP.RESULTS.FITOUT:VOLP')
		time=mdsvalue('dim_of('+etree+'+::TOP.RESULTS.FITOUT:VOLP,0)')
		rmid=mdsvalue('\efit_rmid')	
		mdsclose,etree,shoth
	ENDELSE

	IF keyword_set(tr) THEN BEGIN
		tmp=where(time GE tr[0] AND time LE tr[1])
		IF tmp[0] EQ -1 THEN RETURN,-1
		time=time[tmp]
		volp=volp[tmp,*]
		rmid=rmid[tmp,*]
	ENDIF
	n_time=n(time)+1
	n_grid=n(rhovec)+1

	volterm=fltarr(n_grid,n_time)
	FOR i=0,n_time-1 DO BEGIN
		efitrho=(rmid[i,*]-rmid[i,0])/(max(rmid[i,*])-rmid[i,0])
		dvdrho=deriv(efitrho,volp[i,*])
		dv2drho2=deriv(efitrho,dvdrho)
		IF keyword_set(circular) THEN volterm[*,i]=1.0/rhovec ELSE $
			volterm[*,i]=interpol(dv2drho2,efitrho,rhovec)/interpol(dvdrho,efitrho,rhovec)	;interpolate from efit grid to GENTRAN grid
	ENDFOR
	
	output=volterm
	IF keyword_set(debug) THEN stop
	RETURN,output
END

;7/19/11 if lyman data is not present then neutral density set to 10^-4*ne
;7/5/12 modified to use t,dt to specify multiple time points
PRO lyman2neut,shot,t,dt,temp,dens,lyrho,neut,noly=noly
	IF keyword_set(noly) $
	 THEN status=0 $	;force to use the no-lyman condition if data is corrupted
	 ElSE BEGIN
		;get lyman-a data
		mdsopen,'spectroscopy',shot
			em=mdsvalue('\SPECTROSCOPY::TOP.BOLOMETER.RESULTS.DIODE.LYMID:EMISS',status=status,/quiet)
			r=mdsvalue('dim_of(\SPECTROSCOPY::TOP.BOLOMETER.RESULTS.DIODE.LYMID:EMISS,0)',/quiet)
			tau=mdsvalue('dim_of(\SPECTROSCOPY::TOP.BOLOMETER.RESULTS.DIODE.LYMID:EMISS,1)',/quiet)
		mdsclose,'spectroscopy',shot

		;get efit_data
		mdsopen,'analysis',shot
			rmid=mdsvalue('\efit_rmid')
			efit_t=mdsvalue('dim_of(\efit_rmid)')
		mdsclose,'analysis',shot
	ENDELSE

	ntime=n(t)+1
	IF status THEN BEGIN
		nrho=n(r)+1
		neut=fltarr(nrho,ntime)
		lyrho=fltarr(nrho,ntime)
		FOR i=0,ntime-1 DO BEGIN
			i1=ipt(efit_t,t[i])
			i2=ipt(efit_t,t[i]+dt[i])
			ro=min(sum_array(rmid[i1:i2,*],/j)/(i2-i1+1.0))
			a=max(sum_array(rmid[i1:i2,*],/j)/(i2-i1+1.0))-ro

			ilyrho=(r-ro)/a
			i1=ipt(tau,t[i])
			i2=ipt(tau,t[i]+dt[i])
			lyem=sum_array(em[*,i1:i2],/i)/(i2-i1+1.0)
			temp_int=interpol(temp.temp[*,i],temp.rho[*,i],ilyrho)
			dens_int=interpol(dens.dens[*,i],dens.rho[*,i],ilyrho)
			neut[*,i]=neutral_density(dens_int*1.0e-20,temp_int,lyem)*1.0e20
			lyrho[*,i]=ilyrho
		ENDFOR
	ENDIF ELSE BEGIN
		ilyrho=make(0.87,0.92,20)
		neut=fltarr(20,ntime)
		lyrho=fltarr(20,ntime)
		FOR i=0,ntime-1 DO BEGIN
			lyrho[*,i]=ilyrho
			dens_int=interpol(dens.dens[*,i],dens.rho[*,i],ilyrho)
			neut[*,i]=dens_int*1.0e-4					;fix to small edge neutral density
		ENDFOR
	ENDELSE
END

;+
;NAME:
;	CALC_NEUTSTR
;
;PURPOSE:
;	This function computs the neutral density structure using experimental
;	hydrogenic emission data and recombination equilibrium.
;
;CALLING SEQUENCE:
;	result=CALC_NEUTSTR(shot,t1,t2,temp,dens)
;
;INPUTS:
;	shot	LONG 	shot number
;	t	FLTARR	start of time bins [sec]
;	dt	FLTARR	duration of time bins [sec]
;	temp	STRUC	structure with electron temperature data (see GENTRAN)
;	dens	STRUC	structure with electron density data (see GENTRAN)
;
;OPTIONAL INPUTS:
;	eno	FLTARR	[nrho] of the uncertainty in the neutral density.  If eno is FLOAT it is assume
;			to be the fractional error.	
;
;KEYWORD PARAMETERS:
;	debug	/debug stops the code before the RETURN statement
;
;OUTPUTS:
;	result	STRUC	structure of the neutral density data
;		*.dens	[nrho]	of the neutral density [m^-3]
;		*.rho	[nrho]	of the r/a values
;		*.err	[nrho]	of the absolute error [m^-3]
;		*.bal
;		*.balrho
;		*.ly
;		*.lyrho
;
;MODIFICATION HISTORY:
;	Written by:	M.L. Reinke - Spring 2009
;	6/6/12		M.L. Reinke - modified path to hydrogen ion/rec data to use /usr/local/cmod/idl/atomic_physics/
;       7/5/12		M.L. Reinke - modified to use t and dt instead of t1, t2 in order to get multiple time slices
;
;-

FUNCTION calc_neutstr,shot,t,dt,temp_str,dens_str,eno=eno,debug=debug,noly=noly
	;ion=read_scd_file('/home/mlreinke/idl/impurities/data/adas/scd96_h.dat')
	;rec=read_acd_file('/home/mlreinke/idl/impurities/data/adas/acd96_h.dat')
	ion=read_scd_file('/usr/local/cmod/idl/atomic_physics/adas/scd96_h.dat')
	rec=read_acd_file('/usr/local/cmod/idl/atomic_physics/adas/acd96_h.dat')
	ntime=n(t)+1
	nrho=n(temp_str.rho[*,0])+1	
	neut=fltarr(nrho,ntime)
	balneut=fltarr(nrho,ntime)
	rho=fltarr(nrho,ntime)
	err=fltarr(nrho,ntime)
	lyman2neut,shot,t,dt,temp_str,dens_str,lyrho,lyneut,noly=noly

	FOR i=0,ntime-1 DO BEGIN
		rho[*,i]=temp_str.rho[*,i]
		dens=interpol(dens_str.dens[*,i],dens_str.rho[*,i],rho[*,i])		;this shouldn't be necessary
		temp=temp_str.temp[*,i]

		te_rec_reform=interp_vec_reform(rec.temp,temp)
		ne_rec_reform=interp_vec_reform(rec.dens,dens)
		te_ion_reform=interp_vec_reform(ion.temp,temp)
		ne_ion_reform=interp_vec_reform(ion.dens,dens)
		rec_int=interpolate(rec.rates[1,*,*],te_rec_reform,ne_rec_reform)
		ion_int=interpolate(ion.rates[0,*,*],te_ion_reform,ne_ion_reform)

		balrho=rho[*,i]
		balneut[*,i]=dens*rec_int/ion_int
	
		yc=alog10(balneut[where(balrho LT 0.2),i])
		xc=balrho[where(balrho LT 0.2)]
		lyfit=where(lyrho[*,i] LT 1.0 AND finite(lyneut[*,i]) EQ 1 AND lyneut[*,i] GT 0)
		tmp=min(lyrho[where(finite(lyneut[*,i]) EQ 1 AND lyneut[*,i] GT 0)], lyfit)
		ye=alog10(lyneut[lyfit,i])
		xe=lyrho[lyfit,i]
		xfit=[xc,xe]
		yfit=[yc,ye]
		coefs=poly_fit(xfit,yfit,2)

		tmp=where(rho[*,i] LT 0.975)
		neut[tmp,i]=10^(coefs[0]+coefs[1]*rho[tmp]+coefs[2]*rho[tmp,i]^2)
	
		tmp=where(rho[*,i] GE 0.975)
		lygood=where(finite(lyneut[*,i]) EQ 1 AND lyneut[*,i] GT 0)
		neut[tmp,i]=interpol(lyneut[lygood,i],lyrho[lygood,i],rho[tmp,i])

		IF NOT keyword_set(eno) THEN eno=0.0
		IF n(eno) EQ n(neut[*,i]) THEN err[*,i]=eno ELSE err[*,i]=neut[*,i]*eno[0]
	ENDFOR
	output={dens:neut,rho:rho,err:err,bal:balneut,ly:lyneut,lyrho:lyrho}
	IF keyword_set(debug) THEN stop
	RETURN,output
END

;+
;NAME:
;	WRITE_IONREC_COEFS_TABLE
;
;PURPOSE:
;	This procedure computes the ion/rec rate coefficients and the kinetic profiles on the radial grid
;	defined for the simulation and saves it to a temporary file.
;
;CALLING SEQUENCE:
;	WRITE_IONREC_COEFS_TABLE,rho,ion,rec,cxr,temp,dens,neut,a
;
;INPUTS:
;	rho	FLTARR	[nrho] of the r/a points used in the simulation
;	ion	STRUC	ionization rate structure (see READ_ION_DATA)
;	rec	STRUC	recombination rate structure (see READ_REC_DATA)
;	cxr	STRUC	charge exchange rate structure (see READ_REC_DATA)
;	temp	STRUC	electron temperature structure (see GENTRAN)
;	dens	STRUC	electron density structure (see GENTRAN)
;	neut	STRUC	neutral density structure (see GENTRAN)
;	a	FLOAT	minor radius
;
;OPTIONAL INPUTS:
;	path	path to save IDL savesets computed by this function DEFAULT: '/tmp/ionrec_coefs.dat'
;
;OUTPUTS:
;	The ionrec_coefs structure is saved to the file defined by PATH.  The format is for [0,*] to be the neutral
;	species and [z,*] to be the full stripped.  Ionization rates for [z,*]=0 and recombination rates for [0,*]=0.
;		ionrec_coefs	STRUC	interpolated rate and kinetic profile data to be used in GENTRAN
;				*.i	[z+1,nrho] ionization rate coeff for each charge state [m^3/s]
;				*.di	[z+1,nrho] uncertainty in the ionization rate coeff [m^3/s]
;				*.r	[z+1,nrho] recombination rate coeff for each charge state [m^3/s]
;				*.dr	[z+1,nrho] uncertainty in the recombination rate coeff [m^3/s]
;				*.c	[z+1,nrho] charge exchange rate coeff for each charge state [m^3/s]
;				*.dc	[z+1,nrho] uncertainty in the charge exchange rate coeff [m^3/s]
;				*.d	[nrho] electron density [m^-3]
;				*.derr	[nrho] uncertainty in the electron density [m^-3]
;				*.t	[nrho] electron temperature [eV]
;				*.terr	[nrho] uncertainty in the electron temperature [eV]
;				*.n	[nrho] neutral density [m^-3]
;				*.nerr	[nrho] uncertainty in the neutral density [m^-3]
;
;MODIFICATION HISTORY:
;	Written by:	ML Reinke - Spring 2009
;	6/7/12		ML Reinke - modified the path to save to /tmp/
;
;-

PRO write_ionrec_coefs_table,rho,ion,rec,cxr,temp,dens,neut,a,path=path
	;IF NOT keyword_set(path) THEN path='/home/mlreinke/idl/genie/atomic_physics/ionrec_coefs.dat'
	IF NOT keyword_set(path) THEN path='/tmp/ionrec_coefs.dat'
	
	ff=1.0
	fq=15
	ion.rates[fq,*,*]*=ff
	rec.rates[fq+1,*,*]/=ff
	;stop
		

	dens_int=interpol(dens.dens,dens.rho,rho)
	derr_int=interpol(dens.err,dens.rho,rho)
	temp_int=interpol(temp.temp,temp.rho,rho)
	terr_int=interpol(temp.err,temp.rho,rho)
	z=rec.z
	nrho=n(rho)+1

	ion_int=fltarr(z+1,nrho)
	dion_int=fltarr(z+1,nrho)
	rec_int=fltarr(z+1,nrho)
	drec_int=fltarr(z+1,nrho)
	te_rec_reform=interp_vec_reform(rec.temp,temp_int)
	ne_rec_reform=interp_vec_reform(rec.dens,dens_int)
	te_ion_reform=interp_vec_reform(ion.temp,temp_int)
	ne_ion_reform=interp_vec_reform(ion.dens,dens_int)

	;check for ne/te values that are off scale (too small)
	tmp=where(te_rec_reform EQ -1)
	IF tmp[0] NE -1 THEN te_rec_reform[tmp]=0.0
	tmp=where(ne_rec_reform EQ -1)
	IF tmp[0] NE -1 THEN ne_rec_reform[tmp]=0.0
	tmp=where(te_ion_reform EQ -1)
	IF tmp[0] NE -1 THEN te_ion_reform[tmp]=0.0
	tmp=where(ne_ion_reform EQ -1)
	IF tmp[0] NE -1 THEN ne_ion_reform[tmp]=0.0
	
	FOR i=0,z DO BEGIN
		rec_int[i,*]=interpolate(rec.rates[i,*,*],te_rec_reform,ne_rec_reform)
		ion_int[i,*]=interpolate(ion.rates[i,*,*],te_ion_reform,ne_ion_reform)

		kmin=floor(min(ne_rec_reform))
		kmax=floor(max(ne_rec_reform))+1
		IF total(rec.rates[i,*,kmin]-rec.rates[i,*,kmax]) EQ 0 THEN BEGIN	
			drec_int[i,*]=interpol(deriv(rec.temp,rec.rates[i,*,kmin]),rec.temp,temp_int)
		ENDIF ELSE BEGIN
			IF i EQ 0 THEN print, 'ERROR - no recombination rate data'
		ENDELSE
		kmin=floor(min(ne_ion_reform))
		kmax=floor(max(ne_ion_reform))+1
		IF total(ion.rates[i,*,kmin]-ion.rates[i,*,kmax]) EQ 0 THEN BEGIN	
			dion_int[i,*]=interpol(deriv(ion.temp,ion.rates[i,*,kmin]),ion.temp,temp_int)
		ENDIF ELSE BEGIN
			IF i EQ 0 THEN print, 'ERROR - no ionization rate data'
		ENDELSE
	ENDFOR
	IF size(cxr,/type) EQ 8 THEN BEGIN
		cxr_int=fltarr(z+1,nrho)
		dcxr_int=fltarr(z+1,nrho)
		neut_int=interpol(neut.dens,neut.rho,rho)
		nerr_int=interpol(neut.err,neut.rho,rho)
		te_cxr_reform=interp_vec_reform(rec.temp,temp_int)
		ne_cxr_reform=interp_vec_reform(rec.dens,dens_int)

		;check for ne/te values that are off scale (too small)
		tmp=where(te_cxr_reform EQ -1)
		IF tmp[0] NE -1 THEN te_cxr_reform[tmp]=0.0
		tmp=where(ne_cxr_reform EQ -1)
		IF tmp[0] NE -1 THEN ne_cxr_reform[tmp]=0.0

		FOR i=0,z DO BEGIN
			cxr_int[i,*]=interpolate(cxr.rates[i,*,*],te_cxr_reform,ne_cxr_reform)
			kmin=floor(min(ne_cxr_reform))
			kmax=floor(max(ne_cxr_reform))+1
			IF total(cxr.rates[i,*,kmin]-cxr.rates[i,*,kmax]) EQ 0 THEN BEGIN	
				dcxr_int[i,*]=interpol(deriv(cxr.temp,cxr.rates[i,*,kmin]),cxr.temp,temp_int)
			ENDIF ELSE BEGIN
				IF i EQ 0 THEN print, 'ERROR - no charge exchange rate data'
			ENDELSE
		ENDFOR
	ENDIF ELSE BEGIN
		neut_int=fltarr(nrho)
		nerr_int=fltarr(nrho)
		cxr_int=fltarr(z+1,nrho)
		dcxr_int=fltarr(z+1,nrho)
	ENDELSE
	ionrec_coefs={i:ion_int,di:dion_int,r:rec_int,dr:drec_int,c:cxr_int,dc:dcxr_int,d:dens_int,derr:derr_int,t:temp_int,terr:terr_int,n:neut_int,nerr:nerr_int}

	save,rho,ionrec_coefs,z,a,filename=path
END

;+
;NAME:
;	WRITE_TRANSPORT_COEFS_TABLE
;	
;PURPOSE:
;	This procedure composes the transport coefficient structure and saves it to a temporary file
;
;CALLING SEQUENCE:
;	WRITE_TRANSPORT_COEFS_TABLE,rho,x,dxdrho,d2xdrho2,diff,conv,volterm,a,r
;
;INPUTS:
;	rho		FLTARR	[nrho] of the r/a 
;	x		FLTARR	[nrho]
;	dxdrho		FLTARR	[nrho]
;	d2xdrho2	FLTARR	[nrho]
;	diff		STRUC   
;	conv		STRUC
;	volterm		FLTARR	[nrho]
;	a		FLOAT	minor radius [m]
;-

PRO write_transport_coefs_table,rho,x,dxdrho,d2xdrho2,diff,conv,volterm,a,path=path,ntot_theory=ntot_theory
	;IF NOT keyword_set(path) THEN path='/home/mlreinke/idl/genie/atomic_physics/transport_coefs.dat'
	IF NOT keyword_set(path) THEN path='/tmp/transport_coefs.dat'

	diff_int=interpol(diff.diff,diff.rho,rho)
	diff_dr_int=deriv(rho,diff_int)
	conv_int=interpol(conv.conv,conv.rho,rho)
	conv_dr_int=deriv(rho,conv_int)

	trans_coefs={a:diff_int/a^2*dxdrho^2, b:(diff_int/a^2*volterm+diff_dr_int/a^2-conv_int/a)*dxdrho+diff_int/a^2*d2xdrho2,$
			c:-1.0*(conv_int*volterm/a+conv_dr_int/a)}
	
	ntot_theory=fltarr(n(rho)+1)
	FOR i=0,n(rho)-1 DO ntot_theory[i]=exp(-a*int_tabulated(rho[i:*],conv_int[i:*]/diff_int[i:*]))
	ntot_theory/=ntot_theory[0]
	;stop
	save,rho,x,trans_coefs,filename=path
END

;+
;NAME:
;	FORM_TRANSPORT_MATRIX
;
;PURPOSE:
;	This function computes the transport matrix from the interpolated ion/rec rates and the
;	defined transport (diff/conv) profiles.
;
;CALLING SEQUENCE:
;	result=FORM_TRANSPORT_MATRIX()
;
;OPTIONAL INPUTS:
;	trans_path	STRING	path to save file where the transport coefficients are stored DEFAULT: '/tmp/transport_coefs.dat'
;	ionrec_path	STRING	path to save file where the ion/rec coefficients are stored DEFAULT: '/tmp/ionrec_coefs.dat'
;
;OUTPUTS:
;	result	DBLARR	[z*nrho,z*nrho]
;
;MODIFICATION HISTORY:
;	Written by:	M.L. Reinke - Spring 2009
;	6/7/12		M.L. Reinke - modified to use /tmp for temporary files
;
;-

FUNCTION form_transport_matrix,trans_path=trans_path,ionrec_path=ionrec_path,debug=debug,double=double	
	;IF NOT keyword_set(trans_path) THEN trans_path='/home/mlreinke/idl/genie/atomic_physics/transport_coefs.dat'
	;IF NOT keyword_set(ionrec_path) THEN ionrec_path='/home/mlreinke/idl/genie/atomic_physics/ionrec_coefs.dat'
	IF NOT keyword_set(trans_path) THEN trans_path='/tmp/transport_coefs.dat'
	IF NOT keyword_set(ionrec_path) THEN ionrec_path='/tmp/ionrec_coefs.dat'
	restore,trans_path
	restore,ionrec_path

	rho=x				;partial derivatives are computed in "pixel" units with the scaling handled in a,b,c
	nrho=n(rho)+1
	dens=ionrec_coefs.d
	neut=ionrec_coefs.n

	alpha=dblarr(z*nrho,z*nrho)
	a=trans_coefs.a
	b=trans_coefs.b
	c=trans_coefs.c
	rec=ionrec_coefs.r
	ion=ionrec_coefs.i
	cxr=ionrec_coefs.c
	FOR q=1,z DO BEGIN
		e=-1.0*(rec[q,*]+ion[q,*])
		IF q NE 1 THEN f=ion[q-1,*] ELSE f=fltarr(nrho)
		IF q NE z THEN d=rec[q+1,*] ELSE d=fltarr(nrho)
		FOR k=0,nrho-1 DO BEGIN
			;fill spatial data
			CASE k OF
				0 : BEGIN		;core dn/dx=0, v=0, dv/dx=0
					;drho=rho[k]-rho[k+1]
					drho=rho[k+1]-rho[k]
					alpha[(q-1)*nrho+k,(q-1)*nrho+k]+=-4.0*a[k]/drho^2
					alpha[(q-1)*nrho+(k+1),(q-1)*nrho+k]+=4.0*a[k]/drho^2		
				END
				nrho-1 : BEGIN		;edge
					drho=rho[k]-rho[k-1]
					alpha[(q-1)*nrho+k,(q-1)*nrho+k]+=c[k]-2.0*a[k]/drho^2
				END
				ELSE : BEGIN
					drho=(rho[k+1]-rho[k-1])/2.0
					alpha[(q-1)*nrho+k,(q-1)*nrho+k]+=c[k]-2.0*a[k]/drho^2
					alpha[(q-1)*nrho+(k-1),(q-1)*nrho+k]+=a[k]/drho^2-b[k]/(2.0*drho)
					alpha[(q-1)*nrho+(k+1),(q-1)*nrho+k]+=a[k]/drho^2+b[k]/(2.0*drho)
				END
			ENDCASE
			;fill atomic physics data
			alpha[(q-1)*nrho+k,(q-1)*nrho+k]+=e[k]*dens[k]-cxr[q,k]*neut[k]
			IF q NE z THEN alpha[(q-1+1)*nrho+k,(q-1)*nrho+k]+=d[k]*dens[k]+cxr[q+1,k]*neut[k]
			IF q NE 1 THEN alpha[(q-1-1)*nrho+k,(q-1)*nrho+k]+=f[k]*dens[k]
		ENDFOR
	ENDFOR

	output=alpha
	IF keyword_set(debug) THEN stop
	RETURN,output
END

;+
;NAME:
;	CALC_CSDEN_ERROR
;
;PURPOSE:
;	This function calculates the uncertainty in the charge state density profile
;	due to uncertainties in electron temperature and electron and neutral density.
;
;CALLING SEQUENCE:
;	result=CALC_CSDEN_ERROR(csden)
;
;INPUTS:
;	csden	FLTARR	[z,nrho] of the charge state density profiles 
;
;OPTIONAL INPUTS:
;	trans_path	STRING	path to save file where the transport coefficients are stored DEFAULT: '/tmp/transport_coefs.dat'
;	ionrec_path	STRING	path to save file where the ion/rec coefficients are stored DEFAULT: '/tmp/ionrec_coefs.dat'
;
;KEYWORD PARAMETERS:
;	double		/double inverts the error matrix using double precession (sent to LA_INVERT)
;	debu		/debug stops the code before the RETURN statement
;
;OUTPUTS:
;	result	FLTARR	[z,nrho] of the uncertainty in the csden values
;
;PROCEDURE:
;	See M.L. Reinke APiP 2009 talk/paper
;
;MODIFICATION HISTORY:
;	Written by:	M.L. Reinke - Spring 2009
;	6/7/12		M.L. Reinke - modified paths to usr /tmp/
;
;-

FUNCTION calc_csden_error,csden,trans_path=trans_path,ionrec_path=ionrec_path,debug=debug,double=double
	;IF NOT keyword_set(trans_path) THEN trans_path='/home/mlreinke/idl/genie/atomic_physics/transport_coefs.dat'
	;IF NOT keyword_set(ionrec_path) THEN ionrec_path='/home/mlreinke/idl/genie/atomic_physics/ionrec_coefs.dat'
	IF NOT keyword_set(trans_path) THEN trans_path='/tmp/transport_coefs.dat'
	IF NOT keyword_set(ionrec_path) THEN ionrec_path='/tmp/ionrec_coefs.dat'

	restore,trans_path
	restore,ionrec_path
	rho=x
	nrho=n(rho)+1
	dens=ionrec_coefs.d
	neut=ionrec_coefs.n

	matrix=dblarr(z*nrho,z*nrho)
	tevector=dblarr(z*nrho)
	nevector=dblarr(z*nrho)
	novector=dblarr(z*nrho)
	a=trans_coefs.a
	b=trans_coefs.b
	c=trans_coefs.c
	rec=ionrec_coefs.r
	ion=ionrec_coefs.i
	cxr=ionrec_coefs.c
	drec=ionrec_coefs.dr
	dion=ionrec_coefs.di
	dcxr=ionrec_coefs.dc

	FOR q=1,z DO BEGIN
		IF q NE 1 THEN BEGIN
			izm=ion[q-1,*] 
			nzm=csden[q-1,*]
			dizm=dion[q-1,*]
		ENDIF ELSE BEGIN
			izm=fltarr(nrho)
			nzm=fltarr(nrho)
			dizm=fltarr(nrho)
		ENDELSE
		IF q NE z THEN BEGIN
			nzp=csden[q+1,*]
			rzp=rec[q+1,*]
			drzp=drec[q+1,*]
			czp=cxr[q+1,*]
			dczp=dcxr[q+1,*]
		ENDIF ELSE BEGIN
			rzp=fltarr(nrho)
			nzp=fltarr(nrho)
			drzp=fltarr(nrho)
			czp=fltarr(nrho)
			dczp=fltarr(nrho)
		ENDELSE
		iz=ion[q,*]
		diz=dion[q,*]
		rz=rec[q,*]
		drz=drec[q,*]
		cz=cxr[q,*]
		dcz=dcxr[q,*]
		nz=csden[q,*]
		FOR k=0,nrho-1 DO BEGIN
			CASE k OF
				0 : drho=rho[k+1]-rho[k]
				nrho-1 : drho=rho[k]-rho[k-1]
				ELSE : 	drho=0.5*(rho[k+1]-rho[k-1])
			ENDCASE
			alpha=a[k]/drho^2
			IF k NE 0 THEN BEGIN
				nzkm=csden[q,k-1]
				beta=b[k]/(2.0*drho)
				gamma=c[k]
			ENDIF ELSE BEGIN
 				nzkm=0.0
				beta=0.0
				gamma=0.0
			ENDELSE
			IF k NE nrho-1 THEN nzkp=csden[q,k+1] ELSE nzkp=0.0
	
			matrix[(q-1)*nrho+k,(q-1)*nrho+k]+=-2.0*alpha+gamma-dens[k]*(iz[k]+rz[k])-neut[k]*cz[k]
			IF k NE 0 THEN matrix[(q-1)*nrho+(k-1),(q-1)*nrho+k]+=alpha-beta
			IF k NE nrho-1 THEN matrix[(q-1)*nrho+(k+1),(q-1)*nrho+k]+=alpha+beta
			IF q NE z THEN matrix[(q-1+1)*nrho+k,(q-1)*nrho+k]+=dens[k]*rzp[k]+neut[k]*czp[k]
			IF q NE 1 THEN matrix[(q-1-1)*nrho+k,(q-1)*nrho+k]+=dens[k]*izm[k]
			tevector[(q-1)*nrho+k]=-1.0*dens[k]*drzp[k]*nzp[k]-dens[k]*dizm[k]*nzm[k]+dens[k]*nz[k]*(diz[k]+drz[k])+neut[k]*nz[k]*dcz[k]-neut[k]*nzp[k]*dczp[k]
			nevector[(q-1)*nrho+k]=-1.0*rzp[k]*nzp[k]-izm[k]*nzm[k]+nz[k]*(iz[k]+rz[k])
			novector[(q-1)*nrho+k]=cz[k]*nz[k]-czp[k]*nzp[k]

			;IF k EQ 0 THEN stop
		ENDFOR
	ENDFOR
	ete=ionrec_coefs.terr
	ene=ionrec_coefs.derr
	eno=ionrec_coefs.nerr
	sens_te=reform(la_invert(matrix,double=double)##tevector)
	sigte=fltarr(z*nrho)
	FOR i=0,z-1 DO sigte[nrho*i:nrho*(i+1)-1]=ete
	sens_ne=reform(la_invert(matrix,double=double)##nevector)
	signe=fltarr(z*nrho)
	FOR i=0,z-1 DO signe[nrho*i:nrho*(i+1)-1]=ene
	sens_no=reform(la_invert(matrix,double=double)##novector)
	signo=fltarr(z*nrho)
	FOR i=0,z-1 DO signo[nrho*i:nrho*(i+1)-1]=eno
	errvec_tot=sqrt(sens_te^2*sigte^2+sens_ne^2*signe^2+sens_no^2*signo^2)
	cserr=csden*0.0
	FOR i=0,z-1 DO cserr[i+1,*]=errvec_tot[i*nrho:(i+1)*nrho-1]

	output=cserr
	;stop
	IF keyword_set(debug) THEN stop
	RETURN,output
END

;+
;NAME:
;	FORM_SOURCE_VECTOR
;
;PURPOSE:
;	This function computes the radial profile of the impurity source
;
;CALLING SEQUENCE:
;	result=FORM_SOURCE_VECTOR()
;
;OPTIONAL INPUTS:
;	ionrec_path	STRING	path to save file where the ion/rec coefficients are stored DEFAULT: '/tmp/ionrec_coefs.dat'
;	t_neut		FLOAT	temperature of the neutral species [eV] DEFAULT=1.0
;	ptsource	FLOAT	of the r/a value to place a point source
;
;KEYWORD PARAMETERS:
;	double		/double forms the source vector as a double preceision array (DEFAULT: single)
;	debug		/debug stops the code before the RETURN statement
;
;OPTIONAL OUTPUTS:
;	no		FLTARR	[nrho] of the impurity neutral density profile normalized to its peak.
;
;PROCEDURE:
;	If the ptsource optional input is not used, then the radial impurity neutral density profile is computed by 
;	launching neutrals from the edge boundary at the thermal velocity of t_neut.
;
;MODIFICATION HISTORY:
;	Written by:	M.L. Reinke - Spring 2009	
;	6/7/12		M.L. Reinke - modified path to use /tmp/ionrec_coefs.dat 
;
;-

FUNCTION form_source_vector,ionrec_path=ionrec_path,debug=debug,double=double,t_neut=t_neut,no=no,ptsource=ptsource
	;IF NOT keyword_set(ionrec_path) THEN ionrec_path='/home/mlreinke/idl/genie/atomic_physics/ionrec_coefs.dat'
	IF NOT keyword_set(ionrec_path) THEN ionrec_path='/tmp/ionrec_coefs.dat'
	restore,ionrec_path
	nrho=n(rho)+1
	IF NOT keyword_set(t_neut) THEN t_neut=1.0			;energy of neutral at edge [eV]
	mass=read_atomic_mass(read_atomic_name(z))
	e=1.60e-19			;conversion for eV -> J
	mconv=1.66e-27			;conversion for amu -> kg
	v_th_neut=sqrt(2.0*t_neut*e/(mass*mconv))
	
	no=dblarr(nrho)
	IF keyword_set(ptsource) THEN BEGIN
		no+=1.0e-20
		no[ipt(rho,ptsource)]=1.0
	ENDIF ELSE BEGIN
		integrand=reform(ionrec_coefs.i[0,*]*ionrec_coefs.d)/v_th_neut*a
		FOR i=0,nrho-2 DO no[i]=1.0/rho[i]*exp(-1.0*double(int_tabulated(rho[i:*],integrand[i:*])))
		no[nrho-1]=1
		no[0]=0
	ENDELSE
	q=no*reform(ionrec_coefs.i[0,*]*ionrec_coefs.d)
	q=q

	IF keyword_set(double) THEN BEGIN
		beta=dblarr(z*nrho) 
		beta[0:nrho-1]=-1.0*(q)
	ENDIF ELSE BEGIN
		beta=fltarr(z*nrho)
		beta[0:nrho-1]=-1.0*float(q)
	ENDELSE
	output=beta
	IF keyword_set(debug) THEN stop
	RETURN,output
END

;+
;NAME:
;	GENTRAN
;	
;PURPOSE:
;	Solves the time-independent impurity transport equation given temperature and denisty data
;	as well transport coefficients.
;
;CALLING SEQUENCE:
;	result=GENTRAN(shot,t,z,diff,conv,temp,dens,neut)
;
;INPUTS:
;	shot:	LONG	shot number to use for EFIT data
;	t:	FLT	of the time point to use for EFIT [sec]
;	diff:	STRUC	diffusion structure
;		*.diff	[nrho]	of the assumed diffusivity [m^2/s]
;		*.rho	[nrho] 	of the r/a values
;	conv:	STRUC	convective structure
;		*.diff	[nrho]	of the assumed convective velocity [m/s]
;		*.rho	[nrho] 	of the r/a values
;	temp:	STRUC	electron temperature structure
;		*.temp	[nrho]	of the electron temperature [eV]
;		*.rho	[nrho]	of the r/a values
;		*.err	[nrho]	of the absolute error [eV]
;	dens:	STRUC	electron density structure
;		*.dens	[nrho]	of the electron density [m^-3]
;		*.rho	[nrho]	of the r/a values
;		*.err	[nrho]	of the absolute error [m^3]
;	neut:	STRUC	neutral density structure
;		*.dens	[nrho]	of the neutral density [m^-3]
;		*.rho	[nrho]	of the r/a values
;		*.err	[nrho]	of the absolute error [m^-3]
;
;OPTIONAL INPUTS:
;	rhomax:	FLT	maximum r/a values DEFAULT: (0.905-ro)/a = limiter radius
;	nrho:	INT	number of rho points for which to calculate solution DEFAULT: 50
;	k:	INT	meshing control, k=1 is even mesh k>1 increases density in edge DEFAULT: 2.0
;	rff:	FLT	"recombination fudge factor" to multiply all the recombination rates by DEFAULT: 1.0
;	iff:	FLT	"ionization fudge factor" to multiply all the ionization rates by DEFAULT: 1.0
;
;KEYWORD PARAMETERS:
;	adas:	/adas sends the command to load ADAS rates
;	double:	/double uses double precision variables 
;	verb:	/verb prints status output to console
;	debug:	/debug stops the code just prior to the RETURN statment
;	nodel:	/nodel will skip the step where the temporary transport/ion_rec savesets are removed
;
;OUTPUTS:
;	result:	STRUC	contains the charge state density information and more
;		*.csden		FLTARR	[z+1,nrho] of the charge state density [i,*] of ionization stage i
;		*.cserr		FLTARR	[z+1,nrho] of the error in the charge state density calculated from temp.err, dens.err and neut.err
;		*.rho		FLTARR	[nrho] of the r/a values of *.csden
;		*.rmaj		FLTARR	[nrho] of the major radius values of *.csden [m]
;		*.source	FLTARR	[nrho] of the normalized impurity source used
;		*.ntot_th	FLTARR	[nrho] of the theoretical total impurity density from integrating the conv/diff profile
;		*.irpath	STR	of the path to the interpolated ion/rec data (output of WRITE_IONREC_COEFS_TABLE)
;		*.tpath		STR	of the path to the interpolated transport data (output of WRITE_TRANSPORT_COEFS_TABLE)
;		*.vol		FLTARR	[nrho] of the volterm calculated using CALC_VOLTERM_PROFILE
;
;MODIFICATION HISTORY:
;	Written by: 	ML Reinke 2/09 (adapated from math in Dux 10/29)
;
;-

FUNCTION gentran,shot,t,z,diff,conv,temp,dens,neut,rhomax=rhomax,nrho=nrho,k=k,rff=rff,iff=iff,adas=adas,debug=debug,verb=verb,double=double,$
		ptsource=ptsource,t_neut=t_neut,nodel=nodel

	;get efit_data
	mdsopen,'analysis',shot
	rmid=mdsvalue('\efit_rmid')
	efit_t=mdsvalue('dim_of(\efit_rmid)')
	mdsclose,'analysis',shot
	efit_i=ipt(efit_t,t)
	ro=rmid[efit_i,0]
	a=max(rmid[efit_i,*])-ro
	minor_rad=a
	major_rad=ro

	IF NOT keyword_set(rff) THEN rff=1.0
	IF NOT keyword_set(iff) THEN iff=1.0
	IF NOT keyword_set(nrho) THEN nrho=50
	rhomax=(0.905-ro)/a
	IF NOT keyword_set(k) THEN k=2.0
	x=make(0.0,rhomax^(1.0/k),nrho)
	IF keyword_set(double) THEN rhovec=reverse(rhomax-double(x)^k) ELSE rhovec=reverse(rhomax-x^k)
	;x=reverse(x)
	;rhovec=x^(1.0/k)

	dxdrho=deriv(rhovec,x)
	d2xdrho2=deriv(rhovec,dxdrho)
	;stop
	volterm=calc_volterm_profile(shot,rhovec,tr=[1,1]*efit_t[efit_i],/circ)


	;load recombination and ionization structures
	IF keyword_set(verb) THEN print, 'Loading ION/REC Data'
	rec=read_rec_data(z,adas=adas,cxr=cxr)
	rec.rates*=rff
	ion=read_ion_data(z,adas=adas)
	ion.rates*=iff

	IF keyword_set(verb) THEN print, 'Interpolating ION/REC Data'
	write_ionrec_coefs_table,rhovec,ion,rec,cxr,temp,dens,neut,minor_rad,path=path_ionrec
	IF keyword_set(verb) THEN print, 'Interpolating Transport Data'
        write_transport_coefs_table,rhovec,x,dxdrho,d2xdrho2,diff,conv,volterm,minor_rad,path=path_trans,ntot_theory=ntot_theory
	
	restore, path_ionrec
	restore, path_trans

	IF keyword_set(verb) THEN print, 'Filling Transport Matrix'
	A=form_transport_matrix(double=double)
	IF keyword_set(verb) THEN print, 'Filling Source Vector'
	b=form_source_vector(double=double,no=no,ptsource=ptsource,t_neut=t_neut)
	IF keyword_set(verb) THEN print, 'Inverting Matrix'
	csvec=la_invert(A,double=double)##b
	IF keyword_set(verb) THEN print, 'Matrix Inverted'
	csden=fltarr(z+1,nrho)
	csden[0,*]=no
	FOR i=0,z-1 DO csden[i+1,*]=csvec[i*nrho:(i+1)*nrho-1]
	ctot=total(csden[*,0])
	csden/=ctot
	
	cserr=calc_csden_error(csden,trans_path=trans_path,ionrec_path=ionrec_path,debug=debug,/double)
	;cserr=-1
	source=b[0:nrho-1]
	source/=min(source)
	output={csden:csden,cserr:cserr,rho:rhovec,rmaj:rho*minor_rad+major_rad,source:source,ntot_th:ntot_theory,irpath:path_ionrec,tpath:path_trans,vol:volterm}
	IF keyword_set(debug) THEN stop
	IF NOT keyword_set(nodel) THEN BEGIN
		spawn, 'rm '+path_ionrec
		spawn, 'rm '+path_trans
	ENDIF
	RETURN,output		
END

;legacy save file format
PRO gentran_te_apipfit,path,temp_str,ph=ph,ave=ave,ete=ete
	IF keyword_set(ete) THEN err_min=30.0/1.0e3 ELSE err_min=0.0
	restore,path
	IF NOT keyword_set(ph) THEN ph=0.0
	ipt=ipt(phase,ph)
	IF ipt[0] EQ -1 THEN BEGIN
		IF ph LT min(phase) THEN ipt=0
		IF ph GT max(phase) THEN ipt=n(ph)
	ENDIF
	temp=temp[ipt,*]
	IF keyword_set(ave) THEN temp=temp_ave
	IF NOT keyword_set(ete) THEN ete=0.0
	IF n(ete) EQ n(temp) THEN err=ete ELSE err=temp*ete[0]
;	tmp=where(err LT err_min)
;	IF tmp[0] NE -1 THEN err[tmp]=err_min
			
	temp_str={temp:temp*1.0e3,rho:rho,err:err*1.0e3}
END

;legacy save file format
PRO gentran_ne_apipfit,path,dens_str,ene=ene
	IF keyword_set(ene) THEN err_min=3.0e18/1.0e20 ELSE err_min=0.0
	restore, path
	IF NOT keyword_set(ene) THEN ene=0.0
	IF n(ene) EQ n(dens) THEN err=ene ELSE err=dens*ene[0]
	tmp=where(err LT err_min)
	IF tmp[0] NE -1 THEN err[tmp]=err_min
	dens_str={dens:dens*1.0e20, rho:rho,err:err*1.0e20}
END

;+
;NAME:
;	GENTRAN_TENE_WIDGETFITS
;
;PURPOSE:
;	This procedure loads a fiTS save file and creates the temperature and
;	density structures for use with GENTRAN
;
;CALLING SEQUENCE:
;	GENTRAN_TENE_WIDGETFITS,path,t1,t2,temp,dens
;INPUTS:
;	path	STRING	path to the fiTS save file
;	t1	FLOAT	of the lower time bound [sec]	
;	t2	FLOAT	of the upper time bound [sec]
;
;OPTIONAL INPUTS:
;	ete	FLTARR	[ntemp] of the electron temperature error [eV].  If set as a FLOAT then
;			it is assumed to specify the fractional error of the entire profile.
;	ene	FLTARR	[ndens] of the electron density error [m^-3].  If set as a FLOAT then
;			it is assumed to specify the fractional error of the entire profile.
;
;OUTPUTS:
;	dens	STRUC	of the electron density data in the GENTRAN structure formation
;	temp	STRUC	of the electron temperature data in the GENTRAN structure formation
;
;MODIFICATION HISTORY:
;	Written by:	M.L. Reinke - Spring 2009
;
;-

PRO gentran_tene_widgetfits,path,t1,t2,temp_str,dens_str,ete=ete,ene=ene
	restore,path

	;get efit_data
	mdsopen,'analysis',shot_number
	rmid=mdsvalue('\efit_rmid')
	efit_t=mdsvalue('dim_of(\efit_rmid)')
	mdsclose,'analysis',shot_number
	i1=ipt(efit_t,t1)
	i2=ipt(efit_t,t2)
	ro=min(sum_array(rmid[i1:i2,*],/j)/(i2-i1+1.0))
	a=max(sum_array(rmid[i1:i2,*],/j)/(i2-i1+1.0))-ro
	rlim=0.905

	i1=ipt(ne_fit.time,t1)
	i2=ipt(ne_fit.time,t2)
	dens=reform(sum_array(ne_fit.combined_fit_ne[i1:i2,*],/j)/(i2-i1+1.0)) 	;m^-3
	rmaj=reform(sum_array(ne_fit.rmajor[i1:i2,*],/j)/(i2-i1+1.0))
	tmp=where(rmaj LT rlim)	
	IF NOT keyword_set(ene) THEN ene=0.0
	IF n(ene) EQ n(dens[tmp]) THEN err=ene ELSE err=dens[tmp]*ene[0]
	dens_str={dens:dens[tmp],rho:(rmaj[tmp]-ro)/a,err:err}

	i1=ipt(te_fit.time,t1)
	i2=ipt(te_fit.time,t2)
	temp=reform(sum_array(te_fit.te_comb_fit[i1:i2,*],/j)/(i2-i1+1.0))*1.0e3 	;eV
	rmaj=reform(sum_array(te_fit.rmajor[i1:i2,*],/j)/(i2-i1+1.0))
	tmp=where(rmaj LT rlim)	
	IF NOT keyword_set(ete) THEN ete=0.0
	IF n(ete) EQ n(temp[tmp]) THEN err=ete ELSE err=temp[tmp]*ete[0]

	temp_str={temp:temp[tmp], rho:(rmaj[tmp]-ro)/a,err:err}
END

;+
;NAME:
;	GENTRAN_TENE_QUICKFITS
;
;PURPOSE:
;	This procedure runs qfit.pro and creates the temperature and density 
;	structures for use with GENTRAN
;
;CALLING SEQUENCE:
;	GENTRAN_TENE_QUICKFITS,shot,t1,t2,temp,dens
;INPUTS:
;	shot	LONG	shot number
;	t	FLTARR	[ntime] of the start time of the bins [sec]
;	dt	FLTARR	[ntime] of the length of the time bins [sec]
;
;OPTIONAL INPUTS:
;	ete	FLTARR	[ntemp] of the electron temperature error [eV].  If set as a FLOAT then
;			it is assumed to specify the fractional error of the entire profile.
;	ene	FLTARR	[ndens] of the electron density error [m^-3].  If set as a FLOAT then
;			it is assumed to specify the fractional error of the entire profile.
;
;OUTPUTS:
;	dens	STRUC	of the electron density data in the GENTRAN structure formation
;	temp	STRUC	of the electron temperature data in the GENTRAN structure formation
;
;MODIFICATION HISTORY:
;	Written by:	M.L. Reinke - 6/7/12
;	7/3/12		M.L. Reinke - modified to use the t and dt FLTARRs for specifying multiple
;				      time points.
;
;-

PRO gentran_tene_quickfits_ng,shot,t,dt,temp_str,dens_str,ete=ete,ene=ene,data=data
	if keyword_set(data) $
	 then begin ; use passed data
	 	dens=data.dens
		temp=data.temp
		rmaj=data.rmaj
		time=data.time
		rho =transpose(data.Rrho)

		rmid=data.rmid
		efit_t=data.efit_t
	 endif $
	 else begin ; load needed data
	 	;qfit,shot,dens,temp,rmaj,time,noets=noets
	 	qfit_ng,shot,dens,temp,rmaj,time,data=data ; ANJ switched to qfit_ng for more robust core fits
		rho =transpose(data.Rrho)

		mdsopen,'analysis',shot
		 rmid=mdsvalue('\efit_rmid')
		 efit_t=mdsvalue('dim_of(\efit_rmid)')
		mdsclose,'analysis',shot
	 endelse

	ntime=n(t)+1		 ;number of time points requested
	nrho=n(dens[*,0])+1	 ;number of radial points in qfit
	xrho =fltarr(nrho,ntime)
	xdens=fltarr(nrho,ntime)
	xderr=fltarr(nrho,ntime)
	xtemp=fltarr(nrho,ntime)
	xterr=fltarr(nrho,ntime)

	IF NOT keyword_set(ene) THEN ene=0.0
	IF NOT keyword_set(ete) THEN ete=0.0
	FOR i=0,ntime-1 DO BEGIN
		i1=ipt(time,t[i])
		i2=ipt(time,t[i]+dt[i])
		xrho[*,i]=make(min(rho[*,i1:i2]),max(rho[*,i1:i2]),nrho)
		for j=i1,i2 do begin
			xdens[*,i] += interpol(dens[*,j]*1.e20,rho[*,j],xrho[*,i]) /(i2-i1+1.0) ; m^-3
			xtemp[*,i] += interpol(temp[*,j]*1.e3 ,rho[*,j],xrho[*,i]) /(i2-i1+1.0) ; eV
		endfor
		IF n(ene) EQ nrho-1 THEN xderr[*,i]=ene ELSE xderr[*,i]=xdens[*,i]*ene[0]
		IF n(ete) EQ nrho-1 THEN xterr[*,i]=ete ELSE xterr[*,i]=xtemp[*,i]*ete[0]
	ENDFOR
	dens_str={dens:xdens,rho:xrho,err:xderr}
	temp_str={temp:xtemp, rho:xrho,err:xterr}
END

;+
;NAME:
;	GENTRAN_PROFILES
;	
;PURPOSE:
;	This upper-level procedure runs GENTRAN for a given shot and impurity averaged
;	over a time window.  Links to GENRAD utitlies are included to compute various
;	radiation profiles.  
;
;CALLING SEQUENCE:
;	GENTRAN_PROFILES,shot,z,t1,t2
;	
;INPUTS:
;	shot	LONG	shot number
;	z	INT	impurity element number
;	t1	FLOAT	lower time point for binning [sec]
;	t2	FLOAT	upper time point for binning [sec]
;
;OPTIONAL INPUTS:
;	k	FLOAT	control of non-uniform radial gridding
;	nrho	INT	number of radial grid points to compute 
;	exp	INT	selects from a list of different types of impurity confinement (see PROCEDURE)
;	nz	FLOAT	of the absolute impurity density to be used in radiation comparisons
;	zeff	FLTARR	[nrad] of the zeff profile to be used when computing continuum (contant Zeff use FLOAT)
;	tht	INT	specifies THACO tree number for HIREXSR comparisons
;	ete	FLTARR	[nrad] specify the Te uncertainty (fractional error use FLOAT)
;	ene	FLTARR	[nrad] specify the ne uncertainty (fractional error use FLOAT)
;	eno	FLTARR	[nrad] specify the n0 uncertainty (fractional error use FLOAT)
;	t_neut	FLOAT	temperature of the impurity neutral species [eV] (see GENTRAN)
;	dff	FLOAT	"fudge-factor" multiplied into the diffusion profile DEFAULT=1.0
;	doff	FLOAT	"offset" added to the diffusion profile DEFAULT=0.001 [m^2/s]
;	cff	FLOAT	"fudge-factor" multiplied into the convection profile DEFAULT=1.0
;	xff	FLOAT	"fudge-factor" multiplied into the C-X rates DEFAULT=1.0, 
;			-1.0 turns off charge exchange by zero'ing the n0 profile
;	rff	FLOAT	"fudge-factor" multiplied into REC rates (all charge states) DEFAULT=1.0
;	iff	FLOAT	"fudge-factor" multiplied into ION rates (all charge states) DEFAULT=1.0
;	tpt	FLOAT
;	xph	FLOAT
;	tph	FLOAT
;
;KEYWORD PARAMETERS:
;	qfit		/qfit runs GENIE/qfit.pro "quick fits" for Te and ne profiles 
;	apip		/apip (legacy) - not for general use
;	adas		/adas loads ADAS atomic physics data
;	ptsource	/ptsource assumes a point impurity source at r/a=1
;	plot		/plot outputs the radial n_e,Te,n_0,n_z0,D,v and csden profiles to the active device
;	debug		/debug stops the code at the end of GENTRAN_PROFILES
;	ave
;
;	GENRAD MODELING KEYWORDS	
;	axuv		/axuv will launch AXUV_GENRAD_PROFILES using the simulation output
;	foil		/foil will launch AXUV_GENRAD_PROFILES using the simulation output
;	hirexsr		/hirexsr will launch HIREXSR_GENRAD_THACO_PROFILES using the simulation output
;	mcp		/mcp will launch MCP_GENRAD_PROFILES using the simulation output
;	xtomo		/xtomo will launch XTOMO_GENRAD_PROFILES using the simulation output
;
;OPTIONAL OUTPUTS:
;	csden	FLTARR	[z+1,nr] of the normalized charge state density profile
;	rhovec	FLTARR	[nr] of the r/a values for csden
;	data	STRUC	input data used in the simulation and results
;		*.shot	LONG	shot number
;		*.time	FLOAT	average time point [sec] 0.5*(t1+t2)
;		*.temp	FLTARR	[nr] electron temperature [keV]
;		*.terr	FLTARR	[nr] of the unc. in temp [keV]
;		*.dens	FLTARR	[nr] of the electron density [m^3]
;		*.derr	FLTARR	[nr] of the unc. in the dens [m^3]
;		*.neut	FLTARR	[nr] of the neutral density [m^3]
;		*.nerr	FLTARR	[nr] of the unc. in neut [m^3]
;		*.rho	FLTARR	[nr] of the r/a values
;		*.rmaj	FLTARR	[nr] of the outboard midplane major radius values [m]
;		*.csden	FLTARR	[z+1,nr] of the normalized charge state density profiles
;		*.cserr	FLTARR	[z+1,nr[ of the unc. in csden
;	outa	STRUC	measured and modeled AXUV emissivity data	
;	outf	STRUC	measured and modeled resistive bolometer emissivity data
;	outh	STRUC	measured and modeled HIREXSR emissivity data
;	outm	STRUC	measured and modeled McPherson emissivity data
;	outx	STRUC	measured and modeled XTOMO emissivity data
;	
;MODIFICATION HISTORY:
;	Written by	M.L. Reinke - Spring 2009
;	6/7/12		M.L. Reinke - modified the density loading
;	6/12/12		M.L. Reinke - added documentation, removed dead I/O flags and added the
;				      AXUV data output (outa)
;-

PRO gentran_profiles,shot,z,t1,t2,fq=fq,verb=verb,exp=exp,k=k,nrho=nrho,cff=cff,dff=dff,qfit=qfit,apip=apip,adas=adas,xph=xph,ave=ave,hirexsr=hirexsr,$
		tph=tph,doff=doff,tpt=t,foil=foil,nz=nz,rff=rff,iff=iff,xff=xff,mcp=mcp,zeff=zeff,axuv=axuv,xtomo=xtomo,data=data,profdata=profdata,outm=outm,outh=outh,csden=csden,$
		rhovec=rhovec,ete=ete,ene=ene,eno=eno,outa=outa,outf=outf,outx=outx,ptsource=ptsource,t_neut=t_neut,debug=debug,plot=plot,tht=tht,noly=noly
	IF NOT keyword_set(dff) THEN dff=1.0
	IF NOT keyword_set(cff) THEN cff=1.0
	IF NOT keyword_set(exp) THEN exp=1
	IF !d.name EQ 'PS' THEN ps=1 ELSE ps=0
	IF NOT keyword_set(nz) THEN nz=1.0e18
	IF NOT keyword_set(doff) THEN doff=0.001
	plotwin=10

	pcase=0
	IF keyword_set(qfit) THEN pcase=1
	IF keyword_set(apip) THEN pcase=2
	;load temperature and density profile
	CASE pcase OF
		0 : BEGIN	;default is to use fiTS savefile from the user's home directory
			path='/home/'+logname()+'/fits/fits_'+num2str(shot,1)+'.save'
			gentran_tene_widgetfits,path,t1,t2,temp_str,dens_str,ete=ete,ene=ene
		END

		1 : BEGIN	;use the "quick-fits"
			gentran_tene_quickfits_ng,shot,t1,t2-t1,temp_str,dens_str,ete=ete,ene=ene,data=profdata
		END

		2 : BEGIN	;legacy - use the saveset data from the 2009 APiP talks
			path='/home/mlreinke/presentations/APiP2009/tefit_'+num2str(shot,1)+'.dat'
			gentran_te_apipfit,path,temp_str,ph=tph,ave=ave,ete=ete
			path='/home/mlreinke/presentations/APiP2009/nefit_'+num2str(shot,1)+'.dat'
			gentran_ne_apipfit,path,dens_str,ene=ene
		END
	ENDCASE

	rhomax=max(temp_str.rho) > max(dens_str.rho)
	neut_str=calc_neutstr(shot,t1,t2-t1,temp_str,dens_str,eno=eno,noly=noly)
	
	;openwin,plotwin
	;plot,neut_str.rho,neut_str.dens,/ylog,xtit='r/a',ytit='Neutral Density [m!u-3!n]',chars=1.2
	;oplot,neut_str.balrho,neut_str.bal,color=100
	;oplot,neut_str.lyrho,neut_str.ly,color=200
	
	IF NOT keyword_set(xff) THEN xff=1.0
	IF xff EQ -1.0 THEN xff=0.0
	neut_str.dens=neut_str.dens*xff
	neut_str.err=neut_str.err*xff

	;setup impurity transport structures diff [m^2/s], conv [m/s]
	CASE exp OF
		1 : BEGIN
			diff=[1.0,1.0,1.0]
			diff_rho=[0,0.5,1.0]*rhomax
			conv=[0.0,-0.0,-0.0]
			conv_rho=[0,0.5,1.0]*rhomax
		END	

		2 : BEGIN
			diff=[1.0,1.0,1.0]
			diff_rho=[0,0.5,1.0]*rhomax
			conv=[0.0,-0.5,-1.0]
			conv_rho=[0,0.5,1.0]*rhomax
		END
		
		3 : BEGIN
			diff_rho=[0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0]*rhomax
			diff=diff_rho^2+doff/dff
			conv=[0.0,-0.0,-0.0]
			conv_rho=[0,0.5,1.0]*rhomax
		END

		4 : BEGIN
			diff_rho=[0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0]*rhomax
			diff=diff_rho^2*dff+doff/dff
			conv=[0.0,-0.5,-1.0]
			conv_rho=[0,0.5,1.0]*rhomax
		END
		5 : BEGIN
			diff=[1.0,1.0,1.0,0.10]
			diff_rho=[0,0.5,0.95,1.0]*rhomax
			conv=[0.0,-0.0,-0.0]
			conv_rho=[0,0.5,1.0]*rhomax
		END	
		6 : BEGIN
			diff=[1.0,1.0,1.0,0.10]
			diff_rho=[0,0.5,0.95,1.0]*rhomax
			conv=[0.0,0.5,-1.0,0.0]
			conv_rho=[0.0,0.4,1.0,rhomax]
		END		
		7 : BEGIN
			diff_rho=[0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0]*rhomax
			diff=diff_rho^4+doff/dff
			conv=[0.0,0.15,-1.0,0.0,0.0]*cff
			conv_rho=[0.0,0.35,0.9,1.0,rhomax]

		END		
		8 : BEGIN
			diff=[1.0,1.0,1.0,1.0,0.10]
			diff_rho=[0,0.5,0.95,1.0,rhomax]
			conv=[0.0,0.0,1.0,-1.0,0.0]
			conv_rho=[0.0,0.2,0.4,1.0,rhomax]
		END	
		9 : BEGIN
			diff=[0.3,0.5,0.75,1.0,1.5,0.10]
			diff_rho=[0,0.25,0.35,0.65,1.0,rhomax]
			conv=[0.0,-6.0,9.0,-0.0,-4.5,0.0]
			conv_rho=[0.0,0.175,0.35,0.6,1.0,rhomax]
		END
	
		10 : BEGIN
			diff=[1.0,1.0,1.0,1.0,1.0]
			diff_rho=[0,0.25,0.35,1.0,rhomax]
			conv=[0.0,0.75,1.5,-1.5,-2.5]
			conv_rho=[0.0,0.1,0.4,0.6,rhomax]
		END
		11 : BEGIN
			a=2.0
			b=0.50
			c=0.15
			d=1.75
			diff=[1.0,1.0,1.0,1.0,1.0]
			diff_rho=[0,0.25,0.35,1.0,rhomax]
			conv_rho=dens_str.rho
			conv=-conv_rho*a+exp(-(conv_rho-b)^2/c^2)*d
		END
		12 : BEGIN
			diff=[0.05,0.05,0.5,0.5,0.5]
			diff_rho=[0,0.3,0.4,1.0,rhomax]
			conv=[0.0,-0.0,-0.0]
			conv_rho=[0,0.5,1.0]*rhomax	
		END
		13 : BEGIN
			xo=0.65
			x1=30.0
			x2=1.75
			diff=[0.05,0.05,0.5,0.5,0.5]
			diff_rho=[temp_str.rho]
			diff=atan((diff_rho-xo)*x1)+x2
			conv=[0.0,-0.0,-0.0]
			conv_rho=[0,0.5,1.0]*rhomax	
		END
	END
	diff*=dff
	conv*=cff
	diff_str={diff:diff,rho:diff_rho}
	conv_str={conv:conv,rho:conv_rho}
	
	IF keyword_set(xph) THEN BEGIN
		hirexsr_sawtooth_phase,shot,t1,t2,phase,time,plot=plot,dph=0.65
		tpt=time[ipt(phase,xph)]
		print, xph
		print, tpt
	ENDIF ELSE tpt=(t1+t2)/2.0
	
	out=gentran(shot,tpt,z,diff_str,conv_str,temp_str,dens_str,neut_str,rhomax=rhomax,nrho=nrho,k=k,rff=rff,iff=iff,adas=adas,/double,ptsource=ptsource,t_neut=t_neut,verb=verb,/nodel)

	csden=out.csden
	cserr=out.cserr
	rhovec=out.rho
	ntot_theory=out.ntot_th
	source=out.source

	restore, out.irpath
	spawn, 'rm '+out.irpath
	restore, out.tpath
	spawn, 'rm '+out.tpath
	temp=ionrec_coefs.t
	temperr=ionrec_coefs.terr
	dens=ionrec_coefs.d
	denserr=ionrec_coefs.derr	
	neut=ionrec_coefs.n
	neuterr=ionrec_coefs.nerr	
	
	IF keyword_set(fq) THEN BEGIN
		frac_abund=gentran_fq(z,temp=temp,adas=adas)
		csden=frac_abund.fq
		cserr=csden*0.0
        ENDIF

	data={shot:shot,time:0.5*(t1+t2),temp:temp,terr:temperr,dens:dens,derr:denserr,neut:neut,nerr:neuterr,rho:rhovec,rmaj:out.rmaj,csden:csden,cserr:cserr}
	IF keyword_set(plot) THEN BEGIN
		;setup profile plotting window
		plotwin+=1
		plotstyle=1
		IF keyword_set(ps) THEN BEGIN
			CASE plotstyle OF 
				1 : BEGIN
					!p.multi=[0,0,4,0,0]
					xsize=7.0
					ysize=7.0*1150/850.0
					ls=1.4

				END
				2 : BEGIN
					!p.multi=[0,2,2,0,0]
					xsize=7.5
					ysize=3.5
					ls=0.5
	
				END
				3 : BEGIN
					!p.multi=0
					ls=0.9
					xsize=7.5
					ysize=3.5
					;pos=[0.1,0.135,0.9,0.95]
				END
			ENDCASE

		ENDIF ELSE BEGIN
			xsize=850.0
			ysize=1150.0
			ls=2.0
			!p.multi=[0,0,4,0,0]
		ENDELSE
		IF NOT keyword_set(ps) THEN BEGIN
			device, window_state=var
			IF var[plotwin] EQ 0 THEN window,plotwin,xsize=xsize,ysize=ysize,xpos=1610,ypos=670,title='output profiles,'+num2str(plotwin) $
				ELSE wset,plotwin
		ENDIF ELSE BEGIN
			d_old=!d
			device, xsize=xsize, ysize=ysize, /inches
		ENDELSE	

		xr=[0,max(rhovec)*1.03]
		qhigh=z
		qlow=0

		;plot the kinetic profile properties
		maxd=max(dens/1.0e20)
		maxt=max(temp/1.0e3)
		pltff=1.1
		IF plotstyle NE 3 THEN tit='SHOT: '+num2str(shot,1)+' '+num2str(t1,dp=3)+' < t < '+num2str(t2,dp=3) ELSE tit =' '
		plot,[0],[0],xr=xr,/xsty,xmargin=[10,8],ystyle=5,chars=1.2*ls,yr=[0,maxd*pltff],tit=tit,xtit='r/a',pos=pos
		oplot, [0.05,0.13],[1.0,1.0]*0.21*maxd
		xyouts,0.15,0.2*maxd,'n!lZ,TOT!n GENTRAN',chars=0.8*ls
		oplot, [0.05,0.13],[1.0,1.0]*0.41*maxd,linestyle=2.0
		xyouts,0.15,0.4*maxd,'n!lZ,TOT!n THEORY',chars=0.8*ls
		oploterror,rhovec,dens/1.0e20,denserr/1.0e20,color=100,psym=-8,symsize=0.5*ls,errcolor=100
		oplot,rhovec,ntot_theory/max(ntot_theory)*maxd,color=0,linestyle=2.0
		oplot,rhovec,sum_array(csden,/j)*maxd/max(sum_array(csden,/j))
		axis,color=100,yaxis=0,ytit='n!le!n [10!u20!n m!u-3!n]',chars=1.2*ls,yr=[0,maxd*pltff],/ysty
		oploterror,rhovec,temp/1.0e3*maxd/maxt,temperr/1.0e3*maxd/maxt,color=200,psym=-8,symsize=0.5*ls,errcolor=200
		axis,color=200,yaxis=1,ytit='T!le!n [keV]',chars=1.2*ls,yr=[0,maxd*pltff]*maxt/maxd,/ysty
	
		;plot, the neutral density and impurity source
		minn=min(neut)
		maxn=max(neut)
		maxs=max(source)
		plot,[0],[0],xr=xr,/xsty,xmargin=[10,8],ystyle=5,chars=1.2*ls,yr=[minn,maxn],/ylog,xtit='r/a',pos=pos
		oploterror,rhovec,neut,neuterr,color=30,psym=-8,symsize=0.5*ls,errcolor=30
		axis,color=30,yaxis=0,ytit='n!lo!n [m!u-3!n]',chars=1.2*ls,/ylog
		oplot,rhovec,source*maxn
		axis,yaxis=1,ytit='Normalized Source',chars=1.2*ls,yr=[minn,maxn]*maxs/maxn,/ysty,/ylog

		;plot the transport properties
		maxd=max(diff_str.diff) 
		maxv=max(conv_str.conv) > 1.0
		minv=min(conv_str.conv) < (-1.0)
		plot,[0],[0],xr=xr,/xsty,xmargin=[10,8],ystyle=5,yr=[0,maxd*1.05],chars=1.2*ls,xtit='r/a',pos=pos
		oplot,diff_str.rho,diff_str.diff,color=100,psym=-8
		;oplot,rhovec,source*maxd
		axis,color=100,yaxis=0,ytit='D [m!u2!n/s]',chars=1.2*ls,yr=[0,maxd*1.05],/ysty
		oplot,conv_str.rho,(conv_str.conv-minv*1.025)*maxd*1.05/(maxv*1.05-minv*1.05),color=200,psym=-8
		axis,color=200,yaxis=1,ytit='v [m/s]',chars=1.2*ls,yr=[0,maxd*1.05]*(maxv*1.05-minv*1.05)/(maxd*1.05)+minv*1.025,/ysty
		;xyouts,0.925,0.1*maxd,'GENTRAN source',orient=90,chars=0.6*ls
		oplot,xr,([0,0]-minv*1.025)*maxd*1.05/(maxv*1.05-minv*1.05),linestyle=1,color=200

		;plot the CS profiles
		plot,[0,0],[0,0],xr=xr,yr=[0,max(csden[qlow:qhigh,*])*1.1],xtit='r/a',ytit='n!lq!n/n!lz!n',tit=tit,/xsty,/ysty,chars=1.2*ls,xmargin=[10,8],pos=pos
		lstylelist=[0,2,3]
		l_cntr=0
		FOR i=qlow,qhigh DO BEGIN
			IF l_cntr EQ n(lstylelist)+1 THEN l_cntr=0
			oploterror, rhovec, csden[i,*],cserr[i,*],linestyle=lstylelist[l_cntr]
			mloc=maxloc(reform(csden[i,*]))
			xyouts,rhovec[mloc],csden[i,mloc]*1.01,num2str(i,1),color=100,charsize=0.8*ls
			l_cntr+=1
	        ENDFOR
		!p.multi=0
	ENDIF

	IF keyword_set(hirexsr) THEN  hirexsr_genrad_thaco_profiles,shot,t1,t2,csden,cserr,temp,temperr,dens,denserr,rhovec,plotwin=plotwin,tht=tht
	IF keyword_set(xtomo) THEN  xtomo_genrad_profiles,shot,t1,t2,csden,cserr,temp,temperr,dens,denserr,rhovec,plotwin=plotwin,t=t,nz=nz,zeff=zeff,out=outx
	IF keyword_set(foil) THEN  foil_genrad_profiles,shot,t1,t2,csden,cserr,temp,temperr,dens,denserr,neut,neuterr,rhovec,plotwin=plotwin,nz=nz,zeff=zeff,out=outf,adas=adas
	IF keyword_set(mcp) THEN mcp_genrad_profiles,shot,t1,t2,csden,cserr,temp,temperr,dens,denserr,rhovec,qmin=qmin,qmax=qmax,w=w,nz=nz,plotwin=plotwin,zeff=zeff,out=outm
	IF keyword_set(axuv) THEN  axuv_genrad_profiles,shot,t1,t2,csden,cserr,temp,temperr,dens,denserr,neut,neuterr,rhovec,plotwin=plotwin,nz=nz,t=t,zeff=zeff,out=outa,foil=outf

	IF keyword_set(ps) AND keyword_set(plot) THEN device, xsize=float(d_old.x_size)/d_old.x_px_cm,ysize=float(d_old.y_size)/d_old.y_px_cm
	IF keyword_set(debug) THEN stop
END

;gentran_profiles,1080130025,18,1.0,1.3,k=3,dff=0.5,/fit,tph=0.922,/hirexsr,/back,xph=0.92,tpt=0.995,/xtomo,/axuv,/foil
;gentran_profiles,1080130025,18,1.0,1.3,k=3,dff=0.6,/fit,tph=0.922,/hirexsr,/back,xph=0.92,tpt=0.995,/xtomo,/axuv,/foil,/mcp,exp=3,doff=0.1
;gentran_profiles,1080116014,36,1.0,1.3,k=3,dff=0.5,/fit,/ave,/xtomo,/axuv,/foil,nz=1.75e17
;gentran_profiles,1080520004,18,1.0,1.2,k=3,dff=0.5,/fit,/ave,nz=3.07e17,zeff=1.6,ete=0.05,/foil,/hirexsr,/xtomo,/back,/foil,/axuv

;gentran_profiles,1080130025,18,1.0,1.3,k=3,dff=0.5,/fit,tph=0.922,xph=0.92,tpt=0.995,/axuv,/foil,zeff=3.0,ete=0.1,ene=0.05,/hirexsr,/back,/xtomo,nz=1.0e18
;gentran_profiles,1080130025,42,1.0,1.3,k=3,dff=0.5,/fit,tph=0.922
;gentran_profiles,1080130025,18,1.0,1.3,k=3,dff=0.5/fit,tph=0.922,xph=0.92,tpt=0.995,/axuv,/foil,zeff=3.0,ete=0.1,ene=0.05,/hirexsr,/back,/xtomo,nz=1.0e18,exp=12
;gentran_profiles,1080520004,18,1.0,1.2,k=3,dff=0.5,/fit,/ave,nz=3.07e17,zeff=1.6,ete=0.1,ene=0.05,/foil,/hirexsr,/xtomo,/back,/axuv
;gentran_profiles,1080130025,18,1.0,1.3,k=3,dff=0.5,/fit,tph=0.922,xph=0.92,tpt=0.995,/foil,zeff=3.0,ete=0.1,ene=0.05,/hirex,/back,nz=1.0e18,/adas

;writes the z-independent GENTRAN data to tree to allow parallel GENTRAN_WRITE2TREE
PRO gentran_prof2tree,shot,t,dt,gyro=gyro,transp=transp,fits=fits,noly=noly
	pcase=0
	IF keyword_set(fits) THEN pcase=1		;load using fiTS
	IF keyword_set(gyro) THEN pcase=2		;load from GYRO tree
	IF keyword_set(transp) THEN pcase=3		;load from TRANSP tree
	;load temperature and density profile
	CASE pcase OF
		0 : BEGIN	;use the "quick-fits" as default
			gentran_tene_quickfits_ng,shot,t,dt,temp,dens,ete=ete,ene=ene
		END

		1 : BEGIN	;use fiTS savefile from the user's home directory or specified directory
			IF size(fits,/type) NE 7 THEN path='/home/'+logname()+'/fits/fits_'+num2str(shot,1)+'.save' ELSE path=fits
			gentran_tene_widgetfits,path,t1,t2,temp_str,dens_str,ete=ete,ene=ene
		END

		2 : BEGIN	

                END

		3 : BEGIN

		END
	ENDCASE
	neut=calc_neutstr(shot,t,dt,temp,dens,eno=eno,noly=noly)

	mdsopen,'spectroscopy',shot
    	mdsput,'\SPECTROSCOPY::TOP.IMPSPEC.PROF:DENS',$
        	   'build_signal(build_with_units($1, "m^-3"),*,'+$
                		'build_with_units($2,"r/a"),'+$
                        	'build_with_units($3,"seconds"),'+$
	                     	'build_with_units($4,"seconds"),'+$
                          	'build_with_units($5,"m^-3"))',$
                   dens.dens,dens.rho,t,dt,dens.err
	mdsopen,'spectroscopy',shot
    	mdsput,'\SPECTROSCOPY::TOP.IMPSPEC.PROF:TEMP',$
        	   'build_signal(build_with_units($1, "eV"),*,'+$
                		'build_with_units($2,"r/a"),'+$
                        	'build_with_units($3,"seconds"),'+$
	                     	'build_with_units($4,"seconds"),'+$
                          	'build_with_units($5,"eV"))',$
                   temp.temp,temp.rho,t,dt,temp.err
    	mdsput,'\SPECTROSCOPY::TOP.IMPSPEC.PROF:NEUT',$
        	   'build_signal(build_with_units($1, "m^-3"),*,'+$
                		'build_with_units($2,"r/a"),'+$
                        	'build_with_units($3,"seconds"),'+$
	                     	'build_with_units($4,"seconds"),'+$
                          	'build_with_units($5,"m^-3"))',$
                   neut.dens,neut.rho,t,dt,neut.err
	mdsclose,'spectroscopy',shot
END

;loads profiles into profiles structures formatted for use with GENTRAN
PRO gentran_load_prof,shot,t,dt,dens,temp,neut
	mdsopen,'spectroscopy',shot
	;load electron density profile
	path='\SPECTROSCOPY::TOP.IMPSPEC.PROF:DENS'
	dens=mdsvalue('_sig='+path)
	rho=mdsvalue('dim_of(_sig,0)')
	t=mdsvalue('dim_of(_sig,1)')
	dt=mdsvalue('dim_of(_sig,2)')
	err=mdsvalue('dim_of(_sig,3)')
	dens={dens:dens,rho:rho,err:err}

	;load electron temperature profile
	path='\SPECTROSCOPY::TOP.IMPSPEC.PROF:TEMP'
	temp=mdsvalue('_sig='+path)
	rho=mdsvalue('dim_of(_sig,0)')
	t=mdsvalue('dim_of(_sig,1)')
	dt=mdsvalue('dim_of(_sig,2)')
	err=mdsvalue('dim_of(_sig,3)')
	temp={temp:temp,rho:rho,err:err}

	;load neutral density profile
	path='\SPECTROSCOPY::TOP.IMPSPEC.PROF:NEUT'
	neut=mdsvalue('_sig='+path)
	rho=mdsvalue('dim_of(_sig,0)')
	t=mdsvalue('dim_of(_sig,1)')
	dt=mdsvalue('dim_of(_sig,2)')
	err=mdsvalue('dim_of(_sig,3)')
	neut={dens:neut,rho:rho,err:err}			
END

PRO gentran_write2tree,shot,z,diff=diff,conv=conv
	gentran_load_prof,shot,t,dt,dens,temp,neut
	ntime=n(t)+1

	;set diffusion structure
	IF NOT keyword_set(diff) THEN BEGIN				;set default 0.5 m^2/s constant diff for all times
		idiff=[1.0,1.0,1.0]*0.5
		irho=[0,0.5,1.0]
		npts=n(idiff)+1
		diff=fltarr(npts,ntime)
		rho=fltarr(npts,ntime)
		FOR i=0,ntime-1 DO BEGIN
			rhomax=max(dens.rho[*,i])
			diff[*,i]=idiff
			rho[*,i]=irho*rhomax
		ENDFOR
        ENDIF ELSE BEGIN 
		ndiff=n(diff.diff[0,*])+1
		IF ndiff NE ntime AND ndiff EQ 1 THEN BEGIN		;if a single diff profile input, assume constant
			idiff=diff.diff
			irho=diff.rho
			npts=n(idiff)+1
			diff=fltarr(npts,ntime)
			rho=fltarr(npts,ntime)
			FOR i=0,ntime-1 DO BEGIN
				diff[*,i]=idiff
				rho[*,i]=irho
			ENDFOR
                ENDIF ELSE BEGIN					;if multiple diff profiles, must be same ntime as PROF data
			print, 'number of DIFF time slices different then PROF'
			RETURN
		ENDELSE
        ENDELSE
	diff={diff:diff,rho:rho}
	;define convective velocity structure
	IF NOT keyword_set(conv) THEN BEGIN				;set default 0.0 m/s constant conv for all times
		iconv=[0.0,-0.0,-0.0]
		irho=[0,0.5,1.0]
		npts=n(iconv)+1
		conv=fltarr(npts,ntime)
		rho=fltarr(npts,ntime)
		FOR i=0,ntime-1 DO BEGIN
			rhomax=max(dens.rho[*,i])
			conv[*,i]=iconv
			rho[*,i]=irho*rhomax
		ENDFOR
        ENDIF ELSE BEGIN 
		nconv=n(conv.conv[0,*])+1
		IF nconv NE ntime AND nconv EQ 1 THEN BEGIN		;if a single conv profile input, assume constant
			iconv=conv.conv
			irho=conv.rho
			npts=n(iconv)+1
			conv=fltarr(npts,ntime)
			rho=fltarr(npts,ntime)
			FOR i=0,ntime-1 DO BEGIN
				conv[*,i]=iconv
				rho[*,i]=irho
			ENDFOR
                ENDIF ELSE BEGIN					;if multiple conv profiles, must be same ntime as PROF data
			print, 'number of CONV time slices different then PROF'
			RETURN
		ENDELSE
	ENDELSE
	conv={conv:conv,rho:rho}
	

	FOR i=0,ntime-1 DO BEGIN
		print, 'SHOT: '+num2str(shot,1)+' Z='+num2str(z,1)+': inverting time point '+num2str(i,1)+' of '+num2str(ntime-1,1)
		idiff={diff:diff.diff[*,i],rho:diff.rho[*,i]}
		iconv={conv:conv.conv[*,i],rho:conv.rho[*,i]}
		itemp={temp:temp.temp[*,i],rho:temp.rho[*,i],err:temp.err[*,i]}
		idens={dens:dens.dens[*,i],rho:dens.rho[*,i],err:dens.err[*,i]}
		ineut={dens:neut.dens[*,i],rho:neut.rho[*,i],err:neut.err[*,i]}
		
		out=gentran(shot,t[i]+0.5*dt[i],z,idiff,iconv,itemp,idens,ineut,rhomax=rhomax,nrho=nrho,k=k,rff=rff,iff=iff,adas=adas,/double,ptsource=ptsource,t_neut=t_neut,verb=verb,/nodel)
		IF i EQ 0 THEN BEGIN
			csden=fltarr(z+1,nrho,ntime)
			cserr=fltarr(z+1,nrho,ntime)
			csrho=fltarr(nrho,ntime)
		ENDIF
		csden[*,*,i]=out.csden
		cserr[*,*,i]=out.cserr
		csrho[*,i]=out.rho
        ENDFOR
	stop
END

