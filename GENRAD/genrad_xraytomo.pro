;+
;NAME:
;	XTOMO_CALC_CSPLC
;
;PURPOSE:
;	This function computes/loads the charge-state resolved photon
;	emission coefficients as observed by the C-Mod XTOMO system (50
;	micron Be) for different high-Z impurities.
;
;CALLING SEQUENCE:
;	result=XTOMO_CALC_CSPLC(z)
;
;INPUTS:
;	z	INT atomic number of the element (currently only 18,36)
;
;KEYWORD PARAMETERS:
;	load	/load loads an IDL saveset 
;	adas	/load loads the data computed using ADAS 
;	debug	/debug stops the analysis before the end
;	verb	/verb is passed to the PLC creation functions
;
;PROCEDURE:
;	The 50 micron Be filter transmission is computed from data using:
;		http://henke.lbl.gov/optical_constants/filter2.html
;	loaded via READ_XRAY_DATA.
;
;	If /load is not specificed then ADAS_PEC_CSPLC or
;	HULLAC_CS_PLC are called to compute from available atomic phyiscs
;	data and save the file in ~/atomic_physics/usr/+logname()+/ which
;	assumes the user has the atomic_physics Mercurial database
;	installed.  Please talk w/ M.L. Reinke before establishing new
;	default radiation modeling data.
;
;MODIFICATY HISTORY:
;	Written by:	M.L. Reinke - Spring 2009
;	6/8/12		M.L. Reinke - modified to use new /usr/local/cmod/idl/atomic_physics/ data 
;	6/12/12		M.L. Reinke - fixed error in path to SXR filter data
;
;-

FUNCTION xtomo_calc_csplc,z,debug=debug,load=load,verb=verb,adas=adas

	t_e=10^(make(1.0,4.0,500))
	n_e=[1.0e20]
	CASE z OF
		18 : BEGIN
			;path_save='/home/mlreinke/idl/impurities/data/hullac_kbf/ar/ar_xraytomo_csplc.dat'
			;IF keyword_set(adas) THEN path_save='/home/mlreinke/idl/impurities/data/adas/ar_xraytomo_csplc.dat'
			;path_filter='/home/mlreinke/idl/impurities/data/sxr_filters/be50.dat'
			path_load='/usr/local/cmod/idl/atomic_physics/hullac_kbf/ar/ar_xraytomo_csplc.dat'
			IF keyword_set(adas) THEN path_load='/usr/local/cmod/idl/atomic_physics/adas/ar_xraytomo_csplc.dat'
			path_save='/home/'+logname()+'atomic_physics/usr/'+logname()+'/ar_xraytomo_csplc.dat'
			path_filter='/usr/local/cmod/idl/atomic_physics/sxr_filters/be50.dat'
		END
		36 : BEGIN
			;path_save='/home/mlreinke/idl/impurities/data/hullac_kbf/kr/kr_xraytomo_csplc.dat'
			;IF keyword_set(adas) THEN path_save='/home/mlreinke/idl/impurities/data/adas/kr_xraytomo_csplc.dat'
			;path_filter='/home/mlreinke/idl/impurities/data/sxr_filters/be50.dat'
			path_load='/usr/local/cmod/idl/atomic_physics/hullac_kbf/kr/kr_xraytomo_csplc.dat'
			IF keyword_set(adas) THEN path_load='/usr/local/cmod/idl/atomic_physics/adas/kr_xraytomo_csplc.dat'
			path_save='/home/'+logname()+'atomic_physics/usr/'+logname()+'/kr_xraytomo_csplc.dat'
			path_filter='/usr/local/cmod/idl/atomic_physics/sxr_filters/be50.dat'
		END
		ELSE : RETURN,-1
		
	ENDCASE
	IF NOT keyword_set(load) THEN BEGIN
		filter=read_xray_data(path_filter)
		IF keyword_set(adas) THEN csplc=adas_pec_csplc(z,t_e,n_e,filter=filter,verb=verb) ELSE csplc=hullac_cs_plc(t_e,z,filter=filter,n_e=n_e[0],verb=verb)
		save,csplc,t_e,n_e,filename=path_save
	ENDIF ELSE restore,path_load

	output={csplc:csplc,temp:t_e,dens:n_e}
	IF keyword_set(debug) THEN stop
	RETURN,output
END

;+
;NAME:
;	XTOMO_GENRAD_BREM_EMISS
;
;PURPOSE:
;	This function calculate the continuum emissivity seen by the
;	x-ray tomography system and line-integrated brightness for selected arrays
;	
;CALLING SEQUENCE:
;	result=XTOMO_GENRAD_BREM_EMISS(shot,t,rmaj,temp,temperr,dens,denserr)
;
;INPUTS:
;	shot	LONG	of the shot number
;	t	FLOAT	of the time point to use equilibrium reconstruction
;	rmaj	FLTARR	[nr] of the outboard midplane major radius [m]
;	temp	FLTARR	[nr] electron temperature [keV]
;	temperr	FLTARR	[nr] of the unc. in temp [keV]
;	dens	FLTARR	[nr] of the electron density [m^-3]
;	denserr	FLTARR	[nr] of the unc. in the dens [m^-3]	
;
;OPTIONAL INPUTS:
;	n_s	INT	number of points along line of sight to compute integral over DEFAULT set in LINE_BR
;	zeff	FLTARR	[nr] of the Zeff profile DEFAULT: 1.0 (constant,pure)
;	pos	FLTARR	[4,nch] of the pos DEFAULT shot specific array1, array3 pos arrays from XRAY_TOMO_POS
;	
;KEYWORD PARAMETERS:
;	array5	/array5 will include the 5th (FY12+) array along with
;		array1 and array3 when computing the brightness
;	debug	/debug will stop the code before the RETURN statement
;
;OUTPUTS:
;	result	STRUC	of the computed emissivity and brightness data
;		*.rmaj	[nr] duplication of the input major radius scale
;		*.emiss	[nr] of the emissivity [W/m^3] on the major radius grid
;		*.emerr	[nr] of the uncertainty in the emiss [W/m^3]
;		*.pos	[4,nch] of the POS vectors for the selected [DEFAULT: array1, array3]	
;		*.br	[nch] of the brightess seen by the XTOMO system due to this emissivity [W/m^2]
;		*.brerr	[nch] of the uncertanty in the brightness [W/m^2]
;
;PROCEDURE:
;	The Mewe gaunt factor parameterization is used and the Hutchinson formula (W/m^3/eV) w/o 
;	recombination edges is used for the continuum.
;
;MODIFICATION HISTORY:
;	Written by:	M.L. Reinke - Spring 2009
;	6/2012		M.L. Reinke - modified the filter path to use /usr/local/cmod/idl/atomic_physics
;				      added the ability to compute brightness for array5
;       6/2012          ctyler      - Corrected log_gff to multiply
;                                     temp[i] by 1000 to turn KeV to eV
;-	

FUNCTION xtomo_genrad_brem_emiss,shot,t,rmaj,temp,temperr,dens,denserr,n_s=n_s,debug=debug,zeff=zeff,pos=pos,array5=array5
	IF NOT keyword_set(zeff) THEN zeff=1.0
	;path_filter='/home/mlreinke/idl/impurities/data/sxr_filters/be50.dat'
	path_filter='/usr/local/cmod/idl/atomic_physics/sxr_filters/be50.dat'
	filter=read_xray_data(path_filter)
	eph=filter.e
	lam=ev2ang(eph)
	npts=n(rmaj)+1		;number of radial points
	emiss=fltarr(npts)
	emerr=fltarr(npts)
	h=6.626068e-34	;m^2 kg/s
	c=2.9979e8	;m/s
	e=1.60218e-19	;J/eV
	kb=8.617e-5 	;eV/K
	
	FOR i=0,npts-1 DO BEGIN
		log_gff=0.355*lam^(-0.06)*alog10(lam)+0.3*lam^(-0.066)*alog10(temp[i]*1000/kb*1.0e-6/100.0)+0.0043
		specem=5.03e-14/h*e*4.0*!pi*10^(log_gff)*zeff*dens[i]/1.0e20*dens[i]/1.0e20/sqrt(temp[i])*exp(-eph/temp[i])	;watts/m^3/eV from Hutch
		emiss[i]=int_tabulated(eph,filter.tr*specem)	;W/m^3
		emerr[i]=sqrt((2.0*emiss[i]/dens[i])^2*denserr[i]^2+(0.5*emiss[i]/temp[i])^2*temperr[i]^2)
	ENDFOR
	
	IF NOT keyword_set(pos) THEN BEGIN
		pos1=xray_tomo_pos(array=1,/reform,shot=shot)
		pos3=xray_tomo_pos(array=3,shot=shot)
		IF keyword_set(array5) THEN BEGIN
			pos5=xray_tomo_pos(array=5,shot=shot)
			pos=[[pos1],[pos3],[pos5]]
                ENDIF ELSE pos=[[pos1],[pos3]]
	ENDIF
	br=line_br(pos,emiss,rmaj,[t],shot,t,n_s=n_s,plots=debug)
	brerr=line_br(pos,emerr,rmaj,[t],shot,t,n_s=n_s,plots=debug)

	output={rmaj:rmaj,em:emiss,emerr:emerr,pos:pos,br:br,brerr:brerr}
	IF keyword_set(debug) THEN stop
	RETURN,output
END

;+
;NAME:
;	XTOMO_GENRAD_EMISS
;
;PURPOSE:
;	This function computes the emissivity and brightness profiles seen by the 
;	soft x-ray tomography arrays for a specified impurity 
;
;CALLING SEQUENCE:
;	result=XTOMO_GENRAD_EMISS(shot,t,rmaj,csden,cserr,temp,temperr,dens,denserr)
;
;INPUTS:
;	shot	LONG	of the shot number
;	t	FLOAT	of the time point to use equilibrium reconstruction
;	rmaj	FLTARR	[nr] of the outboard midplane major radius [m]
;	csden	FLTARR	[z+1,nr] of the charge state density profiles [m^-3]
;	cserr	FLTARR	[z+1,nr[ of the unc. in csden [m^-3]
;	temp	FLTARR	[nr] electron temperature [keV]
;	temperr	FLTARR	[nr] of the unc. in temp [keV]
;	dens	FLTARR	[nr] of the electron density [m^-3]
;	denserr	FLTARR	[nr] of the unc. in the dens [m^-3]	
;
;OPTIONAL INPUTS:
;	n_s	INT	number of points along line of sight to compute integral over DEFAULT set in LINE_BR
;	pos	FLTARR	[4,nch] of the pos DEFAULT shot specific array1, array3 pos arrays from XRAY_TOMO_POS
;
;KEYWORD PARAMETERS:
;	array5	/array5 will include the 5th (FY12+) array along with
;		array1 and array3 when computing the brightness
;	adas	/adas will use ADAS rates in XTOMO_CALC_CSPLC rather than HULLAC data
;	debug	/debug will stop the code before the RETURN statement
;
;OUTPUTS:
;	result	STRUC	of the computed emissivity and brightness data
;		*.rmaj	[nr] duplication of the input major radius scale
;		*.emiss	[nr] of the emissivity [W/m^3] on the major radius grid
;		*.emerr	[nr] of the uncertainty in the emiss [W/m^3]
;		*.pos	[4,nch] of the POS vectors for the selected [DEFAULT: array1, array3]	
;		*.br	[nch] of the brightess seen by the XTOMO system due to this emissivity [W/m^2]
;		*.brerr	[nch] of the uncertanty in the brightness [W/m^2]
;
;MODIFICATION HISTORY:
;	Written by:	M.L. Reinke - Spring 2009
;	6/12/2012	M.L. Reinke - modified the filter path to use /usr/local/cmod/idl/atomic_physics
;				      added the ability to compute brightness for array5
;
;-

FUNCTION xtomo_genrad_emiss,shot,t,rmaj,csden,cserr,temp,temperr,dens,denserr,n_s=n_s,debug=debug,adas=adas,pos=pos,array5=array5

	z=n(csden[*,0])
	nr=n(rmaj)+1
	csplc=xtomo_calc_csplc(z,/load,adas=adas)
	emiss=fltarr(nr)
	emerr=fltarr(nr)
	FOR i=0,z DO BEGIN
		csint=interpol(csplc.csplc[i,*],csplc.temp,temp)
		dcsint=interpol(deriv(csplc.temp,csplc.csplc[i,*]),csplc.temp,temp)
		emiss+=csint*dens*csden[i,*]
		emerr+=(dens*csden[i,*]*dcsint)^2*temperr^2+(csint*dens)^2*(cserr[i,*])^2+(csint*csden[i,*])^2*denserr^2
	ENDFOR
	emerr=sqrt(emerr)

	IF NOT keyword_set(pos) THEN BEGIN
		pos1=xray_tomo_pos(array=1,/reform,shot=shot)
		pos3=xray_tomo_pos(array=3,shot=shot)
		IF keyword_set(array5) THEN BEGIN
			pos5=xray_tomo_pos(array=5,shot=shot)
			pos=[[pos1],[pos3],[pos5]]
                ENDIF ELSE pos=[[pos1],[pos3]]
	ENDIF	

	br=line_br(pos,emiss,rmaj,[t],shot,t,n_s=n_s,plots=debug)
	brerr=line_br(pos,emerr,rmaj,[t],shot,t,n_s=n_s,plots=debug)
	
	output={rmaj:rmaj,em:emiss,emerr:emerr,pos:pos,br:br,brerr:brerr}
	IF keyword_set(debug) THEN stop
	RETURN,output
END

;+
;NAME:
;	XTOMO_GENRAD_PROFILES
;
;PURPOSE:
;	This procedure compares measured and modeled emissivity and
;	brightness profiles for the soft x-ray tomography diagnostics
;
;CALLING SEQUENCE:
;	XTOMO_GENRAD_PROFILES,shot,t1,t2,csden,cserr,temp,temperr,dens,denserr,rhovec
;
;INPUTS:
;	shot	LONG	shot number
;	t1	FLOAT	lower time point of averaging window [sec]
;	t2	FLOAT	upper time point of averaging window [sec]
;	csden	FLTARR	[z+1,nr] of the normalized charge state density profiles
;	cserr	FLTARR	[z+1,nr[ of the unc. in csden
;	temp	FLTARR	[nr] electron temperature [keV]
;	temperr	FLTARR	[nr] of the unc. in temp [keV]
;	dens	FLTARR	[nr] of the electron density [m^-3]
;	denserr	FLTARR	[nr] of the unc. in the dens [m^-3]
;	rhovec	FLTARR	[nr] r/a 
;
;OPTIONAL INPUTS:
;	t	FLOAT	of the nearest time point to use rather then averaging over t1 < t < t2
;	nz	FLOAT	of the absolute impurity denisty DEFAULT=10^17 m^-3
;	zeff	FLOAT	of the zeff NOT including the specified impurity DEFAULT=1.0
;
;OUTPUTS:
;	GRAPHICAL to plotwin+1 of the measured and modeled emissivity and brightness profiles
;
;OPTIONAL OUTPUTS:
;	out	STRUC		structure containing the measured and modeled XTOMO data
;		*.xem		experimental emissivity profile [kW/m^3]
;		*.xemerr 	unc. in experimental emissivity [kW/m^3]
;		*.em		modeled emissivity profile [kW/m^3]	
;		*.emerr		unc. in modeled emissivity [kW/m^3]
;		*.rho		normalized minor radius r/a 
;
;MODIFICATION HISTORY:
;	Written by:	M.L. Reinke - Spring 2009
;
;-

PRO xtomo_genrad_profiles,shot,t1,t2,csden,cserr,temp,temperr,dens,denserr,rhovec,plotwin=plotwin,t=t,nz=nz,zeff=zeff,out=out
	IF !d.name EQ 'PS' THEN ps=1 ELSE ps=0
	IF NOT keyword_set(nz) THEN nz=1.0e17
	IF NOT keyword_set(plotwin) THEN plotwin=0

	;get efit_data
	mdsopen,'analysis',shot
	rmid=mdsvalue('\efit_rmid')
	efit_t=mdsvalue('dim_of(\efit_rmid)')
	mdsclose,'analysis',shot
	i1=ipt(efit_t,t1)
	i2=ipt(efit_t,t2)
	ro=min(sum_array(rmid[i1:i2,*],/j)/(i2-i1+1.0))
	a=max(sum_array(rmid[i1:i2,*],/j)/(i2-i1+1.0))-ro

	mdsopen,'xtomo',shot		
	xem=mdsvalue('\XTOMO::TOP.RESULTS.CORE:EMISS',/quiet,status=status)
	xrho=mdsvalue('dim_of(\XTOMO::TOP.RESULTS.CORE:EMISS,0)',/quiet)
	xt=mdsvalue('dim_of(\XTOMO::TOP.RESULTS.CORE:EMISS,1)',/quiet)
	xbrchk=mdsvalue('\XTOMO::TOP.RESULTS.CORE:BRCHK',/quiet)
	xbr=mdsvalue('\XTOMO::TOP.RESULTS.CORE:BRIGHT',/quiet,status=status)
	xgood=mdsvalue('\XTOMO::TOP.RESULTS.CORE.CONFIG:GOOD',/quiet)
	mdsclose,'xtomo',shot
	IF keyword_set(t) THEN BEGIN
		ixtomo=ipt(xt,t)
		xem=xem[*,ixtomo]
		xbr=xbr[*,ixtomo]
		xbrchk=xbrchk[*,ixtomo]	
	ENDIF ELSE BEGIN
		i1=ipt(xt,t1)
		i2=ipt(xt,t2)
		xem=sum_array(xem[*,i1:i2],/i)/(i2-i1+1.0)
		xbr=sum_array(xbr[*,i1:i2],/i)/(i2-i1+1.0)
		xbrchk=sum_array(xbrchk[*,i1:i2],/i)/(i2-i1+1.0)
	ENDELSE	
	xtmp=where(xgood EQ 1)
	xem*=1.0e-3
	xemerr=fltarr(n(xem)+1)+10.0		;from up/down error
	xbr*=1.0e-3
	xbrchk*=1.0e-3

	IF NOT keyword_set(t) THEN tpt=0.5*(t1+t2) ELSE tpt=t

	line_emiss=xtomo_genrad_emiss(shot,tpt,rhovec*a+ro,csden*nz,cserr*nz,temp,temperr,dens,denserr,n_s=n_s,debug=debug)
	brem_emiss=xtomo_genrad_brem_emiss(shot,tpt,rhovec*a+ro,temp,temperr,dens,denserr,n_s=n_s,zeff=zeff)
	emiss_line=line_emiss.em*1.0e-3
	emiss_brem=brem_emiss.em*1.0e-3
	emerr_line=line_emiss.emerr*1.0e-3
	emerr_brem=brem_emiss.emerr*1.0e-3
	br_line=line_emiss.br*1.0e-3
	br_brem=brem_emiss.br*1.0e-3
	brerr_line=line_emiss.brerr*1.0e-3
	brerr_brem=brem_emiss.brerr*1.0e-3
	brtot=br_line+br_brem
	brerr_tot=sqrt(brerr_line^2+brerr_brem^2)
	emtot=emiss_line+emiss_brem
	emerr_tot=sqrt(emerr_line^2+emerr_brem^2)
		
	out={xem:xem,xrho:xrho,xemerr:xemerr,em:emtot,rho:rhovec,emerr:emerr_tot}
	chnorm=[12,13]		

	plotwin+=1
	IF keyword_set(ps) THEN BEGIN
		xsize=7.0
		ysize=7.0*800/1400.0
		ls=1.5
	ENDIF ELSE BEGIN
		xsize=1400.0
		ysize=800.0
		ls=2.0
	ENDELSE
	IF NOT keyword_set(ps) THEN BEGIN
		device, window_state=var
		IF var[plotwin] EQ 0 THEN window,plotwin,xsize=xsize,ysize=ysize,xpos=1610,ypos=670,title='output profiles,'+num2str(plotwin) $
			ELSE wset,plotwin
	ENDIF ELSE BEGIN
		d_old=!d
		device, xsize=xsize, ysize=ysize, /inches
	ENDELSE
	pos1=[0.075,0.3,0.6,0.95]
	pos2=[0.075,0.1,0.6,0.3]
	pos3=[0.7,0.3,0.975,0.975]
	pos4=[0.7,0.1,0.975,0.3]

	pos1=[0.075,0.1,0.6,0.95]
	
	;emissivity comparison
	maxpt=max(emtot) > max(xem)
	tit='SXR Tomography'
	lab=num2str(shot,1)
	IF keyword_set(t) THEN lab=lab+' t='+num2str(t,dp=2) ELSE lab=lab+' '+num2str(t1,dp=2)+' < t < '+num2str(t2,dp=2)
	nztit='n!lz!n = '+num2str(nz,dp=2)+' [m!u-3!n]'
	plot, [0],[0],pos=pos1,yr=[0,maxpt*1.05],ytit='Emissivity [kW/m!u3!n]',/ysty,xr=[0,max(rhovec)],/xsty,chars=0.5*ls,$
		tit=tit,xtit='r/a'
	oploterror,rhovec,emtot,emerr_tot,color=200,errcolor=200
	oploterror,rhovec,emiss_line,emerr_line,color=30,errcolor=30
	oploterror,rhovec,emiss_brem,emerr_brem,color=30,linestyle=2.0,errcolor=30
	oploterror,xrho,xem,xemerr
	xyouts,0.75,maxpt*0.9,nztit,chars=0.5*ls

	;emissivity residual
;	residual=(xem-interpol(emtot,rhovec,xrho))*1.0e-3
;	max=max(residual) > 0
;	min=min(residual) < 0
;	plot, [0],[0],pos=pos2,/noerase,yr=[min,max]*1.05,xtit='r/a',ytit=n2g('Delta'),xr=[0,max(rhovec)],/xsty,/ysty,chars=0.5*ls
;	oplot,xrho,residual,psym=8,symsize=0.5*ls
;	oplot,[0,1],[0,0],linestyle=2.0

	;brightness
	ch=indgen(n(brtot)+1)
	tmp=where(xgood EQ 1)
	maxpt=max(xbr[tmp])
	plot, [0],[0],pos=pos3,/noerase,yr=[0,maxpt*1.05],/ysty,xr=[0,max(ch)],/xsty,ytit='Brightness [kW/m!u2!n]',chars=0.5*ls
	makesym,10
	oplot,ch[tmp],xbr,psym=8,symsize=0.5*ls
	oplot,ch[tmp],xbrchk,color=100
	oploterror,ch,brtot,brerr_tot,color=200,errcolor=200

	;brightness residual
	fit_residual=xbr-xbrchk
	chk_residual=xbr-brtot[tmp]
	max=max(fit_residual) > 0 > max(chk_residual)
	min=min(fit_residual) < 0 < min(chk_residual)
	plot, [0],[0],pos=pos4,/noerase,xr=[0,max(ch)],/xsty,yr=[min,max]*1.05,/ysty,xtit='CH #',ytit=n2g('Delta'),chars=0.5*ls
	oplot,ch[tmp],fit_residual,psym=8,color=100,symsize=0.5*ls
	oploterror,ch[tmp],chk_residual,brerr_tot[tmp],psym=8,color=200,errcolor=200,symsize=0.5*ls
	oplot,[0,100],[0,0],linestyle=2.0
	xyouts,78,min*0.5,lab,orient=90,chars=0.3*ls
	IF keyword_set(ps) THEN device, xsize=float(d_old.x_size)/d_old.x_px_cm,ysize=float(d_old.y_size)/d_old.y_px_cm
END

