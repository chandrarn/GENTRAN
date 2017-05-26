;+
;NAME:
;	MAKE_CMOD_PP
;
;PURPOSE:
;	This procedure is used to generate the spatial profile input
;	file for use with the STRAHL impurity transport code
;
;CALLING SEQUENCE:
;	MAKE_CMOD_PP,shot,time,dt (must specify either /fits or /qfit)
;
;INPUTS:
;	shot	LONG	shot number
;	time	FLTARR	of the time points (lower bound)
;	dt	FLTARR	of the (boxcar) averaging interval
;	*** This makes the center time points time+dt/2.0 ***
;
;OPTIONAL INPUTS:
;	ln_sol	FLOAT	the n_e scale length [cm] for the SOL (DEFAULT uses specified radial profile)
;	lt_sol	FLOAT	the T_e scale length [cm] for the SOL (DEFAULT uses specified radial profile)
;	rho	FLTARR	of the r/a grid points DEFAULT is 100 pts (evenly spaced) between 0 < r/a 1.0 
;	nfrac	FLOAT	sets the neutral density to be nfrac*elec. density
;	trshot	LONG	loads the neutral density profile from TRANSP run under trshot
;	fte	FLOAT	'fudge factor' multipled into the Te profile, used for sensitivity analysis
;	fne	FLOAT	'fudge factor' multipled into the ne profile, used for sensitivity analysis
;	fn0	FLOAT	'fudge factor' multipled into the n0 profile, used for sensitivity analysis
;	temp	STRUC	of the time-evolving electron temperature
;		*.temp	FLTARR 	[xrho,xtime] of the electron temperature [eV]
;		*.rho	FLTARR 	[xrho] of the radial grid [norm. pol. flux]
;		*.time	FLTARR	[xtime] of the time points [sec]
; 	dens	STRUC	of the time-evolving electron temperature
;		*.temp	FLTARR 	[yrho,ytime] of the electron density
;		*.rho	FLTARR 	[yrho] of the radial grid [norm. pol. flux]
;		*.time	FLTARR	[ytime] of the time points [sec]
; 	
;KEYWORD PARAMETERS
;	qfit	/qfit will use the "quick_fits" program in /usr/local/cmod/idl/ to compute the T_e and n_e profiles
;	gpfit	/gpfit will load using GPFIT results for T_e and n_e from /home/USERNAME/gpfit/gpfit_XX_shot_0
;		STRING will look for gpfit_XX_SHOT_0 in path given by gpfit STRING
;	fits	/fits will look in the usrs ~/fits/ directory for fits_SHOTNUMBER.save
;	
;OUTPUTS:
;	The file "ppSHOTNUMBER.0" is written in the ~/strahl/cmod/ directory
;
;OPTIONAL OUTPUTS:
;	filepath	STRING of the path of the output file 
;
;PROCEDURE:
;	This procedure uses the GENTRAN codes GENTRAN_TENE_WIDGETFITS (/fits), GENTRAN_TENE_GPFIT or GENTRAN_TENE_QUICKFITS (/qfit) to produce smooth
;	te & ne profiles for use in STRAHL.  The ln_sol and lt_sol will toggle the 'interp' and 'interpa' flags as per instructions in the
;	STRAHL manual.  GENTRAN code GENTRAN_NEUT_TRANSP specifies the neutral density profile if given TRSHOT
;
;MODIFICATION HISTORY:
;	Written by:	M.L. Reinke - January 2013 (adapted from N.T. Howards's WRITE_STRAHL_FILES.PRO)
;	M.L. Reinke	12/2/2013 - added TRSHOT and fXX optional inputs to specify neutral density profile from TRANSP
;	M.L. Reinke	9/14/2014 - added the /gpfit keyword
;       M.L. Reinke	10/12/2014 - modified use of gpfit keyword to allow it carry path information
;	M.L. Reinke	11/7/14 - added the ability to specify TEMP and DENS structures.
;
;-

PRO make_cmod_pp,shot,time,dt,ln_sol=ln_sol,lt_sol=lt_sol,nfrac=nfrac,trshot=trshot,qfit=qfit,gpfit=gpfit,fte=fte,fne=fne,fn0=fn0,fits=fits,rho=rho,center=center,$
		filepath=filepath,temp=etemp,dens=edens
	
	;set 'fudge factors' for sensitivity analysis
	IF NOT keyword_set(fte) THEN fte=1.0
	IF NOT keyword_set(fne) THEN fne=1.0
	IF NOT keyword_set(fn0) THEN fn0=1.0

	pcase=-1
	IF keyword_set(fits) THEN pcase=0
	IF keyword_set(qfit) THEN pcase=1
	IF keyword_set(gpfit) THEN pcase=2
	IF keyword_set(edens) AND keyword_set(etemp) THEN pcase=3
	ncase=0
	IF keyword_set(trshot) THEN ncase=1
	IF pcase EQ -1 THEN BEGIN
		print, 'please select profile loading (/qfit, /fits)'
		RETURN
	ENDIF
	IF NOT keyword_set(rho) THEN rho=make(0.0,1.0,100)	;default the rho grid to be 100 points evenly spaceed for 0 < rhopol < 1 
	IF NOT keyword_set(nfrac) THEN nfrac=1.0e-8		;set neutral density below level of participating
	CASE pcase OF
		0 : BEGIN	;default is to use fiTS savefile from the user's home directory
			IF size(fits,/type) EQ 7 THEN path=fits ELSE path='/home/'+logname()+'/fits/fits_'+num2str(shot,1)+'.save'
			gentran_tene_widgetfits,path,time,dt,temp,dens,ete=ete,ene=ene,/psin,rho=rho^2,center=center	;interpolate data onto psin=(0->1)^2		
		END

		1 : BEGIN	;use the "quick-fits"
			gentran_tene_quickfits,shot,time,dt,temp,dens,ete=ete,ene=ene,/psin,rho=rho^2,center=center	;interpolate data onto psin=(0->1)^2		
		END

		2 : BEGIN	;use the "gaussian process fits"
			IF size(gpfit,/type) EQ 7 THEN gp_path=gpfit ELSE gp_path=0
			gentran_tene_gpfit,shot,time,dt,temp,dens,ete=ete,ene=ene,/psin,rho=rho^2,center=center,path=gp_path	;interpolate data onto psin=(0->1)^2		
		END
		3 : BEGIN	;density and temperature profiles set, interpolate on the time/space grid
			nrho=n(rho)+1
			ntime=n(time)+1
			temp=fltarr(nrho,ntime)
			dens=fltarr(nrho,ntime)	
			FOR i=0,ntime-1 DO BEGIN
				tmp=where(etemp.time GE time[i] AND etemp.time LT time[i]+dt[i])
				IF n(tmp) EQ 0 THEN temp[*,i]=interpol(etemp.temp[*,tmp[0]],etemp.rho,rho^2) ELSE $
					temp[*,i]=interpol(sum_array(etemp.temp[*,tmp],/i)/(n(tmp)+1.0),etemp.rho,rho^2)	;interpolate data onto psin=(0->1)^2
				tmp=where(edens.time GE time[i] AND edens.time LT time[i]+dt[i])
				IF n(tmp) EQ 0 THEN dens[*,i]=interpol(edens.dens[*,tmp[0]],edens.rho,rho^2) ELSE $
					dens[*,i]=interpol(sum_array(edens.dens[*,tmp],/i)/(n(tmp)+1.0),edens.rho,rho^2)	;interpolate data onto psin=(0->1)^2
			ENDFOR
			temp={temp:temp,rho:rho^2,err:temp*0.0,time:time,dt:dt}
			dens={dens:dens,rho:rho^2,err:dens*0.0,time:time,dt:dt}
		END
	ENDCASE


	CASE ncase OF
		0 : BEGIN	;default uses nfrac multiplied into elec. density
			neut=dens
			neut.dens=neut.dens*nfrac
                END

		1 : BEGIN	;loads the neutral density from TRANSP given in trshot
			gentran_neut_transp,shot,trshot,time,dt,neut,en0=en0,/psin,rho=rho^2,center=center	;interpolate data onto psin=(0->1)^2
		END
	ENDCASE
	ntimes=n(time)+1
	times=time+dt/2.0
	npoints=n(rho)+1
	normdens=dens.dens*fne/1.0e6				;convert to 1/cm^3
	normdens=[normdens[0,*],normdens]
	FOR i=0,ntimes-1 DO normdens[1:*,i]/=normdens[0,i]	;convert to input as normalized profile
	normtemp=temp.temp*fte
	normtemp=[normtemp[0,*],normtemp]
	FOR i=0,ntimes-1 DO normtemp[1:*,i]/=normtemp[0,i]	;convert to input as normalized profile
	normneut=neut.dens*fn0/1.0e6				;convert to 1/cm^3
	normneut=[normneut[0,*],normneut]
	FOR i=0,ntimes-1 DO normneut[1:*,i]/=normneut[0,i]	;convert to input as normalized profile


	IF keyword_set(ln_sol) THEN nefunc='interp' ELSE nefunc='interpa'
	IF keyword_set(lt_sol) THEN tefunc='interp' ELSE tefunc='interpa'

	filepath=strcompress('/home/'+logname()+'/strahl/cmod/pp'+string(shot)+'.0',/remove_all)
	openw,1,filepath
	printf,1,''
	printf,1,'cv    time-vector'		;electron density
	printf,1,int(ntimes)
	printf,1,times
	printf,1,''
	printf,1,'cv    Ne-function'
	printf,1,"      '"+nefunc+"'"
	printf,1,''
	printf,1,''
	printf,1,'cv    x-coordinate'
	printf,1,"      'poloidal rho'"
	printf,1,''
	printf,1,''
	printf,1,'cv    # of interpolation points'
	printf,1,'      ',npoints
	printf,1,''
	printf,1,''
	printf,1,'cv    x-grid for ne-interpolation'
	printf,1,sqrt(dens.rho[*,0])
	printf,1,''
	printf,1,''
	printf,1,'cv    DATA'
	printf,1,normdens
	printf,1,''
	IF keyword_set(ln_sol) THEN BEGIN
		IF n(ln_sol) NE ntimes-1 THEN decay=fltarr(ntimes)+ln_sol[0] ELSE decay=ln_sol
		printf,1,''
		printf,1,'cv    Ne decay length'
		printf,1,decay
		printf,1,''
	ENDIF
	printf,1,''
	printf,1,'cv    time-vector'		;electron temperature
	printf,1,int(ntimes)
	printf,1,times
	printf,1,''
	printf,1,''
	printf,1,'cv    Te-function'
	printf,1,"      '"+tefunc+"'"
	printf,1,''
	printf,1,''
	printf,1,'cv    x-coordinate'
	printf,1,"      'poloidal rho'"
	printf,1,''
	printf,1,''
	printf,1,'cv    # of interpolation points'
	printf,1,'      ',npoints
	printf,1,''
	printf,1,''
	printf,1,'cv    x-grid for Te-interpolation'
	printf,1,sqrt(temp.rho[*,0])
	printf,1,''
	printf,1,''
	printf,1,'cv    DATA'
	printf,1,normtemp
	printf,1,''
	IF keyword_set(lt_sol) THEN BEGIN
		IF n(lt_sol) NE ntimes-1 THEN decay=fltarr(ntimes)+lt_sol[0] ELSE decay=lt_sol
		printf,1,''
		printf,1,'cv    Te decay length'
		printf,1,decay
		printf,1,''
	ENDIF		
	printf,1,''
	printf,1,'cv    time-vector'		;set Ti=Te
	printf,1,0
;printf,1,times
;printf,1,''
;printf,1,''
;printf,1,'cv    Ti-function'
;printf,1,"      'interp'"
;printf,1,''
;printf,1,''
;printf,1,'cv    x-coordinate'
;printf,1,"      'poloidal rho'"
;printf,1,''
;printf,1,''
;printf,1,'cv    # of interpolation points'
;printf,1,'      ',npoints
;printf,1,''
;printf,1,''
;printf,1,'cv    x-grid for Ti-interpolation'
;printf,1,xgrid
;printf,1,''
;printf,1,''
;printf,1,'cv    DATA'
;printf,1,
;printf,1,''
;printf,1,''
;printf,1,'cv    Ti decay length'
;printf,1,decay
	printf,1,''
	printf,1,'cv    time-vector'		;set neutral density to 10^-5 of ne for now
	printf,1,int(ntimes)
	printf,1,times
	printf,1,''
	printf,1,'cv    N0-function'
	printf,1,"      'interpa'"
	printf,1,''
	printf,1,''
	printf,1,'cv    x-coordinate'
	printf,1,"      'poloidal rho'"
	printf,1,''
	printf,1,''
	printf,1,'cv    # of interpolation points'
	printf,1,'      ',npoints
	printf,1,''
	printf,1,''
	printf,1,'cv    x-grid for n0-interpolation'
	printf,1,sqrt(neut.rho[*,0])
	printf,1,''
	printf,1,''
	printf,1,'cv    DATA'
	printf,1,normneut

	close,1
	free_lun,1

END
