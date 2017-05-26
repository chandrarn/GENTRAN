;+
;NAME:
;	RUN_CMOD_STRAHL
;
;PURPOSE:
;	This is a higher level procedure which generates input files for the STRAHL impurity transport code, 
;	executes the code and then outputs selected data from the netCDF file.
;
;CALLING SEQUENCE:
;	RUN_CMOD_STRAHL,shot,z,ta,tb,csden,rho,time
;
;INPUTS:
;	shot	LONG	shot number
;	z	INT	atomic number of impurity to simulate
;	ta	FLOAT	start time of simulation (see optional inputs)
;	tb	FLOAT	stop time of simulation (see optional inputs)
;	
;OPTIONAL INPUTS:
;	ta	FLTARR	[ntau] start times of time bins
;	tb	FLTARR  [ntau] delta times of time bins
;	MAKE_CMOD_GRID Options: /gps (/ps) and /gplots (/plots)
;	MAKE_CMOD_PP Options: ln_sol, lt_sol, fits, gpfit, temp, dens (/qfits is default)
;	MAKE_CMOD_PARAM Options: k,nrho,source,tau,fz,exp,dff,cff,doff,diff,conv,saw
;	GENTRAN_PLOT Options: xr,qhigh,qlow,time (tplot)
;
;	plot	INT	use /plot to plot only the radial profiles using GENTRAN_PLOT and use plot=q to plot both GENTRAN_PLOT
;			and a GENPLT of that time-evolving csden[*,q,*] profile
;	path	STRING	location of STRAHL output file, if set as input skips running of STRAHL and reads the output file
;
;KEYWORD PARAMETERS:
;	nogrid		/nogrid skips MAKE_CMOD_GRID and does not copy the output file from ~/strahl/cmod to ~/strahl/nete
;	nopp		/nopp skips MAKE_CMOD_PP and does not copy the output file from ~/strahl/cmod to ~/strahl/nete
;	noparam		/noparam skips MAKE_CMOD_PARAM and does not copy the output file from ~/strahl/cmod to ~/strahl/param_files
;	nostrahl	/nostrahl skips executing STRAHL as well as the MAKE_CMOD*, useful for viewing already run data
;	debug		/debug stops the code before and after the execution of STRAHL
;	list		/list will list the possible variables that can be retrieved from the STRAHL CDF file
;	quiet		/quiet will suppress terminal notifications 
;
;OUTPUTS:
;	csden	FLTARR	[nrho,z+1,ntime] of the charge state density profiles output from STRAHL
;	rho	FLTARR 	[nrho] of the radial grid in r/a
;	time	FLTARR	[ntime] of the time scale [sec]
;
;OPTIONAL OUTPUTS:
;	path	STRING	of the location of the output file, if set as input leads to different behavior (see INPUTS)
;	term	STRARR	of the terminal output from running STRAHL
;	data	STRUC	of relevant data I/O from the STRAHL run
;		*.shot	LONG	shot number
;		*.time	DBLARR	[ntime] time points [sec]
;		*.temp	DBLARR	[nrho,ntime] of the electron temp. [eV]
;		*.terr	DBLARR	[nrho,ntime] of the unc. in temp. DEFAULT=0.0
;		*.dens	DBLARR	[nrho,ntime] of the electron dens. [m^-3]
;		*.derr	DBLARR	[nrho,ntime] of the unc. in dens. DEFAULT=0.0
;		*.neut	DBLARR	[nrho,ntime] of the neutral dens. [m^-3]
;		*.nerr	DBLARR	[nrho,ntime] of the unc. in neutral dens. DEFAULT=0.0
;		*.psin	DBLARR	[nrho] of the radial grid in norm. pol. flux
;		*.rho	DBLARR	[nrho] of the radial grid in r/a
;		*.rmaj	DBLARR	[nrho] of the radial in LFS midplane major radius [m]
;		*.csden	DBLARR	[nrho,z+1,ntime] of the charge state densities [m^-3]
;		*.cserr	DBLARR	[nrho,z+1,ntime] of the unc. in the charge state densities [m^-3] DEFAULT=0.0
;		*.diff	DBLARR	[nrho,ntime] of the specified diffusion profile [m^2/s]
;		*.conv	DBLARR	[nrho,ntime] of the specified convection profile [m/s]
;		*.dneo	DBLARR	[nrho,ntime] the neoclassical diffusion computed using NEOART [m^2/s]
;		*.vneo	DBLARR	[nrho,ntime] the neoclassical convection computed using NEOART [m/s]
;
;MODIFICATION HISTORY:
;	Written by:	M.L. Reinke - January 2013
;	M.L. Reinke	12/2/2013 - added the TRSHOT and FXX keywords into MAKE_CMOD_PP
;	M.L. Reinke	12/13/2013 - added the PATH keyword to enable read-only operation and output offile location
;	M.L. Reinke 	9/14/2014 - added the GPFIT keyword
;	M.L. Reinke	11/6/2014 - modified how time data is specified to allow for time-evolving kinetic profiles
;
;-

PRO run_cmod_strahl,shot,z,ta,tb,csden,rho,time,nogrid=nogrid,nopp=nopp,noparam=noparam,nostrahl=nostrahl,quiet=quiet,list=list,data=data,$	;main program I/O
			term=term,err=err,debug=debug,plot=plot,xr=xr,qhigh=qhigh,qlow=qlow,tplot=tplot,path=path,$				;main program I/O
			gplot=gplot,gps=gps,$													;grid I/O
			fits=fits,gpfit=gpfit,ln_sol=ln_sol,lt_sol=lt_sol,nfrac=nfrac,trshot=trshot,center=center,fte=fte,fne=fne,fn0=fn0,$	;pp I/O
			temp=temp,dens=dens,$													;pp I/O
			dt=dt,ncyc=ncyc,niter=niter,fiter=fiter,k=k,nrho=nrho,dr=dr,$								;param I/O setup
			tneut=tneut,source=source,tau=tau,fz=fz,exp=exp,dff=dff,cff=cff,doff=doff,diff=diff,conv=conv,saw=saw			;param I/O transport
	
	IF keyword_set(path) THEN nostrahl=1
	IF keyword_set(nostrahl) THEN BEGIN
		nogrid=1
		nopp=1
		noparam=1
	ENDIF
	IF NOT keyword_set(ln_sol) THEN ln_sol=2.0
	IF ln_sol LT 0 THEN ln_sol=0.0
	IF NOT keyword_set(lt_sol) THEN lt_sol=1.0
	IF lt_sol LT 0 THEN lt_sol=0.0
	IF NOT keyword_set(nfrac) THEN nfrac=1.0e-8
	
	shotstr=num2str(shot,1)+'.0'
	elemstr=read_atomic_name(z)
	IF strlen(elemstr) EQ 1 THEN elemstr+='_'
	flxstr=elemstr+'flx'+num2str(shot,1)+'.dat'
	s1=size(ta)
	s2=size(tb)
	IF s1[0] NE 0 AND s2[0] NE 0 THEN BEGIN	
		t0=ta
		deltat=tb
		t1=t0[0]
		t2=last(t0)+last(deltat)
     	ENDIF ELSE BEGIN
		t0=ta
		deltat=tb-ta
		t1=ta
		t2=tb
	ENDELSE

	IF NOT keyword_set(nogrid) THEN BEGIN
		make_cmod_grid,shot,0.5*(t1+t2),filepath=gridpath,plot=gplot,ps=gps
		dest='/home/'+logname()+'/strahl/nete/grid_'+shotstr
		spawn,'cp '+gridpath+' '+dest
		IF NOT keyword_set(quiet) THEN print, 'grid file generated, copied to '+dest
        ENDIF ELSE IF NOT keyword_set(quiet) THEN print, 'skipping MAKE_CMOD_GRID'

	IF NOT keyword_set(nopp) THEN BEGIN
		IF NOT keyword_set(fits) AND NOT keyword_set(gpfit) THEN qfit=1 ELSE qfit=0
		make_cmod_pp,shot,t0,deltat,ln_sol=ln_sol,lt_sol=lt_sol,nfrac=nfrac,trshot=trshot,qfit=qfit,fits=fits,gpfit=gpfit,center=center,filepath=profpath,$
			fte=fte,fne=fne,fn0=fn0,temp=temp,dens=dens
		dest='/home/'+logname()+'/strahl/nete/pp'+shotstr
		spawn,'cp '+profpath+' '+dest
		IF NOT keyword_set(quiet) THEN print, 'profile file generated, copied to '+dest
        ENDIF ELSE IF NOT keyword_set(quiet) THEN print, 'skipping MAKE_CMOD_PP'
	
	IF NOT keyword_set(noparam) THEN BEGIN
		make_cmod_param,shot,z,t1,t2,dt=dt,ncyc=ncyc,ion=ion,tneut=tneut,start=start,niter=niter,fiter=fiter,k=k,nrho=nrho,dr=dr,source=source,tau=tau,fz=fz,exp=exp,$
				dff=dff,cff=cff,doff=doff,diff=diff,conv=conv,saw=saw,filepath=parapath,spath=spath
		IF size(spath,/type) EQ 7 THEN BEGIN
			dest='/home/'+logname()+'/strahl/nete/'+flxstr
			spawn,'cp '+spath+' '+dest
			IF NOT keyword_set(quiet) THEN print, 'source file generated,copied to '+dest
                ENDIF ELSE IF NOT keyword_set(quiet) THEN print,'constant source specified, no source file needed'
		dest='/home/'+logname()+'/strahl/param_files/run_'+shotstr
		spawn,'cp '+parapath+' '+dest
		IF NOT keyword_set(quiet) THEN print, 'para file generated,copied to '+dest
        ENDIF ELSE IF NOT keyword_set(quiet) THEN print, 'skipping MAKE_CMOD_PARAM'

	IF NOT keyword_set(nostrahl) THEN BEGIN
		;Modify strahl.control  correct file
		openw,1,'/home/'+logname()+'/strahl/strahl.control'
		printf,1,'run_'+strtrim(shot,2)+'.0'
		printf,1,'  '+num2str(t2,dp=2)
		printf,1,'E'
		close,1
		free_lun,1

		IF keyword_set(debug) THEN stop
		IF NOT keyword_set(quiet) THEN print, 'running /home/'+logname()+'/strahl/strahl'	

		cd,c=currdir
		cd,'/home/'+logname()+'/strahl/'
		spawn,'./strahl a n',term,err	
		cd,currdir
		IF NOT keyword_set(quiet) THEN print, 'STRAHL run complete, '+num2str(n(term)+1,1)+' lines output to terminal'	
	
		IF keyword_set(debug) THEN stop

		spawn, 'cp /home/'+logname()+'/strahl/result/strahl_result.dat /home/'+logname()+'/strahl/result/'+elemstr+'_'+shotstr+'.dat'
 		IF NOT keyword_set(quiet) THEN print, 'results copied to /home/'+logname()+'/strahl/result/'+elemstr+'_'+shotstr+'.dat'
        ENDIF ELSE if NOT keyword_set(quiet) THEN print, 'skipping running STRAHL'

	IF keyword_set(path) THEN file=path ELSE file=elemstr+'_'+shotstr+'.dat'
	read_strahl,file,'rho_poloidal_grid',rhop
	read_strahl,file,'time',time
	read_strahl,file,'impurity_density',csden
	read_strahl,file,'anomal_diffusion',d
	read_strahl,file,'anomal_drift',v
	read_strahl,file,'classical_diff_coeff',d1
	read_strahl,file,'pfirsch_schlueter_diff_coeff',d2
	read_strahl,file,'banana_plateau_diff_coeff',d3
	dneo=d1+d2+d3
	read_strahl,file,'classical_drift',v1
	read_strahl,file,'pfirsch_schlueter_drift',v2
	read_strahl,file,'banana_plateau_drift',v3
	;read_strahl,file,'classical_drift_t_part',v4
 	;read_strahl,file,'pfirsch_schlueter_drift_t_part',v5
 	;read_strahl,file,'banana_plateau_drift_t_part',v6
 	;read_strahl,file,'ware_pinch',v7
	vneo=v1+v2+v3;+v4+v5+v6+v7
        read_strahl,file,'electron_density',edens
        read_strahl,file,'electron_temperature',etemp
	read_strahl,file,'large_radius_lfs',rmaj			;this is actually RMID in C-Mod terminology
	read_strahl,file,'large_radius',Ro,/gatr
        read_strahl,file,'neutral_hydrogen_density',neut		
	IF NOT keyword_set(neut) THEN neut=0.0*edens 

	a=interpol(rmaj,rhop,1.0)-Ro
	rho=(rmaj-Ro)/a

	data={shot:shot,time:time,temp:etemp,terr:etemp*0.0,dens:edens*1.0e6,derr:edens*0.0,neut:neut*1.0e6,nerr:neut*0.0,$
		psin:rhop^2,rho:rho,rmaj:rmaj/100.0,csden:csden,cserr:csden*0.0,diff:d/1.0e4,conv:v/1.0e2,dneo:dneo/1.0e4,vneo:vneo/1.0e2}
	path='/home/'+logname()+'/strahl/result/'+elemstr+'_'+shotstr+'.dat'
	path=path[0]
	;plot the output data
	IF keyword_set(plot) THEN BEGIN
		gentran_plot,data,xr=xr,qhigh=qhigh,qlow=qlow,time=tplot
		IF plot NE 1 THEN BEGIN
			labels={ilab:'r/a',jlab:'Time [sec]',klab:'Charge State Density',ctit:'Time Evolution of '+elemstr+'!u'+num2str(plot,1)+'+!n',itit:'',jtit:''}
			genplt,reform(data.csden[*,plot,*]),data.rho,data.time,labels=labels,io=io,jo=jo,cct=39,ncntrs=50
		ENDIF
	ENDIF

	;interrogate the result CDF file for available outputs
	IF keyword_set(list) THEN BEGIN
		path='/home/'+logname()+'/strahl/result/'+elemstr+'_'+shotstr+'.dat'
		ncid=ncdf_open(path[0])
		size=ncdf_inquire(ncid)
		FOR i=0,size.nvars-1 DO BEGIN
			var=ncdf_varinq(ncid,long(i))
			print, '  '+var.name
                ENDFOR
		ncdf_close,ncid
	ENDIF
END
