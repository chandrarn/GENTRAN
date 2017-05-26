FUNCTION neuts_to_wmhd,neuts

	l10wmhd=(alog10(neuts)+2.4913)/3.08508
	wmhd=10^l10wmhd
	RETURN,wmhd

END

PRO comp_neut_rate,shots,tr=tr,color=color,titstr=titstr,oplot=oplot,lin=lin
	IF NOT keyword_set(tr) THEN tr=[0.6,1.3]
	IF NOT keyword_set(titstr) THEN titstr=' '
	IF NOT keyword_set(color) THEN color=0
	IF NOT keyword_set(oplot) THEN BEGIN
		neuts = (findgen(400)/2.+1.)*1.e12
		wmhd = neuts_to_wmhd(neuts)
		IF NOT keyword_set(lin) THEN BEGIN
			ylog=1
			yr=[1e12, 4e14]
		ENDIF ELSE BEGIN
			ylog=0
			yr=[0,2.5e14]
		ENDELSE
		plot, wmhd/1.e3, neuts,  yr=yr, xr=[0,300], $
  			ytit='Neut rate', xtit ='Wmhd (kJ)', tit = "Curve calculated from Earl's formula "+titstr, ylog=ylog,/ysty
	ENDIF
	makesym,9
	FOR i=0,n(shots) DO BEGIN
		load_neutrate_data,shots[i],neut,ntime,status=status
		load_wth_data,shots[i],wth,wtime
		tmp=where(wtime GE tr[0] AND wtime LE tr[1])
		IF tmp[0] NE -1 AND status THEN BEGIN
			wth=wth[tmp]/1.0e3
			time=wtime[tmp]
			neut=interpol(neut,ntime,time)
			oplot,wth,neut,psym=8,color=color,symsize=0.5
		ENDIF
	ENDFOR
	makesym,10

END


PRO load_wth_data,shot,y,x
	mdsopen, 'mhd', (shot)
	y=mdsvalue('\ANALYSIS::EFIT_AEQDSK:WPLASM')
	x=mdsvalue('dim_of(\ANALYSIS::EFIT_AEQDSK:WPLASM)')
	mdsclose, 'mhd', (shot)
END

FUNCTION line_gettimes,shot
	mdsopen, 'mhd', (shot)
	t=mdsvalue('dim_of(\MHD::TOP.ANALYSIS.EFIT.RESULTS.G_EQDSK:RMAXIS)',$
		/quiet,status=status)
	IF status THEN output=t ELSE output=-1
	mdsclose, 'mhd', (shot)
	RETURN, t
END

FUNCTION line_getrmid,shot
	mdsopen, 'analysis', shot
	rmid=mdsvalue('\efit_rmid')
	mdsclose, 'analysis', shot
	output=rmid
	RETURN, output
END

PRO load_h2d_data,shot,h2d,time
	mdsopen,'spectroscopy',shot
	h2d=mdsvalue('\SPECTROSCOPY::TOP.HALPH_DALPH.ANALYSIS:H_TO_D_RATIO',status=status)
	time=mdsvalue('dim_of(\SPECTROSCOPY::TOP.HALPH_DALPH.ANALYSIS:H_TO_D_RATIO)',status=status)
	mdsclose,'spectroscopy',shot
	IF NOT status THEN BEGIN
		time=make(0,2.0,100)
		h2d=time*0.0
	ENDIF
END

PRO load_neutrate_data,shot,neuts,time,status=status
	mdsopen,'particles',shot
	neuts=mdsvalue('\PARTICLES::TOP.NEUTRONS.GLOBAL.RESULTS:NEUT_RATE',status=status,/quiet)
	time=mdsvalue('dim_of(\PARTICLES::TOP.NEUTRONS.GLOBAL.RESULTS:NEUT_RATE)',/quiet)
	mdsclose,'particles',shot
END

PRO load_he3neutrate_data,shot,neuts,time,det=det
	IF NOT keyword_set(det) THEN det=4
	CASE abs(det) OF
		1 : path='\PARTICLES::TOP.NEUTRONS.HE_3_BANK.results:he3_nrate_1'
		2 : path='\PARTICLES::TOP.NEUTRONS.HE_3_BANK.results:he3_nrate_2'
		3 : path='\PARTICLES::TOP.NEUTRONS.HE_3_BANK.results:he3_nrate_3'
		4 : path='\PARTICLES::TOP.NEUTRONS.HE_3_BANK.results:he3_nrate_4'
	ENDCASE
	mdsopen,'particles',shot
	neuts=mdsvalue(path)
	time=mdsvalue('dim_of('+path+')')
	mdsclose,'particles',shot
END

PRO load_zave_data,shot,zave,time
	mdsopen,'spectroscopy',shot
	zave=mdsvalue('\SPECTROSCOPY::TOP.Z_METER.ANALYSIS:Z_AVE')
	time=mdsvalue('dim_of(\SPECTROSCOPY::TOP.Z_METER.ANALYSIS:Z_AVE)')
	mdsclose,'spectroscopy',shot
END

PRO load_ti_tree_data,shot,ti,r,t,interp=interp
	mdsopen,'spectroscopy',shot
	ti=mdsvalue('\SPECTROSCOPY::TOP.HIREX_SR:RESULTS:CONSOLIDATED:TI',/quiet,status=status)
	IF NOT status THEN BEGIN
		ti=[-1]
		RETURN
	ENDIF
	rmaj=mdsvalue('\SPECTROSCOPY::TOP.HIREX_SR:RESULTS:CONSOLIDATED:R_MAJOR')
	t=mdsvalue('\SPECTROSCOPY::TOP.HIREX_SR:RESULTS:CONSOLIDATED:TIMES')
	mdsclose,'spectroscopy',shot

	IF keyword_set(interp) THEN BEGIN
		n_r=30
		ti_new=fltarr(n_r,n(t)+1)
		r=make(0.68,0.9,n_r)
		FOR i=0,n(t) DO ti_new[*,i]=interpol(ti[*,i],rmaj[*,i],r)
		ti=ti_new
	ENDIF ELSE r=rmaj
END


PRO load_bsti_tree_data,shot,ti,r,t,redge=redge,tierr=tierr,tht=tht,rho=rho,psin=psin,nr=nr
	IF NOT keyword_set(tht) THEN tht=0
	isrgrid=0
	IF keyword_set(rho) THEN isrgrid=1
	IF keyword_set(psin) THEN isrgrid=2
	IF NOT keyword_set(redge) THEN redge=0.86
	IF keyword_set(tht) THEN rstr='RESULTS'+num2str(tht,1) ELSE rstr='RESULTS'
	path='\SPECTROSCOPY::TOP.HIREXSR.'+rstr+'.BSTI:'
	mdsopen,'spectroscopy',shot
	fit=mdsvalue('_dat='+path+'FIT',/quiet,status=status)
	rho=mdsvalue('dim_of(_dat,0)',/quiet)
	time=mdsvalue('dim_of(_dat,1)',/quiet)
	err=mdsvalue('dim_of(_dat,2)',/quiet)	
	IF NOT status THEN print, 'ERROR - NO BSTI FOR SHOT/THT = '+num2str(shot,1)+'/'+num2str(tht,1)
	inst=mdsvalue('_dat='+path+'INST',/quiet,status=status)
	IF NOT status THEN BEGIN
		print, 'NO INSTRUMENTAL TEMP FOR SHOT/THT = '+num2str(shot,1)+'/'+num2str(tht,1)+', ASSUMING 200 eV'
		inst=fit*0.0+0.2
	ENDIF
	
	mdsopen, 'analysis', shot
	rmid=mdsvalue('\efit_rmid')
	etime=mdsvalue('dim_of(\efit_rmid,0)')
	epsi=mdsvalue('dim_of(\efit_rmid,1)')
	mdsclose, 'analysis', shot

	IF NOT keyword_set(nr) THEN nr=30
	CASE isrgrid OF 
		0 : r=make(0.68,0.89,nr)	;evently spaced in rmaj
		1 : r=make(0,1.0,nr) 		;evenly spaced in r/a
		2 : r=make(0,1.0,nr)^2		;evenly spaced in sqrt(psin)
	ENDCASE
	nt=n(time)+1
	ti_new=fltarr(nr,nt)
	tierr_new=fltarr(nr,nt)
	FOR i=0,nt-1 DO BEGIN
		iefit=ipt(etime,time[i])
		irmid=reform(rmid[iefit,*])
		irho=(irmid-irmid[0])/(last(irmid)-irmid[0])
		tmp=where(rho[*,i] NE -1)
		CASE isrgrid OF
			0 : BEGIN
				rvec=interpol(rmid[iefit,*],irho,rho[tmp,i])
				edgerad=redge
			END
			1 : BEGIN
				rvec=rho[tmp,i]
				edgerad=(redge-irmid[0])/(last(irmid)-irmid[0])
			END
			2 : BEGIN
				rvec=interpol(epsi,irho,rho[tmp,i])
				edgerad=interpol(epsi,rmid[iefit,*],redge)
			END

		ENDCASE
		ti_new[*,i]=interpol(fit[tmp,i]-inst[tmp,i],rvec,r)
		tierr_new[*,i]=interpol(err[tmp,i],rvec,r)
		j=ipt(r,redge)
		ti_new[j:*,i]=ti_new[j,i]*(1.0-(r[j:*]-redge)/(max(r)-redge))
		tierr_new[j:*,i]=0.1
	ENDFOR
	ti=ti_new
	tierr=tierr_new
	t=time

END

PRO load_tinew_tree_data,shot,ti,r,t,z=z,w=w,redge=redge,tierr=tierr,tht=tht,comb=comb,rho=rho,psin=psin,nr=nr
	isrgrid=0
	IF keyword_set(rho) THEN isrgrid=1
	IF keyword_set(psin) THEN isrgrid=2

	IF NOT keyword_set(redge) THEN redge=0.86
	IF NOT keyword_set(z) OR keyword_set(w) THEN z=1
	IF keyword_set(z) THEN detstr='Z' ELSE detstr='W'
	IF NOT keyword_set(tht) THEN tht=0
	IF keyword_set(tht) THEN astr='ANALYSIS'+num2str(tht,1) ELSE astr='ANALYSIS'
	pathstr='\SPECTROSCOPY::TOP.HIREXSR.'+astr+'.HELIKE.PROFILES.'+detstr+':PRO'
	mdsopen,'spectroscopy',shot
	profile=mdsvalue(pathstr)
	psinorm=mdsvalue('dim_of('+pathstr+',0)')
	tau=mdsvalue('dim_of('+pathstr+',1)')
	tinst=mdsvalue('dim_of('+pathstr+',2)')
	pathstr='\SPECTROSCOPY::TOP.HIREXSR.'+astr+'.HELIKE.PROFILES.'+detstr+':PROERR'
	proerr=mdsvalue(pathstr)
	IF keyword_set(comb) THEN BEGIN
		detstr='LYA1'
		pathstr='\SPECTROSCOPY::TOP.HIREXSR.'+astr+'.HLIKE.PROFILES.'+detstr+':PRO'
		hprofile=mdsvalue(pathstr)
		hpsinorm=mdsvalue('dim_of('+pathstr+',0)')
		htau=mdsvalue('dim_of('+pathstr+',1)')
		htinst=mdsvalue('dim_of('+pathstr+',2)')
		pathstr='\SPECTROSCOPY::TOP.HIREXSR.'+astr+'.HLIKE.PROFILES.'+detstr+':PROERR'
		hproerr=mdsvalue(pathstr)
	ENDIF ELSE htau=0.0
	mdsclose,'spectroscopy',shot
	IF total(tau-htau) NE 0 THEN comb=0	;turn off combine if not on the same time base	

	mdsopen, 'analysis', shot
	rmid=mdsvalue('\efit_rmid')
	time=mdsvalue('dim_of(\efit_rmid,0)')
	epsi=mdsvalue('dim_of(\efit_rmid,1)')
	mdsclose, 'analysis', shot

	tmp=where(tau NE -1 AND tau GT time[0] AND tau LT last(time))
	ti=profile[*,tmp,3]
	tierr=proerr[*,tmp,3]
	t=tau[tmp]
	psinorm=psinorm[*,tmp]
	IF keyword_set(comb) THEN BEGIN
		hti=hprofile[*,tmp,3]
		htierr=hproerr[*,tmp,3]
		hpsinorm=hpsinorm[*,tmp]
	ENDIF	

	IF NOT keyword_set(nr) THEN nr=30
	CASE isrgrid OF 
		0 : r=make(0.68,0.89,nr)	;evently spaced in rmaj
		1 : r=make(0,1.0,nr) 		;evenly spaced in r/a
		2 : r=make(0,1.0,nr)^2		;evenly spaced in sqrt(psin)
	ENDCASE
	ti_new=fltarr(nr,n(t)+1)
	tierr_new=fltarr(nr,n(t)+1)
	FOR i=0,n(t) DO BEGIN
		iefit=ipt(time,t[i])
		ipsin=psinorm[where(psinorm[*,i] NE -1),i]
		IF ipsin[0] EQ ipsin[(n(ipsin)+1)/2] THEN ipsin=ipsin[0:(n(ipsin)+1)/2-1]
		npts=n(ipsin)+1
		iti=ti[0:npts-1,i]-tinst
		IF n(ipsin) EQ n(ti[*,i]) THEN itierr=tierr[0:npts-1,i] ELSE itierr=sqrt(tierr[0:npts-1,i]^2+ti[npts:*,i]^2)
		IF keyword_set(comb) THEN BEGIN
			ihpsin=hpsinorm[where(hpsinorm[*,i] NE -1),i]		;assume no /sine term is saved
			ihti=hti[*,i]-htinst
			ihtierr=htierr[*,i]
			etmp=where(ihtierr LT 0.05*ihti)
			IF etmp[0] NE -1 THEN ihtierr[tmp]=0.05*ihti[tmp]
			hetmp=where(ipsin GT comb)
			htmp=where(ihpsin LE comb)
			itemp=[ihti[htmp],iti[hetmp]]
			iterr=[ihtierr[htmp],itierr[hetmp]]
			ipsin=[ihpsin[htmp],ipsin[hetmp]]
		ENDIF ELSE BEGIN
			itemp=iti
			iterr=itierr
			ipsin=ipsin
                ENDELSE
		CASE isrgrid OF
			0 : BEGIN
				rvec=interpol(rmid[iefit,*],epsi,ipsin)
				edgerad=redge
			END
			1 : BEGIN
				Ro=rmid[iefit,0]
				a=last(rmid[iefit,*])-Ro
				rvec=interpol((rmid[iefit,*]-Ro)/a,epsi,ipsin)
				edgerad=(redge-ro)/a
			END
			2 : BEGIN
				rvec=ipsin
				edgerad=interpol(epsi,rmid[iefit,*],redge)
			END

		ENDCASE
		ti_new[*,i]=interpol(itemp,rvec,r)
		tierr_new[*,i]=interpol(iterr,rvec,r)
		j=ipt(r,redge)
		ti_new[j:*,i]=ti_new[j,i]*(1.0-(r[j:*]-redge)/(max(r)-redge))
		tierr_new[j:*,i]=0.1
	ENDFOR
	ti=ti_new
	tierr=tierr_new
	
END	

;+
;NAME:
;	NEUT_RATE
;
;PURPOSE:
;	This function determines the volume-integrated neutron rate given inputs of ion temperature,
;	electron density, hydrogen/deuterium ratio and zeff.  The Bosch and Hale NF v32 pg611 data is used
;	for the <sigma*v>
;
;CALLING SEQUENCE:
;	result=NEUT_RATE(shot,t_i,n_e,h2d,zeff)
;
;INPUTS:
;	shot:		LONG	shot number
;	t_i:		STRUC	ion temperature data
;			*.dat	FLTARR	[r,t] of ion temperature [keV]  **Array size (r,t) can be different for each input**
;			*.rad	FLTARR 	[r] of major radius values [m]
;			*.time	FLTARR	[t] of time points [sec]
;	n_e		STRUC	electron density data
;			*.dat	FLTARR	[r,t] of electron density[10^20 m^-3]	
;			*.rad	FLTARR 	[r] of major radius values [m]
;			*.time	FLTARR	[t] of time points [sec]
;	h2d		STRUC	hydrogen/deuterium fraction
;			*.dat	FLTARR	[r,t] of ratio	
;				**optional to have *.dat be FLTARR [t] and assume constant over plasma radius**
;			*.rad	FLTARR 	[r] of major radius values [m]
;				**if assumeing constant set h2d.rad = [0]
;			*.time	FLTARR	[t] of time points [sec]
;	zeff		STRUC	zeff data
;			*.dat	FLTARR	[r,t] of zeff	
;				**optional to have *.dat be FLTARR [t] and assume constant over plasma radius**
;			*.rad	FLTARR 	[r] of major radius values [m]
;				**if assumeing constant set zeff.rad = [0]
;			*.time	FLTARR	[t] of time points [sec]
;			*.z	INT	of the "average" impurity charge making up the zeff
;
;OPTIONAL INPUTS:
;	time:		FLTARR 	[nt] of the time points at which to calculate the neutron rate DEFAULT: efit_times
;	
;KEYWORD PARAMETERS:
;	plot:		/plot displays for each time point the input radial profiles and the calculated volumetric reaction rate profile
;
;OUTPUTS:
;	result:		FLTARR [nt] of the volume integrated neutron rate [Hz]
;	
;OPTIONAL OUTPUTS:
;	rmaj:		FLTARR [nr,nt] of the efit_rmid values at each time point [m]
;	vol_rate:	FLTARR [nr,nt] of the reaction rate profile at each time point [Hz/m^3]
;
;MODIFICATION HISTORY:
;	Adapted from code by D. Mikkelson
;	Written by:	ML Reinke 1/7/09
;	8/29/11		ML Reinke - modified the calculation of the volume integrated neutron rate.Before used bolo_power_psi1_svd which doesn't return the correct plasma
;					volume for a rate=1.0, now load and calculate directly from EFIT volp
;
;-

FUNCTION neut_rate,shot,t_i,n_e,h2d,zeff,time=time,rmaj=rmaj,vol_rate=vol_rate,debug=debug,plot=plot
	;each is a structure with *.dat,*.rad,*.time with dat[rad,time] and zeff.z
	;ti in units of keV and ne in units of 10^20 m^-3
	;data is interpolated onto EFIT RMID spatial grid at each time point

	mdsopen,'analysis',shot
	efit_volp=mdsvalue('\analysis::efit_fitout:volp',/quiet,status=ok)
	efit_time=mdsvalue('dim_of(\analysis::efit_fitout:volp,0)')
	efit_rmid=mdsvalue('\efit_rmid')
	mdsclose,'analysis',shot

	IF keyword_set(time) THEN time=time[where(time GE min(efit_time) AND time LE max(efit_time))] ELSE time=efit_time
	vol_rate=fltarr(n(efit_rmid[0,*])+1,n(time)+1)
	rmaj=fltarr(n(efit_rmid[0,*])+1,n(time)+1)
	neut_rate=fltarr(n(time)+1)
	;From H.-S. Bosch and G. M. Hale, , Nucl. Fus. 32 (1992) 611, Table VII
	;For D-D neutron reaction  0.2 keV < Ti < 100 keV
	C1=5.4336E-12
	C2=5.85778E-3
	C3=7.68222E-3
	C4=0.
	C5=-2.964E-6
	C6=0.
	C7=0.
	BG=31.397 ; sqrt(keV)
	mrcsq=9.37814E5 ; keV
	
	n_ti=n(t_i.rad)
	n_ne=n(n_e.rad)
	n_h2d=n(h2d.rad)
	n_zeff=n(zeff.rad)

	FOR i=0,n(time) DO BEGIN
		efit_index=ipt(efit_time,time[i])
		rmaj[*,i]=reform(efit_rmid[efit_index,*])
		
		;interpolate data at time point
		ibnd=ibound(t_i.time,time[i])
		IF ibnd[0] NE ibnd[1] THEN ti_i=t_i.dat[*,ibnd[0]]+(time[i]-t_i.time[ibnd[0]])*(t_i.dat[*,ibnd[1]]-t_i.dat[*,ibnd[0]])/(t_i.time[ibnd[1]]-t_i.time[ibnd[0]]) ELSE $
			ti_i=t_i.dat[*,ibnd[0]]

		ibnd=ibound(n_e.time,time[i])
		IF ibnd[0] NE ibnd[1] THEN ne_i=n_e.dat[*,ibnd[0]]+(time[i]-n_e.time[ibnd[0]])*(n_e.dat[*,ibnd[1]]-n_e.dat[*,ibnd[0]])/(n_e.time[ibnd[1]]-n_e.time[ibnd[0]]) ELSE $
			ne_i=n_e.dat[*,ibnd[0]]
			
		IF n_h2d GT 0 THEN BEGIN
			ibnd=ibound(h2d.time,time[i])
			h2d_i=h2d.dat[*,ibnd[0]]+(time[i]-h2d.time[ibnd[0]])*(h2d.dat[*,ibnd[1]]-h2d.dat[*,ibnd[0]])/(h2d.time[ibnd[1]]-h2d.time[ibnd[0]])		
			h2d_rad=h2d.rad
		ENDIF ELSE BEGIN
			IF total(ibound(h2d.time,time[i])) EQ -1 THEN h2d_i=fltarr(n(rmaj)+1) ELSE h2d_i=fltarr(n(rmaj)+1)+interpol(h2d.dat,h2d.time,time[i])
			h2d_rad=rmaj
		ENDELSE
		IF n_zeff GT 0 THEN  BEGIN
			ibnd=ibound(zeff.time,time[i])
			zeff_i=zeff.dat[*,ibnd[0]]+(time[i]-zeff.time[ibnd[0]])*(zeff.dat[*,ibnd[1]]-zeff.dat[*,ibnd[0]])/(zeff.time[ibnd[1]]-zeff.time[ibnd[0]])		
			zeff_rad=zeff.rad
		ENDIF ELSE BEGIN
			zeff_i=fltarr(n(rmaj)+1)+interpol(zeff.dat,zeff.time,time[i])
			zeff_rad=rmaj
		ENDELSE
		;stop

		;calculate reaction rate
		tmp=where(ti_i LE 0.0)
		IF tmp[0] NE -1 THEN ti_i[tmp]=1e-3				;set <= 0 to be 1 eV
 		theta=ti_i/(1 - ti_i*(C2+ti_i*(C4+ti_i*C6))/(1.+ti_i*(C3+ti_i*(C5+ti_i*C7))))
		ksi=(BG^2/(4.*theta))^(1./3.)
   		sigDD=C1*theta*sqrt(ksi/(mrcsq*ti_i^3))*exp(-3.*ksi)*1.0e-6	;<sigma*v> in units of m^3/s

		;interpolate data onto RMID points
		sigDD=interpol(sigDD,t_i.rad,rmaj[*,i])
		ne_i=interpol(ne_i,n_e.rad,rmaj[*,i])
		h2d_i=interpol(h2d_i,h2d_rad,rmaj[*,i])
		zeff_i=interpol(zeff_i,zeff_rad,rmaj[*,i])

		;calculate neutron rate per unit volume
		deut_dens=ne_i*(zeff.z(0)-zeff_i)/(zeff.z(0)-1.0)/(1.0+h2d_i)
		vol_rate_tmp=double(0.5*deut_dens^2*sigDD*1.0e20)
		vol_rate[*,i]=float(1.0e20*vol_rate_tmp)		;reaction rate 1/m^3

		IF keyword_set(plot) THEN BEGIN
			IF !d.name EQ 'PS' THEN ps=1 ELSE ps=0
			IF ps THEN ls=1.0 ELSE ls=1.0
			multi_old=!p.multi
			!p.multi=[0,0,5]
			IF keyword_set(ps) THEN BEGIN
				d_old=!d
				xsize=7.0
				ysize=10.0
				device, xsize=xsize, ysize=ysize, /inches
			ENDIF ELSE openwin,0,ysize=1100
			plot,rmaj[*,i],ne_i,yr=[0,max(ne_i)],chars=1.4,ytit='n!le!n'
			oplot,n_e.rad,n_e.dat[*,ipt(n_e.time,time[i])],color=200
			plot,rmaj[*,i],interpol(ti_i,t_i.rad,rmaj[*,i]),yr=[min(ti_i),max(ti_i)],chars=1.4,ytit='T!lI!n'
			oplot,t_i.rad,t_i.dat[*,ipt(t_i.time,time[i])],color=200
			plot,rmaj[*,i],sigDD,yr=[min(sigDD),max(sigDD)],chars=1.4,ytit='<sigv>'
			plot,rmaj[*,i],h2d_i,yr=[0,max(h2d_i)],chars=1.4,ytit='H2D'
			plot,rmaj[*,i],zeff_i,yr=[0,max(zeff_i)],chars=1.4,ytit='ZEFF'
			!p.multi=multi_old
			IF keyword_set(ps) THEN device, xsize=float(d_old.x_size)/d_old.x_px_cm,ysize=float(d_old.y_size)/d_old.y_px_cm
		ENDIF		
	
		IF keyword_set(debug) OR (keyword_set(plot) AND NOT keyword_set(ps)) THEN stop
		;bolo_power_psi1_svd,shot,rmaj[*,i],time[i],vol_rate[*,i],ptime,neut_rate_i,ok,ipow,rmidp
		;neut_rate[i]=interpol(neut_rate_i,ptime,time[i])
		volp=reform(efit_volp[efit_index,*])
		nr=n(volp)+1
		dv=volp[1:nr-1]-volp[0:nr-2]
		neut_rate[i]=total(dv*0.5*(vol_rate[0:nr-2,i]+vol_rate[0:nr-1,i]))
		
	ENDFOR
	
	output=neut_rate
	IF keyword_set(debug) THEN stop
	RETURN,output
	
END

;+
;NAME:
;	CALC_NEUT_RATE
;
;PURPOSE:
;	This procedure loads/calculates the necessary data and uses NEUT_RATE to calculate the neutron rate
;	and compare it to experiment for a given shot and time(s).
;
;CALLING SEQUENCE:
;	CALC_NEUT_RATE,shot,times,neuts
;
;INPUTS:
;	shot:		LONG	shot number
;	times:		FLTARR	[nt] of the time points to calculate the neutron rate
;	
;OPTIONAL INPUTS:
;	zo:		INT	of the "average" impurity charge DEFAULT: 5
;	fti:		FLT	of the "fudge factor" to multiply the ;Ti profile by DEFAULT: 1.0
;	fne:		FLT	of the "fudge factor" to multiply the exp. electron density profile DEFAULT: 1.0
;	fh:		FLT	of the minority fraction (really it's forcing h2d)
;	fd:		FLT	sets the nD/ne fraction, overriding fh and zeff/zo information
;
;KEYWORD PARAMETERS:
;	pure:		/pure sets zeff=1 and h2d=0 for an upper limit on the neutron rate
;	equilibrium:	/equilibrium sets Ti=Te (from Thomson)
;	plot:		/plot outputs a time dependent plot of the various inputs and the neutron rate compared to experiment
;	ylog:		/ylog displays the neutron rates in log-space
;	neutplot:	/neutplot sends a /plot command to NEUT_RATE for debugging
;	qfit:		/qfit will load the density profile data using QFIT (sent to ZEFF_NEO)
;	zave:		/zave will take zeff from the experimental measure of zave
;	oldti:		/oldti will load Ti profile data from HIREX_SR rather then HIREXSR
;	bsti:		/bsti will load ti profile using LOAD_BSTI_TREE_DATA
;
;OUPUTS:
;	neuts:		FLTARR	[nt] of the neutron rate [1/s]
;
;PROCEDURE:
;	See NEUT_RATE for details on the calculation of the neutron rate.  h2d is loaded from LOAD_H2D_DATA, 
;	zeff is loaded from ZEFF_NEO, Ti is loaded from LOAD_TI_TREE_DATA and ne is loaded from the Thomson data
;	used in ZEFF_NEO
;
;MODIFICATION HISTORY:
;	Written by:	ML Reinke 1/8/09
;	6/8/11		ML Reinke - updated documentation and changed default Ti to HIREXSR nodes.  Use /old to get data from HIREX_SR
;	8/5/11		ML Reinke - added the ability to force the minority fraction (fh) and use the error bars from HIREXSR data
;	9/13/11		ML Reinke - fixed a bug in the HIREXSR Ti loading code that removed the instrumental twice, making predictions drop.
;				    also modified some plotting features
;	7/25/11		M.L. Reinke - added the ability to use the new sbspline Ti profiles using /bsti and specificy the deuteron fraction via fd
;
;-

PRO calc_neut_rate,shot,times,neuts,zo=zo,debug=debug,pure=pure,plot=plot,ylog=ylog,pti=pti,fti=fti,fne=fne,neutplot=neutplot,equilibrium=equilibrium,$
	qfit=qfit,zave=zave,oldti=oldti,bsti=bsti,fh=fh,tht=tht,comb=comb,det=det,fd=fd

	IF NOT keyword_set(fti) THEN fti=1.0
	IF NOT keyword_set(fne) THEN fne=1.0
	IF NOT keyword_set(pti) THEN pti=0.0
	IF NOT keyword_set(zo) THEN zo=5.0
	IF !d.name EQ 'PS' THEN ps=1 ELSE ps=0
	load_h2d_data,shot,h2d,time
	zeff_neo,shot,zeff,timez,ts=ts,efit=efit,qfit=qfit
	IF keyword_set(zave) THEN BEGIN
		load_zave_data,shot,zeff,timez
		npt=n(timez)
		timez=timez[indgen(npt/10)*10]
		zeff=zeff[indgen(npt/10)*10]
	ENDIF
	IF keyword_set(pure) THEN BEGIN 
		zeff/=zeff
		h2d*=0.0
        ENDIF
	IF keyword_set(fd) THEN BEGIN
		zeff=zeff*0.0+1.0+(zo-1.0)*(1.0-fd)
		h2d*=0.0
	ENDIF
	zeff={dat:zeff,rad:[0],time:timez,z:[zo]}
	IF keyword_set(fh) THEN h2d={dat:[time*0.0+fh,fh],rad:[0],time:[time,2.0]} ELSE h2d={dat:h2d,rad:[0],time:time}
	n_e={dat:rotate(ts.ner,4)*1.0e-20*fne,rad:efit.rgrid,time:efit.tgrid}
	;stop
	IF keyword_set(equilibrium) THEN BEGIN
		t_i={dat:rotate(ts.te*1.0e-3,4)*fti,rad:efit.rgrid,time:efit.tgrid,err:0.0*ts.te}
	ENDIF ELSE BEGIN
		fitcase=0
		IF keyword_set(oldti) THEN fitcase=1
		IF keyword_set(bsti) THEN fitcase=2
		CASE fitcase OF 
			0 : load_tinew_tree_data,shot,ti,r,t,redge=redge,tierr=tierr,tht=tht,comb=comb
			1 : BEGIN
				load_ti_tree_data,shot,ti,r,t,/interp
				redge=0.0
				tierr=0.0*ti
                        END
			2 : load_bsti_tree_data,shot,ti,r,t,redge=redge,tierr=tierr,tht=tht
		ENDCASE
		IF ti[0] EQ -1 THEN print, 'No Ti Data Available'
		t_i={dat:ti*fti+pti,rad:r,time:t,err:tierr}
	ENDELSE
	IF total(t_i.err) NE 0 THEN BEGIN											;compute w/ T_i +/- unc. in T_i
		err=1
		t_i.dat+=t_i.err
		neuts_up=neut_rate(shot,t_i,n_e,h2d,zeff,time=times,rmaj=rmaj,vol_rate=vol_rate,debug=debug,plot=neutplot)
		t_i.dat-=2.0*t_i.err
		neuts_low=neut_rate(shot,t_i,n_e,h2d,zeff,time=times,rmaj=rmaj,vol_rate=vol_rate,debug=debug,plot=neutplot)
		t_i.dat+=t_i.err
	ENDIF ELSE err=0
	neuts=neut_rate(shot,t_i,n_e,h2d,zeff,time=times,rmaj=rmaj,vol_rate=vol_rate,debug=debug,plot=neutplot)
	IF keyword_set(plot) THEN BEGIN
		load_neutrate_data,shot,neuts_m,time_m
		load_he3neutrate_data,shot,neuts_he3,time_he3,det=det
		neuts*=1.0e-10
		neuts_m*=1.0e-10
		IF ps THEN ls=1.0 ELSE ls=1.0
		!p.multi=[0,0,5]
		IF keyword_set(ps) THEN BEGIN
			d_old=!d
			xsize=6.0
			ysize=10.0
			device, xsize=xsize, ysize=ysize, /inches
		ENDIF ELSE openwin,0,ysize=1100
		IF NOT keyword_set(xr) THEN xr=[0.4,1.6]
		plot, n_e.time,n_e.dat[0,*],xr=xr,/xst,chars=1.8*ls,tit='SHOT '+num2str(shot,1),ytit='Core n!le!n [10!u20!n m!u-3!n]',yr=[0,max(n_e.dat[0,*])*1.1],/ysty
		maxpt=max(t_i.dat[0,*])*1.1
		plot, t_i.time,t_i.dat[0,*],xr=xr,/xst,chars=1.8*ls,ytit='T!li!n [keV]',yr=[0,maxpt],/ysty
		index=ipt(t_i.rad,0.86)
		oplot, t_i.time,t_i.dat[index,*],linestyle=2.0
		IF err THEN BEGIN
			oploterror, t_i.time,t_i.dat[0,*],t_i.err[0,*]
			oploterror, t_i.time,t_i.dat[index,*],t_i.err[index,*]
		ENDIF
		IF fti NE 1.0 THEN xyouts,0.15,0.3,'ff = '+num2str(fti*100,dp=1)+' [%]'
		oplot, [0.1,0.2],0.82*maxpt*[1.0,1.0]
		xyouts,xr[0]+0.25,0.8*maxpt,'R='+num2str(t_i.rad[0],dp=2)+' [m]',chars=0.9*ls
		oplot, [0.1,0.2],0.67*maxpt*[1.0,1.0],linestyle=2.0
		xyouts,xr[0]+0.25,0.65*maxpt,'R=0.87 [m]',chars=0.9*ls
		IF keyword_set(zave) THEN psym=0 ELSE psym=-5
		tmp=where(zeff.time GE xr[0] AND zeff.time LE xr[1])
		IF tmp[0] NE -1 THEN zeffmax = 4 > max(zeff.dat[tmp]) ELSE zeffmax=4.0
		plot, zeff.time,zeff.dat,xr=xr,/xst,chars=1.8*ls,ytit='Z!leff!n',psym=psym,yr=[0,zeffmax],/ysty
		xyouts,xr[1]-0.25,1.0,'Z='+num2str(zo[0]),chars=1.2*ls
		tmp=where(h2d.time GE times[0] AND h2d.time LE last(times))
		IF tmp[0] NE -1 THEN plot, h2d.time,h2d.dat*100.0,xr=[0,2.0],/xst,chars=1.8*ls,ytit='H/D Ratio [%]',psym=-5,yr=[0,max(h2d.dat[tmp])*1.05]*100.0,/ysty ELSE $
			plot, h2d.time,h2d.dat*100.0,xr=[0,2.0],/xst,chars=1.8*ls,ytit='H/D Ratio [%]',psym=-5
		maxplot=max(neuts) > max(neuts_m[where(time_m GT 0.5 AND time_m LT 1.5)])
		IF keyword_set(ylog) THEN yr=[10.0,maxplot*2.0] ELSE yr=[0,maxplot*1.2]
		plot, times,neuts,xr=xr,/xst,psym=-5,chars=1.8*ls,xtit='Time [sec]',ytit='Neutron Rate 10!u10!n /s',yr=yr,ylog=ylog,/ysty
		IF err THEN oploterror,times,neuts,0.5*(neuts_up-neuts_low)*1.0e-10
		oplot, time_m,neuts_m,color=200
		dt=0.1
		tnorm=make(0.5,1.5,1.0/dt)
		ynorm=fltarr(n(tnorm)+1)
		FOR i=0,n(tnorm) DO BEGIN
			ym=mean(neuts_m[ipt(time_m,tnorm[i]-dt/2.0):ipt(time_m,tnorm[i]+dt/2.0)])
			yhe3=mean(neuts_he3[ipt(time_he3,tnorm[i]-dt/2.0):ipt(time_he3,tnorm[i]+dt/2.0)])
			ynorm[i]=ym/yhe3
		ENDFOR
		tmp=where(time_he3 GE tnorm[0] AND time_he3 LE last(tnorm))
		IF det GT 0 THEN oplot,time_he3[tmp],interpol(ynorm,tnorm,time_he3[tmp])*neuts_he3[tmp],color=150
		neuts_marmar=10^(make(9,15,100))
		wmhd = neuts_to_wmhd(neuts_marmar)
		load_wth_data,shot,wth,wtime
		tmp=where(wth GT 10e3)		;10 kJ floor
		IF tmp[0] NE -1 THEN oplot, wtime[tmp],interpol(neuts_marmar,wmhd,wth[tmp])*1.0e-10,color=100
		!p.multi=0
		stop
		IF keyword_set(ps) THEN device, xsize=float(d_old.x_size)/d_old.x_px_cm,ysize=float(d_old.y_size)/d_old.y_px_cm
	ENDIF
	IF keyword_set(debug) THEN stop
END


	
