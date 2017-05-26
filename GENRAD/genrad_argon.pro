;+
;NAME:
;	AR_XRAY_EMISS
;
;PURPOSE:
;	This function calculates the emissivity profiles for various
;	groups of lines observed by the HIREXSR spectrometer
;
;CALLING SEQUENCE:
;	result=AR_XRAY_EMISS(csden,dens,temp)
;	
;INPUTS:
;	csden	FLTARR	[19,nrho] of the charge state densities [18,*] (fully stripped)
;	dens	FLTARR	[nrho] of the electron density [m^-3]
;	temp	FLTARR	[nrho] of the electron temperature [eV]
;
;OPTIONAL INPUTS:
;	cserr	FLTARR	[19,nrho] of the uncertainty in csden
;	signe	FLTARR	[nrho] of the uncertainty in the electron density [m^-3]
;	sigte	FLTARR	[nrho] of the uncertainty in the electron temperature [keV]
;	wl_roi	FLTARR	[2,nlines] of the wavelength regions (see PROCEDURE for DEFAULTS)
;
;OUTPUTS:
;	result	[nrho,nlines] of the emissivities [ph/s/m^3] integrated over the wavelength ROI's
;
;OPTIONAL OUTPUTS:
;	emerr	[nrho,nlines] of the uncertainty in the emissivity [ph/s/m^3]
;	csemiss	[19,nrho,nlines] of the emission from each charge state making up the emissivity
;
;PROCEDURE:
;	This function calls CALC_GROUP_RATES in read_ark_tables and
;	finds the following emissivity profiles:
;		w+n>=3	[3.9455,3.9600]
;		q+r+a	[3.9780,3.9880]
;		k+j+z	[3.9880,4.0000]
;		lya1+2	[3.7270,3.7385]
;
;MODIFICATION HISTORY:
;	Written by:	M.L. Reinke (2009?)
;
;-

FUNCTION ar_xray_emiss,csden,dens,temp,wl_roi=wl_roi,cserr=cserr,sigte=sigte,signe=signe,emerr=emerr,csemiss=csemiss
	IF NOT keyword_set(cserr) THEN cserr=csden*0.0
	IF NOT keyword_set(signe) THEN signe=0.0
	IF NOT keyword_set(sigte) THEN sigte=0.0	

	z=n(csden[*,0])
	nrho=n(dens)+1
	IF NOT keyword_set(wl_roi) THEN $
		wl_roi=[[3.9455,3.9600],$		;w+n>=3
			[3.9780,3.9880],$		;q,r,a
			[3.9880,4.0000],$		;k+j+z
			[3.7270,3.7385]]		;lya1+lya2	
	nemiss=n(wl_roi[0,*])+1
	
	emiss=fltarr(nrho,nemiss)
	emerr=fltarr(nrho,nemiss)
	csemiss=fltarr(z+1,nrho,nemiss)
	FOR i=0,nemiss-1 DO BEGIN
		calc_group_rates, reform(wl_roi[*,i]), te_rates, rates
		FOR j=0,z DO IF total(rates[j,*]) NE 0 THEN BEGIN
			rates_int=interpol(rates[j,*],te_rates,temp/1.0e3)
			drdte_int=interpol(deriv(te_rates,rates[j,*]),te_rates,temp/1.0e3)
			csemiss[j,*,i]=rates_int*dens*csden[j,*]
			emiss[*,i]+=rates_int*dens*csden[j,*]
			emerr[*,i]+=(dens*csden[j,*]*drdte_int)^2*(sigte/1.0e3)^2+(rates_int*dens)^2*(cserr[j,i])^2+(rates_int*csden[j,*])^2*(signe)^2
		ENDIF
	ENDFOR
	emerr=sqrt(emerr)
	output=emiss
	RETURN,output
	
END

;calculates a mock spectra from the READ_ARK_TABLE files
FUNCTION ar_xray_spec,lamr,csden,temp,nlam=nlam,ti=ti
	IF NOT keyword_set(nlam) THEN nlam=200
	IF NOT keyword_set(ti) THEN ti=temp
	z=18
	ar=reform_ark_data(lamr,/load)
	n_lines=n(ar.num_elec)+1
	linelam=ar.lam
	lam=make(min(linelam),max(linelam),nlam)
	rates=fltarr(z+1,n_lines)
	FOR i=0,n_lines-1 DO BEGIN
		q=z-ar.num_elec[i]
		rates[q+1,i]+=interpol(ar.rec[*,i],ar.temp,temp/1.0e3)
		rates[q,i]+=interpol(ar.exc[*,i],ar.temp,temp/1.0e3)
		rates[q-1,i]+=interpol(ar.ion[*,i],ar.temp,temp/1.0e3)
	ENDFOR


	c=3.0e8 			;speed of light
	e=1.60e-19			;conversion for eV -> J
	mconv=1.66e-27			;conversion for amu -> kg
	mass=read_atomic_mass(z)

	spec_em=fltarr(nlam)
	FOR i=0,n_lines-1 DO BEGIN
		em=total(rates[*,i]*csden)
		conv_factor=(linelam[i]/c)^2*(e/(mass*mconv))
		w=sqrt(ti*conv_factor)
		spec_em+=em/(w*sqrt(2.0*!pi))*exp(-(lam-linelam[i])^2/(2.0*w^2))
	ENDFOR

	output={spec:spec_em,lam:lam}
	RETURN,output
END

;compares photon emission coefficients from
;HULLAC and ADAS
;6/8/12 changed path to /usr/local/cmod/idl/atomic_physics/adas
PRO ar_pec_comp
	
	hull_he=read_brtable2_data(find_path(18,16,'lines'))
	hull_h=read_brtable2_data(find_path(18,17,'lines'))
	hull_li=read_brtable2_data(find_path(18,15,'lines'))
	hull_be=read_brtable2_data(find_path(18,14,'lines'))
	;bp='/home/mlreinke/idl/impurities/data/adas/'
	bp='/usr/local/cmod/idl/atomic_physics/adas/'
 	adas_h=read_pec_file(bp+'transport_llu#ar17.dat')
	adas_he=read_pec_file(bp+'transport_llu#ar16.dat')
	adas_li=read_pec_file(bp+'transport_llu#ar15ic.dat')
	adas_be=read_pec_file(bp+'transport_llu#ar14ic.dat')

	lines=['w','x','y','z']
	FOR i=0,n(lines) DO BEGIN
		read_ar_line_rates,lines[i],ark_te,rec,exc,ion,/load
		IF i EQ 0 THEN BEGIN
			ark_he_exc=fltarr(n(ark_te)+1)
			ark_he_rec=fltarr(n(ark_te)+1)
			ark_he_ion=fltarr(n(ark_te)+1)
		ENDIF
		ark_he_exc+=exc
		ark_he_rec+=rec
		ark_he_ion+=ion
	ENDFOR
	ark_he_te=ark_te
	out=reform_ark_data([3.73,3.74],/load)
	tmp=where(out.num_elec EQ 1)
	ark_h_te=out.temp
	ark_h_exc=sum_array(out.exc[*,tmp],/i)

	openwin,0
	nept=1.0e20
	makesym,10
	ysc=1.0e17
	plot, adas_h.temp,adas_h.pec[0,*,ipt(adas_h.dens,nept)]*ysc,/xlog,yr=[0,1.5],xr=[1.0e2,1.0e4],/ysty,xtit='Electron Temperature [eV]',$
		ytit='<'+n2g('sigma')+'v> [10!u-17!n m!u3!n/s]',chars=1.2
	ipt=ipt(hull_h.dens,nept)
	hullpt=where(hull_h.lam GT 3.7 AND hull_h.lam LT 3.8)
	hull_em=sum_array(hull_h.em[*,*,hullpt],/k)
	oplot, hull_h.temp,hull_em[*,ipt]/hull_h.dens[ipt]*ysc,color=200,psym=-8


	adas_em=adas_he.pec[0,*,ipt(adas_he.dens,nept)]+adas_he.pec[1,*,ipt(adas_he.dens,nept)]+adas_he.pec[2,*,ipt(adas_he.dens,nept)]
	oplot, adas_he.temp,adas_em*ysc,linestyle=2.0
	ipt=ipt(hull_he.dens,nept)
	hullpt=where(hull_he.lam GT 3.9 AND hull_he.lam LT 4.0)
	hull_em=sum_array(hull_he.em[*,*,hullpt],/k)
	oplot, hull_he.temp,hull_em[*,ipt]/hull_he.dens[ipt]*ysc,color=200,psym=-8,linestyle=2.0

	oplot,ark_he_te*1.0e3,ark_he_exc*ysc,color=30,linestyle=2.0
	oplot,ark_h_te*1.0e3,ark_h_exc*ysc,color=30
	;oplot,ark_te*1.0e3,ark_ion,color=30,linestyle=2.0
	;oplot,ark_te*1.0e3,ark_rec,color=30,linestyle=2.0
	
	
	hull_em_he=fltarr(n(hull_he.temp)+1)
	FOR i=0,n(hull_he.lam) DO hull_em_he+=hull_he.em[*,ipt(hull_he.dens,nept),i]*ang2ev(hull_he.lam[i])*1.6e-19/hull_he.dens[ipt(hull_he.dens,nept)]
	hull_em_h=fltarr(n(hull_he.temp)+1)
	FOR i=0,n(hull_h.lam) DO hull_em_h+=hull_h.em[*,ipt(hull_h.dens,nept),i]*ang2ev(hull_h.lam[i])*1.6e-19/hull_h.dens[ipt(hull_h.dens,nept)]
	hull_em_li=fltarr(n(hull_li.temp)+1)
	FOR i=0,n(hull_li.lam) DO hull_em_li+=hull_li.em[*,ipt(hull_li.dens,nept),i]*ang2ev(hull_li.lam[i])*1.6e-19/hull_li.dens[ipt(hull_li.dens,nept)]
	hull_em_be=fltarr(n(hull_be.temp)+1)
	FOR i=0,n(hull_be.lam) DO hull_em_be+=hull_be.em[*,ipt(hull_be.dens,nept),i]*ang2ev(hull_be.lam[i])*1.6e-19/hull_be.dens[ipt(hull_be.dens,nept)]

	adas_em_he=fltarr(n(adas_he.temp)+1)
	FOR i=0,n(adas_he.lam) DO adas_em_he+=adas_he.pec[i,*,ipt(adas_he.dens,nept)]*ang2ev(adas_he.lam[i])*1.6e-19
	adas_em_h=fltarr(n(adas_h.temp)+1)
	FOR i=0,n(adas_h.lam) DO adas_em_h+=adas_h.pec[i,*,ipt(adas_h.dens,nept)]*ang2ev(adas_h.lam[i])*1.6e-19
	adas_em_li=fltarr(n(adas_li.temp)+1)
	FOR i=0,n(adas_li.lam) DO adas_em_li+=adas_li.pec[i,*,ipt(adas_li.dens,nept)]*ang2ev(adas_li.lam[i])*1.6e-19
	adas_em_be=fltarr(n(adas_be.temp)+1)
	FOR i=0,n(adas_be.lam) DO adas_em_be+=adas_be.pec[i,*,ipt(adas_be.dens,nept)]*ang2ev(adas_be.lam[i])*1.6e-19


	openwin,1
	ysc=1.0e33
	plot,[0],[0],/xlog,yr=[0,50.0],xr=[10.0,1.0e4],/ysty,xtit='Electron Temperature [eV]',$
		ytit='PLC [10!u-33!n W*m!u3!n]',chars=1.2
	oplot,hull_he.temp,hull_em_he*ysc,psym=-8,color=200,linestyle=2
	oplot,hull_h.temp,hull_em_h*ysc,psym=-8,color=200
	oplot,hull_li.temp,hull_em_li*ysc,psym=-8,color=200,linestyle=3
	oplot,hull_be.temp,hull_em_be*ysc,psym=-8,color=200,linestyle=4

	oplot,adas_he.temp,adas_em_he*ysc,linestyle=2.0
	oplot,adas_h.temp,adas_em_h*ysc,linestyle=2.0
	oplot,adas_li.temp,adas_em_li*ysc,linestyle=3.0
	oplot,adas_be.temp,adas_em_be*ysc,linestyle=4.0
	
	openwin,2
	plot,[0],[0],/xlog,xr=[10.0,1.0e4],/ysty,xtit='Electron Temperature [eV]',$
			ytit='PLC [10!u-33!n W*m!u3!n]',chars=1.2,yr=[1.0e-34,1.0e-31],/ylog
	;path='/home/mlreinke/idl/impurities/data/adas/plt89_ar.dat'
	path='/usr/local/cmod/idl/atomic_physics/adas/plt89_ar.dat'
	plt=read_pxx_file(path)

	FOR i=17,17 DO BEGIN
		oplot,plt.temp,plt.plc[i,*,ipt(plt.dens,nept)]
		path=find_path(18,i,'lines')
		IF path NE 'NA' THEN BEGIN
			hull=read_brtable2_data(path)
			hull_em=fltarr(n(hull.temp)+1)
			FOR j=0,n(hull.lam) DO hull_em+=hull.em[*,ipt(hull.dens,nept),j]*ang2ev(hull.lam[j])*1.6e-19/hull.dens[ipt(hull.dens,nept)]
			oplot, hull.temp,hull_em,psym=-8,color=200
		ENDIF
	ENDFOR


	stop
END

;compares ADAS and HULLAC based cooling curves
PRO fqcsplc_comp
	z=18
	t_e=10^(make(0.0,4.0,500))
	n_e=[1.0e20]
	
	;path='/home/mlreinke/idl/impurities/data/adas/plt89_ar.dat'
	path='/usr/local/cmod/idl/atomic_physics/adas/plt89_ar.dat'
	plt=read_pxx_file(path)
	ipt=ipt(plt.dens,n_e[0])
	adas_plc=fltarr(z+1,n(t_e)+1)
	FOR i=0,z DO adas_plc[i,*]=interpol(plt.plc[i,*,ipt],plt.temp,t_e)
	adas_temp=t_e

	hull_plc=hullac_cs_plc(t_e,z,filter=filter,n_e=n_e[0],verb=verb,noextrap=noextrap)	
	hull_temp=t_e
	ion=read_ion_data(18,/adas)
 	rec=read_rec_data(18,/adas)
	fq=calc_fracabund(ion,rec,te=adas_temp)
	frac_adas=fq.fq

	ion=read_ion_data(18)
 	rec=read_rec_data(18)
	fq=calc_fracabund(ion,rec,te=hull_temp)
	frac_loch=fq.fq

	fq,z,cs,cs_te
	frac_hull=fltarr(z+1,n(hull_temp)+1)
	FOR i=0,z DO frac_hull[i,*]=interpol(cs[i,*],cs_te,hull_temp)

	openwin,0
	IF keyword_set(ylog) THEN ymin=1.0e-35 ELSE ymin=0.0
	plot,[0],[0],xr=[1,100],/xlog,yr=[ymin,1.5e-31],ylog=ylog
	;FOR i=0,8 DO oplot,adas_temp,adas_plc[i,*]*frac_adas[i,*]
	FOR i=0,8 DO BEGIN
		tmp=where(hull_plc[i,*] NE 0)
		IF tmp[0] NE -1 THEN BEGIN
			oplot,hull_temp[tmp],hull_plc[i,tmp]*frac_hull[i,tmp],color=200 
			oplot,hull_temp[tmp],hull_plc[i,tmp]*frac_loch[i,tmp],color=100
		ENDIF
	ENDFOR

	openwin,1
	plot,[0],[0],xr=[50,5000],/xlog,yr=[ymin,2.5],ylog=ylog,/xsty,xtit='T!le!n [eV]', ytit='Line Radiation Loss [10!u-32!n Wm!u3!n]',chars=1.2
	;FOR i=8,z DO oplot,adas_temp,adas_plc[i,*]*frac_adas[i,*]
	shift=[-1,-1,-1,0,1,1,1,1,1,1]
	FOR i=8,z DO BEGIN
		tmp=where(hull_plc[i,*] NE 0)
		IF tmp[0] NE -1 THEN BEGIN
			oplot,hull_temp[tmp],hull_plc[i,tmp]*frac_hull[i,tmp]*1.0e32,color=200
			oplot,hull_temp[tmp],hull_plc[i,tmp]*frac_loch[i,tmp]*1.0e32,color=100
			maxpt=(max(hull_plc[i,tmp]*frac_hull[i,tmp]) > hull_plc[i,tmp]*frac_loch[i,tmp])*1.0e32
			maxloc=hull_temp[tmp[maxloc(hull_plc[i,tmp]*frac_hull[i,tmp])]]
			xyouts,maxloc,maxpt+0.1,num2rom(i+1)
		ENDIF
	
	ENDFOR

	openwin,2
	plot,[0],[0],xr=[500.0,5.0e3],/xlog,yr=[ymin,3.0e-33],ylog=ylog,/xsty
	;FOR i=z-4,z DO oplot,adas_temp,adas_plc[i,*]*frac_adas[i,*]
	FOR i=z-4,z DO BEGIN
		tmp=where(hull_plc[i,*] NE 0)
		IF tmp[0] NE -1 THEN BEGIN
			oplot,hull_temp[tmp],hull_plc[i,tmp]*frac_hull[i,tmp],color=200
			oplot,hull_temp[tmp],hull_plc[i,tmp]*frac_loch[i,tmp],color=100
		ENDIF
	ENDFOR

	hull_tot=fltarr(n(hull_temp)+1)
	FOR i=0,z DO hull_tot+=hull_plc[i,*]*frac_hull[i,*]
	loch_tot=fltarr(n(hull_temp)+1)
	FOR i=0,z DO loch_tot+=hull_plc[i,*]*frac_loch[i,*]
	adas_tot=fltarr(n(adas_temp)+1)
	FOR i=0,z DO adas_tot+=adas_plc[i,*]*frac_adas[i,*]
	openwin,4
	plot,[0],[0],xr=[50.0,5.0e3],/xlog,yr=[ymin,6.0e-32],ylog=ylog,/xsty

	oplot,adas_temp,adas_tot
	oplot,hull_temp,hull_tot,color=100
	oplot,hull_temp,loch_tot,color=200

	openwin,5
	ymin=1.0e-33
	plot,[0],[0],xr=[1.0,5000.0],/xlog,yr=[0,2.0e-31],/xsty,xtit='T!le!n [eV]',ytit='Line Radiation Loss [Wm!u3!n]',chars=1.2
	oplot,adas_temp,adas_tot
	;oplot,hull_temp,hull_tot,color=100
	oplot,hull_temp,loch_tot,color=200

	stop
END
		
	
