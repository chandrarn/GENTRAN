;6/8/12 - modified to use new /usr/local/cmod/idl/atomic_physics/ data 
FUNCTION foil_calc_line_csplc,z,debug=debug,load=load,verb=verb,adas=adas

	t_e=10^(make(0.0,4.0,500))
	n_e=[1.0e20]
	CASE z OF
		18 : BEGIN
			;path_save='/home/mlreinke/idl/impurities/data/hullac_kbf/ar/ar_foil_csplc.dat'
			;IF keyword_set(adas) THEN path_save='/home/mlreinke/idl/impurities/data/adas/ar_foil_csplc.dat'
			path_load='/usr/local/cmod/idl/atomic_physics/hullac_kbf/ar/ar_foil_csplc.dat'
			IF keyword_set(adas) THEN path_load='/usr/local/cmod/idl/atomic_physics/adas/ar_foil_csplc.dat'
			path_save='/home/'+logname()+'/usr/'+logname()+'/ar_foil_csplc.dat'		
		END
		36 : BEGIN
			;path_save='/home/mlreinke/idl/impurities/data/hullac_kbf/kr/kr_foil_csplc.dat'
			;IF keyword_set(adas) THEN path_save='/home/mlreinke/idl/impurities/data/adas/kr_foil_csplc.dat'
			path_load='/usr/local/cmod/idl/atomic_physics/hullac_kbf/kr/kr_foil_csplc.dat'
			IF keyword_set(adas) THEN path_load='/usr/local/cmod/idl/atomic_physics/adas/kr_foil_csplc.dat'
			path_save='/home/'+logname()+'/usr/'+logname()+'/kr_foil_csplc.dat'
		END
		ELSE : RETURN,-1
		
	ENDCASE
	IF NOT keyword_set(load) THEN BEGIN
		IF keyword_set(adas) THEN csplc=adas_pxx_csplc(z,t_e,n_e,verb=verb) ELSE csplc=hullac_cs_plc(t_e,z,n_e=n_e[0],verb=verb)	
		save,csplc,t_e,n_e,filename=path_save
	ENDIF ELSE restore,path_load

	output={csplc:csplc,temp:t_e,dens:n_e}
	IF keyword_set(debug) THEN stop
	RETURN,output
END

FUNCTION foil_calc_cxr_csplc,z,debug=debug,load=load,verb=verb,adas=adas

	t_e=10^(make(0.0,4.0,500))
	n_e=[1.0e20]
	csplc=adas_pxx_csplc(z,t_e,n_e,verb=verb,/prc)

	output={csplc:csplc,temp:t_e,dens:n_e}
	IF keyword_set(debug) THEN stop
	RETURN,output
END

FUNCTION foil_calc_cont_csplc,z,debug=debug,load=load,verb=verb,adas=adas

	t_e=10^(make(0.0,4.0,500))
	n_e=[1.0e20]
	csplc=adas_pxx_csplc(z,t_e,n_e,verb=verb,/prb)

	output={csplc:csplc,temp:t_e,dens:n_e}
	IF keyword_set(debug) THEN stop
	RETURN,output
END


PRO foil_genrad_em2br,shot,t,emiss,rmaj,br,pos,gpv=gpv

	IF NOT keyword_set(gpv) THEN BEGIN
		pos=foil_midplane_pos(etendue=u)
		pos=pos[*,0:15]
		order=sort(pos[2,*])
		pos=pos[*,order]	;resort
		u=u[0:15]
		u=u[order]
		restore,'/home/mlreinke/idl/genie/data/gpv/foil/genrad_foil_gpv.dat'
		gpv=gpv[0:15,*]
		gpv=gpv[order,*]
		em_grid=grid_profile(ves_cent,emiss,rmaj,[t],shot,tpts=t,/sol)
		br=fltarr(16)
		FOR i=0,n(br) DO BEGIN
			non_zero=where(gpv[i,*] GT 0)
			IF non_zero[0] NE -1 THEN br[i]=total(gpv[i,non_zero]*em_grid[non_zero])/u[i]*4.0*!pi
		ENDFOR

	ENDIF ELSE BEGIN
		pos=foil_midplane_pos()
		pos=pos[*,0:15]
		order=sort(pos[2,*])
		pos=pos[*,order]	;resort
		br=line_br(pos,emiss,rmaj,[t],shot,t,n_s=n_s,plots=debug)
	ENDELSE
END

FUNCTION foil_genrad_line_emiss,shot,t,rmaj,csden,cserr,temp,temperr,dens,denserr,n_s=n_s,debug=debug,adas=adas,gpv=gpv

	z=n(csden[*,0])
	nr=n(rmaj)+1
	csplc=foil_calc_line_csplc(z,/load,adas=adas)
	emiss=fltarr(nr)
	emerr=fltarr(nr)	
	FOR i=0,z DO BEGIN
		csint=interpol(csplc.csplc[i,*],csplc.temp,temp)
		;stop
		dcsint=interpol(deriv(csplc.temp,csplc.csplc[i,*]),csplc.temp,temp)
		emiss+=csint*dens*csden[i,*]
		emerr+=(dens*csden[i,*]*dcsint)^2*temperr^2+(csint*dens)^2*(cserr[i,*])^2+(csint*csden[i,*])^2*denserr^2
	ENDFOR
	emerr=sqrt(emerr)
	foil_genrad_em2br,shot,t,emiss,rmaj,br,pos,/gpv
	foil_genrad_em2br,shot,t,emerr,rmaj,brerr,pos,/gpv
	
	output={rmaj:rmaj,em:emiss,emerr:emerr,pos:pos,br:br,brerr:brerr}
	IF keyword_set(debug) THEN stop
	RETURN,output
END

FUNCTION foil_genrad_cxr_emiss,shot,t,rmaj,csden,cserr,temp,temperr,neut,neuterr,n_s=n_s,debug=debug,adas=adas,gpv=gpv

	z=n(csden[*,0])
	nr=n(rmaj)+1
	csplc=foil_calc_cxr_csplc(z,/load,adas=adas)
	emiss=fltarr(nr)
	emerr=fltarr(nr)		
	FOR i=0,z DO BEGIN
		csint=interpol(csplc.csplc[i,*],csplc.temp,temp)
		dcsint=interpol(deriv(csplc.temp,csplc.csplc[i,*]),csplc.temp,temp)
		emiss+=csint*neut*csden[i,*]
		emerr+=(neut*csden[i,*]*dcsint)^2*temperr^2+(csint*neut)^2*(cserr[i,*])^2+(csint*csden[i,*])^2*neuterr^2
	ENDFOR
	emerr=sqrt(emerr)

	foil_genrad_em2br,shot,t,[emiss,0.0],[rmaj,1.0],br,pos,/gpv
	foil_genrad_em2br,shot,t,[emerr,0.0],[rmaj,1.0],brerr,pos,/gpv
	output={rmaj:rmaj,em:emiss,emerr:emerr,pos:pos,br:br,brerr:brerr}

	IF keyword_set(debug) THEN stop
	RETURN,output
END

FUNCTION foil_genrad_cont_emiss,shot,t,rmaj,csden,cserr,temp,temperr,dens,denserr,n_s=n_s,debug=debug,adas=adas,gpv=gpv
	z=n(csden[*,0])
	nr=n(rmaj)+1
	csplc=foil_calc_cont_csplc(z,/load,adas=adas)
	emiss=fltarr(nr)
	emerr=fltarr(nr)		
	FOR i=0,z DO BEGIN
		csint=interpol(csplc.csplc[i,*],csplc.temp,temp)
		dcsint=interpol(deriv(csplc.temp,csplc.csplc[i,*]),csplc.temp,temp)
		emiss+=csint*dens*csden[i,*]
		emerr+=(dens*csden[i,*]*dcsint)^2*temperr^2+(csint*dens)^2*(cserr[i,*])^2+(csint*csden[i,*])^2*denserr^2
	ENDFOR
	emerr=sqrt(emerr)

	foil_genrad_em2br,shot,t,[emiss,0.0],[rmaj,1.0],br,pos,/gpv
	foil_genrad_em2br,shot,t,[emerr,0.0],[rmaj,1.0],brerr,pos,/gpv
	output={rmaj:rmaj,em:emiss,emerr:emerr,pos:pos,br:br,brerr:brerr}

	IF keyword_set(debug) THEN stop
	RETURN,output
END

FUNCTION foil_genrad_brem_emiss,shot,t,rmaj,temp,temperr,dens,denserr,n_s=n_s,debug=debug,gpv=gpv,zeff=zeff

	IF NOT keyword_set(zeff) THEN zeff=1.0
	eph=10^(make(0.0,5.0,200))
	lam=ev2ang(eph)
	npts=n(rmaj)+1
	emiss=fltarr(npts)
	emerr=fltarr(npts)
	h=6.626068e-34	;m^2 kg/s
	c=2.9979e8	;m/s
	e=1.60218e-19	;J/eV
	kb=8.617e-5 	;eV/K
	
	FOR i=0,npts-1 DO BEGIN
		log_gff=0.355*lam^(-0.06)*alog10(lam)+0.3*lam^(-0.066)*alog10(temp[i]/kb*1.0e-6/100.0)+0.0043
		specem=5.03e-14/h*e*4.0*!pi*10^(log_gff)*zeff*dens[i]/1.0e20*dens[i]/1.0e20/sqrt(temp[i])*exp(-eph/temp[i])	;watts/m^3/eV from Hutch
		emiss[i]=int_tabulated(eph,specem)	;W/m^3
		emerr[i]=sqrt((2.0*emiss[i]/dens[i])^2*denserr[i]^2+(0.5*emiss[i]/temp[i])^2*temperr[i]^2)
	ENDFOR

	foil_genrad_em2br,shot,t,emiss,rmaj,br,pos,/gpv
	foil_genrad_em2br,shot,t,emerr,rmaj,brerr,pos,/gpv
	
	output={rmaj:rmaj,em:emiss,emerr:emerr,pos:pos,br:br,brerr:brerr}
	IF keyword_set(debug) THEN stop
	RETURN,output
END



;breaking this program to work with a newer version of the bolometry tree 8/7/11
PRO foil_genrad_profiles,shot,t1,t2,csden,cserr,temp,temperr,dens,denserr,neut,neuterr,rhovec,plotwin=plotwin,t=t,nz=nz,zeff=zeff,out=out,adas=adas
	IF !d.name EQ 'PS' THEN ps=1 ELSE ps=0
	IF NOT keyword_set(nz) THEN nz=1.0e17
	IF NOT keyword_set(zeff) THEN zeff=1.0
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
	IF NOT keyword_set(t) THEN t=0.5*(t1+t2)
	mdsopen,'spectroscopy',(shot)
	br=mdsvalue('\SPECTROSCOPY::TOP.BOLOMETER.RESULTS.FOIL:BRIGHT')
	rt=mdsvalue('dim_of(\SPECTROSCOPY::TOP.BOLOMETER.RESULTS.FOIL:BRIGHT,0)')
	bt=mdsvalue('dim_of(\SPECTROSCOPY::TOP.BOLOMETER.RESULTS.FOIL:BRIGHT,1)')
	good=mdsvalue('\SPECTROSCOPY::TOP.BOLOMETER.RESULTS.FOIL:GOOD_CHANS')
	good=good[0:n(rt)]
	mdsclose,'spectroscopy',(shot)

	i1=ipt(bt,t1)
	i2=ipt(bt,t2)
	xbr=sum_array(br[*,i1:i2],/i)/(i2-i1+1.0)
	xrt=rt
	br_err=sqrt((0.05*xbr)^2+(1.0e4)^2)
	redge=max(rt)+0.01
	shellinvert,xbr,xrt,em,r,brchk=brchk,br_err=br_err,good=good,npts=25,eps=0.1,em_err=em_err,redge=redge,fmweight=0.15
	xchkrad=brchk.rad
	xbr*=1.0e-6
	xbrerr=br_err*1.0e-6
	xbrchk=brchk.br*1.0e-6
	xrho=(r-ro)/a
	xem=em*1.0e-6
	xemerr=em_err*1.0e-6
	xgood=good
	xtmp=where(good EQ 1)

	line_emiss=foil_genrad_line_emiss(shot,t,rhovec*a+ro,csden*nz,cserr*nz,temp,temperr,dens,denserr,n_s=n_s,debug=debug,adas=adas)
	cxr_emiss=foil_genrad_cxr_emiss(shot,t,rhovec*a+ro,csden*nz,cserr*nz,temp,temperr,neut,neuterr,n_s=n_s,debug=debug)
	cont_emiss=foil_genrad_cont_emiss(shot,t,rhovec*a+ro,csden*nz,cserr*nz,temp,temperr,dens,denserr,n_s=n_s,debug=debug)
	brem_emiss=foil_genrad_brem_emiss(shot,t,rhovec*a+ro,temp,temperr,dens,denserr,n_s=n_s,debug=debug,zeff=zeff)

	emline=line_emiss.em*1.0e-6
	emerr_line=line_emiss.emerr*1.0e-6
	emcxr=cxr_emiss.em*1.0e-6
	emerr_cxr=cxr_emiss.emerr*1.0e-6
	emcont=cont_emiss.em*1.0e-6
	emerr_cont=cont_emiss.emerr*1.0e-6
	embrem=brem_emiss.em*1.0e-6
	emerr_brem=brem_emiss.emerr*1.0e-6
	emtot=emline+emcxr+emcont+embrem
	emerr_tot=sqrt(emerr_line^2+emerr_cxr^2+emerr_cont^2+emerr_brem^2)
	brtot=(line_emiss.br+cxr_emiss.br+cont_emiss.br+brem_emiss.br)*1.0e-6
	brerr_tot=sqrt(line_emiss.brerr^2+cxr_emiss.brerr^2+cont_emiss.brerr^2+brem_emiss.brerr^2)*1.0e-6
	rt=line_emiss.pos[2,*]

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
		device, xsize=xsize, ysize=ysize, /inches
	ENDELSE
	pos1=[0.075,0.1,0.6,0.95]
	pos3=[0.7,0.1,0.975,0.975]
	
	;emissivity comparison
	maxpt=max(emtot+emerr_tot) > max(xem+xemerr)
	tit='FOIL BOLOMETRY'
	lab=num2str(shot,1)
	IF keyword_set(t) THEN lab=lab+' '+num2str(t1,dp=2)+' < t < '+num2str(t2,dp=2)
	nzstr='n!lz!n = '+num2str(nz,dp=2)+' [m!u-3!n]'
	plot, [0],[0],pos=pos1,yr=[0,maxpt*1.05],ytit='Emissivity [MW/m!u3!n]',/ysty,xr=[min(xrho),max(rhovec) > max(xrho)],/xsty,chars=0.5*ls,$
		tit=tit,xtit='r/a'
	oploterror,rhovec,emtot,emerr_tot,color=200,errcolor=200
	oploterror,rhovec,emline,emerr_line,color=30,errcolor=30
	oploterror,rhovec,emcont,emerr_cont,color=30,linestyle=2.0,errcolor=30
	oploterror,rhovec,emcxr,emerr_cxr,color=30,linestyle=3.0,errcolor=30
	oploterror,rhovec,embrem,emerr_brem,color=30,linestyle=4.0,errcolor=30
	oploterror,xrho,xem,xemerr
	xyouts,0.0,0.85*maxpt,nzstr,chars=0.5*ls

	;brightness
	ch=indgen(n(br)+1)
	tmp=where(xgood EQ 1)
	maxpt=max(xbr[tmp]+xbrerr[tmp]) > max(brtot+brerr_tot)
	plot, [0],[0],pos=pos3,/noerase,yr=[0,maxpt*1.05],/ysty,xr=[0.6,0.91],/xsty,ytit='Brightness [MW/m!u2!n]',chars=0.5*ls,xtit='R!lT!n [m]'
	makesym,10
	oploterror,xrt[tmp],xbr[tmp],xbrerr[tmp],psym=8,symsize=0.5*ls
	oplot,xchkrad,xbrchk,color=100
	oploterror,rt,brtot,brerr_tot,color=200,errcolor=200
	;stop
	xyouts,0.923,0.05*maxpt,lab,chars=0.3*ls,orient=90

	out={rho:rhovec,em:emtot,emerr:emerr_tot,xem:xem,xrho:xrho,xemerr:xemerr}
END

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

FUNCTION axuv_calc_csplc,z,debug=debug,load=load,verb=verb,adas=adas

	t_e=10^(make(0.0,4.0,500))
	n_e=[1.0e20]
	CASE z OF
		18 : BEGIN
			;path_save='/home/mlreinke/idl/impurities/data/hullac_kbf/ar/ar_axuv_csplc.dat'
			;IF keyword_set(adas) THEN path_save='/home/mlreinke/idl/impurities/data/adas/ar_axuv_csplc.dat'
			path_load='/usr/local/cmod/idl/atomic_physics/hullac_kbf/ar/ar_axuv_csplc.dat'
			IF keyword_set(adas) THEN path_load='/usr/local/cmod/idl/atomic_physics/adas/ar_axuv_csplc.dat'
			path_save='/home/'+logname()+'/usr/'+logname()+'/ar_axuv_csplc.dat'
		END
		36 : BEGIN
			;path_save='/home/mlreinke/idl/impurities/data/hullac_kbf/kr/kr_axuv_csplc.dat'
			;IF keyword_set(adas) THEN path_save='/home/mlreinke/idl/impurities/data/adas/kr_axuv_csplc.dat'
			path_load='/usr/local/cmod/idl/atomic_physics/hullac_kbf/kr/kr_axuv_csplc.dat'
			IF keyword_set(adas) THEN path_load='/usr/local/cmod/idl/atomic_physics/adas/kr_axuv_csplc.dat'
			path_save='/home/'+logname()+'/usr/'+logname()+'/kr_axuv_csplc.dat'
		END
		ELSE : RETURN,-1
		
	ENDCASE
	IF NOT keyword_set(load) THEN BEGIN
		filter=axuv_trans()
		IF keyword_set(adas) THEN csplc=adas_pec_csplc(z,t_e,n_e,filter=filter,verb=verb) ELSE csplc=hullac_cs_plc(t_e,z,filter=filter,n_e=n_e[0],verb=verb)
		save,csplc,t_e,n_e,filename=path_save
	ENDIF ELSE restore,path_load

	output={csplc:csplc,temp:t_e,dens:n_e}
	IF keyword_set(debug) THEN stop
	RETURN,output
END

FUNCTION axuv_genrad_emiss,shot,t,rmaj,csden,cserr,temp,temperr,dens,denserr,n_s=n_s,debug=debug

	z=n(csden[*,0])
	nr=n(rmaj)+1
	csplc=axuv_calc_csplc(z,/load)
	emiss=fltarr(nr)
	emerr=fltarr(nr)
	FOR i=0,z DO BEGIN
		csint=interpol(csplc.csplc[i,*],csplc.temp,temp)
		dcsint=interpol(deriv(csplc.temp,csplc.csplc[i,*]),csplc.temp,temp)
		emiss+=csint*dens*csden[i,*]
		emerr+=(dens*csden[i,*]*dcsint)^2*temperr^2+(csint*dens)^2*(cserr[i,*])^2+(csint*csden[i,*])^2*denserr^2
	ENDFOR
	emerr=sqrt(emerr)
	axuv_genrad_em2br,shot,t,emiss,rmaj,br,pos,/gpv
	axuv_genrad_em2br,shot,t,emerr,rmaj,brerr,pos,/gpv
	
	output={rmaj:rmaj,em:emiss,emerr:emerr,pos:pos,br:br,brerr:brerr}
	IF keyword_set(debug) THEN stop
	RETURN,output
END

PRO axuv_genrad_em2br,shot,t,emiss,rmaj,br,pos,gpv=gpv

	IF keyword_set(gpv) THEN BEGIN
		pos=axuv_midplane_pos(etendue=u)
		pos=pos[*,0:21]
		u=u[0:21]
		restore,'/home/mlreinke/idl/genie/data/gpv/foil/genrad_axuv_gpv.dat'
		gpv=gpv[0:21,*]
		em_grid=grid_profile(ves_cent,emiss,rmaj,[t],shot,tpts=t,/sol)
		br=fltarr(22)
		FOR i=0,n(br) DO BEGIN
			non_zero=where(gpv[i,*] GT 0)
			IF non_zero[0] NE -1 THEN br[i]=total(gpv[i,non_zero]*em_grid[non_zero])/u[i]*4.0*!pi
		ENDFOR
	ENDIF ELSE BEGIN
		pos=axuv_midplane_pos()
		pos=pos[*,0:21]
		br=line_br(pos,emiss,rmaj,[t],shot,t,n_s=n_s,plots=debug)
	ENDELSE
END

FUNCTION axuv_genrad_brem_emiss,shot,t,rmaj,temp,temperr,dens,denserr,n_s=n_s,debug=debug,gpv=gpv,zeff=zeff

	filter=axuv_trans()
	IF NOT keyword_set(zeff) THEN zeff=1.0
	eph=filter.e
	lam=ev2ang(eph)
	npts=n(rmaj)+1
	emiss=fltarr(npts)
	emerr=fltarr(npts)
	h=6.626068e-34	;m^2 kg/s
	c=2.9979e8	;m/s
	e=1.60218e-19	;J/eV
	kb=8.617e-5 	;eV/K
	
	FOR i=0,npts-1 DO BEGIN
		log_gff=0.355*lam^(-0.06)*alog10(lam)+0.3*lam^(-0.066)*alog10(temp[i]/kb*1.0e-6)+0.0043
		specem=5.03e-14/h*e*4.0*!pi*10^(log_gff)*zeff*dens[i]/1.0e20*dens[i]/1.0e20/sqrt(temp[i])*exp(-eph/temp[i])	;watts/m^3/eV from Hutch
		emiss[i]=int_tabulated(eph,filter.tr*specem)	;W/m^3
		emerr[i]=sqrt((2.0*emiss[i]/dens[i])^2*denserr[i]^2+(0.5*emiss[i]/temp[i])^2*temperr[i]^2)
	ENDFOR

	axuv_genrad_em2br,shot,t,emiss,rmaj,br,pos,/gpv
	axuv_genrad_em2br,shot,t,emerr,rmaj,brerr,pos,/gpv
	
	output={rmaj:rmaj,em:emiss,emerr:emerr,pos:pos,br:br,brerr:brerr}
	IF keyword_set(debug) THEN stop
	RETURN,output
END

;6/12/12 - added the out optional output
PRO axuv_genrad_profiles,shot,t1,t2,csden,cserr,temp,temperr,dens,denserr,neut,neuterr,rhovec,plotwin=plotwin,t=t,nz=nz,zeff=zeff,out=out,foil=foil
	IF !d.name EQ 'PS' THEN ps=1 ELSE ps=0
	IF NOT keyword_set(nz) THEN nz=1.0e17

	;get efit_data
	mdsopen,'analysis',shot
	rmid=mdsvalue('\efit_rmid')
	efit_t=mdsvalue('dim_of(\efit_rmid)')
	mdsclose,'analysis',shot
	i1=ipt(efit_t,t1)
	i2=ipt(efit_t,t2)
	ro=min(sum_array(rmid[i1:i2,*],/j)/(i2-i1+1.0))
	a=max(sum_array(rmid[i1:i2,*],/j)/(i2-i1+1.0))-ro

	mdsopen,'spectroscopy',(shot)
	em=mdsvalue('\SPECTROSCOPY::TOP.BOLOMETER.RESULTS.DIODE.AXA:EMISS')
	r=mdsvalue('dim_of(\SPECTROSCOPY::TOP.BOLOMETER.RESULTS.DIODE.AXA:EMISS,0)')
	bt=mdsvalue('dim_of(\SPECTROSCOPY::TOP.BOLOMETER.RESULTS.DIODE.AXA:EMISS,1)')
	br=mdsvalue('\SPECTROSCOPY::TOP.BOLOMETER.RESULTS.DIODE.AXA:BRIGHT')
	rt=mdsvalue('dim_of(\SPECTROSCOPY::TOP.BOLOMETER.RESULTS.DIODE.AXA:BRIGHT,0)')
	brchk=mdsvalue('\SPECTROSCOPY::TOP.BOLOMETER.RESULTS.DIODE.AXA:BRCHK')
	good=mdsvalue('\SPECTROSCOPY::TOP.BOLOMETER.DIODE_CALIB.AXA:GOOD')
	s=mdsvalue('\SPECTROSCOPY::TOP.BOLOMETER.DIODE_CALIB.AXA:S')
	chkrad=reverse(rt[where(good EQ 1)])
	mdsclose,'spectroscopy',(shot)
	IF keyword_set(t) THEN BEGIN
		xem=em[*,ipt(bt,t)]
		xbr=br[*,ipt(bt,t)]
		xbrchk=brchk[*,ipt(bt,t)]
	ENDIF ELSE BEGIN
		i1=ipt(bt,t1)
		i2=ipt(bt,t2)
		xem=sum_array(em[*,i1:i2],/i)/(i2-i1+1.0)
		xbr=sum_array(br[*,i1:i2],/i)/(i2-i1+1.0)
		xbrchk=sum_array(brchk[*,i1:i2],/i)/(i2-i1+1.0)
	ENDELSE
	xrt=rt
	xchkrt=chkrad
	xrho=(r-ro)/a
	xgood=good
	xtmp=where(good EQ 1)
	IF NOT keyword_set(t) THEN t=0.5*(t1+t2)
	br_err=calc_diodebr_error(shot,'axa',0.01)*1.0e-6
	br_err=sqrt((0.02*xbr)^2+(br_err)^2)

	shellinvert,xbr,xrt,xem,r,brchk=brchk,br_err=br_err,good=good,npts=30,eps=0.01,em_err=xemerr,redge=0.915
	xchkrt=brchk.rad
	xbr_err=brchk.br_err
	xbrchk=brchk.br
	xrho=(r-ro)/a

	;set to maximum sensitivity
	xem*=s/0.28*1.0e-6
	xemerr*=s/0.28*1.0e-6
	xbr*=s/0.28*1.0e-6
	xbr_err*=s/0.28*1.0e-6
	xbrchk*=s/0.28*1.0e-6
	line_emiss=axuv_genrad_emiss(shot,t,rhovec*a+ro,csden*nz,cserr*nz,temp,temperr,dens,denserr,n_s=n_s,debug=debug)
	brem_emiss=axuv_genrad_brem_emiss(shot,t,rhovec*a+ro,temp,temperr,dens,denserr,n_s=n_s,debug=debug,zeff=zeff)

	emline=line_emiss.em*1.0e-6
	emerr_line=line_emiss.emerr*1.0e-6
	embrem=brem_emiss.em*1.0e-6
	emerr_brem=brem_emiss.emerr*1.0e-6
	emtot=emline+embrem
	emerr_tot=sqrt(emerr_line^2+emerr_brem^2)
	brtot=(line_emiss.br+brem_emiss.br)*1.0e-6
	brerr_tot=sqrt(line_emiss.brerr^2+brem_emiss.brerr^2)*1.0e-6

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
		device, xsize=xsize, ysize=ysize, /inches
	ENDELSE
	pos1=[0.075,0.1,0.6,0.95]
	pos3=[0.7,0.1,0.975,0.975]
	
	;emissivity comparison
	maxpt=max(emtot+emerr_tot) > max(xem+xemerr)
	tit='AXUV DIODE'
	lab=num2str(shot,1)
	IF keyword_set(t) THEN lab=lab+' t='+num2str(t,dp=2) ELSE lab=lab+' '+num2str(t1,dp=2)+' < t < '+num2str(t2,dp=2)
	nzstr='n!lz!n = '+num2str(nz,dp=2)+' [m!u-3!n]'
	plot, [0],[0],pos=pos1,yr=[0,maxpt*1.05],ytit='Emissivity [MW/m!u3!n]',/ysty,xr=[min(xrho),max(rhovec)],/xsty,chars=0.5*ls,$
		tit=tit,xtit='r/a'
	oploterror,rhovec,emtot,emerr_tot,color=200,errcolor=200
	oploterror,rhovec,emline,emerr_line,color=30,errcolor=30
	oploterror,rhovec,embrem,emerr_brem,color=30,linestyle=4.0,errcolor=30
	oploterror,xrho,xem,xemerr
	IF keyword_set(foil) THEN oplot, foil.xrho,foil.xem,linestyle=2.0
	xyouts,0.0,0.85*maxpt,nzstr,chars=0.5*ls

	;brightness
	ch=indgen(n(br)+1)
	tmp=where(xgood EQ 1)
	maxpt=max(xbr[tmp]+xbr_err[tmp]) > max(brtot+brerr_tot)
	plot, [0],[0],pos=pos3,/noerase,yr=[0,maxpt*1.05],/ysty,xr=[0.6,0.91],/xsty,ytit='Brightness [MW/m!u2!n]',chars=0.5*ls,xtit='R!lT!n [m]'
	makesym,10
	oploterror,xrt[tmp],xbr[tmp],xbr_err[tmp],psym=8,symsize=0.5*ls
	oplot,xchkrt,xbrchk,color=100
	oploterror,rt,brtot,brerr_tot,color=200,errcolor=200
	xyouts,0.923,0.05*maxpt,lab,chars=0.3*ls,orient=90

	out={rho:rhovec,em:emtot,emerr:emerr_tot,xem:xem,xrho:xrho,xemerr:xemerr}

	IF keyword_set(foil) THEN BEGIN
		plotwin+=1
		IF keyword_set(ps) THEN BEGIN
			xsize=7.0
			ysize=7.0*900/1100.0
			ls=2.0
		ENDIF ELSE BEGIN
			xsize=1100.0
			ysize=900.0
			ls=2.0
		ENDELSE
		IF NOT keyword_set(ps) THEN BEGIN
			device, window_state=var
			IF var[plotwin] EQ 0 THEN window,plotwin,xsize=xsize,ysize=ysize,xpos=1610,ypos=670,title='output profiles,'+num2str(plotwin) $
				ELSE wset,plotwin
		ENDIF ELSE BEGIN
			device, xsize=xsize, ysize=ysize, /inches
		ENDELSE

		maxpt=max(foil.em+foil.emerr) > max(foil.xem+foil.xemerr)
		lab=num2str(shot,1)
		IF keyword_set(t) THEN lab=lab+' t='+num2str(t,dp=2) ELSE lab=lab+' '+num2str(t1,dp=2)+' < t < '+num2str(t2,dp=2)
	
		plot, [0],[0],yr=[0,maxpt*1.05],ytit='Emissivity [MW/m!u3!n]',/ysty,xr=[min(xrho),max(rhovec)],/xsty,chars=0.5*ls,xtit='r/a'
		oploterror,rhovec,emtot,emerr_tot,color=200,errcolor=200,linestyle=2.0
		oploterror,xrho,xem,xemerr,linestyle=2.0
		oploterror,foil.xrho,foil.xem,foil.xemerr
		oploterror,foil.rho,foil.em,foil.emerr,color=200,errcolor=200
		oplot, [0.0,0.15],maxpt*0.8*[1.0,1.0],linestyle=2.0
		oplot, [0.0,0.15],maxpt*0.7*[1.0,1.0]
		xyouts,0.2,0.68*maxpt,'FOIL',chars=0.8*ls
		xyouts,0.2,0.78*maxpt,'AXUV',chars=0.8*ls
	ENDIF
	
END
