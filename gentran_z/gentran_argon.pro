PRO calib_fit_func,x,a,f
	f=a[0]+a[1]/x
END

PRO xeus_argon_calib,plot=plot,verb=verb,xtomo=xtomo,foil=foil,fitplot=fitplot,dff=dff,pcol=pcol
	IF NOT keyword_set(dff) THEN dff=0.5
	shot=1110125006
	t1=[0.75,1.00,1.25,1.43]
	t2=[0.80,1.05,1.30,1.48]
	CASE dff OF
		1.0 : nz=[6.79,5.28,3.49,2.37]*1.0e17			;densities derived from XTOMO brightness
		0.5 : nz=[7.15,5.35,3.53,2.57]*1.0e17
		0.1 : nz=[7.87,5.81,3.76,2.76]*1.0e17
		ELSE : RETURN
	ENDCASE
	;nz=[7.3,5.3,3.7,2.63]*1.0e17					;densities derived from core FOIL emissivity
	lam0=10
	lam1=40
	nlam=3000
	xeus_pos=[2.561,0.210,0.179,6.50*!pi/180.0]
	genpos_pos_reform,xeus_pos,[0.44,1.0,0.6,-0.6]
	dlam=0.04
	lam=make(lam0,lam1,nlam)
	ntime=n(t1)+1
	charge=findgen(19)
	zeff_bck=1.15
	FOR j=0,ntime-1 DO BEGIN
		efit=0
		rmid=0
		data=0
		gentran_profiles,shot,18,t1[j],t2[j],/pt,dff=dff,data=data,plot=plot
		print, 'csden profile calculated '+num2str(shot,1)+' '+num2str(t1[j],dp=2)+' < t < '+num2str(t2[j],dp=2)
		pec=adas_pec_spec(18,data.temp,data.dens,[lam0,lam1])
	
		nlines=n(pec.lam)+1
		nrad=n(pec.temp)+1
		IF j EQ 0 THEN BEGIN
			emiss=fltarr(nrad,nlines,ntime)
			bright=fltarr(nlines,ntime)
			spec=fltarr(nlam,ntime)
			zeff=fltarr(nrad,ntime)
			br_lilike=fltarr(ntime)
		ENDIF
		FOR i=0,nlines-1 DO emiss[*,i,j]=data.csden[pec.q[i],*]*data.dens*pec.sigv[i,*]
		print, 'emissivities calculated'
	
		FOR i=0,nlines-1 DO bright[i,j]=line_br(xeus_pos,emiss[*,i,j],data.rmaj,[data.time],data.shot,data.time,plots=plot,efit=efit,rmid=rmid,verb=verb)
		print, 'line-integrated brightness calculated'

		FOR i=0,nlines-1 DO spec[*,j]+=bright[i,j]/(dlam*sqrt(2.0*!pi))*exp(-(lam-pec.lam[i])^2/(2*dlam^2))
		print, 'line-integrated spectra formed'
	
		;use the xtomo to derive an impurity density	
		IF keyword_set(xtomo) THEN xtomo_genrad_profiles,shot,t1[j],t2[j],data.csden,data.cserr,data.temp,data.terr,data.dens,data.derr,data.rho,$
			plotwin=plotwin,t=t,nz=nz[j],zeff=zeff_bck,out=outx

		IF keyword_set(foil) THEN foil_genrad_profiles,shot,t1[j],t2[j],data.csden,data.cserr,data.temp,data.terr,data.dens,data.derr,$
			data.neut,data.nerr,data.rho,plotwin=plotwin,t=t,nz=nz[j],zeff=zeff_bck,out=outf
		IF keyword_set(foil) THEN stop

		FOR i=0,nrad-1 DO zeff[i,j]=total(data.csden[*,j]*charge^2*nz[j]/data.dens[i])
		br_lilike[j]=total(bright[where(pec.lam EQ 23.5),j])					;ph/s/m^2/nz
	ENDFOR
	tau=0.5*(t1+t2)
	zeff_neo,shot,zeff_neo,zeff_time,dt=0.05
	openwin,0
	plot,zeff_time,zeff_neo,xr=[0,2.0],yr=[0,3.2],xtit='Time [sec]',ytit='Z!lEFF!n'
	oplot,tau,zeff[0,*]+zeff_bck,psym=-8,color=100
	xyouts,0.2,1.5,'ZEFF_NEO'
	xyouts,0.2,1.0,'ZEFF_GENTRAN = 1.15+'+n2g('Delta')+'ZEFF_AR',color=100


	mdsopen,'spectroscopy',shot
	cnts=mdsvalue('\SPECTROSCOPY::TOP.XEUS.RAW:DATA')
	mdsclose,'spectroscopy',shot
	brxeus=sum_array(cnts,/j)
	bl=mean(brxeus[950:*])	
	brxeus-=bl
	brxeus[0]=0

	lampt=[18.97,21.6,48.59,40.98,38.87]
	xpts=[1026.0,948.9,321.2,475.9,522.3]

	dt=1.954/(899-25)
	cnts/=dt					;cnts/frame -> cnts/sec
	xeus_time=(indgen(1000)-25)*dt
	coefs=poly_fit(xpts,lampt,2)
	xplot=findgen(1340)
	xeus_lam=xplot^2*coefs[2]+xplot*coefs[1]+coefs[0]
	IF keyword_set(fitplot) THEN BEGIN
		openwin,0
		plot,[0],[0],xr=[0,1340],yr=[10,75],xtit='Pixel #', ytit=n2g('lambda')+'[A\ng]',chars=1.2,/xsty,/ysty
		oplot, xpts,lampt,psym=8
		oplot,xplot,xeus_lam,color=200
	ENDIF
	
	nxtime=n(xeus_time)+1
	br_lilike_xeus=fltarr(nxtime)
	tmp=where(xeus_lam GE 23.2 AND xeus_lam LE 23.8)
	xfit=xeus_lam[tmp]
	start=ipt(xeus_time,0.5)
	stop=ipt(xeus_time,1.5)
	FOR i=start,stop DO BEGIN
		yfit=cnts[tmp,i]
		a=[max(yfit)-500,23.5,0.1]
		ain=a
		fit=gaussfit(xfit,yfit,a,nterms=5)
		IF total(a-ain) NE 0 THEN br_lilike_xeus[i]=a[0]*(a[2]*sqrt(2.0*!pi))					;cnts/s integrated over line
	ENDFOR
	br_meas=fltarr(ntime)
	FOR i=0,ntime-1 DO br_meas[i]=mean(br_lilike_xeus[where(xeus_time GE t1[i] AND xeus_time LE t2[i])])
	print, 'measured xeus intensity found'
	!p.multi=[0,0,2]
	openwin,1
	plot,xeus_time,br_lilike_xeus,xr=[0,2.0],yr=[0,max(br_lilike_xeus)*1.05],xtit='Time [sec]',ytit='[AU/s]',tit='23.5 Ang XEUS Intensity'
	plot,tau,br_lilike,xr=[0,2.0],yr=[0,max(br_lilike)*1.05],xtit='Time [sec]',ytit='[ph/s/m!u2!n/n!lz!n]',psym=-8,tit='23.5 Ang GENRAD Intensity'
	!p.multi=0
	const=br_lilike*nz/br_meas

	openwin,2
	tit=num2str(shot,1)+' dff='+num2str(dff,dp=1)
	plot,nz*1.0e-17,const*1.0e-14,psym=8,xtit='Ar Density x10!u17!n from XTOMO',ytit='XEUS Const x10!u14!n [ph/m!u2!n/AU]',xr=[0,20],yr=[0.0,3.5]
	xfit=nz*1.0e-17
	yfit=const*1.0e-14
	a=[2.5,2.0]
	out=curvefit(xfit,yfit,weights,a,function_name='calib_fit_func',/noderiv)
	xpl=make(0.1,25,100)
	oplot,xpl,a[0]+a[1]/xpl,color=200
	xyouts,2,1,'XEUS Calib @ 23.5 Anstroms: '+num2str(a[0],dp=1)+'x10!u14!n [ph/m!u2!n/AU]',chars=1.2,color=200
	
END

PRO hirexsr_load_thaco_bright,shot,group,br,pos,tau

	CASE group OF
		1 : BEGIN	;w+n>=3
			pospath='\SPECTROSCOPY::TOP.HIREXSR.ANALYSIS.HELIKE.MOMENTS.W:POS'
			cpath='\SPECTROSCOPY::TOP.HIREXSR.ANALYSIS.HELIKE.FITS.WN3:COEFS'
			labels='\SPECTROSCOPY::TOP.HIREXSR.ANALYSIS.HELIKE.FITS.WN3:LABELS'
			lo=['w','wn5','wn4','wn3']
		END
		2 : BEGIN ;xy
			pospath='\SPECTROSCOPY::TOP.HIREXSR.ANALYSIS.HELIKE.MOMENTS.X:POS'
			cpath='\SPECTROSCOPY::TOP.HIREXSR.ANALYSIS.HELIKE.FITS.XY:COEFS'
			labels='\SPECTROSCOPY::TOP.HIREXSR.ANALYSIS.HELIKE.FITS.XY:LABELS'
			lo=['n','x','y','st','yn4','yn3']
		END
		3 : BEGIN ;qra
			pospath='\SPECTROSCOPY::TOP.HIREXSR.ANALYSIS.HELIKE.MOMENTS.Z:POS'
			cpath='\SPECTROSCOPY::TOP.HIREXSR.ANALYSIS.HELIKE.FITS.ZJK:COEFS'
			labels='\SPECTROSCOPY::TOP.HIREXSR.ANALYSIS.HELIKE.FITS.ZJK:LABELS'
			lo=['q','r','a']
		END

		4 : BEGIN ;zjk
			pospath='\SPECTROSCOPY::TOP.HIREXSR.ANALYSIS.HELIKE.MOMENTS.Z:POS'
			cpath='\SPECTROSCOPY::TOP.HIREXSR.ANALYSIS.HELIKE.FITS.ZJK:COEFS'
			labels='\SPECTROSCOPY::TOP.HIREXSR.ANALYSIS.HELIKE.FITS.ZJK:LABELS'
			lo=['z','j','k']
		END
		5 : BEGIN	;lya1,2
			pospath='\SPECTROSCOPY::TOP.HIREXSR.ANALYSIS.HLIKE.MOMENTS.LYA1:POS'
			cpath='\SPECTROSCOPY::TOP.HIREXSR.ANALYSIS.HLIKE.FITS.LYA:COEFS'
			labels='\SPECTROSCOPY::TOP.HIREXSR.ANALYSIS.HLIKE.FITS.LYA:LABELS'
			lo=['lyas1','lya1','lyas2','lya2','lyas3','lyas4']
		END
		6 : BEGIN	;lyaTJ
			pospath='\SPECTROSCOPY::TOP.HIREXSR.ANALYSIS.HLIKE.MOMENTS.J:POS'
			cpath='\SPECTROSCOPY::TOP.HIREXSR.ANALYSIS.HLIKE.FITS.TJ:COEFS'
			labels='\SPECTROSCOPY::TOP.HIREXSR.ANALYSIS.HLIKE.FITS.TJ:LABELS'
			lo=['T','K','Q','B','R','A','J']
		END
		ELSE : RETURN
	ENDCASE

	mdsopen,'spectroscopy',shot
	pos=mdsvalue(pospath)
	coefs=mdsvalue(cpath)
	time=mdsvalue('dim_of('+cpath+',1)')
	labels=mdsvalue(labels)
	mdsclose,'spectroscopy',shot
	
	;select time range of interest
	tmp=where(time NE -1)
	tau=time[tmp]
	coefs=coefs[*,tmp,*]
	
	;select spatial coords of interest (assume only one CHMAP)
	tmp=where(pos[0,*,0] NE -1)
	pos=pos[*,tmp]
	coefs=coefs[tmp,*,*]

	ntime=n(tau)+1
	nch=n(pos[0,*])+1
	
	br=fltarr(nch,ntime)
	FOR i=0,n(lo) DO BEGIN
		index=where(labels EQ lo[i])
		FOR j=0,nch-1 DO BEGIN
			FOR k=0,ntime-1 DO br[j,k]+=coefs[j,k,index*3]*coefs[j,k,index*3+2]*sqrt(2.0*!pi)
		ENDFOR
	ENDFOR

END


PRO hirexsr_argon_calib,group=group,prof=prof,dff=dff
	IF NOT keyword_set(group) THEN group=4
	CASE group OF 
		1 : BEGIN ;wn3
			wl_roi=[3.9455,3.9600]
			lab='wn3'
			ymax=10.0
		END
		2 : BEGIN ;xy
			wl_roi=[3.9600,3.9750]
			lab='xy'
			ymax=3.5
		END
		3 : BEGIN ;qra
			wl_roi=[3.9780,3.9875]
			lab='qra'
			ymax=2.0
		END
		4 : BEGIN ;zjk
			wl_roi=[3.9875,4.0000]
			lab='zjk'
			ymax=2.0
		END
		5 : BEGIN ;lya
			wl_roi=[3.725,3.746]
			lab='Lya1,2'
			ymax=3.5
		END
		6 : BEGIN ;tj
			wl_roi=[3.750,3.780]
			lab='TJ'
			ymax=3.5
		END
		ELSE : RETURN
	ENDCASE
	IF NOT keyword_set(dff) THEN dff=0.5
	shot=1110125006
	t1=[0.75,1.00,1.25,1.43]
	t2=[0.80,1.05,1.30,1.48]
	CASE dff OF
		1.0 : nz=[6.79,5.28,3.49,2.37]*1.0e17			;densities derived from XTOMO brightness
		0.5 : nz=[7.15,5.35,3.53,2.57]*1.0e17
		0.1 : nz=[7.87,5.81,3.76,2.76]*1.0e17
		ELSE : RETURN
	ENDCASE
	;nz=[7.3,5.3,3.7,2.63]*1.0e17					;densities derived from core FOIL emissivity

	ntime=n(t1)+1
	charge=findgen(19)
	zeff_bck=1.15
	hirexsr_load_thaco_bright,shot,group,br,pos,thirex
	print, 'measured HIREXSR intensity found'
	IF group EQ 5 OR group EQ 6 THEN normch=3 ELSE normch=24
	FOR j=0,ntime-1 DO BEGIN
		efit=0
		rmid=0
		data=0
		gentran_profiles,shot,18,t1[j],t2[j],/pt,dff=dff,data=data,plot=plot
		print, 'csden profile calculated '+num2str(shot,1)+' '+num2str(t1[j],dp=2)+' < t < '+num2str(t2[j],dp=2)
		iemiss=ar_xray_emiss(data.csden,data.dens,data.temp,wl_roi=wl_roi,cserr=data.cserr,sigte=data.terr,signe=data.derr,emerr=emerr,csemiss=csem)
		nrad=n(iemiss)+1
		nch=n(pos[0,*])+1
		IF j EQ 0 THEN BEGIN
			emiss=fltarr(nrad,ntime)
			bright=fltarr(nch,ntime)
			brcore=fltarr(ntime)
			zeff=fltarr(nrad,ntime)
		ENDIF
		emiss[*,j]=iemiss
		print, 'emissivities calculated'
	
		IF keyword_set(prof) THEN bright[*,j]=line_br(pos,emiss[*,j],data.rmaj,[data.time],data.shot,data.time,plots=plot,verb=verb)
		brcore[j]=line_br(pos[*,normch],emiss[*,j],data.rmaj,[data.time],data.shot,data.time,plots=plot,verb=verb)
		print, 'line-integrated brightness calculated'
	
		;use the xtomo to derive an impurity density	
		IF keyword_set(xtomo) THEN xtomo_genrad_profiles,shot,t1[j],t2[j],data.csden,data.cserr,data.temp,data.terr,data.dens,data.derr,data.rho,$
			plotwin=plotwin,t=t,nz=nz[j],zeff=zeff_bck,out=outx

		;use foil bolometry to derive an impurity density
		IF keyword_set(foil) THEN foil_genrad_profiles,shot,t1[j],t2[j],data.csden,data.cserr,data.temp,data.terr,data.dens,data.derr,$
			data.neut,data.nerr,data.rho,plotwin=plotwin,t=t,nz=nz[j],zeff=zeff_bck,out=outf

		FOR i=0,nrad-1 DO zeff[i,j]=total(data.csden[*,j]*charge^2*nz[j]/data.dens[i])
	ENDFOR
	tau=0.5*(t1+t2)
	zeff_neo,shot,zeff_neo,zeff_time,dt=0.05
	openwin,0
	plot,zeff_time,zeff_neo,xr=[0,2.0],yr=[0,3.2],xtit='Time [sec]',ytit='Z!lEFF!n'
	oplot,tau,zeff[0,*]+zeff_bck,psym=-8,color=100
	xyouts,0.2,1.5,'ZEFF_NEO'
	xyouts,0.2,1.0,'ZEFF_GENTRAN = 1.15+'+n2g('Delta')+'ZEFF_AR',color=100
	

	br_meas=fltarr(ntime)
	br_prof_meas=fltarr(nch,ntime)
	FOR i=0,ntime-1 DO BEGIN
		tmp=where(thirex GE t1[i] AND thirex LE t2[i])
		IF tmp[0] EQ -1 THEN tmp=[ipt(thirex,0.5*(t1[i]+t2[i]))]
		br_meas[i]=mean(br[normch,tmp])
		FOR j=0,nch-1 DO br_prof_meas[j,i]=mean(br[j,tmp])
	ENDFOR
	!p.multi=[0,0,2]
	openwin,1
	plot,thirex,br[normch,*],xr=[0,2.0],yr=[0,max(br[normch,*])*1.05],xtit='Time [sec]',ytit='10!u17!n [ph/s/m!u2!n]',tit='HIREXSR '+lab+' Intensity'
	plot,tau,brcore,xr=[0,2.0],yr=[0,max(brcore)*1.05],xtit='Time [sec]',ytit='[ph/s/m!u2!n/n!lz!n]',psym=-8,tit='GENRAD '+lab+' Intensity'
	!p.multi=0
	const=brcore*nz/(br_meas*1.0e17)

	openwin,2
	tit=num2str(shot,1)+' dff='+num2str(dff,dp=1)
	plot,nz*1.0e-17,const*1.0e-2,psym=8,xtit='Ar Density x10!u17!n from XTOMO',ytit='HIREXSR '+lab+' Const x10!u2!n',xr=[0,20],yr=[0.0,ymax]
	xfit=nz*1.0e-17
	yfit=const*1.0e-2
	a=[2.5,2.0]
	out=curvefit(xfit,yfit,weights,a,function_name='calib_fit_func',/noderiv)
	xpl=make(0.1,25,100)
	oplot,xpl,a[0]+a[1]/xpl,color=200
	xyouts,2,0.25*ymax,'HIREXSR '+lab+' Calib: '+num2str(a[0],dp=2)+'x10!u2!n (REAL/MEASURED)',chars=1.2,color=200
	
	IF keyword_set(prof) THEN BEGIN
		color=[0,30,100,200]
		!p.multi=[0,2,2]
		openwin,3
		ch=indgen(nch)+1
		FOR i=0,n(color) DO BEGIN

			IF group EQ 5 OR group EQ 6 THEN ymax=1.05 ELSE ymax=max(br_prof_meas[10:40,i]/br_prof_meas[normch,i])
			plot,[0],[0],xr=[0,nch+1],yr=[0,ymax],/xsty,/ysty,xtit='CH #', ytit='Normalized Intensity',$
				tit=lab+' '+num2str(shot,1)+' t='+num2str(0.5*(t1[i]+t2[i]),dp=3)
			makesym,10
			oplot,ch,br_prof_meas[*,i]/br_prof_meas[normch,i],psym=-8,color=color[i]	
			makesym,9
			oplot,ch,bright[*,i]/bright[normch,i],psym=-8,color=color[i],linestyle=2.0
		ENDFOR
		!p.multi=0
		makesym,10
	ENDIF

END

PRO hirexsr_argon_density,shot,t1,t2,dff=dff,plot=plot
	IF n(t1) NE n(t2) THEN RETURN
	ntime=n(t1)+1
	IF NOT keyword_set(dff) THEN dff=0.5
	IF NOT keyword_set(group) THEN group=4
	CASE group OF 
		1 : BEGIN ;wn3
			wl_roi=[3.9455,3.9600]
			lab='wn3'
			eta=4.98e2
		END
		2 : BEGIN ;xy
			wl_roi=[3.9600,3.9750]
			lab='xy'
			eta=1.72e2
		END
		3 : BEGIN ;qra
			wl_roi=[3.9780,3.9875]
			lab='qra'
			eta=1.21e2
		END
		4 : BEGIN ;zjk
			wl_roi=[3.9875,4.0000]
			lab='zjk'
			eta=1.1e2
		END
		5 : BEGIN ;lya
			wl_roi=[3.725,3.746]
			lab='Lya1,2'
			eta=2.06e2
		END
		6 : BEGIN ;tj
			wl_roi=[3.750,3.780]
			lab='TJ'
			eta=1.92e2
		END
		ELSE : RETURN
	ENDCASE
	eta*=1.0e17
	hirexsr_load_thaco_bright,shot,group,br,pos,tau
	nch=n(pos[0,*])+1
	normch=int(nch/2)

	FOR j=0,ntime-1 DO BEGIN
		efit=0
		rmid=0
		data=0

		gentran_profiles,shot,18,t1[j],t2[j],/pt,dff=dff,data=data,plot=plot
		print, 'csden profile calculated '+num2str(shot,1)+' '+num2str(t1[j],dp=2)+' < t < '+num2str(t2[j],dp=2)
		iemiss=ar_xray_emiss(data.csden,data.dens,data.temp,wl_roi=wl_roi,cserr=data.cserr,sigte=data.terr,signe=data.derr,emerr=emerr,csemiss=csem)
		nrad=n(iemiss)+1
		IF j EQ 0 THEN BEGIN
			emiss=fltarr(nrad,ntime)
			bright=fltarr(nch,ntime)
			brcore=fltarr(ntime)
			zeff=fltarr(nrad,ntime)
		ENDIF
		emiss[*,j]=iemiss
		print, 'emissivities calculated'

		IF keyword_set(prof) THEN bright[*,j]=line_br(pos,emiss[*,j],data.rmaj,[data.time],data.shot,data.time,plots=plot,verb=verb)
		brcore[j]=line_br(pos[*,normch],emiss[*,j],data.rmaj,[data.time],data.shot,data.time,plots=plot,verb=verb)
		print, 'line-integrated brightness calculated'
	
	ENDFOR

	nz=fltarr(ntime)
	FOR j=0,ntime-1 DO BEGIN
		tmp=where(tau GE t1[j] AND tau LE t2[j]) 
		IF tmp[0] NE -1 THEN ave_br=mean(br[normch,tmp]) ELSE ave_br=br[normch,ipt(0.5*(t1[j]+t2[j]),tau)]
		nz[j]=eta*ave_br/brcore[j]
	ENDFOR



	stop
END

