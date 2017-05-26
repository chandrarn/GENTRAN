FUNCTION read_dielsattauto_data,save=save,load=data,path=path

	IF NOT keyword_set(path) THEN path='/usr/local/cmod/idl/atomic_physics/'
	data_path=path+'DielSattAuto.dat'
	save_path='DielSattAuto.sav'

	IF keyword_set(load) THEN BEGIN
		restore, save_path
		RETURN,data
	ENDIF

	npts=1000
	nsat=intarr(npts)
	wave=fltarr(npts)
	Es=lonarr(npts)
	Ar=fltarr(npts)
	Aa=fltarr(npts)
	F2=fltarr(npts)
	gs=intarr(npts)
	g0=intarr(npts)
	Config=strarr(npts)
	
	openr,lun,data_path,/get_lun
	line=''	
	FOR i=0,79 DO readf,lun,line
	readf,lun,line
	numlev=int(line)
	cntr=0
	isat=2
	WHILE isat LT 7 DO BEGIN
		readf,lun,line
		IF total(strsplit(line,'$',/extract) EQ line) EQ 0 THEN BEGIN
			isat=int(last(strsplit(line,'=',/extract)))
			readf,lun,line
			header=strsplit(line,'!',/extract)
			nlines=int(header[0])
			FOR i=0,nlines-1 DO BEGIN
				readf,lun,line
				data=strsplit(line, '	',/extract)
				nsat[cntr]=isat
				wave[cntr]=float(data[0])
				Es[cntr]=long(data[1])
				Ar[cntr]=float(data[2])*1.0e13
				Aa[cntr]=float(data[3])*1.0e13
				F2[cntr]=float(data[4])*1.0e13
				gs[cntr]=int(data[5])
				c0=last(strsplit(data[6],'-',/extract))
				index=strpos(c0,'/')
				frac=strsplit(strmid(c0,index-1,3),'/',/extract)
				g0[cntr]=float(frac[0])/float(frac[1])*2.0+1
				config[cntr]=data[6]
				cntr+=1
			ENDFOR
                ENDIF
		isat+=1
	ENDWHILE
	tmp=where(wave NE 0)
	nsat=nsat[tmp]
	wave=wave[tmp]
	Es=Es[tmp]
	Ar=Ar[tmp]
	Aa=Aa[tmp]
	F2=F2[tmp]
	gs=gs[tmp]
	g0=g0[tmp]
	config=config[tmp]
	data={nsat:nsat,wave:wave,Es:Es,Ar:Ar,Aa:Aa,F2:F2,gs:gs,g0:g0,config:config}

	close,lun
	free_lun, lun
		
	IF keyword_set(save) THEN save,data,filename='/home/'+logname()+'/atomic_physics/'+save_path
	RETURN,data
END

PRO marchuk_w_rates,etemp,exc,rec
	;these are manually entered from Appendix C and D from Marchuk's thesis
	etemp=make(0.2,4.0,39)*1.0e3
	exc=[4.51e-6,7.24e-4,8.99e-3,4.03e-2,0.109,0.222,0.377,0.569,0.791,1.03,1.29,1.56,1.84,2.21,2.40,2.67,2.95,3.21,3.48,$
		3.73,3.98,4.22,4.46,4.69,4.91,5.12,5.33,5.53,5.73,5.91,6.1,6.27,6.44,6.61,6.77,6.93,7.08,7.22,7.36]*1.0e-18 	;m^3/s
	etrec=make(0.2,4.0,20)*1.0e3
	rec=[8.78,4.95,4.04,3.82,3.75,3.83,3.63,3.52,3.33,3.15,2.98,2.88,2.72,2.58,2.44,2.35,2.22,2.11,2.01,1.93]*1.0e-19	;m^3/s
	rec=10^(interpol(alog10(rec),etrec,etemp))	
END	

PRO marchuk_sat_rates,etemp,rec,data
	data=READ_DIELSATTAUTO_DATA()
	IF NOT keyword_set(etemp) THEN etemp=double(10^make(0.2,4,1000))
	a0=5.29177e-11	;bohr radius in meters
	ry=13.605	;rydberg in eV
	ntemp=n(etemp)+1
	ntrans=n(data.nsat)+1
	rec=dblarr(ntemp,ntrans)

	FOR i=0,ntrans-1 DO rec[*,i]=4.0*!pi^1.5*a0^3/(etemp/ry)^1.5*data.f2[i]*exp(-1.0*cm2ev(data.es[i])/etemp)
END

PRO plot_dielsatt_rates
	data=READ_DIELSATTAUTO_DATA()
	IF NOT keyword_set(etemp) THEN etemp=double(10^make(0.2,4,1000))
	a0=5.29177e-11	;bohr radius in meters
	ry=13.605	;rydberg in eV
	c=2.998e8 			;speed of light
	e=1.602e-19			;conversion for eV -> J
	mconv=1.661e-27			;conversion for amu -> kg
	mass=39.95			;atomic mass of argon

	ntemp=n(etemp)+1
	ntrans=n(data.nsat)+1
	rates=dblarr(ntemp,ntrans)

	FOR i=0,ntrans-1 DO rates[*,i]=4.0*!pi^1.5*a0^3/(etemp/ry)^1.5*data.f2[i]*exp(-1.0*cm2ev(data.es[i])/etemp)

	etk=etemp/8.617e-5
	rates*=1.0e6
	yr=[0,1.2e-12]
	openwin,0
	plot,[0],[0],xr=[4.0e6,5.0e7],/xsty,/xlog,yr=yr,/ysty,xtit='T!ie!n [K]',ytit=n2g('alpha')+'!iDR!n [cm!u3!n/s]',chars=1.2
	tmp=where(data.nsat EQ 2)
	oplot,etk,sum_array(rates[*,tmp],/i),color=200
	tmp=where(data.nsat EQ 3)
	oplot,etk,sum_array(rates[*,tmp],/i),color=30
	tmp=where(data.nsat EQ 4)
	oplot,etk,sum_array(rates[*,tmp],/i),color=100
	oplot,etk,sum_array(rates,/i),color=120
	
	openwin,1
	marchuk_w_rates,wtemp,exc,rec
	rates*=1.0e-6
	plot,[0],[0],xr=[0.5,3.0],yr=[0.001,1.0],/xsty,/ysty,/ylog,xtit='T!ie!n [keV]',ytit='Ratio',chars=1.2
	oplot,etemp*1.0e-3,rates[*,9]/interpol(exc,wtemp,etemp),color=100
	;oplot,etemp*1.0e-3,rates[*,11]/interpol(exc,wtemp,etemp),color=200
	ark=reform_ark_data([3.94,3.96],/load)
	read_ar_line_rates,'k',kte,krec,kexc,kion,/load
	read_ar_line_rates,'j',jte,jrec,jexc,jion,/load
	read_ar_line_rates,'w',wte,wrec,wexc,wion,/load
	oplot,kte,krec/wexc,color=100,linestyle=2.0
	;oplot,jte,jrec/wexc,color=200,linestyle=2.0

	tmp=where(data.nsat EQ 3 AND data.wave LT 3.96)
	oplot,etemp*1.0e-3,sum_array(rates[*,tmp],/i)/interpol(exc,wtemp,etemp),color=120
	atmp=where(ark.type EQ '30201') 
	oplot,ark.temp,sum_array(ark.rec[*,atmp],/i)/ark.exc[*,0],linestyle=2.0,color=120
	tmp=where(data.nsat EQ 4 AND data.wave LT 3.96)
	oplot,etemp*1.0e-3,sum_array(rates[*,tmp],/i)/interpol(exc,wtemp,etemp),color=30
	atmp=where(ark.type EQ '40201') 
	oplot,ark.temp,sum_array(ark.rec[*,atmp],/i)/ark.exc[*,0],linestyle=2.0,color=30
	tmp=where(data.nsat GE 5)
	;oplot,etemp*1.0e-3,sum_array(rates[*,tmp],/i)/interpol(exc,wtemp,etemp),color=0
	atmp=where(ark.type NE '30201' AND ark.type NE '40201' AND ark.type NE '201') 
	;oplot,ark.temp,sum_array(ark.rec[*,atmp],/i)/ark.exc[*,0],linestyle=2.0,color=0

	openwin,2
	plot,[0],[0],xr=[0.3,4.0],yr=[0.01,1.0],/ylog,/xsty,/ysty,xtit='Elec. Temperature [keV]',$
		ytit='Emissivity Ratio '+n2g('epsilon')+'!in=3!n/('+n2g('epsilon')+'!iw!n+'+n2g('epsilon')+'!in>=3!n)',chars=1.3
	tmpa=where(data.nsat EQ 3 AND data.wave GT 3.953 AND data.wave LT 3.96)
	tmpb=where(data.nsat GE 3 AND data.wave LT 3.96)
	oplot,etemp*1.0e-3,sum_array(rates[*,tmpa],/i)/(sum_array(rates[*,tmpb],/i)+interpol(exc,wtemp,etemp)),linestyle=3.0
	ark=reform_ark_data([3.94,3.96],/load)
	atmp=where(ark.type EQ '30201' AND ark.lam GT 3.953) 
	btmp=where(ark.type NE '201')
	oplot,ark.temp,sum_array(ark.rec[*,atmp],/i)/(sum_array(ark.rec[*,btmp],/i)+ark.exc[*,0]),linestyle=2.0
	ratio=10^make(alog10(0.02),alog10(0.25),1000)
	oplot,0.1645*ratio^(-0.7642),ratio,color=200
	oplot,[2.2,2.8],[0.8,0.8],linestyle=2.0
	oplot,[2.2,2.8],[0.4,0.4],linestyle=3.0
	oplot,[2.2,2.8],[0.2,0.2],linestyle=0.0,color=200
	xyouts,2.9,0.76,'FAC'
	xyouts,2.9,0.38,'AUTOSTRUCTURE'
	xyouts,2.9,0.188,'EXPERIMENT',color=200
	openwin,3
	itemp=ark.temp[ipt(ark.temp,1.0)]*1.0e3
	nlam=1000
	lam=make(3.945,3.96,nlam)
	spec=fltarr(nlam,3,2)		;w,n=3,n=4,n >=5

	index=ipt(etemp,itemp)
	FOR i=0,n(data.wave) DO BEGIN
		conv_factor=(data.wave[i]/c)*sqrt(e/(mass*mconv))		;conversion factor for Ti to Gaussian width in angstroms
		sigma=sqrt(itemp)*conv_factor
		ispec=rates[index,i]/(sqrt(2.0*!pi)*sigma)*exp(-(lam-data.wave[i])^2/(2.0*sigma^2))
		IF data.nsat[i] EQ 3 THEN spec[*,0,0]+=ispec
		IF data.nsat[i] EQ 4 THEN spec[*,1,0]+=ispec
		IF data.nsat[i] GE 5 THEN spec[*,2,0]+=ispec
     	ENDFOR
	
	index=ipt(ark.temp,itemp/1.0e3)
	FOR i=0,n(ark.lam) DO BEGIN
		conv_factor=(ark.lam[i]/c)*sqrt(e/(mass*mconv))		;conversion factor for Ti to Gaussian width in angstroms
		sigma=sqrt(itemp)*conv_factor
		ispec=ark.rec[index,i]/(sqrt(2.0*!pi)*sigma)*exp(-(lam-ark.lam[i])^2/(2.0*sigma^2))
		IF ark.type[i] EQ '30201' THEN spec[*,0,1]+=ispec
		IF ark.type[i] EQ '40201' THEN spec[*,1,1]+=ispec
		IF ark.type[i] NE '30201' AND ark.type[i] NE '40201' AND ark.type[i] NE '201' THEN spec[*,2,1]+=ispec
   
        ENDFOR

	yr=[0,max(sum_array(spec[*,*,1],/i))]*1.02e17
	plot, [0],[0],xr=[3.945,3.96],/xsty,yr=yr,/ysty,$
		xtit='Wavelength [Ang]',ytit='Norm. Spec. Emissivity '+n2g('epsilon')+'!l'+n2g('lambda')+'!n/n!iAr!nn!ie!n [10!u17!n m!u3!n/s/Ang]',chars=1.0
	oplot,lam,spec[*,0,0]*1.0e17,color=120,linestyle=3.0
	oplot,lam,spec[*,1,0]*1.0e17,color=30,linestyle=3.0
	oplot,lam,spec[*,2,0]*1.0e17,color=150,linestyle=3.0
	oplot,lam,sum_array(spec[*,*,0],/i)*1.0e17,linestyle=3.0,thick=4

	oplot,lam,spec[*,0,1]*1.0e17,color=120,linestyle=2.0
	oplot,lam,spec[*,1,1]*1.0e17,color=30,linestyle=2.0
	oplot,lam,spec[*,2,1]*1.0e17,color=150,linestyle=2.0
	oplot,lam,sum_array(spec[*,*,1],/i)*1.0e17,linestyle=2.0,thick=4

	

	stop
END

