FUNCTION pecwl_check_files,header
	path='/home/mlreinke/atomic_physics/adas/adf15/pec_wl/'
	spawn, 'ls '+path+header+'_*',output,error
	IF error EQ '' THEN RETURN, output ELSE RETURN,-1
END

FUNCTION pecwl_compile_pec,header,debug=debug,verb=verb,load=load
	files=pecwl_check_files(header)
	IF size(files,/type) NE 7 THEN RETURN,-1
	name=strsplit(files[0],'_wl0',/extract,/regex)
	save_path=name[0]+'.sav'
	IF keyword_set(load) THEN BEGIN
		restore,save_path
		RETURN,output
	ENDIF
	nfiles=n(files)+1
	nlam=fltarr(nfiles)
	FOR i=0,nfiles-1 DO BEGIN				;check for total # of lines and ensure gridding is the same
		ipec=read_pec_file(files[i])
		IF keyword_set(verb) THEN help, ipec.lam
		nlam[i]=n(ipec.lam)+1
		IF i EQ 0 THEN BEGIN
			ntemp=n(ipec.temp)+1
			ndens=n(ipec.dens)+1
			dens=ipec.dens
			temp=ipec.temp
			z=ipec.z
			q=ipec.q		
                ENDIF ELSE BEGIN
			IF total(dens-ipec.dens) NE 0 OR ndens NE n(ipec.dens)+1 THEN BEGIN
				print, 'density grids do not match'
				RETURN,-2
                        ENDIF
			IF total(temp-ipec.temp) NE 0 OR ntemp NE n(ipec.temp)+1 THEN BEGIN
				print, 'density grids do not match'
				RETURN,-3
                        ENDIF
		ENDELSE
		IF i EQ nfiles-1 THEN path=ipec.path		;will store path as final, indicating all WLXX < file
        ENDFOR

	pec=fltarr(total(nlam),ntemp,ndens)
	lam=fltarr(total(nlam))
	type=strarr(total(nlam))
	FOR i=0,nfiles-1 DO BEGIN
		ipec=read_pec_file(files[i])	
		IF i EQ 0 THEN pec[0:nlam[i]-1,*,*]=ipec.pec ELSE pec[total(nlam[0:i-1]):total(nlam[0:i])-1,*,*]=ipec.pec
		IF i EQ 0 THEN lam[0:nlam[i]-1]=ipec.lam ELSE lam[total(nlam[0:i-1]):total(nlam[0:i])-1]=ipec.lam
		IF i EQ 0 THEN type[0:nlam[i]-1]=ipec.type ELSE type[total(nlam[0:i-1]):total(nlam[0:i])-1]=ipec.type
        ENDFOR
	output={pec:pec,lam:lam,temp:temp,dens:dens,z:z,q:q,type:type,path:path}
	save,output,file=save_path
	IF keyword_set(debug) THEN stop
	RETURN,output
END

FUNCTION pec2csplc,pec,filter=filter,recom=recom
	IF keyword_set(filter) THEN BEGIN
		trans=interpol(filter.tr,alog10(filter.e),alog10(ang2ev(pec.lam)))
        ENDIF ELSE BEGIN
		trans=pec.lam*0.0+1.0
	ENDELSE
	IF keyword_set(recom) THEN type='RECOM' ELSE type='EXCIT'

	csplc=reform(pec.pec[0,*,*])*0.0
	FOR i=0,n(pec.lam) DO IF pec.type[i] EQ type THEN csplc+=pec.pec[i,*,*]*ang2ev(pec.lam[i])*trans[i]*1.60217657e-19
	RETURN,csplc
END

PRO mo_pecwl2plt,type=type,load=load
	IF NOT keyword_set(type) THEN type='ls'
	lines=[0,2,3]
	fq=gentran_fq(42)
	openwin,0
	plot,fq.temp,fq.fq[42,*],/xlog,yr=[0,1.0],/xsty
	for i=0,42 DO oplot, fq.temp,fq.fq[i,*],linestyle=lines[i mod 3]
	plc=plc(fq.temp,42)
	z=42
	dpt=1.0e20
	FOR i=0,z-1 DO BEGIN
		print, z-1-i
		header=type+'#mo'+strtrim(z-1-i,1)
		pec=pecwl_compile_pec(header,load=load)
		IF i EQ 0 THEN BEGIN
			index=ipt(pec.dens,dpt)
			pectemp=pec.temp
			ntemp=n(pec.temp)+1
			csplc=fltarr(ntemp,z+1)		;0+ is [*,0], 18+ is [*,18]
		ENDIF
		IF size(pec,/type) EQ 8 THEN BEGIN
			ics=pec2csplc(pec)
			csplc[*,z-1-i]=ics[*,index]
		ENDIF
	ENDFOR
	tit=type+'#mo'
	pecplc=fltarr(ntemp)
	FOR i=0,z DO pecplc+=interpol(fq.fq[i,*],fq.temp,pectemp)*csplc[*,i]
	
	plt=read_pxx_file('/home/mlreinke/atomic_physics/adas/plt89_mo.dat')
	prb=read_pxx_file('/home/mlreinke/atomic_physics/adas/prb89_mo.dat')
	pltplc=plt.temp*0.0
	prbplc=prb.temp*0.0
	index=ipt(plt.dens,dpt)
	FOR i=0,z DO pltplc+=interpol(fq.fq[i,*],fq.temp,plt.temp)*plt.plc[i,*,index]
	index=ipt(prb.dens,dpt)
	FOR i=0,z DO prbplc+=interpol(fq.fq[i,*],fq.temp,prb.temp)*prb.plc[i,*,index]

	openwin,1
	plot, fq.temp,plc,/xlog,xr=[10.0,2.0e4],yr=[1.0e-34,4.0e-31],ylog=0,/xsty,ytit='W*m!u3!n',xtit='Elec. Temp. [eV]',chars=1.2,tit=tit
	oplot, pectemp,pecplc,psym=8
	oplot,plt.temp,pltplc,color=200
	oplot,plt.temp,pltplc+interpol(prbplc,prb.temp,plt.temp),color=200,linestyle=2.0
	stop

	openwin,2
	plot, fq.temp,plc*1.0e13,/xlog,xr=[10.0,2.0e4],yr=[1.0e-24,1.0e-17],/ylog,/xsty,ytit='egs/s cm!u3!n',xtit='Elec. Temp. [eV]',chars=1.2,tit=tit
	oplot, pectemp,pecplc*1.0e13,psym=8
	oplot,plt.temp,pltplc*1.0e13,color=200
	
	openwin,3
	plot, [0],[0],/xlog,xr=[10.0,2.0e4],yr=[1.0e-35,4.0e-32],/ylog,/xsty,ytit='W*m!u3!n',xtit='Elec. Temp. [eV]',chars=1.2,tit=tit
	oplot, pectemp,pecplc
	FOR i=0,z DO BEGIN
		oplot, pectemp,csplc[*,i]*interpol(fq.fq[i,*],fq.temp,pectemp),linestyle=3.0
	ENDFOR

	stop		

END

PRO mofilt_pec,nofilter=nofilter
	path_filter='/usr/local/cmod/idl/atomic_physics/sxr_filters/be50.dat'
	IF NOT keyword_set(nofilter) THEN BEGIN
		filter=read_xray_data(path_filter)
		ytit='Filtered Power Loss [10!u-32!n W m^3]'
        ENDIF ELSE ytit='Power Loss [10!u-32!n W m^3]'

	type=['ls']
	color=[0]
	fq=gentran_fq(42)

	z=42
	dpt=1.0e20
	openwin,0
	!p.multi=0
	FOR j=0,n(type) DO BEGIN
		tit=type[j]+'#mo'
		plot,[0],[0],xr=[100,20000],/xlog,/xsty,yr=[0,1.3],xtit='E. Temp [eV]',ytit=ytit,tit=tit,chars=1.2
		FOR i=0,16 DO BEGIN
			print, z-1-i
			header=type[j]+'#mo'+strtrim(z-1-i,1)
			pec=pecwl_compile_pec(header,/load)
			IF i EQ 0 THEN BEGIN
				index=ipt(pec.dens,dpt)
				pectemp=pec.temp
				ntemp=n(pec.temp)+1
				csplc=fltarr(ntemp,z+1,4)		;0+ is [*,0], 18+ is [*,18]
			ENDIF
			IF size(pec,/type) EQ 8 THEN BEGIN
				ics=pec2csplc(pec,filter=filter)
				csplc[*,z-1-i,j]=ics[*,index]
				pow=interpol(fq.fq[z-1-i,*],fq.temp,pec.temp)*ics[*,index]*1.0e32
				oplot,pec.temp,pow,color=color[j]
				k=maxloc(pow)
				xyouts,pec.temp[k],pow[k]*1.02,'+'+num2str(z-1-i,1),color=color[j]
                        ENDIF

		ENDFOR

		
	ENDFOR

	stop		

	!p.multi=0
END

PRO arfilt_pec,nofilter=nofilter
	path_filter='/usr/local/cmod/idl/atomic_physics/sxr_filters/be50.dat'
	IF NOT keyword_set(nofilter) THEN BEGIN
		filter=read_xray_data(path_filter)
		ytit='Filtered Power Loss [10!u-33!n W m!u3!n]'
        ENDIF ELSE ytit='Power Loss [10!u-33!n W m!u3!n]'

	IF !d.name EQ 'PS' THEN ps=1 ELSE ps=0
	plotwin=0
	IF keyword_set(ps) THEN BEGIN
		xsize=8.0
		ysize=3.0
		ls=1.3
	ENDIF ELSE BEGIN
		xsize=1100.0
		ysize=1100.0*4/7
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

	type=['ls','ls','ls']
	noarf=[0,0,1]
	mlr=[0,1,0]
	color=[30,100,200]
	color2=[70,120,150]
	fq=gentran_fq(18)

	z=18
	dpt=1.0e20
	openwin,0
	!p.multi=[0,4,0]
	FOR j=0,n(type) DO BEGIN
		IF noarf[j] THEN tit=type[j]+'#ar' ELSE tit='arf40_'+type[j]+'#ar'
		IF mlr[j] THEN tit='mlr13_'+type[j]+'#ar'
		plot,[0],[0],xr=[100,10000],/xlog,/xsty,yr=[0,3.0],xtit='E. Temp [eV]',ytit=ytit,tit=tit,chars=1.2*ls
		FOR i=0,3 DO BEGIN
			print, z-1-i
			IF noarf[j] THEN header=type[j]+'#ar'+strtrim(z-1-i,1) ELSE header='arf40_'+type[j]+'#ar'+strtrim(z-1-i,1)
			IF mlr[j] THEN header='mlr13_'+type[j]+'#ar'+strtrim(z-1-i,1)
			pec=pecwl_compile_pec(header,/load)
			IF i EQ 0 THEN BEGIN
				index=ipt(pec.dens,dpt)
				pectemp=pec.temp
				ntemp=n(pec.temp)+1
			ENDIF
			IF size(pec,/type) EQ 8 THEN BEGIN
				ics=pec2csplc(pec,filter=filter)
				pow=interpol(fq.fq[z-1-i,*],fq.temp,pec.temp)*ics[*,index]*1.0e33
				oplot,pec.temp,pow,color=color[j]
				k=maxloc(pow)
				xyouts,pec.temp[k],pow[k]*1.02,'+'+num2str(z-1-i,1),color=color[j],chars=0.5*ls
				ics=pec2csplc(pec,filter=filter,/recom)
				pow=interpol(fq.fq[z-i,*],fq.temp,pec.temp)*ics[*,index]*1.0e33
				oplot,pec.temp,pow,color=color2[j],linestyle=3.0
				xyouts,pec.temp[k],pow[k]*1.02,'+'+num2str(z-i,1),color=color2[j],chars=0.5*ls

                        ENDIF

		ENDFOR

		
	ENDFOR
	plot,[0],[0],xr=[100,10000],/xlog,/xsty,yr=[0,3.0],xtit='E. Temp [eV]',ytit=ytit,tit='HULLAC',chars=1.2*ls
	IF keyword_set(nofilter) THEN BEGIN
		t_e=10^(make(1.0,4.0,500))
		n_e=t_e*0.0+1.0e20
		csplc=hullac_cs_plc(t_e,18,n_e=n_e)
        ENDIF ELSE restore,'/usr/local/cmod/idl/atomic_physics/hullac_kbf/ar/ar_xraytomo_csplc.dat'
	FOR i=0,3 DO BEGIN
		pow=interpol(fq.fq[z-1-i,*],fq.temp,t_e)*csplc[z-1-i,*]*1.0e33
		oplot,t_e,pow,color=0
		k=maxloc(pow)
		xyouts,t_e[k],pow[k]*1.02,'+'+num2str(z-1-i,1),color=0,chars=0.5*ls
 	ENDFOR

	stop		

	!p.multi=0
END

PRO ar_pecwl2plt,type=type,noarf=noarf,load=load,mlr=mlr
	IF NOT keyword_set(type) THEN type='ls'
	lines=[0,2,3]
	fq=gentran_fq(18)
	openwin,0
	plot,fq.temp,fq.fq[18,*],/xlog,yr=[0,1.0]
	for i=0,18 DO oplot, fq.temp,fq.fq[i,*],linestyle=lines[i mod 3]
	plc=plc(fq.temp,18)
	z=18
	dpt=1.0e20
	FOR i=0,z-1 DO BEGIN
		print, z-1-i
		IF keyword_set(noarf) THEN header=type+'#ar'+strtrim(z-1-i,1) ELSE header='arf40_'+type+'#ar'+strtrim(z-1-i,1)
		IF keyword_set(mlr) THEN header='mlr13_'+type+'#ar'+strtrim(z-1-i,1)
		pec=pecwl_compile_pec(header,load=load)
		IF i EQ 0 THEN BEGIN
			index=ipt(pec.dens,dpt)
			pectemp=pec.temp
			ntemp=n(pec.temp)+1
			csplc=fltarr(ntemp,z+1)		;0+ is [*,0], 18+ is [*,18]
		ENDIF
		IF size(pec,/type) EQ 8 THEN BEGIN
			ics=pec2csplc(pec)
			csplc[*,z-1-i]=ics[*,index]
			ics=pec2csplc(pec,/recom)
			csplc[*,z-i]+=ics[*,index]
		ENDIF
	ENDFOR
	IF keyword_set(noarf) THEN tit=type+'#ar' ELSE tit='arf40_'+type+'#ar'
	IF keyword_set(mlr) THEN tit='mlr13_'+type+'#ar'
	pecplc=fltarr(ntemp)
	FOR i=0,18 DO pecplc+=interpol(fq.fq[i,*],fq.temp,pectemp)*csplc[*,i]
	
	plt=read_pxx_file('/home/mlreinke/atomic_physics/adas/plt89_ar.dat')
	prb=read_pxx_file('/home/mlreinke/atomic_physics/adas/prb89_ar.dat')
	pltplc=plt.temp*0.0
	prbplc=plt.temp*0.0
	index=ipt(plt.dens,dpt)
	FOR i=0,z DO pltplc+=interpol(fq.fq[i,*],fq.temp,plt.temp)*plt.plc[i,*,index]
	FOR i=0,z DO prbplc+=interpol(fq.fq[i,*],fq.temp,plt.temp)*prb.plc[i,*,index]

	plt40=read_pxx_file('/home/mlreinke/atomic_physics/adas/plt40_ar.dat')
	plt40plc=plt.temp*0.0
	index=ipt(plt40.dens,dpt)
	FOR i=0,z DO plt40plc+=interpol(fq.fq[i,*],fq.temp,plt40.temp)*plt40.plc[i,*,index]


	;hcsplc=hullac_cs_plc(pectemp,18)
	;hplc=fltarr(ntemp)
	;FOR i=0,z DO hplc+=interpol(fq.fq[i,*],fq.temp,pectemp)*hcsplc[i,*]
	

	openwin,1
	plot, fq.temp,plc,/xlog,xr=[2.0,2.0e4],yr=[1.0e-34,4.0e-31],/ylog,/xsty,ytit='W*m!u3!n',xtit='Elec. Temp. [eV]',chars=1.2,tit=tit
	oplot, pectemp,pecplc,psym=8
	oplot,plt.temp,pltplc+prbplc,color=200
	;oplot,pectemp,hplc,color=100
	oplot,plt40.temp,plt40plc+interpol(prbplc,plt.temp,plt40.temp),color=30
	FOR i=0,z DO BEGIN
		;oplot, pectemp,csplc[*,i]*interpol(fq.fq[i,*],fq.temp,pectemp)
        ENDFOR

	openwin,2
	plot, fq.temp,plc*1.0e13,/xlog,xr=[2.0,2.0e4],yr=[1.0e-24,1.0e-17],/ylog,/xsty,ytit='egs/s cm!u3!n',xtit='Elec. Temp. [eV]',chars=1.2,tit=tit
	oplot, pectemp,pecplc*1.0e13,psym=8
	oplot,plt.temp,pltplc*1.0e13,color=200
	;oplot,pectemp,hplc*1.0e13,color=100
	oplot,plt40.temp,plt40plc*1.0e13,color=30
	FOR i=10,z DO BEGIN
		;oplot, pectemp,csplc[*,i]*interpol(fq.fq[i,*],fq.temp,pectemp)*1.0e13,psym=5
		;oplot, plt.temp,interpol(fq.fq[i,*],fq.temp,plt.temp)*plt.plc[i,*,index]*1.0e13,color=200,linestyle=2.0
	ENDFOR
	
	openwin,3
	plot, [0],[0],/xlog,xr=[2.0,2.0e4],yr=[1.0e-35,4.0e-32],/ylog,/xsty,ytit='W*m!u3!n',xtit='Elec. Temp. [eV]',chars=1.2,tit=tit
	oplot, pectemp,pecplc
	oplot,plt40.temp,plt40plc,color=30
	FOR i=0,z DO BEGIN
		oplot, pectemp,csplc[*,i]*interpol(fq.fq[i,*],fq.temp,pectemp),linestyle=3.0
		oplot,plt40.temp,interpol(fq.fq[i,*],fq.temp,plt40.temp)*plt40.plc[i,*,index],color=30,linestyle=2.0
	ENDFOR

	stop		

END
