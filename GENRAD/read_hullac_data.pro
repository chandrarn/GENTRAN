FUNCTION find_path,z,q,type
	;index_path='/home/mlreinke/idl/impurities/data/hullac_kbf/index'
	index_path='/usr/local/cmod/idl/atomic_physics/hullac_kbf/index'
	openr,lun,index_path,/get_lun
	line=strarr(1)
	readf,lun,line
	match='z='+num2str(z,1)+','+type
	WHILE eof(lun) NE 1 AND strmatch(line, match, /fold_case) EQ 0 DO BEGIN
		readf,lun,line
	ENDWHILE
	FOR i=0,q DO readf,lun,line
	tmp=strsplit(line,',',/extract)
	path=tmp[1]

	close,lun
	free_lun,lun
	output=path
	RETURN,output
END

FUNCTION read_brtable2_data,path,debug=debug,load=load,verb=verb
	file_path=path+'brtable2'
	save_path=path+'brtable_idl.dat'
	IF keyword_set(verb) THEN print, 'Loading: '+file_path
	IF keyword_set(load) THEN BEGIN
		restore, save_path
		RETURN,output
	ENDIF
	;default sizes
	maxlevels=100
	ndens=6
	dens=fltarr(ndens)

	openr,lun,file_path,/get_lun
	line=strarr(1)
	readf,lun,line
	
	;find temperature
	match='temperature'
	tmp=strsplit(line,' ',/extract)
	WHILE total(strmatch(tmp,match,/fold_case)) EQ 0 DO BEGIN
		readf,lun,line
		tmp=strsplit(line,' ',/extract)
	ENDWHILE
	temp=float(tmp[1:n(tmp)])
	ntemp=n(temp)+1
	readf,lun,line

	;find levels
	tmp=strsplit(line,' ',/extract)
	match='level'
	WHILE total(strmatch(tmp,match,/fold_case)) EQ 0 DO BEGIN
		readf,lun,line
		tmp=strsplit(line,' ',/extract)
	ENDWHILE
	level=intarr(maxlevels)
	name=strarr(maxlevels)
	jlevel=fltarr(maxlevels)
	readf,lun,line
	tmp=strsplit(line,' ',/extract)
	match='transitions'
	cntr=0
	WHILE total(strmatch(tmp,match,/fold_case)) EQ 0 DO BEGIN
	;	level[cntr]=int(tmp[0])
	;	name[cntr]=tmp[1]
	;	jlevel[cntr]=(flt(tmp[2])-1)/2.0
		readf,lun,line
		tmp=strsplit(line,' ',/extract)
		cntr+=1
	ENDWHILE
	ntrans=long(tmp[2])
	IF keyword_set(verb) THEN print, 'NTRANS: '+num2str(ntrans,1)	

	;find emissivity
	emiss=fltarr(ntemp,ndens,ntrans)
	lam=fltarr(ntrans)
	low=intarr(ntrans)
	up=intarr(ntrans)

	FOR i=0L,ntrans-1 DO BEGIN
		FOR k=0,1 DO readf,lun,line
		tmp=strsplit(line,' ',/extract)
		lam[i]=float(tmp[n(tmp)-1])
		tmp=tmp[3]
		tmp=strsplit(tmp,'<-',/extract)
		low[i]=int(tmp[0])
		up[i]=int(tmp[1])
		FOR k=0,1 DO readf,lun,line
		FOR j=0,ndens-1 DO BEGIN
			readf,lun,line
			tmp=strsplit(line,'|',/extract)	
			IF i EQ 0 THEN dens[j]=float(tmp[0]*1.0e6)		;density in m^-3
			emiss[*,j,i]=float(strsplit(tmp[1],' ',/extract))	;ph/s/iondens (has n_e multiplied in)
		ENDFOR
	ENDFOR		

	close,lun
	free_lun,lun

	output={em:emiss,dens:dens,temp:temp,lam:lam}
	IF keyword_set(debug) THEN stop
	save,output,filename=save_path
	RETURN,output
END

;interpolate as per eq (7-9) in AD&NDT argon paper
FUNCTION interp_emdat,n_e,t_e,emdat,debug=debug,noextrap=noextrap
	n_pts=n(n_e)+1
	IF n_pts NE n(t_e)+1 THEN RETURN,-1
	em_int=fltarr(n_pts)
	dens=emdat.dens
	n_dens=n(dens)+1
	temp=emdat.temp
	n_temp=n(temp)+1
	FOR i=0,n_pts-1 DO BEGIN
		dpt=ipt(dens,n_e[i])
		IF dpt[0] EQ -1 AND n_e[i] GT max(dens) THEN dpt=n_dens-1
		IF dpt[0] EQ -1 AND n_e[i] LT min(dens) THEN dpt=[0]
		em=reform(emdat.em[*,dpt])/dens[dpt]
;		tbnd=ibound(temp,t_e[i])
;		IF tbnd[0] EQ -1 AND t_e[i] GT max(temp) THEN em_int[i]=em[n_temp-2]*(1.0-(t_e[i]-temp[n_temp-2])$
;			/(temp[n_temp-1]-temp[n_temp-2]))+em[n_temp-1]*(t_e[i]-temp[n_temp-2])/(temp[n_temp-1]-temp[n_temp-2])
;		IF tbnd[0] EQ -1 AND t_e[i] LT min(temp) THEN em_int[i]=em[0]*t_e[i]/temp[0]
;		IF tbnd[0] NE -1 THEN em_int[i]=em[tbnd[0]]*(1.0-(t_e[i]-temp[tbnd[0]])$
;			/(temp[tbnd[1]]-temp[tbnd[0]]))+em[tbnd[1]]*(t_e[i]-temp[tbnd[0]])/(temp[tbnd[1]]-temp[tbnd[0]])
		IF NOT keyword_set(noextrap) THEN IF t_e[i] LT min(temp) THEN em_int[i]=em[0]*t_e[i]/temp[0]
		IF NOT keyword_set(noextrap) THEN IF t_e[i] GT max(temp) THEN em_int[i]=interpol(em,alog10(temp),alog10(t_e[i]))
		IF t_e[i] GE min(temp) AND t_e[i] LE max(temp) THEN em_int[i]=interpol(em,temp,t_e[i])
	ENDFOR

	output=em_int
	IF keyword_set(debug) THEN stop
	RETURN,output
END
				

FUNCTION calc_line_plc,t_e,z,n_e=n_e,plot=plot,filter=filter
	ntemp=n(t_e)+1
	IF NOT keyword_set(n_e) THEN n_e=1.0e20
	n_e=fltarr(ntemp)+n_e
	plc=dblarr(ntemp)
	fq,z,cs,cs_te
	frac=fltarr(z+1,ntemp)
	FOR i=0,z DO frac[i,*]=interpol(cs[i,*],cs_te,t_e)			

	FOR i=0,z DO BEGIN
		print, i
		path=find_path(z,i,'lines')
		IF path NE 'NA' THEN BEGIN
			emdat=read_brtable2_data(path,/load)
			ph_energy=ang2ev(emdat.lam)	
			IF keyword_set(filter) THEN trans=interpol(filter.tr,filter.e,ph_energy) ELSE trans=fltarr(n(ph_energy)+1)+1.0
			emiss=fltarr(n(emdat.temp)+1,n(emdat.dens)+1)
			FOR j=0L,n(ph_energy) DO emiss+=emdat.em[*,*,j]*ph_energy[j]*trans[j]*1.6e-19
			emdat_sum={em:emiss,temp:emdat.temp,dens:emdat.dens}
			plc+=interp_emdat(n_e,t_e,emdat_sum)*frac[i,*]
		ENDIF
	ENDFOR
	IF keyword_set(plot) THEN BEGIN
		plcplt,zz=[z]
		oplot,t_e,plc,psym=8
	ENDIF
	output=plc
	RETURN,output
END
	
FUNCTION hullac_cs_plc,t_e,z,n_e=n_e,filter=filter,verb=verb,raw=raw,noextrap=noextrap
	IF keyword_set(raw) THEN BEGIN
		ntemp=6
		t_e=fltarr(z+1,ntemp)
	ENDIF ELSE ntemp=n(t_e)+1
	IF NOT keyword_set(n_e) THEN n_e=1.0e20
	n_e=fltarr(ntemp)+n_e
	plc=dblarr(z+1,ntemp)

	FOR i=0,z DO BEGIN
		path=find_path(z,i,'lines')
		IF path NE 'NA' THEN BEGIN
			IF keyword_set(verb) THEN print, 'Reading: '+path
			emdat=read_brtable2_data(path,/load)
			ph_energy=ang2ev(emdat.lam)	
			IF keyword_set(filter) THEN trans=interpol(filter.tr,filter.e,ph_energy) ELSE trans=fltarr(n(ph_energy)+1)+1.0
			emiss=fltarr(n(emdat.temp)+1,n(emdat.dens)+1)
			FOR j=0L,n(ph_energy) DO emiss+=emdat.em[*,*,j]*ph_energy[j]*trans[j]*1.6e-19
			IF keyword_set(raw) THEN BEGIN
				t_e[i,*]=emdat.temp
				dpt=ipt(emdat.dens,n_e[0])
				plc[i,*]=emiss[*,dpt]/emdat.dens[dpt]
			ENDIF ELSE BEGIN
				emdat_sum={em:emiss,temp:emdat.temp,dens:emdat.dens}
				plc[i,*]=interp_emdat(n_e,t_e,emdat_sum,noextrap=noextrap)
			ENDELSE
		ENDIF
	ENDFOR
	output=plc
	RETURN,output
END

FUNCTION calc_spec,z,t_e,dlam,n_e=n_e,emin=emin,qmin=qmin,qmax=qmax,terr=terr
	IF NOT keyword_set(n_e) THEN n_e=1.0e20
	n_temp=n(t_e)+1
	spec=fltarr(1000000,n_temp)
	err=fltarr(1000000,n_temp)
	lam=fltarr(1000000)
	q=fltarr(1000000)
	cntr=0L
	IF NOT keyword_set(qmax) THEN qmax=z
	IF NOT keyword_set(qmin) THEN qmin=0
	IF NOT keyword_set(emin) THEN emin=0
	FOR i=qmin,qmax DO BEGIN
		path=find_path(z,i,'lines')
		IF path NE 'NA' THEN BEGIN
			emdat=read_brtable2_data(path,/load)
			tmp=where(emdat.lam GE dlam[0] AND emdat.lam LE dlam[1] AND emdat.em GE emin)
			IF tmp[0] NE -1 THEN BEGIN; AND max(emdat.temp) GE t_e/3.0 AND min(emdat.temp) LE t_e*3.0 THEN BEGIN
				print, i
				ntmp=n(tmp)+1
				ndens=n(emdat.dens)+1
				temp=emdat.temp
				em=reform((emdat.em[*,ndens-1,tmp]/emdat.dens[ndens-1]+emdat.em[*,ndens-2,tmp]/emdat.dens[ndens-2])/2.0	)
				lam[cntr:cntr+ntmp-1]=emdat.lam[tmp]
				FOR j=0L,ntmp-1 DO BEGIN
					spec[cntr+j,*]=interpol(em[*,j],temp,t_e)
					IF keyword_set(terr) THEN err[cntr+j,*]=interpol(deriv(temp,em[*,j]),temp,t_e)^2*terr^2
				ENDFOR
				q[cntr:cntr+ntmp-1]=fltarr(ntmp)+i
				FOR k=0L,n(tmp) DO BEGIN
					check=where(lam[cntr:cntr+ntmp-1] EQ emdat.lam[tmp[k]])
					IF n(check) NE 0 THEN BEGIN
						spec[cntr+check[0],*]+=spec[cntr+check[1],*]
						err[cntr+check[0],*]+=err[cntr+check[1],*]
						spec[cntr+check[1],*]=0
						err[cntr+check[1],*]=0
					ENDIF
					
				ENDFOR
						
				cntr+=ntmp
			ENDIF
		ENDIF
	ENDFOR

	;check for intersection
	tmp=where(lam NE 0)
	output={sigv:spec[tmp,*],lam:lam[tmp],q:q[tmp],temp:t_e,z:z,err:sqrt(err[tmp,*])}
	RETURN,output
END

PRO plot_spec,z,temp,dlam,qmin=qmin,qmax=qmax,labmax=labmax,emin=emin
	IF NOT keyword_set(labmax) THEN labmax = 0.2

	ntemp=n(temp)+1

	!p.multi=[0,0,ntemp]

	FOR i=0,ntemp-1 DO BEGIN
		out=calc_spec(z,temp[i],dlam,qmin=qmin,qmax=qmax,emin=emin)
		print, max(out.sigv)
		out.sigv/=max(out.sigv)	
		plot,out.lam,out.sigv,/nodata,xr=dlam,/xsty,yr=[0,1.2],/ysty,xtit='Wavelength [Ang]',ytit='Normalized Emissivity',$
			tit='Z = '+num2str(z,1)+' T!le!n='+num2str(temp*1.0e-3,dp=2)+' [keV]',chars=1.2
		FOR j=0L,n(out.lam) DO BEGIN
			oplot,out.lam[j]*[1,1],[0.0,out.sigv[j]],color=200
			IF out.sigv[j] GT labmax THEN xyouts, out.lam[j]+(dlam[1]-dlam[0])*0.005,out.sigv[j]+0.03,num2rom(out.q[j]+1),orient=90,color=100
		ENDFOR
	ENDFOR
	!p.multi=0
END
