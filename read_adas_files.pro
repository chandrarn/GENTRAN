;6/6/12 - modified to use version controlled atomic physics directory
PRO adas_path,z,ion,rec,cxr
	;base_path='/home/mlreinke/idl/impurities/data/adas/'
	base_path='/usr/local/cmod/idl/atomic_physics/adas/'
	CASE z OF 
		1 : BEGIN
			rpath='acd96_h.dat'
			ipath='scd96_h.dat'
			cpath='NA'
		END

		2 : BEGIN
			rpath='acd96_he.dat'
			ipath='scd96_he.dat'
			cpath='NA'
		END

		7 : BEGIN
			rpath='acd96_n.dat'
			ipath='scd96_n.dat'
			cpath='NA'
		END

		10 : BEGIN
			rpath='acd96_ne.dat'
			ipath='scd96_ne.dat'
			cpath='NA'
		END

		18 : BEGIN
			rpath='acd89_ar.dat'
			ipath='scd89_ar.dat'
			cpath='ccd89_ar.dat'
		END

		20 : BEGIN
			rpath='acd85_ca.dat'
			ipath='scd85_ca.dat'
			cpath='NA'
		END

		36 : BEGIN
			rpath='acd89_kr.dat'
			ipath='scd89_kr.dat'
			cpath='ccd89_kr.dat'
		END	

		42 : BEGIN
			rpath='acd89_mo.dat'
			ipath='scd89_mo.dat'
			cpath='ccd89_mo.dat'
		END
		54 : BEGIN
			rpath='acd89_xe.dat'
			ipath='scd89_xe.dat'
			cpath='ccd89_xe.dat'
		END
		74 : BEGIN
			rpath='acd50_w.dat'
			ipath='scd50_w.dat'
			cpath='NA'
		END
		ELSE : BEGIN
			ipath='NA'
			rpath='NA'
			cpath='NA'
		END
	ENDCASE
	IF ipath NE 'NA' THEN ion=base_path+ipath
	IF rpath NE 'NA' THEN rec=base_path+rpath
	IF cpath NE 'NA' THEN cxr=base_path+cpath ELSE cxr=cpath
END


FUNCTION read_scd_file,path,debug=debug
	nrow=8

	openr,lun,path,/get_lun
	line=strarr(1)
	readf,lun,line
	tmp=strsplit(line,' ',/extract)
	z=int(tmp[0])
	ndens=int(tmp[1])	;number of density
	ntemp=int(tmp[2])	;number of temperature
	qmin=int(tmp[3])-1	;lowest charge state
	qmax=int(tmp[4])	;highest charge state
	readf,lun,line	
	rates=dblarr(z+1,ntemp,ndens)

	;read density
	nread_dens=ndens/nrow
	IF ndens MOD nrow NE 0 THEN nread_dens+=1
	FOR i=0,nread_dens-1 DO BEGIN
		readf,lun,line
		IF i EQ 0 THEN superline=line ELSE superline=superline+' '+line
	ENDFOR
	tmp=strsplit(superline,' ',/extract)
	dens=10^(float(tmp))*1.0e6 	;density in m^-3

	;read temperature
	nread_temp=ntemp/nrow
	IF ntemp MOD nrow NE 0 THEN nread_temp+=1
	FOR i=0,nread_temp-1 DO BEGIN
		readf,lun,line
		IF i EQ 0 THEN superline=line ELSE superline=superline+' '+line
	ENDFOR
	tmp=strsplit(superline,' ',/extract)
	temp=10^(float(tmp)) 		;electron temperature in eV


	;fill rates data
	FOR k=qmin,qmax-1 DO BEGIN
		readf,lun,line
		tmp=strtrim( strsplit(line,'/',/extract), 2 )
		whrZ1=where(strmid(tmp,0,2) eq 'Z1') ; select the Z1 tag
		tmp=strsplit(tmp[whrZ1], '=',/extract)
		q=int(tmp[1])-1
		FOR j=0,ntemp-1 DO BEGIN
			FOR i=0,nread_dens-1 DO BEGIN
				readf,lun,line
				IF i EQ 0 THEN superline=line ELSE superline=superline+' '+line
			ENDFOR
			tmp=strsplit(superline,' ',/extract)
			rates[q,j,*]=10^(float(tmp))*1.0e-6	;rates in m^-3/s
		ENDFOR
	ENDFOR

	close,lun
	free_lun,lun

	output={rates:rates,temp:temp,dens:dens,z:z,path:path}
	IF keyword_set(debug) THEN stop
	RETURN,output
END

FUNCTION read_acd_file,path,debug=debug
	nrow=8

	openr,lun,path,/get_lun
	line=strarr(1)
	readf,lun,line
	tmp=strsplit(line,' ',/extract)
	z=int(tmp[0])
	ndens=int(tmp[1])	;number of density
	ntemp=int(tmp[2])	;number of temperature
	qmin=int(tmp[3])-1	;lowest charge state
	qmax=int(tmp[4])	;highest charge state
	readf,lun,line	
	rates=dblarr(z+1,ntemp,ndens)

	;read density
	nread_dens=ndens/nrow
	IF ndens MOD nrow NE 0 THEN nread_dens+=1
	FOR i=0,nread_dens-1 DO BEGIN
		readf,lun,line
		IF i EQ 0 THEN superline=line ELSE superline=superline+' '+line
	ENDFOR
	tmp=strsplit(superline,' ',/extract)
	dens=10^(float(tmp))*1.0e6 	;density in m^-3

	;read temperature
	nread_temp=ntemp/nrow
	IF ntemp MOD nrow NE 0 THEN nread_temp+=1
	FOR i=0,nread_temp-1 DO BEGIN
		readf,lun,line
		IF i EQ 0 THEN superline=line ELSE superline=superline+' '+line
	ENDFOR
	tmp=strsplit(superline,' ',/extract)
	temp=10^(float(tmp)) 		;electron temperature in eV


	;fill rates data
	FOR k=qmin,qmax-1 DO BEGIN
		readf,lun,line
		tmp=strtrim( strsplit(line,'/',/extract), 2 )
		whrZ1=where(strmid(tmp,0,2) eq 'Z1') ; select the Z1 tag
		tmp=strsplit(tmp[whrZ1], '=',/extract)
		q=int(tmp[1])-1
		FOR j=0,ntemp-1 DO BEGIN
			FOR i=0,nread_dens-1 DO BEGIN
				readf,lun,line
				IF i EQ 0 THEN superline=line ELSE superline=superline+' '+line
			ENDFOR
			tmp=strsplit(superline,' ',/extract)
			rates[q+1,j,*]=10^(float(tmp))*1.0e-6	;rates in m^-3/s
		ENDFOR
	ENDFOR

	close,lun
	free_lun,lun
	
	output={rates:rates,temp:temp,dens:dens,z:z,path:path}
	IF keyword_set(debug) THEN stop
	RETURN,output
END

FUNCTION read_pxx_file,path,debug=debug,type=type
	nrow=8
	
	openr,lun,path,/get_lun
	line=strarr(1)
	readf,lun,line
	tmp=strsplit(line,' ',/extract)
	z=int(tmp[0])
	ndens=int(tmp[1])	;number of density
	ntemp=int(tmp[2])	;number of temperature
	qmin=int(tmp[3])-1	;lowest charge state
	qmax=int(tmp[4])	;highest charge state
	readf,lun,line	
	plc=dblarr(z+1,ntemp,ndens)
	IF NOT keyword_set(type) THEN type='-'
 	wave=fltarr(z+1)

	;read density
	nread_dens=ndens/nrow
	IF ndens MOD nrow NE 0 THEN nread_dens+=1
	FOR i=0,nread_dens-1 DO BEGIN
		readf,lun,line
		IF i EQ 0 THEN superline=line ELSE superline=superline+' '+line
	ENDFOR
	tmp=strsplit(superline,' ',/extract)
	dens=10^(float(tmp))*1.0e6 	;density in m^-3

	;read temperature
	nread_temp=ntemp/nrow
	IF ntemp MOD nrow NE 0 THEN nread_temp+=1
	FOR i=0,nread_temp-1 DO BEGIN
		readf,lun,line
		IF i EQ 0 THEN superline=line ELSE superline=superline+' '+line
	ENDFOR
	tmp=strsplit(superline,' ',/extract)
	temp=10^(float(tmp)) 		;electron temperature in eV


	;fill rates data
	FOR k=qmin,qmax-1 DO BEGIN
		readf,lun,line
		tmp=strsplit(line,'/',/extract)
		IF strlowcase(type) EQ 'pls' THEN wave[k]=float(tmp[3])
		tmp=strsplit(line, 'Z1=',/extract,/regex)
		IF tmp[0] NE line THEN BEGIN	;most common usage
			tmp=strsplit(tmp[1],/extract)
			q=int(tmp[0])-1
                ENDIF ELSE BEGIN		;this handles a specific case of SXRAUG data
			tmp=strsplit(line,'*',/extract)
			q=int(tmp[0])
		ENDELSE
		FOR j=0,ntemp-1 DO BEGIN
			FOR i=0,nread_dens-1 DO BEGIN
				readf,lun,line
				IF i EQ 0 THEN superline=line ELSE superline=superline+' '+line
			ENDFOR
			tmp=strsplit(superline,' ',/extract)
			plc[q,j,*]=10^(float(tmp))*1.0e-6	;rates in MW/m^-3
		ENDFOR
	ENDFOR

	close,lun
	free_lun,lun
	
	output={plc:plc,temp:temp,dens:dens,z:z,path:path,wave:wave,type:type}
	IF keyword_set(debug) THEN stop
	RETURN,output
END

FUNCTION read_sxb_file,path,debug=debug,type=type
	nrow=8
	
	openr,lun,path,/get_lun
	;get header information
	line=strarr(1)
	readf,lun,line
	tmp=strsplit(line,'/',/extract)
	ntrans=int(tmp[0])
	tmp=strsplit(tmp[1],' ',/extract)
	IF NOT strmatch(tmp[1], '+') THEN BEGIN
		tmp0=strsplit(tmp[0],'+',/extract)
		zstr=tmp0[0]
		q=int(tmp[1])
	ENDIF ELSE BEGIN
		zstr=tmp[0]
		q=int(tmp[2])
	ENDELSE
		

	readf,lun,line
	tmp=strsplit(line,' ',/extract)
	ndens=int(tmp[2])	;number of density		(assume all transitions are the same)
	ntemp=int(tmp[3])	;number of temperature		(assume all transitions are the same)
	sxb=dblarr(ntrans,ntemp,ndens)
	lam=fltarr(ntrans)
	lam[0]=float(strsplit(tmp[0],'A',/extract))					;first wavelength

	;read density
	nread_dens=ndens/nrow
	IF ndens MOD nrow NE 0 THEN nread_dens+=1
	FOR i=0,nread_dens-1 DO BEGIN
		readf,lun,line
		IF i EQ 0 THEN superline=line ELSE superline=superline+' '+line
	ENDFOR
	tmp=strsplit(superline,' ',/extract)
	dens=float(tmp)*1.0e6 	;density in m^-3

	;read temperature
	nread_temp=ntemp/nrow
	IF ntemp MOD nrow NE 0 THEN nread_temp+=1
	FOR i=0,nread_temp-1 DO BEGIN
		readf,lun,line
		IF i EQ 0 THEN superline=line ELSE superline=superline+' '+line
	ENDFOR
	tmp=strsplit(superline,' ',/extract)
	temp=float(tmp) 		;electron temperature in eV


	;fill rates data
	FOR k=0,ntrans-1 DO BEGIN		
		FOR j=0,ndens-1 DO BEGIN
			FOR i=0,nread_temp-1 DO BEGIN
				readf,lun,line
				IF i EQ 0 THEN superline=line ELSE superline=superline+' '+line
			ENDFOR
			tmp=strsplit(superline,' ',/extract)
			sxb[k,*,j]=double(tmp)	;rates in ionization/photon
		ENDFOR
		IF k NE ntrans-1 THEN BEGIN
			readf,lun,line
			tmp=strsplit(line,' ',/extract)
			lam[k+1]=float(tmp[0])
			FOR i=0,(nread_temp+nread_dens-1) DO readf,lun,line
		ENDIF
	ENDFOR

	close,lun
	free_lun,lun
	
    	output={sxb:sxb,lam:lam,temp:temp,dens:dens,z:zstr,q:q,path:path}
	IF keyword_set(debug) THEN stop
	RETURN,output
END

;7/24/13 - fixed the header read which was failing on single character
;          charge states (i.e. AR+ 9 vs. AR+10)
;9/16/14 - modified header read to read pec#40 format (should be backcompat)

FUNCTION read_pec_file,path,debug=debug
	nrow=8
	
	openr,lun,path,/get_lun
	;get header information
	line=strarr(1)
	readf,lun,line
	tmp=strsplit(line,'/',/extract)
	ntrans=int(tmp[0])
	IF strmatch(line, '*+*') THEN tmp=strsplit(tmp[1],'+',/extract)		;two options for header
	IF strmatch(line, '*:*') THEN tmp=strsplit(tmp[1],':',/extract)
	zstr=tmp[0]
	tmp=strsplit(tmp[1],' ',/extract)
	q=int(tmp[0])		

	readf,lun,line
	WHILE strmatch(line, '*isel*') EQ 0 AND strmatch(line, '*ISEL*') EQ 0 DO readf,lun,line
	tmp=strsplit(line,' ',/extract)
	IF tmp[1] EQ 'A' THEN BEGIN	;two options for line header
		ndens=int(tmp[2])	;number of density		(assume all transitions are the same)
		ntemp=int(tmp[3])	;number of temperature		(assume all transitions are the same)
        ENDIF ELSE BEGIN
		ndens=int(tmp[1])	;number of density		(assume all transitions are the same)
		ntemp=int(tmp[2])	;number of temperature		(assume all transitions are the same)
  	ENDELSE
	pec=dblarr(ntrans,ntemp,ndens)
	lam=fltarr(ntrans)
	lam[0]=float(tmp[0])					;first wavelength
	type=strarr(ntrans)
	IF strmatch(line, '*EXCIT*') OR strmatch (line, '*excit*') THEN type[0]='EXCIT' ELSE type[0]='RECOM'
	
	;read density
	nread_dens=ndens/nrow
	IF ndens MOD nrow NE 0 THEN nread_dens+=1
	FOR i=0,nread_dens-1 DO BEGIN
		readf,lun,line
		IF i EQ 0 THEN superline=line ELSE superline=superline+' '+line
	ENDFOR
	tmp=strsplit(superline,' ',/extract)
	dens=float(tmp)*1.0e6 	;density in m^-3

	;read temperature
	nread_temp=ntemp/nrow
	IF ntemp MOD nrow NE 0 THEN nread_temp+=1
	FOR i=0,nread_temp-1 DO BEGIN
		readf,lun,line
		IF i EQ 0 THEN superline=line ELSE superline=superline+' '+line
	ENDFOR
	tmp=strsplit(superline,' ',/extract)
	temp=float(tmp) 		;electron temperature in eV


	;fill rates data
	FOR k=0,ntrans-1 DO BEGIN		
		FOR j=0,ndens-1 DO BEGIN
			FOR i=0,nread_temp-1 DO BEGIN
				readf,lun,line
				IF i EQ 0 THEN superline=line ELSE superline=superline+' '+line
			ENDFOR
			tmp=strsplit(superline,' ',/extract)
			pec[k,*,j]=double(tmp)*1.0e-6	;rates in m^-3/s
		ENDFOR
		IF k NE ntrans-1 THEN BEGIN
			readf,lun,line
			tmp=strsplit(line,' ',/extract)
			lam[k+1]=float(tmp[0])
			IF strmatch(line, '*EXCIT*') OR strmatch(line, '*excit*') THEN type[k+1]='EXCIT' ELSE type[k+1]='RECOM'
			FOR i=0,(nread_temp+nread_dens-1) DO readf,lun,line
		ENDIF
	ENDFOR

	close,lun
	free_lun,lun
	
    	output={pec:pec,lam:lam,temp:temp,dens:dens,z:zstr,q:q,path:path,type:type}
	IF keyword_set(debug) THEN stop
	RETURN,output
END
;6/6/12 - modified to use version controlled atomic physics directory
FUNCTION adas_pec_path,z,q
	CASE z OF
		 5 : BEGIN
			;bp='/home/mlreinke/idl/impurities/data/adas/'
			bp='/usr/local/cmod/idl/atomic_physics/adas/'
			CASE q OF
				0 : pec_path = bp+'pec93#b_llu#b0.dat'
				1 : pec_path =bp+'pec93#b_llu#b0.dat'
				2 : pec_path =bp+'pec93#b_llu#b0.dat'
				3 : pec_path =bp+'pec93#b_llu#b0.dat'
				4 : pec_path =bp+'pec93#b_llr#b4.dat'
				ELSE : pec_path=''
			ENDCASE
		END
		18 : BEGIN
			;bp='/home/mlreinke/idl/impurities/data/adas/'
			bp='/usr/local/cmod/idl/atomic_physics/adas/'
			CASE q OF
				17 : pec_path=bp+'transport_llu#ar17.dat'
				16 : pec_path=bp+'transport_llu#ar16.dat'
				15 : pec_path=bp+'transport_llu#ar15ic.dat'
				14 : pec_path=bp+'transport_llu#ar14ic.dat'
				ELSE : pec_path=''
			ENDCASE
		END
		74 : BEGIN
			;bp='/home/mlreinke/idl/impurities/data/adas/'
			bp='/usr/local/cmod/idl/atomic_physics/adas/'
			CASE q OF
				44 : pec_path=bp+'W_44_adf15.dat'
				ELSE : pec_path=''
			ENDCASE
		END
		ELSE : pec_path=''
	ENDCASE
	RETURN,pec_path
 END

FUNCTION adas_pec_spec,z,temp,dens,dlam
	pec_path=strarr(z+1)
	FOR i=0,z DO pec_path[i]=adas_pec_path(z,i)
	n_temp=n(temp)+1
	spec=fltarr(1000,n_temp)
	lam=fltarr(1000)
	q=fltarr(1000)
	cntr=0
	ivec=alog10(temp)
	jvec=alog10(dens)

	FOR i=0,z DO IF pec_path[i] NE '' THEN BEGIN
		out=read_pec_file(pec_path[i])
		IF keyword_set(verb) THEN BEGIN
			print, 'PEC from: '+pec_path[i]
			print, num2str(min(out.temp),dp=2)+' < Te [eV] < '+num2str(max(out.temp),dp=2)
			print, num2str(min(out.dens),dp=2)+' < ne [m^-3] < '+num2str(max(out.dens),dp=2)
		ENDIF
		tmp=where(out.lam GE dlam[0] AND out.lam LE dlam[1])
		IF tmp[0] NE -1 THEN FOR j=0,n(tmp) DO BEGIN
			pect=alog10(out.temp)
			pecd=alog10(out.dens)
			ipts=interp_vec_reform(pect,ivec)
			jpts=interp_vec_reform(pecd,jvec)
			spec[cntr,*]=interpolate(reform(out.pec[tmp[j],*,*]),ipts,jpts)
			q[cntr]=out.q
			lam[cntr]=out.lam[tmp[j]]
			cntr+=1
		ENDFOR
	ENDIF
	tmp=where(lam NE 0)
	output={sigv:spec[tmp,*], lam:lam[tmp],q:q[tmp],temp:temp,z:z}
	RETURN,output
END

;6/6/12 - modified to use version controlled atomic physics directory
FUNCTION adas_pec_csplc,z,temp,dens,filter=filter,verb=verb,debug=debug

	pec_path=strarr(z+1)
	CASE z OF
		 5 : BEGIN
			;bp='/home/mlreinke/idl/impurities/data/adas/'
			bp='/usr/local/cmod/idl/atomic_physics/adas/'
			FOR i=0,3 DO pec_path[i]=bp+'pec93#b_llu#b'+num2str(i,1)+'.dat'
			pec_path[4]=bp+'pec93#b_llr#b4.dat'
		END
		18 : BEGIN
			;bp='/home/mlreinke/idl/impurities/data/adas/'
			bp='/usr/local/cmod/idl/atomic_physics/adas/'
			pec_path[17]=bp+'transport_llu#ar17.dat'
			pec_path[16]=bp+'transport_llu#ar16.dat'
			pec_path[15]=bp+'transport_llu#ar15ic.dat'
			pec_path[14]=bp+'transport_llu#ar14ic.dat'
		END

		ELSE : RETURN,-1
	ENDCASE
	ntemp=n(temp)+1
	ndens=n(dens)+1
	ivec=alog10(temp)
	jvec=alog10(dens)
	plc=dblarr(z+1,ntemp,ndens)

	FOR i=0,z DO IF pec_path[i] NE '' THEN BEGIN
		out=read_pec_file(pec_path[i])
		IF keyword_set(verb) THEN BEGIN
			print, 'PEC from: '+pec_path[i]
			print, num2str(min(out.temp),dp=2)+' < Te [eV] < '+num2str(max(out.temp),dp=2)
			print, num2str(min(out.dens),dp=2)+' < ne [m^-3] < '+num2str(max(out.dens),dp=2)
		ENDIF
		pect=alog10(out.temp)
		pecd=alog10(out.dens)
		ipts=interp_vec_reform(pect,ivec)
		jpts=interp_vec_reform(pecd,jvec)
		nlam=n(out.lam)
		ph_energy=ang2ev(out.lam)	
		IF keyword_set(filter) THEN trans=interpol(filter.tr,filter.e,ph_energy) ELSE trans=fltarr(nlam)+1.0
		FOR j=0,nlam-1 DO plc[i,*,*]+=interpolate(reform(out.pec[j,*,*]),ipts,jpts,/grid)*trans[j]*ph_energy[j]*1.6e-19
	ENDIF
	output=plc
	IF keyword_set(debug) THEN stop
	RETURN,output

END

;6/6/12 - modified to use version controlled atomic physics directory
FUNCTION adas_pxx_csplc,z,temp,dens,filter=filter,verb=verb,debug=debug,prb=prb,prc=prc


	CASE z OF
		 2 : BEGIN
			bp='/usr/local/cmod/idl/atomic_physics/adas/'
			pltpath=bp+'plt96_he.dat'
		 END
		 5 : BEGIN
			;bp='/home/mlreinke/idl/impurities/data/adas/'
			bp='/usr/local/cmod/idl/atomic_physics/adas/'
			pltpath=bp+'plt89_b.dat'
			prbpath=bp+'prb89_b.dat'
			
		END
		 7 : BEGIN
			bp='/usr/local/cmod/idl/atomic_physics/adas/'
			pltpath=bp+'plt96_n.dat'
		 END
		 9 : BEGIN
			;bp='/home/mlreinke/idl/impurities/data/adas/'
			bp='/usr/local/cmod/idl/atomic_physics/adas/'
			pltpath=bp+'plt89_f.dat'
			prbpath=bp+'prb89_f.dat'
			
                 END
		 10 : BEGIN
			bp='/usr/local/cmod/idl/atomic_physics/adas/'
			pltpath=bp+'plt96_ne.dat'
		 END
		 18 : BEGIN
			;bp='/home/mlreinke/idl/impurities/data/adas/'
			bp='/usr/local/cmod/idl/atomic_physics/adas/'
			pltpath=bp+'plt89_ar.dat'
			prbpath=bp+'prb89_ar.dat'
			prcpath=bp+'prc89_ar.dat'
	
			
		END
		 36 : BEGIN
			;bp='/home/mlreinke/idl/impurities/data/adas/'
			bp='/usr/local/cmod/idl/atomic_physics/adas/'
			pltpath=bp+'plt89_kr.dat'
			prbpath=bp+'prb89_kr.dat'
			prcpath=bp+'prc89_kr.dat'
			
		END
		 42 : BEGIN
			;bp='/home/mlreinke/idl/impurities/data/adas/'
			bp='/usr/local/cmod/idl/atomic_physics/adas/'
			pltpath=bp+'plt89_mo.dat'
			prbpath=bp+'prb89_mo.dat'
			prcpath=bp+'prc89_mo.dat'
		
		END


		ELSE : RETURN,-1
	ENDCASE
	path=pltpath
	IF keyword_set(prb) THEN path=prbpath
	IF keyword_set(prc) THEN path=prcpath

	ntemp=n(temp)+1
	ndens=n(dens)+1
	ivec=alog10(temp)
	jvec=alog10(dens)
	plc=dblarr(z+1,ntemp,ndens)
	out=read_pxx_file(path)
	pxxt=alog10(out.temp)
	pxxd=alog10(out.dens)
	ipts=interp_vec_reform(pxxt,ivec)
	jpts=interp_vec_reform(pxxd,jvec)
	FOR i=0,z DO plc[i,*,*]=interpolate(reform(out.plc[i,*,*]),ipts,jpts,/grid)
	
	output=plc
	IF keyword_set(debug) THEN stop
	RETURN,output

END
