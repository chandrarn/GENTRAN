PRO znlike_emiss,csden,temp,dens,wave,chg,emiss

	path='/home/mlreinke/idl/impurities/data/adas/W_44_adf15.dat'
	pec=read_pec_file(path)
	nlines=10			;take the first 10 which are really the brightest onces
	ntemp=n(temp)+1
	emiss=fltarr(nlines,ntemp)
	dindex=ipt(dens[0],pec.dens)	;assume density sensitivity is weak
	lte_pec=alog10(pec.temp)
	lte=alog10(temp)
	FOR i=0,nlines-1 DO emiss[i,*]=interpol(reform(pec.pec[i,*,dindex]),lte_pec,lte)*dens*csden[44,*]
	wave=pec.lam[0:nlines-1]
	chg=intarr(nlines)+44
END

FUNCTION mcp_pos,shot,jack=jack
	IF NOT keyword_set(jack) THEN BEGIN
		mdsopen, 'spectroscopy',(shot)
		jack=mdsvalue('\SPECTROSCOPY::TOP.VUV:JACK_POS')
	        mdsclose, 'spectroscopy',(shot)
	ENDIF

	ang=ACOS((72.375*72.375+87.5*87.5-(.9536*jack+27.457)*(.9536*jack+27.457))/(2.0*72.375*87.5))-.280392
	pivot=[3.738,.3315]
	;adjust z0 so that r0 = 1.1.  This will make it so that the default n_s will work
	zo=pivot[1]+tan(-1*ang)*(pivot[0]-1.1)

	output=[1.1,zo,0.0,ang]
	RETURN, output
END

PRO load_mcpline_data,shot,line,br,t,verb=verb

	mdsopen,'spectroscopy',(shot)
	id=mdsvalue('\spectroscopy::top.vuv.analysis.time_hist:line_ids')
	fws=mdsvalue('\spectroscopy::top.vuv.analysis.time_hist:fwhms')
	backs=mdsvalue('\spectroscopy::top.vuv.analysis.time_hist:back_wave')
	bri=mdsvalue('\spectroscopy::top.vuv.analysis.time_hist:b_v_time')
	t2=mdsvalue('dim_of(\spectroscopy::top.vuv.analysis.time_hist:b_v_time)')
	w2=mdsvalue('dim_of(\spectroscopy::top.vuv.analysis.time_hist:b_v_time,1)')
	id2=mdsvalue('dim_of(\spectroscopy::top.vuv.analysis.time_hist:b_v_time,2)')
	mdsclose, 'spectroscopy', (shot)

	;set values less than zero to zero
	badpts=where(bri lt -50)
	IF badpts[0] NE -1 THEN bri[badpts]=0

	IF keyword_set(verb) THEN BEGIN
		print,'Time histories for '+num2str(n(id2)+1)+' lines in the tree'
		print,'index     line ID       exp. wavelength '
		FOR i=0,n(id) DO print,i,id(i),w2(i)
	ENDIF

	ID_match=[97.17,114.98,116.0,118.213,121.985,127.81,134.880,132.85,129.871,112.8,99.12,232.0,234.2,$
		255.74,265.7,237.35,111.1,158.18,178.93,220.2,242.5]
	t=t2
	br_index=where(w2 EQ ID_match[line-1])
	IF br_index[0] NE -1 THEN BEGIN
		br=bri[br_index,*]
		t=t2
		IF keyword_set(verb) THEN print, id(br_index), w2(br_index)
	ENDIF ELSE BEGIN
		IF keyword_set(verb) THEN print, 'Line not taken for that shot'
		br=-1
		t=-1
	ENDELSE
END

PRO mcp_tungsten_density,shot,t1,t2,dff=dff,plot=plot,xtomo=xtomo
	IF n(t1) NE n(t2) THEN RETURN
	IF NOT keyword_set(dff) THEN dff=0.25

	;setup view
	mcp_pos=mcp_pos(shot)
	genpos_pos_reform,mcp_pos,[0.44,1.0,0.6,-0.6]

	ntime=n(t1)+1
	charge=findgen(74)

	FOR j=0,ntime-1 DO BEGIN
		efit=0
		rmid=0
		data=0

		gentran_profiles,shot,74,t1[j],t2[j],/pt,dff=dff,data=data,plot=plot
		print, 'csden profile calculated '+num2str(shot,1)+' '+num2str(t1[j],dp=2)+' < t < '+num2str(t2[j],dp=2)

		znlike_emiss,data.csden,data.temp,data.dens,wave,chg,iemiss
		print, 'emissivities calculated'

		nlines=n(wave)+1
		nrad=n(data.rho)+1
		IF j EQ 0 THEN BEGIN
			emiss=fltarr(nlines,nrad,ntime)
			bright=fltarr(nlines,ntime)
		ENDIF
		emiss[*,*,j]=iemiss

		FOR i=0,nlines-1 DO BEGIN
			bright[i,j]=line_br(mcp_pos,reform(emiss[i,*,j]),data.rmaj,[data.time],data.shot,data.time,plots=plot,efit=efit,rmid=rmid,verb=verb)
			IF keyword_set(plot) then stop
		ENDFOR
		print, 'line-integrated brightness calculated'
	ENDFOR
	br_znlike=reform(bright[1,*])
	eta_belike=1.05e19				;mcpherson calibration constant from Be-like molybdenum
		
	load_mcpline_data,shot,8,br_znlike_mcp,tau	;load Fe line brightness which is degenerate with 132 Ang Zn-like W line

	nz_znlike=fltarr(ntime)
	FOR j=0,ntime-1 DO BEGIN
		tmp=where(tau GE t1[j] AND tau LE t2[j]) 
		IF tmp[0] NE -1 THEN ave_br=mean(br_znlike_mcp[tmp]) ELSE ave_br=br_znlike_mcp[ipt(0.5*(t1[j]+t2[j]),tau)]
		nz_znlike[j]=eta_belike*ave_br/br_znlike[j]
	ENDFOR



	stop
END