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

PRO check_lines
	print, 'Lines currently on file'
	print, 'ID	NAME		WL (ang)'
	print, '-----------------------------------'
	print, '1	B V (2nd)	48.585 x 2'
	print, '2	Ti XII		114.98'
	print, '3	Mo XXXI		115.999'
	print, '4	Ti XVI		118.21'
	print, '5	Ti XIV		121.985'
	print, '6	Mo XXXII	127.81'
	print, '7       F VII		134.880'
	print, '8	Fe XXIII 	132.85'
	print, '9	O VI 		129.785'
	print, '10	F VII		112.8'
	print, '11	N VII (4th)	24.78 x 4'
	print, '12	Mo XXXI (2nd)	115.99 x 2'
	print, '13	Cu XVIII 	234.2'
	print, '14	Mo XXXII (2nd)  127.87 x 2'
	print, '15 	Fe XXIII (2nd)	132.85 x 2'
	print, '16	Cu XIX 		237.35'
	print, '17 	Cu XXIII	111.07'
	print, ''
	print, 'Use ID in load_mcpline_data,shot,ID,br,t'
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

PRO mcp_moly_density,shot,t1,t2,dff=dff,plot=plot,xtomo=xtomo
	IF n(t1) NE n(t2) THEN RETURN
	IF NOT keyword_set(dff) THEN dff=0.5

	;setup view
	mcp_pos=mcp_pos(shot)
	genpos_pos_reform,mcp_pos,[0.44,1.0,0.6,-0.6]

	lam0=90
	lam1=150
	nlam=3000
	dlam=0.27
	lam=make(lam0,lam1,nlam)
	ntime=n(t1)+1
	charge=findgen(42)

	FOR j=0,ntime-1 DO BEGIN
		efit=0
		rmid=0
		data=0

		gentran_profiles,shot,42,t1[j],t2[j],/pt,dff=dff,data=data,plot=plot
		print, 'csden profile calculated '+num2str(shot,1)+' '+num2str(t1[j],dp=2)+' < t < '+num2str(t2[j],dp=2)

		sxr_vuv_molybdenum,data.csden,data.temp,data.dens,wave,chg,iemiss
		print, 'emissivities calculated'

		nlines=n(wave)+1
		nrad=n(data.rho)+1
		IF j EQ 0 THEN BEGIN
			emiss=fltarr(nlines,nrad,ntime)
			bright=fltarr(nlines,ntime)
			spec=fltarr(nlam,ntime)
			zeff=fltarr(nrad,ntime)
		ENDIF
		emiss[*,*,j]=iemiss

		FOR i=0,nlines-1 DO BEGIN
			bright[i,j]=line_br(mcp_pos,reform(emiss[i,*,j]),data.rmaj,[data.time],data.shot,data.time,plots=plot,efit=efit,rmid=rmid,verb=verb)
			IF keyword_set(plot) then stop
		ENDFOR
		print, 'line-integrated brightness calculated'

		FOR i=0,nlines-1 DO spec[*,j]+=bright[i,j]/(dlam*sqrt(2.0*!pi))*exp(-(lam-wave[i])^2/(2*dlam^2))
		print, 'line-integrated spectra formed'
		FOR i=0,nrad-1 DO zeff[i,j]=total(data.csden[*,j]*charge^2/data.dens[i])

		IF keyword_set(xtomo) THEN xtomo_genrad_profiles,shot,t1[j],t2[j],data.csden,data.cserr,data.temp,data.terr,data.dens,data.derr,data.rho,$
			plotwin=plotwin,t=t,nz=nz_belike[0],zeff=zeff_bck,out=outx
	ENDFOR
	br_lilike=reform(bright[11,*])
	br_belike=reform(bright[3,*])
	
	eta_lilike=0.65e19
	eta_belike=1.05e19
	
	load_mcpline_data,shot,3,br_belike_mcp,tau
	load_mcpline_data,shot,6,br_lilike_mcp,tau

	nz_lilike=fltarr(ntime)
	FOR j=0,ntime-1 DO BEGIN
		tmp=where(tau GE t1[j] AND tau LE t2[j]) 
		IF tmp[0] NE -1 THEN ave_br=mean(br_lilike_mcp[tmp]) ELSE ave_br=br_lilike_mcp[ipt(0.5*(t1[j]+t2[j]),tau)]
		nz_lilike[j]=eta_lilike*ave_br/br_lilike[j]
	ENDFOR

	nz_belike=fltarr(ntime)
	FOR j=0,ntime-1 DO BEGIN
		tmp=where(tau GE t1[j] AND tau LE t2[j]) 
		IF tmp[0] NE -1 THEN ave_br=mean(br_belike_mcp[tmp]) ELSE ave_br=br_belike_mcp[ipt(0.5*(t1[j]+t2[j]),tau)]
		nz_belike[j]=eta_belike*ave_br/br_belike[j]
	ENDFOR


	stop
END
