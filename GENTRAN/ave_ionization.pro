;;;;;;;;;;;;;;
;TUNGSTEN
;;;;;;;;;;;;;;

FUNCTION w_avez,te,options
	
	;INPUT a vector of electron temperatures [eV]
	;OUTPUT a vector of the average ionizations states at those Te values

	;this data compiled using Post and Jensen data Atomic Data and Nuclear Data Tables Vol(20) pg 397
	;

	t_pts= [[100.0,200.0],$
		[200.1,2.0e3],$
		[2.001e3,2.0e4],$
		[2.001e4,1.0e5]]
	
	M=     [[ 1.555776e2,  9.075264e2,   2.344661e3,   3.050606e3,   1.995100e3,   5.250867e2],$
		[ 2.560687e1,  2.746452e1,   1.880356e1,   2.109782e1,  -1.779688e1,  -5.671642e1],$
		[-8.081505e1,  7.807107e2,  -1.793494e3,   1.856495e3,  -8.288591e2,   1.207945e2],$
		[ 6.270438e3, -1.827260e4,   2.132807e4,  -1.235512e4,   3.557784e3,  -4.078162e2]]

	avez=fltarr(n(te)+1)	;create	
	
	FOR i=0,n(te) DO BEGIN
		j=where(t_pts[0,*] LE round(te[i]) AND t_pts[1,*] GE round(te[i]))
		IF (j NE -1) THEN BEGIN
			Mvec=M[*,j]		
			FOR k=0,5 DO avez[i]+=(alog10(te[i]/1000.0))^k*Mvec[k]
		ENDIF ELSE BEGIN
			avez[i]=-2
		ENDELSE
	ENDFOR

	OUTPUT=(avez)

	RETURN, output		

END

;;;;;;;;;;;;;;
;MOLYBDENUM
;;;;;;;;;;;;;;

FUNCTION mo_avez,te,options
	
	;INPUT a vector of electron temperatures [eV]
	;OUTPUT a vector of the average ionizations states at those Te values

	;this data compiled using Post and Jensen data Atomic Data and Nuclear Data Tables Vol(20) pg 397
	;

	t_pts= [[60.0,200.0],$
		[200.1,2.0e3],$
		[2.001e3,2.0e4],$
		[2.001e4,1.0e5]]
	
	M=     [[ 4.992398e1,  8.506396e1,  -7.973177e1,  -3.870386e2,  -3.718835e2,  -1.120589e2],$
		[ 2.449023e1,  2.684258e1,   1.299348e1,  -4.145381e1,  -1.496935e2,  -1.385436e2],$
		[ 3.412437e1,  4.468081e0,  -1.301926e2,   3.624758e2,  -3.300305e2,   9.866466e1],$
		[ 1.367133e2, -2.287314e2,   2.005241e2,  -7.899527e1,   1.288693e1,  -4.473430e-1]]

	avez=fltarr(n(te)+1)	;create	
	
	FOR i=0,n(te) DO BEGIN
		j=where(t_pts[0,*] LE round(te[i]) AND t_pts[1,*] GE round(te[i]))
		IF (j NE -1) THEN BEGIN
			Mvec=M[*,j]		
			FOR k=0,5 DO avez[i]+=(alog10(te[i]/1000.0))^k*Mvec[k]
		ENDIF ELSE BEGIN
			avez[i]=-2
		ENDELSE
	ENDFOR

	OUTPUT=(avez)

	RETURN, output	

END

;;;;;;;;;;;;;;
;ARGON
;;;;;;;;;;;;;;

FUNCTION ar_avez,te,options
	
	;INPUT a vector of electron temperatures [eV]
	;OUTPUT a vector of the average ionizations states at those Te values

	;this data compiled using Post and Jensen data Atomic Data and Nuclear Data Tables Vol(20) pg 397
	;

	t_pts= [[30.0,200.0],$
		[200.1,2.0e3],$
		[2.001e3,2.0e4],$
		[2.001e4,1.0e5]]
	
	M=     [[-6.351877e1, -4.145188e2,  -8.502200e2,  -8.074783e2,  -3.621262e2,  -6.187830e1],$
		[ 1.591067e1, -7.886774e-1,  2.874539e0,   3.361189e1,  -3.306891e1,  -7.162801e1],$
		[ 1.296383e1,  1.833252e1,  -2.834795e1,   2.267165e1,  -9.219744e0,   1.507133e0],$
		[-7.890005e1,  2.939922e2,  -3.551084e2,   2.133685e2,  -6.376566e1,   7.582572e0]]

	avez=fltarr(n(te)+1)	;create	
	
	FOR i=0,n(te) DO BEGIN
		j=where(t_pts[0,*] LE round(te[i]) AND t_pts[1,*] GE round(te[i]))
		IF (j NE -1) THEN BEGIN
			Mvec=M[*,j]		
			FOR k=0,5 DO avez[i]+=(alog10(te[i]/1000.0))^k*Mvec[k]
		ENDIF ELSE BEGIN
			avez[i]=-2
		ENDELSE
	ENDFOR

	OUTPUT=(avez)

	RETURN, output	

END


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;MAIN PROGRAMS FOR USING AVERAGE IONIAZATION CURVES
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;


FUNCTION avez,te, z, q=quiet
	;automates the process of filling plc for multiple elements
	
	;INPUT te: m-element vector of electron temperatures in [eV]
	;INPUT z:  n-element vector of element atomic number 
	;OUTPUT plc: m x n array of average ionizations with z not on file set to 0, te out of range set to -1

	;OPTIONS
	;  /q = quiet mode, ERRORS still shown, but updates supressed

	zonfile=[18,42,74]
	
	;find # of,if any, matches to z on file
	cntr=0
	FOR i=0,n(zonfile) DO IF (where(z EQ zonfile[i]) NE -1) THEN cntr+=1
	IF (cntr EQ 0 ) THEN BEGIN
		print, 'ERROR: No z-values matched those on file'
		RETURN,-1
	ENDIF ELSE BEGIN
		IF (keyword_set(quiet)) THEN quiet=1  ELSE print, 'Matches found, filling '+num2str(cntr)+' of '+num2str(int(n(z)+1),1)+' z-values requested.'
		avez=fltarr(n(te)+1,n(z)+1)		;init plc array
		FOR i=0,n(z) DO BEGIN	
			;argon
			IF (z[i] EQ 18) THEN BEGIN
				ar=ar_avez(te)
				avez[*,i]=ar
			ENDIF
			;molybdenum
			IF (z[i] EQ 42) THEN BEGIN
				mo=mo_avez(te)
				avez[*,i]=mo
			ENDIF
			;tungsten
			IF (z[i] EQ 74) THEN BEGIN
				w=w_avez(te)
				avez[*,i]=w
			ENDIF
		
		ENDFOR
			
		output=avez
		RETURN, output
	ENDELSE
END


;;;;;;;;;;;;;;;;;;;;;;;;;;;
;GENERAL UTILITIES
;;;;;;;;;;;;;;;;;;;;;;;;;;;


PRO checkAVEZ
	;script to show which elements are on file
	z=findgen(92)+1
	te=100.0
	avez=avez(te,z, /q)
	
	print, 'For 1 < Z < 92 the following are on file'
	for i=0,n(z) DO IF (avez[0,i] NE 0) THEN print, '  '+num2str(i+1)+' '+num2elem_name(i+1)
END

FUNCTION chkAVEZ,z			;a binary version of checkZ using a function call
	te=1.0e3
	avez=avez(te,z, /q)
	
	IF (avez GT 0) THEN output=1 ELSE output=0

	RETURN,output
END


;;;;;;;;;;;;;;;;;;;;;;;;
;PLOTTING AVERAGE IONIZATION CURVES
;;;;;;;;;;;;;;;;;;;;;;;;

PRO avezplt,temp=te_r,zrange=z_r,zz=zval,max=maxpt,min=minpt, ps=plot2ps
	
	;setup variable space
	IF (keyword_set(te_r)) THEN te=findgen(5001)/5000.0*(te_r[1]-te_r[0])+te_r[0] ELSE te=findgen(5001)/5000.0*(1.0e4-10.0)+10.0
	IF (keyword_set(z_r)) THEN z=findgen(z_r[1]-z_r[0]+1)+z_r[0] ELSE z=findgen(92)+1
	IF (keyword_set(zval)) THEN z=zval

	;run plc to fill array
	avez=avez(te,z,/q)

	;fix out of range values to -1
	FOR i=0,n(z) DO BEGIN
		IF (avez[0,i] NE 0) THEN BEGIN
			bad=where(avez[*,i] LE -1.99)
			IF (bad[0] NE -1) THEN FOR j=0,n(bad) DO avez[bad[j],i]=-1
		ENDIF
	ENDFOR
	
	;define color map for the good channels
	n_good=n(where(avez[0,*] NE 0))+1
	colormap=intarr(n_good)
	n_colors=!d.table_size-30	;those last few colors suck
	FOR i=0,n_good-1 DO BEGIN
		colormap[i]=i*(n_colors)/(n_good)
	ENDFOR


	;setup ps plot, if chosen
	IF (keyword_set(plot2ps)) THEN BEGIN
		current_device=!d.name
		psplot
	ENDIF ELSE BEGIN
		wreset,10
	ENDELSE

	IF (keyword_set(maxpt)) THEN maxplot=maxpt ELSE maxplot=max(avez)
	IF (keyword_set(minpt)) THEN minplot=minpt ELSE minplot=0

	;define plot
	plot,[1],[1],/xlog,title='Average Ionization State for Various Impurities, Z', ytitle='<Z>', $
		xtitle='Te [eV]', yrange=[minplot,maxplot],xrange=[te[0], te[n(te)]],chars=1.3
	
	;run through array and plot valid entries
	cntr=0
	FOR i=0, n(z) DO BEGIN
		IF (avez[0,i] NE 0) THEN BEGIN
			oplot,te,avez[*,i],color=colormap[cntr]	
			xyouts,.16,.13+.03*cntr,/norm,'Z='+num2str(floor(z[i]))+' '+num2elem(floor(z[i])),color=colormap[cntr]
			cntr+=1
		ENDIF
	ENDFOR
		

	;close ps plot if chosen
	IF (keyword_set(plot2ps)) THEN BEGIN
		device,/close
		set_plot,current_device
		spawn, 'pwd',out_string
		IF (strmatch(out_string, '*mlreinke*') EQ 1) THEN BEGIN
			xwplot					
			spawn, 'cp idl.ps /home/mlreinke/idl/plots/avezplt.ps'
		ENDIF
	ENDIF
	
	
END



	
