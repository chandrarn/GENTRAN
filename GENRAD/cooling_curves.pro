;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;Cooling Curves for various Impurities
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

;you will need to compile mlr_functions before compiling this code


;;;;;;;;;;;;;;
;Tungsten
;;;;;;;;;;;;;;

FUNCTION w_plc,te,options
	
	;INPUT a vector of electron temperatures [eV]
	;OUTPUT a vector of power loss coeffiencts at those Te values [W m^3]
	;to compare to emissivty take ne*ni*plc [W/m^3]

	;this data compiled using Post and Jensen data Atomic Data and Nuclear Data Tables Vol(20) pg 397
	;

	;NOTE THESE COOLING CURVES WERE GENERATED USING A CORONAL EQUILIBRIUM CHARGE STATE DISTRIBUTION

	;Post data in erg* cm^3/s using 1 erg*cm^3/s = 1e-13 W*m^3 to transform

	t_pts= [[100.0,200.0],$
		[200.1,2.0e3],$
		[2.001e3,2.0e4],$
		[2.001e4,1.0e5]]
	
	M=     [[ 5.340828e0,  1.560876e2,   4.171704e2,   5.502576e2,   3.567583e2,   9.042786e1],$
		[-1.723894e1,  5.423752e-2, -1.221070e0,   4.411812e-1, -4.485821e0,  -7.836137e0],$
		[-1.474880e1, -1.439542e1,   2.105855e1,  -4.394746e0,  -1.106006e1,   5.616985e0],$
		[-2.624260e2,  7.125586e2,  -8.250168e2,   4.742407e2,  -1.355175e2,   1.541889e1]]

	plc=fltarr(n(te)+1)	;create	
	
	FOR i=0,n(te) DO BEGIN
		j=where(t_pts[0,*] LE round(te[i]) AND t_pts[1,*] GE round(te[i]))
		IF (j NE -1) THEN BEGIN
			Mvec=M[*,j]		
			FOR k=0,5 DO plc[i]+=(alog10(te[i]/1000.0))^k*Mvec[k]
		ENDIF ELSE BEGIN
			;print, num2str(te[i])+' out of range, please enter '+num2str(t_pts[0])+' < Te < '+num2str(t_pts[n(t_pts)])
			;print, 'Filling array value with 1'
			plc[i]=13
		ENDELSE
	ENDFOR

	OUTPUT=10^(plc)*1.0e-13

	RETURN, output		;in units of W*m^3

END

PRO wplt,te_r,max=maxpt,color=color, ps=plot2ps
	;INPUT [te_o, te_1] for the plotting interval
	;NOTE: set plotting window before running
	
	te=findgen(5001)/5000.0*(te_r[1]-te_r[0])+te_r[0]
	plc=w_plc(te)

	IF keyword_set(maxpt) THEN maxplot=maxpt ELSE maxplot=max(plc)
	IF keyword_set(color) THEN cplt=color ELSE cplt=1	
	
	IF (keyword_set(plot2ps)) THEN BEGIN
		current_device=!d.name
		set_plot,'ps'
	ENDIF	
	
	plot,[0],[0],/xlog,xtitle='Te [eV]',ytitle='PLC [W*m^3]', title='Tungsten Power Loss Coefficient - Post  ADaNDT 20(397)',chars=1.3,yrange=[0,maxplot],xrange=[te_r[0],te_r[1]]
	oplot,te,plc,color=cplt

	IF (keyword_set(plot2ps)) THEN BEGIN
		device,/close
		spawn, 'cp idl.ps /home/mlreinke/idl/plots/W_plc.ps'
		set_plot,current_device
		xplot
	ENDIF

END

;;;;;;;;;;;;;;
;Xenon
;;;;;;;;;;;;;;

FUNCTION xe_plc,te,options
	
	;INPUT a vector of electron temperatures [eV]
	;OUTPUT a vector of power loss coeffiencts at those Te values [W m^3]
	;to compare to emissivty take ne*ni*plc [W/m^3]

	;this data compiled using Post and Jensen data Atomic Data and Nuclear Data Tables Vol(20) pg 397
	;

	;NOTE THESE COOLING CURVES WERE GENERATED USING A CORONAL EQUILIBRIUM CHARGE STATE DISTRIBUTION

	;Post data in erg* cm^3/s using 1 erg*cm^3/s = 1e-13 W*m^3 to transform

	t_pts= [[80.0,200.0],$
		[200.1,2.0e3],$
		[2.001e3,2.0e4],$
		[2.001e4,1.0e5]]
	
	M=     [[-2.027267e1, -1.530175e1,  -3.074557e1,  -3.155124e1,  -1.600739e1,  -3.091098e0],$
		[-1.778249e1,  2.776323e-1,  1.901048e0,  -5.727093e0,  -6.456918e0,   2.998205e0],$
		[-2.445709e1,  5.504901e1,  -1.625266e2,   2.174682e2,  -1.363026e2,   3.242958e1],$
		[-3.693018e1,  6.802254e1,  -9.562685e1,   6.445815e1,  -2.107990e1,   2.700386e0]]

	plc=fltarr(n(te)+1)	;create	
	
	FOR i=0,n(te) DO BEGIN
		j=where(t_pts[0,*] LE round(te[i]) AND t_pts[1,*] GE round(te[i]))
		IF (j NE -1) THEN BEGIN
			Mvec=M[*,j]		
			FOR k=0,5 DO plc[i]+=(alog10(te[i]/1000.0))^k*Mvec[k]
		ENDIF ELSE BEGIN
			;print, num2str(te[i])+' out of range, please enter '+num2str(t_pts[0])+' < Te < '+num2str(t_pts[n(t_pts)])
			;print, 'Filling array value with 1'
			plc[i]=13
		ENDELSE
	ENDFOR

	OUTPUT=10^(plc)*1.0e-13

	RETURN, output		;in units of W*m^3

END

PRO xeplt,te_r,max=maxpt,color=color, ps=plot2ps
	;INPUT [te_o, te_1] for the plotting interval
	;NOTE: set plotting window before running
	
	te=findgen(5001)/5000.0*(te_r[1]-te_r[0])+te_r[0]
	plc=xe_plc(te)

	IF keyword_set(maxpt) THEN maxplot=maxpt ELSE maxplot=max(plc)
	IF keyword_set(color) THEN cplt=color ELSE cplt=1	
	
	IF (keyword_set(plot2ps)) THEN BEGIN
		current_device=!d.name
		set_plot,'ps'
	ENDIF	
	
	plot,[0],[0],/xlog,xtitle='Te [eV]',ytitle='PLC [W*m^3]', title='Xenon Power Loss Coefficient - Post  ADaNDT 20(397)',chars=1.3,yrange=[0,maxplot],xrange=[te_r[0],te_r[1]]
	oplot,te,plc,color=cplt

	IF (keyword_set(plot2ps)) THEN BEGIN
		device,/close
		spawn, 'cp idl.ps /home/mlreinke/idl/plots/Xe_plc.ps'
		set_plot,current_device
		xplot
	ENDIF

END

;;;;;;;;;;;;;;
;MOLYBDENUM
;;;;;;;;;;;;;;

FUNCTION mo_plc,te,options
	
	;INPUT a vector of electron temperatures [eV]
	;OUTPUT a vector of power loss coeffiencts at those Te values [W m^3]
	;to compare to emissivty take ne*ni*plc [W/m^3]

	;this data compiled using Fournier data from Nuclear Fusion Vol(38) pg 639
	;(see Vol(37) pg 825 for original paper on the methods used)

	;NOTE THESE COOLING CURVES WERE GENERATED USING A CORONAL EQUILIBRIUM CHARGE STATE DISTRIBUTION

	;Fournier data in erg* cm^3/s use 1 erg*cm^3/s = 1e-13 W*m^3 to transform

	t_pts= [[6.0,25.0],$
		[25.1,50.0],$
		[50.1,105.0],$
		[105.1,1.0e3],$
		[1.001e3,2.0e3],$
		[2.001e3,4.6e3],$
		[4.601e3,1.0e4],$
		[1.0001e4,9.0e4]]
	
	M=     [[4.9293e-18,-1.9057e-18, 2.7930e-19, -1.8465e-20, 5.8204e-22, -7.0863e-24],$
		[1.9549e-18, -2.3368e-19, 1.8209e-20, -5.6108e-22, 8.6529e-24,-5.4087e-26],$
		[-1.8936e-16, 1.2754e-17, -3.2845e-19, 4.1157e-21, -2.5272e-23, 6.1089e-26],$
		[2.9673e-19, 6.1339e-21, -1.4524e-23, 2.4232e-26, -2.469e-29, 9.7442e-33],$
		[3.5572e-18, -5.1609e-21, 4.5234e-24, -2.6200e-27, 9.1977e-31, -1.3673e-34],$
		[-2.1404e-20, 8.0170e-22, -2.0621e-25, -3.3560e-29, 1.6380e-32, -1.4314e-36],$
		[-1.1072e-18, 9.9453e-22, -2.4313e-25, 2.7879e-29, -1.5434e-33, 3.3294e-38],$
		[4.4309e-19, -2.0372e-23, 6.7661e-28, -1.1843e-32, 1.0234e-37, -3.4369e-43]]

	plc=fltarr(n(te)+1)	;create	
	
	FOR i=0,n(te) DO BEGIN
		j=where(t_pts[0,*] LE round(te[i]) AND t_pts[1,*] GE round(te[i]))
		IF (j NE -1) THEN BEGIN
			Mvec=M[*,j]		
			FOR k=0,5 DO plc[i]+=(te[i])^k*Mvec[k]
		ENDIF ELSE BEGIN
			;print, num2str(te[i])+' out of range, please enter '+num2str(t_pts[0])+' < Te < '+num2str(t_pts[n(t_pts)])
			;print, 'Filling array value at 1'
			plc[i]=1e13
		ENDELSE
	ENDFOR

	OUTPUT=plc*1.0e-13

	RETURN, output		;in units of W*m^3

END

PRO moplt,te_r, max=maxpt,color=color, ps=plot2ps
	;INPUT [te_o, te_1] for the plotting interval
	;NOTE: set plotting window before running
	
	te=findgen(5001)/5000.0*(te_r[1]-te_r[0])+te_r[0]
	plc=mo_plc(te)

	IF keyword_set(maxpt) THEN maxplot=maxpt ELSE maxplot=max(plc)
	IF keyword_set(color) THEN cplt=color ELSE cplt=1

	IF (keyword_set(plot2ps)) THEN BEGIN
		current_device=!d.name
		set_plot,'ps'
	ENDIF

	plot,[0],[0],/xlog,xtitle='Te [eV]',ytitle='PLC [W*m^3]', title='Molybdenm Power Loss Coefficient - Fournier NucFus 38(639)',chars=1.3,yrange=[0,maxplot],xrange=[te_r[0],te_r[1]]
	oplot,te,plc,color=cplt

	IF (keyword_set(plot2ps)) THEN BEGIN
		device,/close
		spawn, 'cp idl.ps /home/mlreinke/idl/plots/Mo_plc.ps'
		set_plot,current_device
	ENDIF
		
END


;;;;;;;;;;;;;;;;;;;;;
;MOLY - POST & Jensen
;;;;;;;;;;;;;;;;;;;;;

FUNCTION moPJ_plc,te,options
	
	;INPUT a vector of electron temperatures [eV]
	;OUTPUT a vector of power loss coeffiencts at those Te values [W m^3]
	;to compare to emissivty take ne*ni*plc [W/m^3]

	;this data compiled using Post and Jensen data Atomic Data and Nuclear Data Tables Vol(20) pg 397
	;

	;NOTE THESE COOLING CURVES WERE GENERATED USING A CORONAL EQUILIBRIUM CHARGE STATE DISTRIBUTION

	;Post data in erg* cm^3/s using 1 erg*cm^3/s = 1e-13 W*m^3 to transform

	t_pts= [[60.0,200.0],$
		[200.1,2.0e3],$
		[2.001e3,2.0e4],$
		[2.001e4,1.0e5]]
	
	M=     [[-1.391054e2, -6.493335e2,  -1.365838e3,  -1.406464e3,  -7.086213e2,  -1.400571e2],$
		[-1.772591e1, -1.058217e0,  -3.583172e0,   1.660089e0,   8.595372e0,   4.532909e0],$
		[-1.385096e1, -3.678452e1,   1.140587e2,  -1.635634e2,   1.076260e2,  -2.642488e1],$
		[ 3.992683e1, -1.757093e2,   2.074927e2,  -1.214589e2,   3.531804e1,  -4.083832e0]]

	plc=fltarr(n(te)+1)	;create	
	
	FOR i=0,n(te) DO BEGIN
		j=where(t_pts[0,*] LE round(te[i]) AND t_pts[1,*] GE round(te[i]))
		IF (j NE -1) THEN BEGIN
			Mvec=M[*,j]		
			FOR k=0,5 DO plc[i]+=(alog10(te[i]/1000.0))^k*Mvec[k]
		ENDIF ELSE BEGIN
			;print, num2str(te[i])+' out of range, please enter '+num2str(t_pts[0])+' < Te < '+num2str(t_pts[n(t_pts)])
			;print, 'Filling array value with 1'
			plc[i]=13
		ENDELSE
	ENDFOR

	OUTPUT=10^(plc)*1.0e-13

	RETURN, output		;in units of W*m^3

END

;;;;;;;;;;;;;;
;Copper
;;;;;;;;;;;;;;

FUNCTION cu_plc,te,options
	
	;INPUT a vector of electron temperatures [eV]
	;OUTPUT a vector of power loss coeffiencts at those Te values [W m^3]
	;to compare to emissivty take ne*ni*plc [W/m^3]

	;this data compiled using Post and Jensen data Atomic Data and Nuclear Data Tables Vol(20) pg 397
	;

	;NOTE THESE COOLING CURVES WERE GENERATED USING A CORONAL EQUILIBRIUM CHARGE STATE DISTRIBUTION

	;Post data in erg* cm^3/s using 1 erg*cm^3/s = 1e-13 W*m^3 to transform

	t_pts= [[30.0,200.0],$
		[200.1,2.0e3],$
		[2.001e3,2.0e4],$
		[2.001e4,1.0e5]]
	
	M=     [[-1.751314e1,  6.909930e0,   2.271756e1,   2.718317e1,   1.415977e1,   2.801611e0],$
		[-1.830902e1,  2.274986e-1, -1.275432e0,  -1.123179e1,  -6.848201e0,   4.333373e0],$
		[-1.695900e1, -7.326308e0,   2.601280e0,   1.201293e1,  -1.377229e1,   4.259765e0],$
		[-2.648274e1,  2.298584e1,  -2.815227e1,   1.675992e1,  -4.884827e0,   5.635291e-1]]

	plc=fltarr(n(te)+1)	;create	
	
	FOR i=0,n(te) DO BEGIN
		j=where(t_pts[0,*] LE round(te[i]) AND t_pts[1,*] GE round(te[i]))
		IF (j NE -1) THEN BEGIN
			Mvec=M[*,j]		
			FOR k=0,5 DO plc[i]+=(alog10(te[i]/1000.0))^k*Mvec[k]
		ENDIF ELSE BEGIN
			;print, num2str(te[i])+' out of range, please enter '+num2str(t_pts[0])+' < Te < '+num2str(t_pts[n(t_pts)])
			;print, 'Filling array value with 1'
			plc[i]=13
		ENDELSE
	ENDFOR

	OUTPUT=10^(plc)*1.0e-13

	RETURN, output		;in units of W*m^3

END

PRO cuplt,te_r,max=maxpt,color=color, ps=plot2ps
	;INPUT [te_o, te_1] for the plotting interval
	;NOTE: set plotting window before running
	
	te=findgen(5001)/5000.0*(te_r[1]-te_r[0])+te_r[0]
	plc=cu_plc(te)

	IF keyword_set(maxpt) THEN maxplot=maxpt ELSE maxplot=max(plc)
	IF keyword_set(color) THEN cplt=color ELSE cplt=1	
	
	IF (keyword_set(plot2ps)) THEN BEGIN
		current_device=!d.name
		set_plot,'ps'
	ENDIF	
	
	plot,[0],[0],/xlog,xtitle='Te [eV]',ytitle='PLC [W*m^3]', title='Copper Power Loss Coefficient - Post  ADaNDT 20(397)',chars=1.3,yrange=[0,maxplot],xrange=[te_r[0],te_r[1]]
	oplot,te,plc,color=cplt

	IF (keyword_set(plot2ps)) THEN BEGIN
		device,/close
		spawn, 'cp idl.ps /home/mlreinke/idl/plots/Cu_plc.ps'
		set_plot,current_device
		xplot
	ENDIF

END

;;;;;;;;;;;;;;
;Nickel
;;;;;;;;;;;;;;

FUNCTION ni_plc,te,options
	
	;INPUT a vector of electron temperatures [eV]
	;OUTPUT a vector of power loss coeffiencts at those Te values [W m^3]
	;to compare to emissivty take ne*ni*plc [W/m^3]

	;this data compiled using Post and Jensen data Atomic Data and Nuclear Data Tables Vol(20) pg 397
	;

	;NOTE THESE COOLING CURVES WERE GENERATED USING A CORONAL EQUILIBRIUM CHARGE STATE DISTRIBUTION

	;Post data in erg* cm^3/s using 1 erg*cm^3/s = 1e-13 W*m^3 to transform

	t_pts= [[30.0,200.0],$
		[200.1,2.0e3],$
		[2.001e3,2.0e4],$
		[2.001e4,1.0e5]]
	
	M=     [[-1.203248e1,  3.253908e1,   6.790773e1,   6.529924e1,   2.973465e1,   5.271279e0],$
		[-1.830482e1, -3.319243e-3, -3.332313e0,  -1.112798e1,   1.053073e-1,  9.448907e0],$
		[-1.697378e1,  -9.495470e0,   1.109362e1,   4.045904e-2, -6.521934e0,   2.654915e0],$
		[-2.864081e1,  2.999289e1,  -3.726082e1,   2.258060e1,  -6.716598e0,   7.911687e-1]]

	plc=fltarr(n(te)+1)	;create	
	
	FOR i=0,n(te) DO BEGIN
		j=where(t_pts[0,*] LE round(te[i]) AND t_pts[1,*] GE round(te[i]))
		IF (j NE -1) THEN BEGIN
			Mvec=M[*,j]		
			FOR k=0,5 DO plc[i]+=(alog10(te[i]/1000.0))^k*Mvec[k]
		ENDIF ELSE BEGIN
			;print, num2str(te[i])+' out of range, please enter '+num2str(t_pts[0])+' < Te < '+num2str(t_pts[n(t_pts)])
			;print, 'Filling array value with 1'
			plc[i]=13
		ENDELSE
	ENDFOR

	OUTPUT=10^(plc)*1.0e-13

	RETURN, output		;in units of W*m^3

END

PRO niplt,te_r,max=maxpt,color=color, ps=plot2ps
	;INPUT [te_o, te_1] for the plotting interval
	;NOTE: set plotting window before running
	
	te=findgen(5001)/5000.0*(te_r[1]-te_r[0])+te_r[0]
	plc=ni_plc(te)

	IF keyword_set(maxpt) THEN maxplot=maxpt ELSE maxplot=max(plc)
	IF keyword_set(color) THEN cplt=color ELSE cplt=1	
	
	IF (keyword_set(plot2ps)) THEN BEGIN
		current_device=!d.name
		set_plot,'ps'
	ENDIF	
	
	plot,[0],[0],/xlog,xtitle='Te [eV]',ytitle='PLC [W*m^3]', title='Nickel Power Loss Coefficient - Post  ADaNDT 20(397)',chars=1.3,yrange=[0,maxplot],xrange=[te_r[0],te_r[1]]
	oplot,te,plc,color=cplt

	IF (keyword_set(plot2ps)) THEN BEGIN
		device,/close
		spawn, 'cp idl.ps /home/mlreinke/idl/plots/Ni_plc.ps'
		set_plot,current_device
		xplot
	ENDIF

END

;;;;;;;;;;;;;;
;Iron
;;;;;;;;;;;;;;

FUNCTION fe_plc,te,options
	
	;INPUT a vector of electron temperatures [eV]
	;OUTPUT a vector of power loss coeffiencts at those Te values [W m^3]
	;to compare to emissivty take ne*ni*plc [W/m^3]

	;this data compiled using Post and Jensen data Atomic Data and Nuclear Data Tables Vol(20) pg 397
	;

	;NOTE THESE COOLING CURVES WERE GENERATED USING A CORONAL EQUILIBRIUM CHARGE STATE DISTRIBUTION

	;Post data in erg* cm^3/s using 1 erg*cm^3/s = 1e-13 W*m^3 to transform

	t_pts= [[20.0,200.0],$
		[200.1,2.0e3],$
		[2.001e3,2.0e4],$
		[2.001e4,1.0e5]]
	
	M=     [[-2.752599e1, -3.908228e1,  -6.469423e1,  -5.555048e1,  -2.405568e1,  -4.093160e0],$
		[-1.834973e1, -1.252028e0,  -7.533115e0,  -3.289693e0,   2.866739e1,   2.830249e1],$
		[-1.671042e1, -1.646143e1,   3.766238e1,  -3.944080e1,   1.918529e1,  -3.509238e0],$
		[-2.453957e1,  1.795222e1,  -2.356360e1,   1.484503e1,  -4.542323e0,  5.477462e-1]]

	plc=fltarr(n(te)+1)	;create	
	
	FOR i=0,n(te) DO BEGIN
		j=where(t_pts[0,*] LE round(te[i]) AND t_pts[1,*] GE round(te[i]))
		IF (j NE -1) THEN BEGIN
			Mvec=M[*,j]		
			FOR k=0,5 DO plc[i]+=(alog10(te[i]/1000.0))^k*Mvec[k]
		ENDIF ELSE BEGIN
			;print, num2str(te[i])+' out of range, please enter '+num2str(t_pts[0])+' < Te < '+num2str(t_pts[n(t_pts)])
			;print, 'Filling array value with 1'
			plc[i]=13
		ENDELSE
	ENDFOR

	OUTPUT=10^(plc)*1.0e-13

	RETURN, output		;in units of W*m^3

END

PRO feplt,te_r,max=maxpt,color=color, ps=plot2ps
	;INPUT [te_o, te_1] for the plotting interval
	;NOTE: set plotting window before running
	
	te=findgen(5001)/5000.0*(te_r[1]-te_r[0])+te_r[0]
	plc=fe_plc(te)

	IF keyword_set(maxpt) THEN maxplot=maxpt ELSE maxplot=max(plc)
	IF keyword_set(color) THEN cplt=color ELSE cplt=1	
	
	IF (keyword_set(plot2ps)) THEN BEGIN
		current_device=!d.name
		set_plot,'ps'
	ENDIF	
	
	plot,[0],[0],/xlog,xtitle='Te [eV]',ytitle='PLC [W*m^3]', title='Iron Power Loss Coefficient - Post  ADaNDT 20(397)',chars=1.3,yrange=[0,maxplot],xrange=[te_r[0],te_r[1]]
	oplot,te,plc,color=cplt

	IF (keyword_set(plot2ps)) THEN BEGIN
		device,/close
		spawn, 'cp idl.ps /home/mlreinke/idl/plots/Fe_plc.ps'
		set_plot,current_device
		xplot
	ENDIF

END


;;;;;;;;;;;;;;
;Titanium
;;;;;;;;;;;;;;

FUNCTION ti_plc,te,options
	
	;INPUT a vector of electron temperatures [eV]
	;OUTPUT a vector of power loss coeffiencts at those Te values [W m^3]
	;to compare to emissivty take ne*ni*plc [W/m^3]

	;this data compiled using Post and Jensen data Atomic Data and Nuclear Data Tables Vol(20) pg 397
	;

	;NOTE THESE COOLING CURVES WERE GENERATED USING A CORONAL EQUILIBRIUM CHARGE STATE DISTRIBUTION

	;Post data in erg* cm^3/s using 1 erg*cm^3/s = 1e-13 W*m^3 to transform

	t_pts= [[20.0,200.0],$
		[200.1,2.0e3],$
		[2.001e3,2.0e4],$
		[2.001e4,1.0e5]]
	
	M=     [[ 2.391331e1,  1.833595e2,   3.019617e2,   2.376019e2,   9.049792e1,   1.345090e1],$
		[-1.899097e1, -3.403261e0,   1.439830e0,   1.735567e1,   2.804832e-1, -1.943971e1],$
		[-1.929037e1, -3.260377e0,   1.454427e1,  -2.383997e1,   1.642804e1,  -4.084697e0],$
		[-1.341780e1, -1.675967e1,   1.843320e1,  -1.033234e1,   2.960530e0,  -3.423194e-1]]

	plc=fltarr(n(te)+1)	;create	
	
	FOR i=0,n(te) DO BEGIN
		j=where(t_pts[0,*] LE round(te[i]) AND t_pts[1,*] GE round(te[i]))
		IF (j NE -1) THEN BEGIN
			Mvec=M[*,j]		
			FOR k=0,5 DO plc[i]+=(alog10(te[i]/1000.0))^k*Mvec[k]
		ENDIF ELSE BEGIN
			;print, num2str(te[i])+' out of range, please enter '+num2str(t_pts[0])+' < Te < '+num2str(t_pts[n(t_pts)])
			;print, 'Filling array value with 1'
			plc[i]=13
		ENDELSE
	ENDFOR

	OUTPUT=10^(plc)*1.0e-13

	RETURN, output		;in units of W*m^3

END

PRO tiplt,te_r,max=maxpt,color=color, ps=plot2ps
	;INPUT [te_o, te_1] for the plotting interval
	;NOTE: set plotting window before running
	
	te=findgen(5001)/5000.0*(te_r[1]-te_r[0])+te_r[0]
	plc=fe_plc(te)

	IF keyword_set(maxpt) THEN maxplot=maxpt ELSE maxplot=max(plc)
	IF keyword_set(color) THEN cplt=color ELSE cplt=1	
	
	IF (keyword_set(plot2ps)) THEN BEGIN
		current_device=!d.name
		set_plot,'ps'
	ENDIF	
	
	plot,[0],[0],/xlog,xtitle='Te [eV]',ytitle='PLC [W*m^3]', title='Titanium Power Loss Coefficient - Post  ADaNDT 20(397)',chars=1.3,yrange=[0,maxplot],xrange=[te_r[0],te_r[1]]
	oplot,te,plc,color=cplt

	IF (keyword_set(plot2ps)) THEN BEGIN
		device,/close
		spawn, 'cp idl.ps /home/mlreinke/idl/plots/Ti_plc.ps'
		set_plot,current_device
		xplot
	ENDIF

END


;;;;;;;;;;;;;;
;ARGON
;;;;;;;;;;;;;;

FUNCTION ar_plc,te,options
	
	;INPUT a vector of electron temperatures [eV]
	;OUTPUT a vector of power loss coeffiencts at those Te values [W m^3]
	;to compare to emissivty take ne*ni*plc [W/m^3]

	;this data compiled using Fournier data from Atomic Data and Nuclear Data Tables Vol(70) pg 231.
	;2/3/09 updated from the erratum v74 pg 333

	;NOTE THESE COOLING CURVES WERE GENERATED USING A CORONAL EQUILIBRIUM CHARGE STATE DISTRIBUTION

	;Fournier data in erg* cm^3/s use 1 erg*cm^3/s = 1e-13 W*m^3 to transform

	t_pts= [[2.0,20.0],$
		[20.1,65.0],$
		[65.1,200.0],$
		[200.1,500.0],$
		[500.1,900.0],$
		[900.1,5.0e3],$
		[5.001e3,1.9e4]]

	
	M=     [[-2.1649e-20, 1.4578e-17, -3.6783e-15, 1.9365e-12, -1.2132e-10, 2.2273e-9],$
		[-4.0750e-17, 5.2641e-15, -2.3998e-13, 5.1080e-12, -5.1923e-11, 2.0402e-10],$
		[-9.4925e-20, 1.8178e-17, -4.2915e-16, 4.4330e-15, -1.9550e-14, 3.1178e-14],$
		[1.4444e-18, -3.3185e-17, 3.1617e-16, -1.2354e-15, 2.1150e-15, -1.3297e-15],$
		[2.2333e-18, -1.1327e-17, 2.3024e-17, -2.2213e-17, 9.5271e-18, -1.1878e-18],$
		[3.7985e-20,  9.1616e-21, -2.7730e-21,-2.1838e-21, 7.7637e-22, -6.7637e-23],$
		[6.2632e-20, -1.9337e-20, 2.8671e-21, -2.1029e-22, 7.5285e-24, -1.0456e-25]]

	plc=fltarr(n(te)+1)	;create	
	
	FOR i=0,n(te) DO BEGIN
		j=where(t_pts[0,*] LE round(te[i]) AND t_pts[1,*] GE round(te[i]))
		IF (j NE -1) THEN BEGIN
			Mvec=M[*,j]		
			FOR k=0,5 DO plc[i]+=(te[i]/1000.0)^k*Mvec[k]
		ENDIF ELSE BEGIN
			;print, num2str(te[i])+' out of range, please enter '+num2str(t_pts[0])+' < Te < '+num2str(t_pts[n(t_pts)])
			;print, 'Filling array value with 1'
			plc[i]=1e13
		ENDELSE
	ENDFOR

	OUTPUT=plc*1.0e-13

	RETURN, output		;in units of W*m^3

END

PRO arplt,te_r,max=maxpt,color=color,ps=plot2ps
	;INPUT [te_o, te_1] for the plotting interval
	;NOTE: set plotting window before running
	
	te=findgen(5001)/5000.0*(te_r[1]-te_r[0])+te_r[0]
	plc=ar_plc(te)

	IF keyword_set(maxpt) THEN maxplot=maxpt ELSE maxplot=max(plc)
	IF keyword_set(color) THEN cplt=color ELSE cplt=1

	IF (keyword_set(plot2ps)) THEN BEGIN
		current_device=!d.name
		set_plot,'ps'
	ENDIF	

	plot,[0],[0],/xlog,xtitle='Te [eV]',ytitle='PLC [W*m^3]', title='Argon Power Loss Coefficient - Fournier ADaNDT 70(231)',chars=1.3,yrange=[0,maxplot],xrange=[te_r[0],te_r[1]]
	oplot,te,plc,color=cplt

	IF (keyword_set(plot2ps)) THEN BEGIN
		device,/close
		spawn, 'cp idl.ps /home/mlreinke/idl/plots/Ar_plc.ps'
		set_plot,current_device
	ENDIF

END

;;;;;;;;;;;;;;
;KRYPTON
;;;;;;;;;;;;;;

FUNCTION kr_plc,te,options
	
	;INPUT a vector of electron temperatures [eV]
	;OUTPUT a vector of power loss coeffiencts at those Te values [W m^3]
	;to compare to emissivty take ne*ni*plc [W/m^3]

	;this data compiled using Fournier data from Nuclear Fusion Vol(40) pg 847
	;

	;NOTE THESE COOLING CURVES WERE GENERATED USING A CORONAL EQUILIBRIUM CHARGE STATE DISTRIBUTION

	;Fournier data in erg* cm^3/s use 1 erg*cm^3/s = 1e-13 W*m^3 to transform

	t_pts= [[2.0,20.0],$
		[20.1,60.0],$
		[60.1,1.0e3],$
		[1.001e3,4.0e3],$
		[4.001e3,1.0e4],$
		[1.001e4,9.8e4]]
	
	M=     [[-7.161e-20, 7.5601e-17, -1.4173e-14, 3.1602e-12, -1.5418e-10, 1.4161e-9],$
		[2.2367e-17, -2.7255e-15, 1.3254e-13, -3.1678e-12, 3.7230e-11,-1.7206e-10],$
		[-3.5748e-19, 1.4953e-17, -3.3640e-17, 2.6368e-17, -7.3632e-18, 7.3027e-19],$
		[1.7893e-18, -2.6622e-18, 2.2494e-18, -8.7655e-19, 1.5437e-19, -1.0062e-20],$
		[2.8065e-18, -1.6507e-18, 4.1501e-19, -5.2627e-20, 3.3394e-21, -8.4462e-23],$
		[1.3618e-19, -1.6018e-21, -8.8831e-24, 7.9631e-25, -9.4893e-27, 3.5754e-29]]

	plc=fltarr(n(te)+1)	;create	
	
	FOR i=0,n(te) DO BEGIN
		j=where(t_pts[0,*] LE round(te[i]) AND t_pts[1,*] GE round(te[i]))
		IF (j NE -1) THEN BEGIN
			Mvec=M[*,j]		
			FOR k=0,5 DO plc[i]+=(te[i]/1000.0)^k*Mvec[k]
		ENDIF ELSE BEGIN
			;print, num2str(te[i])+' out of range, please enter '+num2str(t_pts[0])+' < Te < '+num2str(t_pts[n(t_pts)])
			;print, 'Filling array value with 1'
			plc[i]=1e13
		ENDELSE
	ENDFOR

	OUTPUT=plc*1.0e-13

	RETURN, output		;in units of W*m^3

END

PRO krplt,te_r,max=maxpt,color=color, ps=plot2ps
	;INPUT [te_o, te_1] for the plotting interval
	;NOTE: set plotting window before running
	
	te=findgen(5001)/5000.0*(te_r[1]-te_r[0])+te_r[0]
	plc=kr_plc(te)

	IF keyword_set(maxpt) THEN maxplot=maxpt ELSE maxplot=max(plc)
	IF keyword_set(color) THEN cplt=color ELSE cplt=1	
	
	IF (keyword_set(plot2ps)) THEN BEGIN
		current_device=!d.name
		set_plot,'ps'
	ENDIF	
	
	plot,[0],[0],/xlog,xtitle='Te [eV]',ytitle='PLC [W*m^3]', title='Krypton Power Loss Coefficient - Fournier NucFus 40(847)',chars=1.3,yrange=[0,maxplot],xrange=[te_r[0],te_r[1]]
	oplot,te,plc,color=cplt

	IF (keyword_set(plot2ps)) THEN BEGIN
		device,/close
		spawn, 'cp idl.ps /home/mlreinke/idl/plots/Kr_plc.ps'
		set_plot,current_device
		xplot
	ENDIF

END

;;;;;;;;;;;;;;
;NEON
;;;;;;;;;;;;;;

FUNCTION ne_plc,te,options
	
	;INPUT a vector of electron temperatures [eV]
	;OUTPUT a vector of power loss coeffiencts at those Te values [W m^3]
	;to compare to emissivty take ne*ni*plc [W/m^3]

	;this data compiled using Post and Jensen data Atomic Data and Nuclear Data Tables Vol(20) pg 397
	;

	;NOTE THESE COOLING CURVES WERE GENERATED USING A CORONAL EQUILIBRIUM CHARGE STATE DISTRIBUTION

	;Post data in erg* cm^3/s using 1 erg*cm^3/s = 1e-13 W*m^3 to transform

	t_pts= [[8.0,20.0],$
		[20.1,200.0],$
		[200.1,2.0e3],$
		[2.001e3,2.0e4],$
		[2.001e4,1.0e5]]
	
	M=     [[ 6.767667e1,  1.613787e2,   1.104767e2,   3.198749e1,   3.143149e0,   0.000000e0],$
		[-6.805016e1, -2.580876e2,  -5.346855e2,  -5.284263e2,  -2.476333e2,  -4.425504e1],$
		[-2.022685e1, -1.141439e0,   1.428329e0,   1.089571e0,  -6.969004e0,  -3.701472e0],$
		[-2.021988e1, -1.100588e0,   8.377123e-1,  1.707851e-1, -3.522400e-1,  1.041226e-1],$
		[-2.064906e1,  8.835572e-2, -2.864026e-1,   4.440143e-1, -1.881675e-1,  2.719145e-2]]

	plc=fltarr(n(te)+1)	;create	
	
	FOR i=0,n(te) DO BEGIN
		j=where(t_pts[0,*] LE round(te[i]) AND t_pts[1,*] GE round(te[i]))
		IF (j NE -1) THEN BEGIN
			Mvec=M[*,j]		
			FOR k=0,5 DO plc[i]+=(alog10(te[i]/1000.0))^k*Mvec[k]
		ENDIF ELSE BEGIN
			;print, num2str(te[i])+' out of range, please enter '+num2str(t_pts[0])+' < Te < '+num2str(t_pts[n(t_pts)])
			;print, 'Filling array value with 1'
			plc[i]=13
		ENDELSE
	ENDFOR

	OUTPUT=10^(plc)*1.0e-13

	RETURN, output		;in units of W*m^3

END

PRO neplt,te_r,max=maxpt,color=color, ps=plot2ps
	;INPUT [te_o, te_1] for the plotting interval
	;NOTE: set plotting window before running
	
	te=findgen(5001)/5000.0*(te_r[1]-te_r[0])+te_r[0]
	plc=ne_plc(te)

	IF keyword_set(maxpt) THEN maxplot=maxpt ELSE maxplot=max(plc)
	IF keyword_set(color) THEN cplt=color ELSE cplt=1	
	
	IF (keyword_set(plot2ps)) THEN BEGIN
		current_device=!d.name
		set_plot,'ps'
	ENDIF	
	
	plot,[0],[0],/xlog,xtitle='Te [eV]',ytitle='PLC [W*m^3]', title='Neon Power Loss Coefficient - Post  ADaNDT 20(397)',chars=1.3,yrange=[0,maxplot],xrange=[te_r[0],te_r[1]]
	oplot,te,plc,color=cplt

	IF (keyword_set(plot2ps)) THEN BEGIN
		device,/close
		spawn, 'cp idl.ps /home/mlreinke/idl/plots/Ne_plc.ps'
		set_plot,current_device
		xplot
	ENDIF

END

;;;;;;;;;;;;;;
;FLOURINE
;;;;;;;;;;;;;;

FUNCTION f_plc,te,options
	
	;INPUT a vector of electron temperatures [eV]
	;OUTPUT a vector of power loss coeffiencts at those Te values [W m^3]
	;to compare to emissivty take ne*ni*plc [W/m^3]

	;this data compiled using Post and Jensen data Atomic Data and Nuclear Data Tables Vol(20) pg 397
	;

	;NOTE THESE COOLING CURVES WERE GENERATED USING A CORONAL EQUILIBRIUM CHARGE STATE DISTRIBUTION

	;Post data in erg* cm^3/s using 1 erg*cm^3/s = 1e-13 W*m^3 to transform

	t_pts= [[6.0,20.0],$
		[20.1,200.0],$
		[200.1,2.0e3],$
		[2.001e3,2.0e4],$
		[2.001e4,1.0e5]]
	
	M=     [[ 9.771138e2,  2.608148e3,  2.730338e3,   1.426168e3,   3.715237e2,   3.861789e1],$
		[-5.855138e1, -1.994465e2, -3.879919e2,  -3.528638e2,  -1.5008052e2, -2.412274e1],$
		[-2.044910e1, -9.224362e-1, 8.079891e-1, -9.429368e-1,  9.734815e-1,  6.783201e0],$
		[-2.043091e1, -1.093509e0,  1.314363e0,  -5.879604e-1,  1.307366e-1, -9.586941e-3],$
		[-2.085060e1,  2.631718e-1,-4.309072e-1,  5.229414e-1, -2.143155e-1,  3.089496e-2]]

	plc=fltarr(n(te)+1)	;create	
	
	FOR i=0,n(te) DO BEGIN
		j=where(t_pts[0,*] LE round(te[i]) AND t_pts[1,*] GE round(te[i]))
		IF (j NE -1) THEN BEGIN
			Mvec=M[*,j]		
			FOR k=0,5 DO plc[i]+=(alog10(te[i]/1000.0))^k*Mvec[k]
		ENDIF ELSE BEGIN
			;print, num2str(te[i])+' out of range, please enter '+num2str(t_pts[0])+' < Te < '+num2str(t_pts[n(t_pts)])
			;print, 'Filling array value with 1'
			plc[i]=13
		ENDELSE
	ENDFOR

	OUTPUT=10^(plc)*1.0e-13

	RETURN, output		;in units of W*m^3

END

PRO fplt,te_r,max=maxpt,color=color, ps=plot2ps
	;INPUT [te_o, te_1] for the plotting interval
	;NOTE: set plotting window before running
	
	te=findgen(5001)/5000.0*(te_r[1]-te_r[0])+te_r[0]
	plc=f_plc(te)

	IF keyword_set(maxpt) THEN maxplot=maxpt ELSE maxplot=max(plc)
	IF keyword_set(color) THEN cplt=color ELSE cplt=1	
	
	IF (keyword_set(plot2ps)) THEN BEGIN
		current_device=!d.name
		set_plot,'ps'
	ENDIF	
	
	plot,[0],[0],/xlog,xtitle='Te [eV]',ytitle='PLC [W*m^3]', title='Flourine Power Loss Coefficient - Post  ADaNDT 20(397)',chars=1.3,yrange=[0,maxplot],xrange=[te_r[0],te_r[1]]
	oplot,te,plc,color=cplt

	IF (keyword_set(plot2ps)) THEN BEGIN
		device,/close
		spawn, 'cp idl.ps /home/mlreinke/idl/plots/F_plc.ps'
		set_plot,current_device
		xplot
	ENDIF

END

;;;;;;;;;;;;;
;NITROGEN
;;;;;;;;;;;;;;

FUNCTION n_plc,te,options
	
	;INPUT a vector of electron temperatures [eV]
	;OUTPUT a vector of power loss coeffiencts at those Te values [W m^3]
	;to compare to emissivty take ne*ni*plc [W/m^3]

	;this data compiled using Post and Jensen data Atomic Data and Nuclear Data Tables Vol(20) pg 397
	;

	;NOTE THESE COOLING CURVES WERE GENERATED USING A CORONAL EQUILIBRIUM CHARGE STATE DISTRIBUTION

	;Post data in erg* cm^3/s using 1 erg*cm^3/s = 1e-13 W*m^3 to transform

	t_pts= [[4.0,20.0],$
		[20.1,200.0],$
		[200.1,2.0e3],$
		[2.001e3,2.0e4],$
		[2.001e4,1.0e5]]
	
	M=     [[-1.967182e2, -2.429049e2, -7.454123e1,   3.126366e1,   2.166881e1,   3.300054e0],$
		[-3.615155e1, -3.943802e1, -5.564129e0,   5.140343e1,   4.369243e1,   1.027448e1],$
		[-2.093912e1, -5.677397e-1, 7.664689e-1, -2.610450e-1,  3.464473e-1,  6.723385e-1],$
		[-2.093039e1, -6.617905e-1, 1.146777e0,  -7.664689e-1,  3.042676e-1,  -6.024562e-2],$
		[-9.452522e0, -3.583144e1,  4.386446e1,  -2.639331e1,   7.890268e0,  -9.366682e-1]]

	plc=fltarr(n(te)+1)	;create	
	
	FOR i=0,n(te) DO BEGIN
		j=where(t_pts[0,*] LE round(te[i]) AND t_pts[1,*] GE round(te[i]))
		IF (j NE -1) THEN BEGIN
			Mvec=M[*,j]		
			FOR k=0,5 DO plc[i]+=(alog10(te[i]/1000.0))^k*Mvec[k]
		ENDIF ELSE BEGIN
			;print, num2str(te[i])+' out of range, please enter '+num2str(t_pts[0])+' < Te < '+num2str(t_pts[n(t_pts)])
			;print, 'Filling array value with 1'
			plc[i]=13
		ENDELSE
	ENDFOR

	OUTPUT=10^(plc)*1.0e-13

	RETURN, output		;in units of W*m^3

END

PRO nplt,te_r,max=maxpt,color=color, ps=plot2ps
	;INPUT [te_o, te_1] for the plotting interval
	;NOTE: set plotting window before running
	
	te=findgen(5001)/5000.0*(te_r[1]-te_r[0])+te_r[0]
	plc=n_plc(te)

	IF keyword_set(maxpt) THEN maxplot=maxpt ELSE maxplot=max(plc)
	IF keyword_set(color) THEN cplt=color ELSE cplt=1	
	
	IF (keyword_set(plot2ps)) THEN BEGIN
		current_device=!d.name
		set_plot,'ps'
	ENDIF	
	
	plot,[0],[0],/xlog,xtitle='Te [eV]',ytitle='PLC [W*m^3]', title='Nitrogen Power Loss Coefficient - Post  ADaNDT 20(397)',chars=1.3,yrange=[0,maxplot],xrange=[te_r[0],te_r[1]]
	oplot,te,plc,color=cplt

	IF (keyword_set(plot2ps)) THEN BEGIN
		device,/close
		spawn, 'cp idl.ps /home/mlreinke/idl/plots/N_plc.ps'
		set_plot,current_device
		xplot
	ENDIF

END


;;;;;;;;;;;;;
;CARBON
;;;;;;;;;;;;;;

FUNCTION c_plc,te,options
	
	;INPUT a vector of electron temperatures [eV]
	;OUTPUT a vector of power loss coeffiencts at those Te values [W m^3]
	;to compare to emissivty take ne*ni*plc [W/m^3]

	;this data compiled using Post and Jensen data Atomic Data and Nuclear Data Tables Vol(20) pg 397
	;

	;NOTE THESE COOLING CURVES WERE GENERATED USING A CORONAL EQUILIBRIUM CHARGE STATE DISTRIBUTION

	;Post data in erg* cm^3/s using 1 erg*cm^3/s = 1e-13 W*m^3 to transform

	t_pts= [[3.0,20.0],$
		[20.1,200.0],$
		[200.1,2.0e3],$
		[2.001e3,2.0e4],$
		[2.001e4,1.0e5]]
	
	M=     [[ 1.965300e3,  4.572039e3,  4.159590e3,   1.871560e3,   4.173889e2,   3.699382e1],$
		[ 7.467599e1,  4.549038e2,  8.372937e2,   7.402515e2,   3.147607e2,   5.164578e1],$
		[-2.120151e1, -3.668933e-1, 7.295099e-1, -1.944827e-1, -1.263576e-1, -1.491027e-1],$
		[-2.121979e1, -2.346986e-1, 4.093794e-1,  7.874548e-2, -1.841379e-1,  5.590744e-2],$
		[-2.476796e1,  9.408181e0, -9.657446e0,   4.999161e0,  -1.237382e0,   1.160610e-1]]

	plc=fltarr(n(te)+1)	;create	
	
	FOR i=0,n(te) DO BEGIN
		j=where(t_pts[0,*] LE round(te[i]) AND t_pts[1,*] GE round(te[i]))
		IF (j NE -1) THEN BEGIN
			Mvec=M[*,j]		
			FOR k=0,5 DO plc[i]+=(alog10(te[i]/1000.0))^k*Mvec[k]
		ENDIF ELSE BEGIN
			;print, num2str(te[i])+' out of range, please enter '+num2str(t_pts[0])+' < Te < '+num2str(t_pts[n(t_pts)])
			;print, 'Filling array value with 1'
			plc[i]=13
		ENDELSE
	ENDFOR

	OUTPUT=10^(plc)*1.0e-13

	RETURN, output		;in units of W*m^3

END

PRO cplt,te_r,max=maxpt,color=color, ps=plot2ps
	;INPUT [te_o, te_1] for the plotting interval
	;NOTE: set plotting window before running
	
	te=findgen(5001)/5000.0*(te_r[1]-te_r[0])+te_r[0]
	plc=c_plc(te)

	IF keyword_set(maxpt) THEN maxplot=maxpt ELSE maxplot=max(plc)
	IF keyword_set(color) THEN cplt=color ELSE cplt=1	
	
	IF (keyword_set(plot2ps)) THEN BEGIN
		current_device=!d.name
		set_plot,'ps'
	ENDIF	
	
	plot,[0],[0],/xlog,xtitle='Te [eV]',ytitle='PLC [W*m^3]', title='Carbon Power Loss Coefficient - Post  ADaNDT 20(397)',chars=1.3,yrange=[0,maxplot],xrange=[te_r[0],te_r[1]]
	oplot,te,plc,color=cplt

	IF (keyword_set(plot2ps)) THEN BEGIN
		device,/close
		spawn, 'cp idl.ps /home/mlreinke/idl/plots/C_plc.ps'
		set_plot,current_device
		xplot
	ENDIF

END


;;;;;;;;;;;;;;
;BORON
;;;;;;;;;;;;;;

FUNCTION b_plc,te,options
	
	;INPUT a vector of electron temperatures [eV]
	;OUTPUT a vector of power loss coeffiencts at those Te values [W m^3]
	;to compare to emissivty take ne*ni*plc [W/m^3]

	;this data compiled using Post and Jensen data Atomic Data and Nuclear Data Tables Vol(20) pg 397
	;

	;NOTE THESE COOLING CURVES WERE GENERATED USING A CORONAL EQUILIBRIUM CHARGE STATE DISTRIBUTION

	;Post data in erg* cm^3/s using 1 erg*cm^3/s = 1e-13 W*m^3 to transform

	t_pts= [[2.0,20.0],$
		[20.1,200.0],$
		[200.1,2.0e3],$
		[2.001e3,2.0e4],$
		[2.001e4,1.0e5]]
	
	M=     [[-1.508695e3, -3.512267e3, -3.286126e3,  -1.520070e3,  -3.470698e2,  -3.127689e1],$
		[-6.370160e1, -2.156758e2, -4.308101e2,  -4.222842e2,  -2.008412e2,  -3.687482e1],$
		[-2.147874e1, -1.565300e-1, 6.181287e-1, -2.477378e-1, -1.060488e-1, -4.537644e-2],$
		[-2.147337e1, -1.829426e-1, 6.678447e-1, -3.864809e-1,  1.165920e-1, -1.400226e-2],$
		[-2.446008e1,  9.264960e0, -1.110684e1,  -6.837004e1,  -2.065106e0,   2.458526e-1]]

	plc=fltarr(n(te)+1)	;create	
	
	FOR i=0,n(te) DO BEGIN
		j=where(t_pts[0,*] LE round(te[i]) AND t_pts[1,*] GE round(te[i]))
		IF (j NE -1) THEN BEGIN
			Mvec=M[*,j]		
			FOR k=0,5 DO plc[i]+=(alog10(te[i]/1000.0))^k*Mvec[k]
		ENDIF ELSE BEGIN
			;print, num2str(te[i])+' out of range, please enter '+num2str(t_pts[0])+' < Te < '+num2str(t_pts[n(t_pts)])
			;print, 'Filling array value with 1'
			plc[i]=13
		ENDELSE
	ENDFOR

	OUTPUT=10^(plc)*1.0e-13

	RETURN, output		;in units of W*m^3

END

PRO bplt,te_r,max=maxpt,color=color, ps=plot2ps
	;INPUT [te_o, te_1] for the plotting interval
	;NOTE: set plotting window before running
	
	te=findgen(5001)/5000.0*(te_r[1]-te_r[0])+te_r[0]
	plc=b_plc(te)

	IF keyword_set(maxpt) THEN maxplot=maxpt ELSE maxplot=max(plc)
	IF keyword_set(color) THEN cplt=color ELSE cplt=1	
	
	IF (keyword_set(plot2ps)) THEN BEGIN
		current_device=!d.name
		set_plot,'ps'
	ENDIF	
	
	plot,[0],[0],/xlog,xtitle='Te [eV]',ytitle='PLC [W*m^3]', title='Boron Power Loss Coefficient - Post  ADaNDT 20(397)',chars=1.3,yrange=[0,maxplot],xrange=[te_r[0],te_r[1]]
	oplot,te,plc,color=cplt

	IF (keyword_set(plot2ps)) THEN BEGIN
		device,/close
		spawn, 'cp idl.ps /home/mlreinke/idl/plots/B_plc.ps'
		set_plot,current_device
		xplot
	ENDIF

END




;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;MAIN PROGRAMS FOR USING COOLING CURVES
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

;+
;NAME:
;	PLC
;
;PURPOSE:
;	This function calculates the radiative plasma loss coefficient [W*m^3] for given element at
;	a given temperature.
;
;CALLING SEQUENCE:
;	result=PLC(te,z)
;
;INPUTS:
;	te:	1D FLTARR of electron temperatures [eV] of length - M
;	Z:	1D INTARR of element numbers of lenght - N
;
;KEYWORD_PARAMETERS:
;	q:	/q suppresses non-error updates to terminal
;	
;OUTPUTS:
;	result:	1D FLTARR (M x N) of plasma loss coefficients [W*m^3]
;		 - z not on file, values are filled with 0
;		 - te out of range, values are filled with -1
;
;RESTRICTIONS:
;	Although 1 < Z < 92 are available in Post & Jensen, not all have been input into the table
;	due to the large amount of spare time that would entail.  Type CHECKZ to view the elements that are
;	on file and if you're drink isn't on the menu ask the bartender (mlreinke@mit.edu) to make it for you.
;
;PROCEDURE:
;	Values are taken from Post & Jensen, Atomic Data and Nuclear Data Tables Vol(20) pg 397 for
;	all data except for Ar, Kr and Mo which are from Fournier's papers (see /home/mlreinke/idl/impurities/data/table_ref)
;	for full references.  To get the P&J curves for where Fournier's available use Z=-Z in the input.
;
;RESTRICTIONS:
;	This needs MLR_FUNCTIONS to use properly.  Compile using @/home/mlreinke/idl/impurities/cc_ini.bat
;
;MODIFICATION HISTORY:
;	Written by:	ML Reinke, Spring 2005
;	9-28-05:	ML Reinke - adding proper documentation to work with HELPME
;	10-5-05:	ML Reinke - added Mo P&J curve
;
;-


FUNCTION plc,te, z, q=quiet
	;automates the process of filling plc for multiple elements
	
	;INPUT te: m-element vector of electron temperatures in [eV]
	;INPUT z:  n-element vector of element atomic number 
	;OUTPUT plc: m x n array of plasma loss coefficents [W*m^3] with z not on file set to 0, te out of range set to -1

	;OPTIONS
	;  /q = quiet mode, ERRORS still shown, but updates supressed

	zonfile=[5,6,7,9,10,18,22,26,29,36,42,54,74,-42]
	
	;find # of,if any, matches to z on file
	cntr=0
	FOR i=0,n(zonfile) DO IF (where(z EQ zonfile[i]) NE -1) THEN cntr+=1
	IF (cntr EQ 0 ) THEN BEGIN
		print, 'ERROR: No z-values matched those on file'
		RETURN,-1
	ENDIF ELSE BEGIN
		IF (keyword_set(quiet)) THEN quiet=1  ELSE print, 'Matches found, filling '+num2str(cntr)+' of '+num2str(int(n(z)+1),1)+' z-values requested.'
		plc=dblarr(n(te)+1,n(z)+1)		;init plc array
		FOR i=0,n(z) DO BEGIN
			;boron
			IF (z[i] EQ 5) THEN BEGIN
				b=b_plc(te)
				plc[*,i]=b
			ENDIF
				
			;carbon
			IF (z[i] EQ 6) THEN BEGIN
				c=c_plc(te)
				plc[*,i]=c
			ENDIF

			;nitrogen
			IF (z[i] EQ 7) THEN BEGIN
				n=n_plc(te)
				plc[*,i]=n
			ENDIF

			;flourine
			IF (z[i] EQ 9) THEN BEGIN
				f=f_plc(te)
				plc[*,i]=f
			ENDIF

			;neon
			IF (z[i] EQ 10) THEN BEGIN
				neon=ne_plc(te)
				plc[*,i]=neon
			ENDIF

			;argon
			IF (z[i] EQ 18) THEN BEGIN
				ar=ar_plc(te)
				plc[*,i]=ar
			ENDIF

			;titanium
			IF (z[i] EQ 22) THEN BEGIN
				ti=ti_plc(te)
				plc[*,i]=ti
			ENDIF
		
			;iron
			IF (z[i] EQ 26) THEN BEGIN
				fe=fe_plc(te)
				plc[*,i]=fe
			ENDIF

			;nickel
			IF (z[i] EQ 28) THEN BEGIN
				ni=ni_plc(te)
				plc[*,i]=ni
			ENDIF
		
			;copper
			IF (z[i] EQ 29) THEN BEGIN
				cu=cu_plc(te)
				plc[*,i]=cu
			ENDIF
	
			;krypton
			IF (z[i] EQ 36) THEN BEGIN
				kr=kr_plc(te)
				plc[*,i]=kr
			ENDIF
		
			;molybdenum
			IF (z[i] EQ 42) THEN BEGIN
				mo=mo_plc(te)
				plc[*,i]=mo
			ENDIF

			IF (z[i] EQ -42) THEN BEGIN
				mo=moPJ_plc(te)
				plc[*,i]=mo
			ENDIF

			;xenon
			IF (z[i] EQ 54) THEN BEGIN
				xe=xe_plc(te)
				plc[*,i]=xe
			ENDIF

			;tungsten
			IF (z[i] EQ 74) THEN BEGIN
				w=w_plc(te)
				plc[*,i]=w
			ENDIF			
		ENDFOR
			
		output=plc
		RETURN, output
	ENDELSE
END


;;;;;;;;;;;;;;;;;;;;;;;;;;;
;GENERAL UTILITIES
;;;;;;;;;;;;;;;;;;;;;;;;;;;


FUNCTION num2elem,n

	elem=strarr(93)
	elem[5]='B'
	elem[6]='C'
	elem[7]='N'
	elem[9]='F'
	elem[10]='Ne'
	elem[18]='Ar'
	elem[20]='Ca'
	elem[22]='Ti'
	elem[26]='Fe'
	elem[28]='Ni'
	elem[29]='Cu'
	elem[36]='Kr'
	elem[42]='Mo'
	elem[54]='Xe'
	elem[74]='W'


	IF (elem[n] EQ '') THEN OUTPUT='XX' ELSE OUTPUT=elem[n]
	
	RETURN, output
END

FUNCTION num2elem_name,n


	elem=strarr(93)
	elem[5]='Boron'
	elem[6]='Carbon'
	elem[7]='Nitrogen'
	elem[8]='Oxygen'
	elem[9]='Fluorine'
	elem[10]='Neon'
	elem[13]='Aluminum'
	elem[17]='Chlorine'
	elem[18]='Argon'
	elem[20]='Calcium'
	elem[21]='Scandium'
	elem[22]='Titanium'
	elem[26]='Iron'
	elem[28]='Nickel'
	elem[29]='Copper'
	elem[36]='Krypton'
	elem[42]='Molybdenum'
	elem[54]='Xenon'
	elem[74]='Tungsten'

	IF (elem[n] EQ '') THEN OUTPUT='XX' ELSE OUTPUT=elem[n]
	
	RETURN, output

END

PRO checkZ
	;script to show which elements are on file
	z=findgen(92)+1
	te=100.0
	plc=plc(te,z, /q)
	
	print, 'For 1 < Z < 92 the following are on file'
	for i=0,n(z) DO IF (plc[0,i] NE 0) THEN print, '  '+num2str(i+1)+' '+num2elem_name(i+1)
END

FUNCTION chkZ,z			;a binary version of checkZ using a function call
	te=1.0e3
	plc=plc(te,z, /q)
	
	IF (plc GT 0) THEN output=1 ELSE output=0

	RETURN,output
END

;;;;;;;;;;;;;;;;;;;;;;;;
;PLOTTING COOLING CURVES
;;;;;;;;;;;;;;;;;;;;;;;;

;+
;NAME:
;	PLCPLT
;	
;PURPOSE:
;	This procedure provides an easy interface for plotting multiple cooling-curves using
;	data from the function PLC
;
;CALLING_SEQUENCE:
;	PLCPLT
;
;OPTIONAL INPUTS:
;	temp:		temp=[te0,te1] sets the range [always log-log plot]
;	zrange:		zrange=[z0,z1] selects a range of Z's and plots all available
;	zz:		zz=[z0,z1...zN] selects specific Z's to plot
;	max:		max=FLT	value of the plot maximum
;	min:		min=FLT value of the plot minimum
;
;KEYWORD_PARAMETERS:
;	ps:		/ps to plot to postscript (output to /home/<username>/idl/plots/plcplt.ps)
;
;MODIFICATION HISTORY:
;	Written by:	ML Reinke, Spring 2005
;	9-28-05:	ML Reinke - adding proper documentation to work with HELPME
;
;-


PRO plcplt,temp=te_r,zrange=z_r,zz=zval,max=maxpt,min=minpt, ps=plot2ps,win=win,label=label
	
	;setup variable space
	IF (keyword_set(te_r)) THEN te=findgen(5001)/5000.0*(te_r[1]-te_r[0])+te_r[0] ELSE te=findgen(5001)/5000.0*(1.0e4-10.0)+10.0
	IF (keyword_set(z_r)) THEN z=findgen(z_r[1]-z_r[0]+1)+z_r[0] ELSE z=findgen(92)+1
	IF (keyword_set(zval)) THEN z=zval

	;run plc to fill array
	plc=plc(te,z,/q)

	;fix out of range values to -1
	FOR i=0,n(z) DO BEGIN
		IF (plc[0,i] NE 0) THEN BEGIN
			bad=where(plc[*,i] GT 0.99)
			IF (bad[0] NE -1) THEN FOR j=0,n(bad) DO plc[bad[j],i]=-1
		ENDIF
	ENDFOR
	
	;define color map for the good channels
	n_good=n(where(plc[0,*] NE 0))+1
	colormap=intarr(n_good)
	n_colors=!d.table_size-30	;those last few colors suck
	FOR i=0,n_good-1 DO BEGIN
		colormap[i]=i*(n_colors)/(n_good)
	ENDFOR


	;setup ps plot, if chosen
	IF (keyword_set(plot2ps)) THEN BEGIN
		current_device=!d.name
		psplot_3
	ENDIF ELSE BEGIN
		IF keyword_set(win) THEN IF win NE -1 THEN wreset,10
	ENDELSE

	IF (keyword_set(maxpt)) THEN maxplot=maxpt ELSE maxplot=max(plc)
	IF (keyword_set(minpt)) THEN minplot=minpt ELSE minplot=1e-34

	;define plot
	IF keyword_set(plot2ps) THEN tit='Plasma Power Loss Coefficents' ELSE tit='Plasma Power Loss Coefficents for Various Impurities, Z'
	plot,[1],[1],/xlog,/ylog,title=tit, ytitle='PLC [W*m!u3!n]', xtitle='Te [eV]', yrange=[minplot,maxplot],xrange=[te[0], te[n(te)]],chars=1.3,/xsty
	
	;run through array and plot valid entries
	cntr=0
	IF NOT keyword_set(label) THEN label=strarr(n(z)+1)
	IF keyword_set(plot2ps) THEN BEGIN
		xo=0.2
		yo=0.2
		del=0.05
	ENDIF ELSE BEGIN
		xo=0.16
		yo=0.13
		del=0.05
	ENDELSE

	FOR i=0, n(z) DO BEGIN
		IF (plc[0,i] NE 0) THEN BEGIN
			oplot,te,plc[*,i],color=colormap[cntr]	
			xyouts,xo,yo+del*cntr,/norm,'Z='+num2str(abs(z[i]),1)+' '+num2elem(floor(abs(z[i])))+' '+label[i],color=colormap[cntr]
			cntr+=1
		ENDIF
	ENDFOR
		

	;close ps plot if chosen
	IF (keyword_set(plot2ps)) THEN BEGIN
		device,/close
		set_plot,current_device
		spawn, 'pwd',out_string
		spawn, 'cp idl.ps /home/'+logname()+'/idl/plots/plcplt.ps'
		IF (strmatch(out_string, '*mlreinke*') EQ 1) THEN BEGIN
			xwplot					
		ENDIF
	ENDIF
	
	
END 



	
