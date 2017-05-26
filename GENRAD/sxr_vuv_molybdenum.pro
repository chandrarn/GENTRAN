;***7/19/11 - adapted from lines_for_moly.pro which uses COMMMON
;             blocks from MIST setup.  Want to use with GENTRAN outputs
;
;	dens 	[nq,nr] charge state [m-3]
;	te	[nr]	electron temp [eV]
;	nelec	[nr]	electron density [m-3]
;
;12-06-05 added Ne-like emission at 3.74 Ang using MEWE 

;+
;NAME:
;	SXR_VUV_MOLYBDENUM
;
;PURPOSE:
;	This procedure calculates the emissivity of various SXR and VUV lines from
;	highly ionized molybdenum for given Te, ne and molybdenum charge state
;	density profiles
;
;CALLING SEQUENCE:
;	SXR_VUV_MOLYBDENUM,dens,te,nelec,molywave,molychg,molyemiss
;	
;INPUTS:
;	dens:		FLTARR	[nq,nr]	of the density profile for each charge state (index is i=q+)[m^-3]
;	te:		FLTARR	[nr]	of the electron temperature for each radial location [eV]
;	nelec:		FLTARR	[nr]	of the electron density for each radial location [m^-3]
;
;OUTPUTS:
;	molywave	FLTARR	[nlines] of the wavelength of the modeled transitions [Ang]
;	molychg		INTARR	[nlines] of the charge state for each line (index is i=q+)
;	molyemiss	FLTARR	[nlines,nr] of the emissivity for each line at each radial location [ph/s/m^3] 
;
;MODIFICATION HISTORY:
;	7/19/11		M.L. Reinke - adapated from LINES_FOR_MOLY which had hardcoded common blocks
;				for use with MIST utilities.
;	6/8/12		M.L. Reinke - moved moly data from mist director to /usr/local/cmod/idl/atomic_physics/
;-
	
PRO sxr_vuv_molybdenum,dens,te,nelec,molywave,molychg,molyemiss,verb=verb,molyatomic_dat=molyatomic_dat
	dens_in=dens
	temp_in=te
	nelec_in=nelec
	
	dens/=1.0e6 ;[/m^3 -> /cm^3]
	te/=1.0e3 ;[eV -> keV]
	nelec/=1.0e6 ;[/m^3 -> /cm^3]

;common moly, highw,vzone, totdpr
;common emi,molyemis,emis
;common emi2, molywave,molychg
;common pars,te,nelec,rad
;COMMON SHOT,FNT,ELEMENT,jxz
;common den,dens,neuts
;common lemi,lemist,emhead,nlines
 
;========================= Define variables===================

	molyline = 300
        numcoll = 10
	nrad=n(te)+1					;number of radial bins
	
        collte = fltarr(molyline,numcoll)   		; collision rate temperatures
        molybij = fltarr(molyline)  			; moly line branching ratios
        molychg = fltarr(molyline)      		; charge state of line
        molychg[*] = -1.0
        molytext = strarr(molyline)     		; text to descibe each line
        molyarr = intarr(molyline)      		; number of cx and te in array

        molywave = fltarr(molyline) 			; Transition wavelength

        molywave1 = fltarr(molyline)    		; referenced wavelength for 1st cx
        collex = fltarr(molyline,numcoll)   		; collisional excitation (Te, Ion)
        molyn = fltarr(molyline)        		; fraction of moly ions in ground state/total ions of a charge state

        collconst = fltarr(molyline,numcoll)  		; collisional excitation const

	; 2nd set
        molywave2 = fltarr(molyline)     		; reference wave for 2nd cx value
        collex2 = fltarr(molyline,numcoll)		; 2nd cx value for moly charge states
        molyn2 = fltarr(molyline) 			; fraction of moly ions in ground state/total ions of a charge state
        collconst2 = fltarr(molyline,numcoll) 		; const for collisional excitation

	;3rd set
        molywave3 = fltarr(molyline)     		; reference wave for 3rd cx value
        collex3 = fltarr(molyline,numcoll)		; 3rd cx value for moly charge states
        molyn3 = fltarr(molyline) 			; fraction of moly ions in ground state/total ions of a charge state
        collconst3 = fltarr(molyline,numcoll) 		; const for collisional excitation

	; 4th set
        molywave4 = fltarr(molyline)     		; reference wave for 4th cx value
        collex4 = fltarr(molyline,numcoll)		; 4th cx value for moly charge states
        molyn4 = fltarr(molyline) 			; fraction of moly ions in ground state/total ions of a charge state
        collconst4 = fltarr(molyline,numcoll) 		; const for collisional excitation

        moly_cool = fltarr(molyline,nrad) 		; Cooling rate to compare with PJ data
        molyemis = fltarr(molyline,nrad) 		; emissivity per radial (MIST) for each line
        molydpr = fltarr(molyline,nrad)   		; line power per MIST radii for each line
        totdpr = fltarr(nrad) 				; total line power for each MIST radii
                                                                                  

;========================================================
;	Load moly data
;=======================================================

  if keyword_set(molyatomic_dat) then begin
  	nams=tag_names(molyatomic_dat)
	for j=0,n(nams) do st=execute(nams[j]+'=molyatomic_dat.'+nams[j]) ; unpack structure into namespace...
  endif else begin

	ss = ''
	i = 0

	;openr,2,'/usr/local/cmod/codes/spectroscopy/mist/molyatomicplm1_exp_wave.dat'
	openr,2,'/usr/local/cmod/idl/atomic_physics/molyatomicplm1_exp_wave.dat'
	REPEAT BEGIN
		readf,2,ss
		IF ( strpos(ss,';') EQ 0 OR strpos(ss,'@') EQ 0 ) THEN goto, skip
		IF keyword_set(verb) THEN print, 'reading excit. data for line ', ss
		molywave[i] = float(ss)
		readf,2,ss
		IF keyword_set(verb) THEN print, ss
		molytext[i] = ss
		readf,2,ss
		IF keyword_set(verb) THEN print, ss
		molyarr[i] = fix(ss)
		FOR j = 0,molyarr[i]-1 DO BEGIN
			readf,2,ss
			IF keyword_set(verb) THEN print, ss
			collte[i,j] = float(ss)	
		ENDFOR
		readf,2,ss
		molywave1(i) = float(ss)
		FOR j = 0,molyarr(i)-1 DO BEGIN
        		readf,2,ss
        		collex[i,j] = float(ss)
			collconst[i,j] = collex[i,j]*sqrt(collte[i,j])*exp(12400.0/(molywave1[i]*collte[i,j]))
			IF keyword_set(verb) THEN print, collconst(i,j)
		ENDFOR
		readf,2,ss
		molyn[i] = float(ss)    ; read density of  state
		readf,2,ss
		molychg[i] = float(ss)  ; read charge state
		readf,2,ss
		molybij[i] = float(ss)  ; read branching ratios
		i = i + 1
		skip:
	ENDREP UNTIL (strpos(ss,'@') GE 0 OR i GT 399)
	close,2

;==============================================================
;	Load file Number 2
;==============================================================

	ss = ''
	;openr,2,'/usr/local/cmod/codes/spectroscopy/mist/molyatomicplm2_exp_wave.dat'
	openr,2,'/usr/local/cmod/idl/atomic_physics/molyatomicplm2_exp_wave.dat'
	REPEAT BEGIN
		readf,2,ss
		IF ( strpos(ss,';') EQ 0 OR strpos(ss,'@') EQ 0 )THEN goto, skip2
		IF keyword_set(verb) THEN print, ss
		IF keyword_set(verb) THEN print, 'reading excit. data for line ', ss
		molywave[i] = float(ss)
		readf,2,ss
		molytext[i] = ss
		readf,2,ss
		molyarr[i] = fix(ss)
		FOR j = 0,molyarr(i)-1 DO BEGIN
			readf,2,ss
			collte[i,j] = float(ss)	
		ENDFOR
		readf,2,ss
		molywave1[i] = float(ss)
		FOR j = 0,molyarr[i]-1 DO BEGIN
		          readf,2,ss
		          collex[i,j] = float(ss)
		          collconst[i,j] = collex[i,j]*sqrt(collte[i,j])*exp(12400.0/(molywave1[i]*collte[i,j]))
		ENDFOR
		readf,2,ss
		molyn[i]= float(ss) 
		readf,2,ss
		molywave2[i] = float(ss)
		FOR j = 0,molyarr[i]-1 DO BEGIN
		        readf,2,ss
			IF keyword_set(verb) THEN print, ss
		        collex2[i,j] = float(ss)
			collconst2[i,j] = collex[i,j]*sqrt(collte[i,j])*exp(12400.0/(molywave2[i]*collte[i,j]))
			IF keyword_set(verb) THEN print, collconst[i,j]
		ENDFOR
		readf,2,ss
		molyn2[i]=float(ss)
		readf,2,ss
		molychg[i] = float(ss)
		readf,2,ss
		molybij[i] = float(ss)
		i = i + 1
		skip2:
	ENDREP UNTIL (strpos(ss,'@') GE 0 OR i GT 399)
	close,2

;==============================================================
;       Load file Number 3
;       file_4 = [8,14,15,16,17,18,19]
;       contains Mo XXIX, Mo XXVII, Mo XXVIII, Mo XXIV, Mo XXV, Mo XXVI, Mo XXXV
;==============================================================

        ss = ''
	;openr,2,'/usr/local/cmod/codes/spectroscopy/mist/molyatomicplm3_exp_wave.dat'
	openr,2,'/usr/local/cmod/idl/atomic_physics/molyatomicplm3_exp_wave.dat'
        REPEAT BEGIN
        	readf,2,ss
        	IF ( strpos(ss,';') EQ 0 OR strpos(ss,'@') EQ 0 )THEN goto, skip3
	        IF keyword_set(verb) THEN print, ss
       		IF keyword_set(verb) THEN print, i
		IF keyword_set(verb) THEN print, 'reading excit. data for line ', ss
        	molywave[i] = float(ss)
        	readf,2,ss
        	molytext[i] = ss
        	readf,2,ss
	        molyarr[i] = fix(ss)
        	FOR j = 0,molyarr[i]-1 DO BEGIN
                	readf,2,ss
	                collte[i,j] = float(ss)
        	ENDFOR
       		readf,2,ss
	        molywave1[i] = float(ss)
        	FOR j = 0,molyarr[i]-1 DO BEGIN
                	readf,2,ss
	                collex[i,j] = float(ss)
        	        collconst[i,j] = collex[i,j] *sqrt(collte[i,j])*exp(12400.0/(molywave1[i]*collte[i,j]))
        	ENDFOR
        	readf,2,ss
        	molyn[i]= float(ss) ; density of first state feeding line

		;second line
        	readf,2,ss
        	molywave2[i] = float(ss) 
		; wavelength of second line feeding state
	        FOR j = 0,molyarr[i]-1 DO BEGIN
        	        readf,2,ss
                	collex2[i,j] = float(ss)
	                collconst2[i,j] = collex2[i,j]*sqrt(collte[i,j])*exp(12400.0/(molywave2[i]*collte[i,j]))
        	ENDFOR
	        readf,2,ss
        	molyn2[i]=float(ss) ; density of second state feeding line

		;third line
	        readf,2,ss
        	molywave3[i] = float(ss) 
		; wavelength of second line feeding state
		FOR j = 0,molyarr[i]-1 DO BEGIN
                	readf,2,ss
	                collex3[i,j] = float(ss)
        	        collconst3[i,j] = collex3[i,j]*sqrt(collte[i,j])*exp(12400.0/(molywave3[i]*collte[i,j]))
	        ENDFOR
        	readf,2,ss
	        molyn3[i]=float(ss) ; density of second state feeding line

		;Fourth line
	        readf,2,ss
        	molywave4[i] = float(ss) 
		; wavelength of second line feeding state
        	FOR j = 0,molyarr[i]-1 DO BEGIN
                	readf,2,ss
	                collex4[i,j] = float(ss)
        	        collconst4[i,j] = collex4[i,j]*sqrt(collte[i,j])*exp(12400.0/(molywave4[i]*collte[i,j]))
	        ENDFOR
        	readf,2,ss
	        molyn4[i]=float(ss) ; density of second state feeding line
	
        	readf,2,ss
	        molychg[i] = float(ss); charge state
        	readf,2,ss
	        molybij[i] = float(ss); branch ratio of line
        	i = i + 1
	        skip3:
        ENDREP UNTIL strpos(ss,'@') ge 0
        close,2
        skipx2:
	skip_file3:

molyatomic_dat={molywave:molywave, molytext:molytext, molyarr:molyarr, collte:collte, $
		molywave1:molywave1, collex:collex,   collconst:collconst,   molyn:molyn, $ 
		molywave2:molywave2, collex2:collex2, collconst2:collconst2, molyn2:molyn2, $
		molywave3:molywave3, collex3:collex3, collconst3:collconst3, molyn3:molyn3, $
		molywave4:molywave4, collex4:collex4, collconst4:collconst4, molyn4:molyn4, $
		molychg:molychg, molybij:molybij, i:i $
		}
  endelse ; keyword_set(molyatomic_dat)


;=========================================================================
;	DO MOLY CALCLATIONS
;=========================================================================
;
	;Transitional rates and Calculate emissivities
	molyline = i						;this ends up being the #lines that were actually filled
	molycx = fltarr(molyline,nrad)
	molycx2 = fltarr(molyline,nrad)
	molycx3 = fltarr(molyline,nrad)
	molycx4 = fltarr(molyline,nrad)
	molywave = molywave(0:molyline-1)

	FOR j = 0,molyline-1 DO BEGIN  ; do for each line
	        cxtemp = where (collex(j,*) gt 0.0)
        	molyfit1 = poly_fit(collte(j,cxtemp),collconst(j,cxtemp),1)
	        molyfit2 = poly_fit(collte(j,cxtemp),collconst2(j,cxtemp),1)
        	molyfit3 = poly_fit(collte(j,cxtemp),collconst3(j,cxtemp),1)
	        molyfit4 = poly_fit(collte(j,cxtemp),collconst4(j,cxtemp),1)

        	mdens = dens (molychg(j),*)
	        molycxtemp = molyfit1(0) + molyfit1(1)*te(*)*1000.
        	molycxtemp2 = molyfit2(0) + molyfit2(1)*te(*)*1000.
	        molycxtemp3 = molyfit3(0) + molyfit3(1)*te(*)*1000.
        	molycxtemp4 = molyfit4(0) + molyfit4(1)*te(*)*1000.
	        molycx(j,*) =(molycxtemp/sqrt(te(*)*1000.))*exp(-12400/(molywave1(j)*te(*)*1000.))
        	if (molywave2(j) ne 0.0) then molycx2(j,*)=(molycxtemp2/sqrt(te(*)*1000.))*exp(-12400/(molywave2(j)*te(*)*1000.))
        	if (molywave3(j) ne 0.0) then molycx3(j,*)=(molycxtemp3/sqrt(te(*)*1000.))*exp(-12400/(molywave3(j)*te(*)*1000.))
        	if (molywave4(j) ne 0.0) then molycx4(j,*)=(molycxtemp4/sqrt(te(*)*1000.))*exp(-12400/(molywave4(j)*te(*)*1000.))
        	molyemis(j,*) = mdens* molybij(j)* nelec (*) * (molycx(j,*)*molyn(j)+molycx2(j,*)*molyn2(j)$
			+molycx3(j,*)*molyn3(j)+molycx4(j,*)*molyn4(j))

		;Prevents emission fit problems
        	molyemis = molyemis >0
        	moly_cool(j,*) = (1.99E-15 / molywave(j))*molybij(j)*$
                  	(molycx(j,*)*molyn(j) +molycx2(j,*)*molyn2(j) $
                	+molycx3(j,*)*molyn3(j) +molycx4(j,*)*molyn4(j))
        	molydpr(j,*) = molyemis(j,*) * 1.99E-15 / molywave(j)
	ENDFOR
	totdpr = total (molydpr,1)

	;MLR 12-06-05 Adding Ne-like emission from HIREX
	; This uses the MEWE method that is employed in LINES.PRO
	; to calculate the Ne-like line emission that is monitoried by
	; HIREX from Mo XXXIII at 3.74 Ang.

	;**data taken from atdata.dat**
	;------------------------------------
	;42  33  3.740   12   /Ne-like 4d-2p
	; 3315.0  .517   .09  .14  0.0  .28
	;------------------------------------

	emhead=[42,  33,  3.740,   12]
	empar=[0,3315.0,  .517,   .09,  .14,  0.0,  .28]
	ll=0
	mo_nelike=fltarr(nrad)
	j=0 ;first time slice (time indepenendt hard coded)
	FOR i=0,nrad-1 DO BEGIN
		kzone=i
        	pn=nelec[kzone]
        	if (emhead(3,ll) eq 12) then begin
            		zy=empar(1,ll)/(te(kzone)*1000)
            
            		IF (ZY GE 1.5) THEN begin
                		FF=ALOG((ZY+1)/ZY) - (.36+.03*SQRT(ZY+.01))/(ZY+1)/(ZY+1)
	            	endif ELSE begin
        	       		FF=ALOG((ZY+1)/ZY) - (.36+.03/SQRT(ZY+.01))/(ZY+1)/(ZY+1)
            		ENDelse
            
	            	gg=empar(3,ll)+(empar(4,ll)*zy-empar(5,ll)*zy*zy+empar(6,ll))*ff+empar(5,ll)*zy 
        	    
            		ssnalines=1.58e-5/sqrt(te(kzone)*1000)/empar(1,ll)*empar(2,ll)*gg*exp(-zy)
	            	mo_nelike[i]=pn*dens(emhead(1,ll)-1,kzone,j)*ssnalines
        	endif
	ENDFOR

	;add on the ne-like data to the existing infrastructure
	;with a hard coded wavelength and charge state
	molywave=[molywave,3.740] 
	tmp=where(molychg EQ -1)
	molychg[tmp[0]]=32
	molyemis[tmp[0],*]=mo_nelike

	;PROOF THAT molyemis=[mo_nelike,transpose(molyemis)] actually works
	;-------------------
	;a=intarr(4,4)+1
	;b=intarr(4)
	;c=[a,transpose(b)]	;adds as another i-index to the end
	;------------------

	num_lines=n_elements(molywave)
	nlines=num_lines
	lemist = fltarr(nlines,nrad,200)
	emis = molyemis(0:num_lines-1,*)
	lemist(*,*,0) = molyemis(0:num_lines-1,*)

	tmp=where(molychg NE -1)
	molywave=molywave
	molychg=int(molychg[tmp])
	molyemiss=molyemis[tmp,*]*1.0e6		;convert /cm^3 to /m^3
	dens=dens_in
	te=temp_in
	nelec=nelec_in

END

FUNCTION nelike_mo_4d,dens,te,nelec

	nrad=n(te)+1					;number of radial bins
	emhead=[42,  33,  3.740,   12]
	empar=[0,3315.0,  .517,   .09,  .14,  0.0,  .28]
	ll=0
	emiss=fltarr(nrad)
	j=0 ;first time slice (time indepenendt hard coded)
	FOR i=0,nrad-1 DO BEGIN
		kzone=i
        	pn=nelec[kzone]
        	if (emhead(3,ll) eq 12) then begin
            		zy=empar(1,ll)/(te(kzone))
            
            		IF (ZY GE 1.5) THEN begin
                		FF=ALOG((ZY+1)/ZY) - (.36+.03*SQRT(ZY+.01))/(ZY+1)/(ZY+1)
	            	endif ELSE begin
        	       		FF=ALOG((ZY+1)/ZY) - (.36+.03/SQRT(ZY+.01))/(ZY+1)/(ZY+1)
            		ENDelse
            
	            	gg=empar(3,ll)+(empar(4,ll)*zy-empar(5,ll)*zy*zy+empar(6,ll))*ff+empar(5,ll)*zy 
        	    
            		ssnalines=1.58e-5/sqrt(te(kzone))/empar(1,ll)*empar(2,ll)*gg*exp(-zy)
	            	emiss[i]=pn*dens(emhead(1,ll)-1,kzone,j)*ssnalines
        	endif
	ENDFOR
	RETURN,emiss*1.0e-6
END
