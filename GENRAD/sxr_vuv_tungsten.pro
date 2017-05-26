;+
;NAME:
;	SXR_VUV_TUNGSTEN
;
;PURPOSE:
;	This procedure calculates the emissivity of various SXR and VUV lines from
;	highly ionized tungsten for given Te, ne and tungsten charge state
;	density profiles
;
;CALLING SEQUENCE:
;	SXR_VUV_TUNGSTN,dens,te,nelec,wave,chg,emiss
;	
;INPUTS:
;	dens:		FLTARR	[nq,nr]	of the density profile for each charge state (index is i=q+)[m^-3]
;	te:		FLTARR	[nr]	of the electron temperature for each radial location [eV]
;	nelec:		FLTARR	[nr]	of the electron density for each radial location [m^-3]
;
;OUTPUTS:
;	wave		FLTARR	[nlines] of the wavelength of the modeled transitions [Ang]
;	chg		INTARR	[nlines] of the charge state for each line (index is i=q+)
;	emiss		FLTARR	[nlines,nr] of the emissivity for each line at each radial location [ph/s/m^3] 
;
;MODIFICATION HISTORY:
;	Written by:	M.L. Reinke 9/16/14 (adapted from SXR_VUV_MOLYBDENUM)
;
;-
	
PRO sxr_vuv_tungsten,dens,te,nelec,wave,chg,emiss,data=data
	;manually linked lines from Ralenchko J.Phys. B v40 pg3861 to specifictransitions in ADAS files
	wave=[128.52,134.88,131.24,139.14,128.95,129.41,134.95,126.29,128.17,134.81,132.88,134.8,127.12]
	chg=[40,40,41,41,42,42,42,43,43,43,44,44,45]
	pecnum=[47,48,48,49,46,49,48,48,47,49,49,50,50]
	path='/usr/local/cmod/idl/atomic_physics/adas/pec40#w_ic#w'	
	telog=alog10(te)
	nelog=alog10(nelec)

	nrad=n(te)+1
	nlines=n(wave)+1
	emiss=fltarr(nlines,nrad)

	FOR i=0,nlines-1 DO BEGIN
		IF NOT keyword_set(data) THEN BEGIN
			file=path+num2str(chg[i],1)+'.dat'
 			idata=read_pec_file(file)
			IF i EQ 0 THEN exdata=create_struct('line'+num2str(i,1),idata) ELSE exdata=create_struct(exdata,'line'+num2str(i,1),idata)
                ENDIF ELSE idata=data.(i)
		pec=reform(idata.pec[pecnum[i]-1,*,*])
		jvec=interp_vec_reform(alog10(idata.dens),nelog)
		ivec=interp_vec_reform(alog10(idata.temp),telog)
		intpec=interpolate(pec,ivec,jvec,missing=0.0)
		emiss[i,*]=nelec*dens[chg[i],*]*intpec
	ENDFOR
	IF NOT keyword_set(data) THEN data=exdata
END

