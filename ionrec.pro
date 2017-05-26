;+
;NAME:
;	READ_REC_DATA
;
;PURPOSE:
;	This function loads the REC stucture by calling z-specific routines.
;
;CALLING SEQUENCE:
;	result=READ_REC_DATA(z)
;
;INPUTS:
;	z:	INT	the atomic number for the rates of interest
;
;KEYWORD PARAMETERS:
;	adas:	/adas loads the data from adas tables
;
;OUTPUTS:
;	result:	STRUC	containing the recombination rates
;		*.rates	FLTARR 	[z+1,nt,nd] of the recombination rates [i,*,*] for ionization stage [i] for temperature and density [m^3/s]
;		*.temp	FLTARR 	[nt] of the electron temperatures for the rates [eV]
;		*.dens	FLTARR	[nd] of the electron densities for the rates [m^-3]
;		*.z	INT	z input to the function
;		*.path	STR	of the path to the datafile where the rates were loaded
;
;OPTIONAL OUTPUTS:
;	cxr:	STRUC	containing the charge exchange recombination rates [same format as result]
;
;MODIFICATION HISTORY:
;	Written by: 	ML Reinke - 3/09
;	6/7/12		ML Reinke - modified pathing for tungsten ADAS file to /usr/local/cmod/idl/atomic_physics/
;	6/11/12		ML Reinke - modified to force a /load going into READ_LOCH_REC_DATA
;	6/12/12		ML Reinke - removed READ_LOCH_REC_DATA /load
;
;-


FUNCTION read_rec_data,z,load=load,adas=adas,cxr=cxr
	cxr=-1
	IF not keyword_set(adas) THEN BEGIN
	 CASE z OF 
		18 : BEGIN
			adas_path,z,ion,rec,cxr
			rec=read_loch_rec_data()
			cxr=read_acd_file(cxr)
		END
		36 : BEGIN
			adas_path,z,ion,rec,cxr
			rec=read_matt06_rec_data(z)
			cxr=read_acd_file(cxr)
		END
		42 : BEGIN
			adas_path,z,ion,rec,cxr
			rec=read_matt06_rec_data(z)
			cxr=read_acd_file(cxr)
		END
		74 : BEGIN
			;path='/home/mlreinke/idl/impurities/data/adas/acd96_w.dat'
			path='/usr/local/cmod/idl/atomic_physics/adas/acd96_w.dat'
			rec=read_putterich_file(path)
		END
		ELSE : BEGIN
			print,'defaulting to ADAS recombination rates...'
			adas=1
		END
	 ENDCASE
	ENDIF

	IF keyword_set(adas) THEN BEGIN
		adas_path,z,ion,rec,cxr
		IF rec NE 'NA' THEN rec=read_acd_file(rec)
		IF cxr NE 'NA' THEN cxr=read_acd_file(cxr)
		RETURN,rec
	ENDIF

	output=rec
	RETURN,output
END

;+
;NAME:
;	READ_ION_DATA
;
;PURPOSE:
;	This function loads the ION stucture by calling z-specific routines.
;
;CALLING SEQUENCE:
;	result=READ_ION_DATA(z)
;
;INPUTS:
;	z:	INT	the atomic number for the rates of interest
;
;KEYWORD PARAMETERS:
;	adas:	/adas loads the data from adas tables
;
;OUTPUTS:
;	result:	STRUC	containing the ionization rates
;		*.rates	FLTARR 	[z+1,nt,nd] of the ionization rates [i,*,*] for ionization stage [i] for temperature and density [m^3/s]
;		*.temp	FLTARR 	[nt] of the electron temperatures for the rates [eV]
;		*.dens	FLTARR	[nd] of the electron densities for the rates [m^-3]
;		*.z	INT	z input to the function
;		*.path	STR	of the path to the datafile where the rates were loaded
;
;MODIFICATION HISTORY:
;	Written by: 	ML Reinke - 3/09
;	6/7/12		ML Reinke - modified pathing for tungsten ADAS file to /usr/local/cmod/idl/atomic_physics/
;	6/11/12		ML Reinke - modified to force a /load going into READ_LOCH_ION_DATA
;	6/12/12		ML Reinke - removed READ_LOCH_REC_DATA /load
;
;-

FUNCTION read_ion_data,z,load=load,adas=adas
	IF not keyword_set(adas) THEN BEGIN
	 CASE z OF 
		18 : BEGIN
			ion=read_loch_ion_data()
		END
		36 : BEGIN
			ion=read_matt06_ion_data(z)
		END
		42 : BEGIN
			ion=read_matt06_ion_data(z)
		END
		74 : BEGIN
			;path='/home/mlreinke/idl/impurities/data/adas/scd96_w.dat'
			path='/usr/local/cmod/idl/atomic_physics/adas/scd96_w.dat'
			ion=read_scd_file(path)
			;path='/home/mlreinke/idl/impurities/data/adas/scd01_w.dat'
			;ion=read_putterich_file(path)
		END
		ELSE : BEGIN
			print,'defaulting to ADAS ionization rates...'
			adas=1
		END
	 ENDCASE
	ENDIF

	IF keyword_set(adas) THEN BEGIN
		adas_path,z,ion,rec
		IF ion NE 'NA' THEN ion=read_scd_file(ion)
		RETURN,ion
	ENDIF

	output=ion
	RETURN,output
END

PRO write_scd_file,ion
	path='/home/'+logname()+'/scd.dat'
	openw,lun,path,/get_lun
	x=size(ion.rates)
	ntemp=x[2]
	ndens=8
	dinfo=strsplit(systime(),' ', /extract)
	year=last(dinfo)
	date=dinfo[1]+', '+dinfo[2]+' '+year

	header='  '+num2str(ion.z,1)+'   '+num2str(ndens,1)+'   '+num2str(ntemp,1)+'   1   '+num2str(ion.z,1)+'   /'+strupcase(num2elem_name(ion.z))+$
		'     '+ion.path+' '+year
	printf,lun,header
	printf,lun,'-------------------------------------------------------------------------------'
	printf,lun,make(10,16,ndens),format='(8(f9.5,x))'		;in alog10(cm^-3)
	printf,lun,alog10(ion.temp),format='(8(f9.5,x))'		;in alog10(eV)
	
	FOR i=0,ion.z-1 DO BEGIN
		printf,lun,'------------------/ IPRT= 1  / IGRD= 1  /IONIS   / Z1= '+num2str(i+1,1)+'   / DATE= '+date
		FOR j=0,ntemp-1 DO BEGIN
			data=dblarr(ndens)+ion.rates[i,j,0]*1.0e6
			tmp=where(data EQ 0.0)
			IF tmp[0] NE -1 THEN data[tmp]=double(10^(-36.0))
			printf,lun,alog10(data),format='(8(f9.5,x))'
		ENDFOR
	ENDFOR
	printf,lun,'C-------------------------------------------------------------------------------'
	printf,lun,'C'
	printf,lun,'C  Generate ionisation from : '
	printf,lun,'C     .'+ion.path
	printf,lun,'C'
	printf,lun,'C  CODE     : CONVERTED FROM IDL SAVESET'
	printf,lun,'C  PRODUCER : GENTRAN - '+logname()              
	printf,lun,'C  DATE     : '+date
	printf,lun,'C'
	printf,lun,'C-------------------------------------------------------------------------------'


	close,lun
	free_lun,lun
END

PRO write_acd_file,rec
	path='/home/'+logname()+'/acd.dat'
	openw,lun,path,/get_lun
	x=size(rec.rates)
	ntemp=x[2]
	ndens=8
	dinfo=strsplit(systime(),' ', /extract)
	year=last(dinfo)
	date=dinfo[1]+', '+dinfo[2]+' '+year

	header='  '+num2str(rec.z,1)+'   '+num2str(ndens,1)+'   '+num2str(ntemp,1)+'   1   '+num2str(rec.z,1)+'   /'+strupcase(num2elem_name(rec.z))+$
		'     '+rec.path+' '+year
	printf,lun,header
	printf,lun,'-------------------------------------------------------------------------------'
	printf,lun,make(10,16,ndens),format='(8(f9.5,x))'		;in alog10(cm^-3)
	printf,lun,alog10(rec.temp),format='(8(f9.5,x))'		;in alog10(eV)
	
	FOR i=1,rec.z DO BEGIN
		printf,lun,'------------------/ IPRT= 1  / IGRD= 1  /IONIS   / Z1= '+num2str(i,1)+'   / DATE= '+date
		FOR j=0,ntemp-1 DO BEGIN
			data=dblarr(ndens)+rec.rates[i,j,0]*1.0e6
			tmp=where(data EQ 0.0)
			IF tmp[0] NE -1 THEN data[tmp]=double(10^(-36.0))
			printf,lun,alog10(data),format='(8(f9.5,x))'
		ENDFOR
	ENDFOR
	printf,lun,'C-------------------------------------------------------------------------------'
	printf,lun,'C'
	printf,lun,'C  Generate recombination from : '
	printf,lun,'C     .'+rec.path
	printf,lun,'C'
	printf,lun,'C  CODE     : CONVERTED FROM IDL SAVESET'
	printf,lun,'C  PRODUCER : GENTRAN - '+logname()              
	printf,lun,'C  DATE     : '+date
	printf,lun,'C'
	printf,lun,'C-------------------------------------------------------------------------------'


	close,lun
	free_lun,lun
END
