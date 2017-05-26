;+
;NAME:
;	READ_XRAY_DATA
;
;PURPOSE:
;	This function reads soft x-ray transmission tables output from
;	the LBL Henke tables (http://henke.lbl.gov/optical_constants/)
;
;CALLING SEQUENCE:
;	result=READ_XRAY_DATA(path)
;
;INPUT:
;	path	STR of the file path to the transmission data table
;
;OUTPUT:
;	result:	STRUC
;		*.e	FLTARR [n] of the photon energy [eV]
;		*.tr	FTLARR [n] of the transmission at that energy
;		*.elem	STR of the element or compound the filter material
;		*.thick	FLT of the filter thickness [microns]
;
;MODIFICATION HISTORY:
;	Written by: 	M.L. Reinke - back in 07 or 08 maybe?
;
;-

FUNCTION read_xray_data,path,debug=debug

	openr,lun, path,/get_lun
	line=strarr(1)
	
	readf,lun,line
	out=strsplit(line, ' ',/extract)
	symbol=out[0]
	a=strpos(out[2], '=')+1
	thick=float(strmid(out[2],a,(strlen(out[2])-a)))

	readf,lun,line
	
	energy=fltarr(500)
	trans=fltarr(500)

	readf,lun,line
	cntr=0
	WHILE eof(lun) NE 1  DO BEGIN
		out=strsplit(line, ' ',/extract)
		energy[cntr]=float(out[0])
		trans[cntr]=double(out[1])
		cntr+=1	
		readf,lun,line
	ENDWHILE
	tmp=where(energy EQ 0.0)
	IF tmp[0] NE -1 THEN BEGIN
		energy=energy[0:tmp[0]-1]
		trans=trans[0:tmp[0]-1]
	ENDIF

	output={e:energy,tr:trans,elem:symbol,thick:thick}
	
	close,lun
	free_lun, lun
	IF keyword_set(debug) THEN stop
	RETURN,output
END

FUNCTION read_suthgff_data,path,debug=debug
	openr,lun,path,/get_lun
	line=strarr(1)
	npts=3321
	u=fltarr(npts)
	gsq=fltarr(npts)
	gff=fltarr(npts)
	FOR i=0,4 DO readf,lun,line
	FOR i=0,npts-1 DO BEGIN
		readf,lun,line
		data=strsplit(line,' ',/extract)
		u[i]=float(data[0])
		gsq[i]=float(data[1])
		gff[i]=float(data[2])
	ENDFOR

	ufix=gsq[0:80]		;this is due to a mislabeling of the data table
	gsq_fix=u[indgen(41)*81]	;this is due to a mislabeling of the data table
	gff=reform(gff,81,41)
	output={u:ufix,gsq:gsq_fix,gff:gff}
	
	close,lun
	free_lun, lun
	IF keyword_set(debug) THEN stop	

	RETURN,output
END
