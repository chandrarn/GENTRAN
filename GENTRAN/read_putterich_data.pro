FUNCTION read_putterich_file,path,debug=debug
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
		tmp=strsplit(line,'_',/extract)
		tmp=strsplit(tmp[0], '=',/extract)
		q=int(tmp[1])
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