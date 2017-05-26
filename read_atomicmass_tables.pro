;10/15/20 - updated to use the environmental variable
;

FUNCTION read_atomic_name,z, debug=debug, verb=verb
	path=getenv('ATOMIC_PHYSICS_PATH')
	;filename='/home/mlreinke/idl/impurities/data/atomic_mass'
	;filename='/usr/local/cmod/idl/atomic_physics/atomic_mass'
	filename=path+'atomic_mass'
	openr,lun,filename,/get_lun
	ident_str='*Atomic Number*'

	output=''
	stopsearch=-1
	line=strarr(1)	

	readf,lun,line
	
	WHILE eof(lun) NE 1 AND stopsearch EQ -1 DO BEGIN
		IF strmatch(line, ident_str, /fold_case) THEN BEGIN
			a=strpos(line, '=')+1
			number=int(strmid(line,a,(strlen(line)-a)))
			IF number EQ z  THEN BEGIN
				readf,lun,line
				stopsearch=1
				a=strpos(line, '=')+1
				output=strmid(line,a,(strlen(line)-a))
			ENDIF 
		ENDIF
		readf,lun,line
	ENDWHILE
	
	close,lun
	free_lun, lun
	IF keyword_set(debug) THEN stop
	output=strmid(output,1,strlen(output)-1) ;removes a leading space that somehow got in there	
	RETURN, output[0]
END

FUNCTION read_atomic_mass,elem_str,full=full,debug=debug, verb=verb,z=z
	path=getenv('ATOMIC_PHYSICS_PATH')
	x=size(elem_str,/type)
	input=elem_str
	IF x EQ 2 OR x EQ 4 OR x EQ 5 THEN elem_str=read_atomic_name(elem_str)
	z=elem_str
;	IF keyword_set(z) THEN elem_str=read_atomic_name(elem_str)
	
	;filename='/home/mlreinke/idl/impurities/data/atomic_mass'
	;filename='/usr/local/cmod/idl/atomic_physics/atomic_mass'
	filename=path+'atomic_mass'
	openr,lun,filename,/get_lun
	ident_str='*Atomic Symbol*'
	data_str='*= '+elem_str

	output=-1
	stopsearch=-1
	startsearch=-1
	iso_cntr=0
	line=strarr(1)
	mass=fltarr(20)
	isocomp=fltarr(20)	

	readf,lun,line
	
	WHILE eof(lun) NE 1 AND stopsearch EQ -1 DO BEGIN
		IF strmatch(line, ident_str, /fold_case) THEN BEGIN
			IF strmatch(line, data_str, /fold_case) THEN BEGIN
				startsearch=1
				readf,lun,line
				readf,lun,line
				a=strpos(line, '=')+1
				b=strpos(line, '(')
				mass[iso_cntr]=float(strmid(line,a,(b-a)))
				readf,lun,line
				a=strpos(line, '=')+1
				b=strpos(line, '(')
				isocomp[iso_cntr]=float(strmid(line,a,(b-a)))
				IF keyword_set(verb) THEN print, mass[iso_cntr], isocomp[iso_cntr]
				iso_cntr+=1
			ENDIF ELSE IF startsearch EQ 1 THEN stopsearch=1
		ENDIF
		readf,lun,line
	ENDWHILE
	
	IF NOT keyword_set(full) THEN output=total(mass*isocomp/100.0) ELSE $
		output={mass:mass[where(mass NE 0)], frac:isocomp[where(mass NE 0)]/100.0}
	
	close,lun
	free_lun, lun
	elem_str=input
	IF keyword_set(debug) THEN stop
	
	RETURN, output
END

FUNCTION read_atomic_charge,elem_str, debug=debug, verb=verb
	path=getenv('ATOMIC_PHYSICS_PATH')

	;filename='/home/mlreinke/idl/impurities/data/atomic_mass'
	;filename='/usr/local/cmod/idl/atomic_physics/atomic_mass'
	filename=path+'atomic_mass'
	openr,lun,filename,/get_lun
	ident_str='*Atomic Number*'
	data_str='*= '+elem_str

	output=-1
	stopsearch=-1
	line=strarr(1)	

	readf,lun,line
	
	WHILE eof(lun) NE 1 AND stopsearch EQ -1 DO BEGIN
		IF strmatch(line, ident_str, /fold_case) THEN BEGIN
			a=strpos(line, '=')+1
			number=int(strmid(line,a,(strlen(line)-a)))
			readf,lun,line
			IF strmatch(line, data_str, /fold_case) THEN BEGIN
				stopsearch=1
				output=number
			ENDIF 
		ENDIF
		readf,lun,line
	ENDWHILE
	
	close,lun
	free_lun, lun
	IF keyword_set(debug) THEN stop
	
	RETURN, output[0]
END
