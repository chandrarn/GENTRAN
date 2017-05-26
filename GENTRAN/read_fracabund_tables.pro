;+
;NAME:
;	READ_FQ_DATA
;
;PURPOSE:
;	The purpose of this function is load fractional abundance data for 1 <= Z <= 28
;	for Te > 1.0 (to fully stripped).  Data is loaded directly from 'table2.dat'
;	referenced in Mazzotta 'Astronomy & Astrophysics Supplement Series' v133 p403
;	1998 without interpolation.  The closest value Te is used.  Fournier tables are used
;	for Z=18,36,42
;
;CALLING SEQUENCE:
;	result=READ_FQ_DATA(z,te)
;
;INPUTS:
;	z:	INT nuclear charge to identify the element
;	te:	FLT electron temperature in [eV]
;
;KEYWORD PARAMETERS:
;	quiet:	/quiet to supress error messages to terminal
;	debug:	/debug to stop after the line has been identified and decompiled 
;		(use print, line to print the line from the table that was used)
;	
;OUTPUTS:
;	result:	STRUC that holds the data read from tables
;		*.q 	INTARR of ionization states
;		*.fq	FLTARR of fractional abundances
;		*.te	FLT the actual Te that was given in the table [eV]
;
;OPTIONAL OUTPUTS:
;	result=-1	Z(> 28 or < 1 if you're a tool) Te (< 1 eV) is out of range
;	result=-2	Te requested above the upper range of the table, assume fully stripped
;
;RESTRICTIONS:
;	This code requiers MLR_FUNCTIONS and READ_ATOMIC_NAME.  Use @genie_ini.bat to load the
;	correct depenencies.
;
;MODIFICATION HISTORY
;	Written by: 	ML Reinke, 9-03-05, I can't believe I didn't do this sooner
;	6/2/08:		ML Reinke - Added Fournier's Ar,Kr data and changed the formatting of the Mo data.
;	6/8/12:		ML Reinke - moved the data and pathing to /usr/local/cmod/idl/atomic_physics	
;
;-


FUNCTION read_fq_data,z,te,quiet=quiet,debug=debug
	IF z GT 28 OR te LT 0.85 AND z NE 42 AND z NE 36  THEN BEGIN
		IF NOT keyword_set(quiet) THEN print, 'Z (1-28) or Te (> 0.85 eV) out of range '
		RETURN,-1
	ENDIF
        IF z LE 28 AND z NE 18 THEN BEGIN
                ;filename='/home/mlreinke/idl/impurities/data/table2.dat'
		filename='/usr/local/cmod/idl/atomic_physics/table2.dat'
		openr,lun,filename,/get_lun
		zname=read_atomic_name(z)
		ident_str='*'+zname+' *'
		fail_str='*'+read_atomic_name(z+1)+' *'
		stopsearch=0
		log_te_pt_old=0	
		output=0
		line=strarr(1)
		kb=1.38065e-23	;J/K
		e=1.609e-19	;J/eV
		log_te=alog10(te*e/kb)

		readf,lun,line
		WHILE NOT stopsearch AND NOT eof(lun) DO BEGIN
			IF strmatch(line, fail_str, /fold_case) THEN BEGIN
				IF NOT keyword_set(quiet) THEN print, 'Te above upper range, assume fully stripped'
				RETURN,-2
			ENDIF
			IF strmatch(line, ident_str, /fold_case) THEN BEGIN
				log_te_pt=float(strmid(line,3,4))		
				IF log_te_pt GE log_te AND log_te_pt_old LE log_te THEN BEGIN
					IF abs(log_te_pt-log_te) GT abs(log_te_pt_old-log_te) THEN BEGIN
						line=line_old
						log_te_pt=log_te_pt_old
					ENDIF
					val=strsplit(line, '-', /extract)
					tmp=strsplit(line, '-')
					q0=(tmp[1]-1-7)/6		;6 pts in '-x.xxx', 7 in 'ZZ x.xx', 1 for the money	
					q=indgen(n_elements(tmp)-1)+q0[0]
					fq=10^(-1*float(val[1:n(val)]))
					stopsearch=1
					IF keyword_set(debug) THEN stop			
				ENDIF
				log_te_pt_old=log_te_pt
				line_old=line
			ENDIF
			readf,lun,line
		ENDWHILE	

		close,lun
		free_lun,lun
		output={q:q, fq:fq,te:10^(log_te_pt)*kb/e}
	ENDIF

	IF z EQ 18 THEN BEGIN
		;filename='/home/mlreinke/idl/impurities/data/ar_charge_state.sav'
		filename='/usr/local/cmod/idl/atomic_physics/ar_charge_state.sav'
		restore,filename
		q=indgen(19)
		output={q:q,fq:frac,te:te}
	ENDIF

	IF z EQ 36 THEN BEGIN
		;filename='/home/mlreinke/idl/impurities/data/kr_coronal.sav'
		filename='/usr/local/cmod/idl/atomic_physics/kr_coronal.sav'
		restore, filename
		q=indgen(37)
		output={q:q,fq:frac,te:te}
	ENDIF

	IF z EQ 42 THEN BEGIN
		;filename='/home/mlreinke/idl/impurities/data/frac_abund_Mo.sav'
		filename='/usr/local/cmod/idl/atomic_physics/frac_abund_Mo.sav'
		restore, filename
		q=indgen(43)
		output={q:q,fq:frac,te:te*1.0e3}
	ENDIF

	RETURN,output
END


;This will get the vector of Te(Z) that the is used in the data table.
FUNCTION read_te_vector,z,quiet=quiet
	IF z GT 28 THEN BEGIN
		IF NOT keyword_set(quiet) THEN print, 'Z (1-28) out of range '
		RETURN,-1
	ENDIF
	;filename='/home/mlreinke/idl/impurities/data/table2.dat'
	filename='/usr/local/cmod/idl/atomic_physics/table2.dat'
	openr,lun,filename,/get_lun
	zname=read_atomic_name(z)
	ident_str='*'+zname+' *'
	fail_str='*'+read_atomic_name(z+1)+' *'
	stopsearch=0
	cntr=0
	te_vec=fltarr(100)	;oversize and trim down
	line=strarr(1)
	kb=1.38065e-23	;J/K
	e=1.609e-19	;J/eV

	readf,lun,line
	WHILE NOT stopsearch AND NOT eof(lun) DO BEGIN
		IF strmatch(line, ident_str, /fold_case) THEN BEGIN
			log_te_pt=float(strmid(line,3,4))		
			te_vec[cntr]=10^(log_te_pt)*kb/e
			cntr+=1
		ENDIF
		readf,lun,line
		IF strmatch(line, fail_str, /fold_case) THEN stopsearch=1
	ENDWHILE
	close,lun
	free_lun,lun
	output=te_vec[0:cntr-1]
	RETURN,output
END

;+
;NAME:
;	FQ
;
;PURPOSE:
;	This function reads data from the frac abund tables using READ_FQ_DATA and reformats it into easy
;	to manage arrays
;
;CALLING SEQUENCE:
;	FQ,z,frac,te
;
;INPUTS:
;	z	INT	If you can't figure this out, go directly to jail and do not pass go (
;
;OUTPUTS:
;	frac:	FLTARR	[z+1,n_te] of the fractional abundance.  frac[i,*] is for the ionization state so
;			frac[0,*] is neutral and frac[z,*] is fully stripped.
;	te:	FLTARR	[n_te] of the electron temperature values from the table [eV]
;
;OPTIONAL OUTPUTS:
;	If z is outside the range of READ_FQ_DATA then frac and te will be returned as -1
;
;MODIFICATION HISTORY:
;	Written by:	ML Reinke - 8/16/07 (why the hell didn't I write this sooner?)
;	6/2/08:		ML Reinke - made the Fournier tables read out seperately from the Mazzotta data
;					for better Te gridding
;
;-

PRO fq,z,frac,te
	IF z LT 28 AND z NE 18 THEN BEGIN
		te=read_te_vector(z)
		n_te=n(te)+1
		IF te[0] EQ -1 THEN BEGIN
			frac=-1
			RETURN
		ENDIF
		frac=fltarr(z+1,n_te)
		FOR i=0,n_te-1 DO BEGIN
			out=read_fq_data(z,te[i])
			frac[out.q,i]=out.fq
		ENDFOR
	ENDIF ELSE BEGIN
		out=read_fq_data(z,1.0)
		te=out.te
		frac=out.fq
	ENDELSE
END
	


;+
;NAME:
;	FQPLT
;
;PURPOSE:
;	This procedure will use READ_FQ_DATA to create a plot of fractional abundance
;	versus electron temperature for 1 <= Z <= 28.  The resolution of the plots are
;	limited to the Te meshing of 'table2.dat' so these plots aren't prime-time plots
;	but are good for looking at what charge state would dominate in a particular
;	Te range.  Noble gas-like charge states are easily observable.
;
;CALLING SEQUENCE:
;	FQPLT,z
;
;INTPUS:
;	z:	nuclear charge of the element of interest
;
;OPTIONAL INPUTS:
;	tr:	FLTARR te=[te0, te1] sets the te plotting range [DEFAULT [min,max]]
;
;KEYWORD PARAMETERS:
;	ps:	/ps to use the postscript defice (setup using ps_plot_3 in MLR_FUNCTIONS)
;		file is saved to /home/username/idl/plots
;
;OUTPUTS
;	This procedure, when /ps isn't invoked, will take over window 20 
;
;SIDE EFFECTS:
;	Your PS and X-Win devices will be reconfigured as determined by the scripts XWPLOT
;	and PS_PLOT_3 in /home/mlreinke/idl/general/mlr_functions.  Once again, I am the color
;	table tyrant.
;
;RESTRICTIONS:
;	All restrictions of READ_FQ_DATA apply here as well.
;
;MODIFICATION HISTORY:
;	Written by: 	ML Reinke, 9-03-05
;	6/2/08:		ML Reinke - forced a default tr and made data load through the FQ procedure,added Kr
;-


PRO fqplt,z,ps=ps,tr=tr,ylog=ylog
	IF NOT keyword_set(tr) THEN tr=[10,1.0e4]
	IF z GT 28 OR z LT 1 AND z NE 42 AND z NE 36 THEN BEGIN
		print, 'Z out of range'
		RETURN
	ENDIF
	IF keyword_set(ps) THEN psplot_3 ELSE BEGIN $
		device, window_state=var
		IF var[20] EQ 0 THEN BEGIN
			xwplot
			xwplot
			window,20,xsize=1025,ysize=830,xpos=550,ypos=285, title='FracAbund,20'
		ENDIF
	ENDELSE

	;fill data
	fq,z,fq,te_vec
	tr[1]=max(te_vec) < tr[1]
	IF z NE 42 THEN tr[0]=min(te_vec) < tr[0]

	;plot
	lines=[2,3,4]
	xtit='Te [eV]'
	ytit='Fractional Abundance'
	tit='f_q(Te) for '+read_atomic_name(z)+' (Z='+num2str(z,1)+')'
	IF keyword_set(tr) THEN xr=tr ELSE xr=[min(te_vec),max(te_vec)]
	if keyword_set(ylog) then begin
		fqmn=min(fq[where(fq gt 0.0)])
		print,'fqmn='+num2str(fqmn)
		yr=[fqmn,1.1]
	endif else yr=[0,1.1]
	plot,[0,0],[0,0],xr=xr,yr=yr,xtit=xtit,ytit=ytit,tit=tit,/xsty, /ysty,chars=1.3,/xlog,ylog=ylog
	lstylelist=[0,2,3]
	l_cntr=0
	FOR i=0,z DO BEGIN
		IF l_cntr EQ n(lstylelist)+1 THEN l_cntr=0
		oplot, te_vec,fq[i,*],linestyle=lstylelist[l_cntr]
		mloc=maxloc(reform(fq[i,*]))
		xyouts,te_vec[mloc],fq[i,mloc]*1.01,num2str(i,1),color=100,charsize=1.35
		l_cntr+=1
	ENDFOR	

	IF keyword_set(ps) THEN BEGIN
		device, /close
		spawn, 'cp idl.ps /home/'+logname()+'/idl/plots/fqplt.ps'
		xwplot
	ENDIF

END
