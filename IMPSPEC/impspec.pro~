;+
;NAME:
;	READ_ISPEC_PATH
;
;PURPOSE:
;	This function returns the path to the ISPEC file for 
;	a given impurity, Z.
;
;CALLING SEQUENCE:
;	result=READ_ISPEC_PATH(z)
;	
;INPUTS:
;	z	INT	atomic number of impurity 
;
;OUTPUTS:
;	result	STRING	path
;
;PROCEDURE:
;	Currently, all ISPEC files are kept in /home/mlreinke/idl/vuv/
;	during the initial testing and development of IMPSPEC
;
;	If the Z is not on file, then path='none' will be returned
;
;MODIFICATION HISTORY:
;	Written by:	M.L. Reinke - 3/2012
;	6/8/12		M.L. Reinke - modified the ispec paths to usr /usr/local/cmod/idl/GENIE/IMPSPEC/ispec/
;	7/24/14		M.L. Reinke - added path for W.ispec
;	10/19/15	M.L. Reinke - modified to read ispec files from GENIE_PATH
;
;-

FUNCTION read_ispec_path,z
	gpath=getenv('GENIE_PATH')
	IF gpath EQ '' THEN gpath='/usr/local/cmod/idl/GENIE/'	;default to C-Mod path
	CASE z OF 
		5  : spath='IMPSPEC/ispec/B.ispec'
		7  : spath='IMPSPEC/ispec/N.ispec'
		8  : spath='IMPSPEC/ispec/O.ispec'
		9  : spath='IMPSPEC/ispec/F.ispec'
		10 : spath='IMPSPEC/ispec/Ne.ispec'
		18 : spath='IMPSPEC/ispec/Ar.ispec'
		42 : spath='IMPSPEC/ispec/Mo.ispec'
		74 : spath='IMPSPEC/ispec/W.ispec'
		ELSE : spath='none'
	ENDCASE
	path=gpath+spath
	RETURN,path
END

;+
;NAME:
;	READ_ISPEC_FILE
;
;PURPOSE:
;	This function reads an ASCII file in the "ISPEC" format that
;	can be used in the IMPSPEC line-fitting code.
;
;CALLING SEQUENCE:
;	result=READ_ISPEC_FILE(path)
;
;INPUTS:
;	path	STRING	of the path to the ISPEC file
;
;OPTIONAL INPUTS:
;	z	INT	of the atomic number of interest (calls READ_ISPEC_PATH)
;
;OUTPUTS:
;	result	STRUC	"ISPEC" file formatted for use in other IMPSPEC codse
;		*.LINE#		STRUC	information describing how to fit LINE number "#"
;			*.DLAM		FLTARR	[lam0,lam1] to truncate wavelength range
;			*.SPEC		STRING	spectrometer identifier (i.e. 'xeus','loweus',etc)
;			*.INST		FLOAT	maximum line width for (< 0 specifies a lorentzian, > 0 gaussian)
;			*.ILINES	INT	number of spectral lines in the group (index=0 is line of interest)
;			*.LAM		FLTARR	[ilines] center wavelengths
;			*.Z		INTARR	[ilines] atomic number
;			*.LABEL		STRARR	[ilines] line label in spectroscopy notation
;			*.ISO		STRARR	[ilines] isoelectronic element symbol
;		*.ELEM		STRING	elemental symbol
;		*.Z		INT	atomic number
;		*.NLINES	INT	number of lines specified in the file
;
;MODIFICATION HISTORY:
;	Written by:	M.L. Reinke - 3/12
;
;-
PRO readfn,lun,line,str=str ; reads with reads if str is set, readf in not set
 if keyword_set(str) then begin
  line=lun[0]
  if n(lun) ne 0 then lun=lun[1:n(lun)] else lun=''
 endif else readf,lun,line
END
FUNCTION eoffn,lun,str=str
 if keyword_set(str) then return,(lun[0] eq '') and (n(lun) eq 0) else return,eof(lun)
END

FUNCTION read_ispec_file,path,z=z,str=str
; ANJ working, need to conditionally change readf to reads when str keyword is passed (ie path contains file text instead of path to text)
	IF keyword_set(z) THEN $
         if z ne 0 $
          then path=read_ispec_path(z) $
          else str=1
	IF path[0] EQ 'none' THEN RETURN,-1
	IF keyword_set(str) THEN lun=strsplit(path,10B,/extract) ELSE openr,lun,path,/get_lun
	line=strarr(1)

	;HEADER
	;-------------------------
	readfn,lun,line,str=str & tmp=strsplit(line,',',/extract)
	 elem=tmp[0]
	 z=int(tmp[1])
	readfn,lun,line,str=str & tmp=strsplit(line,'=',/extract)
	 nlines=int(tmp[1])
	readfn,lun,line,str=str
	WHILE eoffn(lun,str=str) NE 1 $
              AND strmatch(line,'*#DATA#*',/fold_case) EQ 0 $
              DO readfn,lun,line,str=str
	if keyword_set(str) then data_start=lun else point_lun,-lun,data_start ; -ve lun gets pointer for data_start

	;DATA
	;-------------------------
	FOR i=0,nlines-1 DO BEGIN
		if keyword_set(str) then lun=data_start else point_lun,lun,data_start ; +ve lun sets pointer to data_start return to beginning of data
		WHILE eoffn(lun,str=str) NE 1 $
                      AND strmatch(line,'line='+num2str(i,1),/fold_case) EQ 0 $
                      DO readfn,lun,line,str=str ;search for line=i
		readfn,lun,line,str=str & tmp=strsplit(line,'=',/extract)
		 spec=tmp[1]
		readfn,lun,line,str=str & tmp=strsplit(line,',',/extract)
		 dlam=float(tmp)
		readfn,lun,line,str=str & tmp=strsplit(line,'=',/extract)
		 inst=float(tmp[1])
		readfn,lun,line,str=str & tmp=strsplit(line,'=',/extract)
		 ilines=int(tmp[1])
		 lam=fltarr(ilines)
		 iz=intarr(ilines)
		 label=strarr(ilines)
		 iso=strarr(ilines)
		FOR j=0,ilines-1 DO BEGIN
			readfn,lun,line,str=str & tmp=strsplit(line,',',/extract)
			 lam[j]=float(tmp[0])
			 label[j]=tmp[1]
			 iso[j]=tmp[2]
			 tmp=strsplit(label[j],'  ',/extract)
			 iz[j]=read_atomic_charge(tmp[0])
                ENDFOR
		result=execute('line'+num2str(i,1)+'={dlam:dlam,spec:spec,inst:inst,ilines:ilines,lam:lam,z:iz,label:label,iso:iso}')
	ENDFOR
	ispec={elem:elem,z:z,nlines:nlines}
	FOR i=0,nlines-1 DO BEGIN
		j=nlines-1-i
		addstr="'line"+num2str(j,1)+"',line"+num2str(j,1)
		result=execute("ispec=create_struct("+addstr+",ispec)")
        ENDFOR
	if not keyword_set(str) then begin
	 close,lun
	 free_lun,lun
	endif

	return,ispec
END

;+
;NAME:
;	ISPEC2STRING
;
;PURPOSE:
;	This function converts the ISPEC ASCII file into a single
;	string for storage into the tree.
;
;CALLING SEQUENCE:
;	result=ISPEC2STRING(path)
;
;INPUTS:
;	path	STRING	of the path to the ISPEC file
;
;OPTIONAL INPUTS:
;	z	INT	of the atomic number of interest (calls READ_ISPEC_PATH)
;
;OUTPUTS:
;	result	STRING	of the ISPEC files with string(10B) inserted
;			for carriage returns.
;
;MODIFICATION HISTORY:
;	Written by:	M.L. Reinke - 3/12
;
;-

FUNCTION ispec2string,path,z=z
	IF keyword_set(z) THEN path=read_ispec_path(z)
	openr,lun,path,/get_lun
	line=strarr(1)
	readf,lun,line
	string=line+string(10B)
	WHILE eof(lun) NE 1 DO BEGIN
		readf,lun,line
		string=string+line+string(10B)
	ENDWHILE
	close,lun
	free_lun,lun
	RETURN,string[0]
END

;+
;NAME:
;	READ_ISPEC_TREE
;
;PURPOSE:
;	This function reads data stored in the tree to create
;	the ISPEC structure
;
;CALLING SEQUECNE:
;	result=READ_ISPEC_TREE(shot,z)
;
;INPUTS:
;	shot	LONG	shot number
;	z	INT	atomic number of element
;
;OUTPUTS:
;	result	STRUC	ISPEC structure as described in READ_ISPEC_FILE
;			A value of -1 is returned for an error
;PROCEDURE:
;	This loads the string stored in the tree and creates a ASCII
;	file so that READ_ISPEC_FILE can be used.  [/tmp/tmp.ispec]
;
;MODIFICATION HISTORY:
;	Written by:	M.L. Reinke - 3/12
;	M.L. Reinke 	5/22/12 - added a Z suffix to the tmp.ispec file so parallel z's can be run
;	M.L. Reinke	10/20/15 - modified to use environment variables for tree/nodes
;-

FUNCTION read_ispec_tree,shot,z
	tree=getenv('IMPSPEC_MDS_TREE')
	node=getenv('IMPSPEC_MDS_PATH')

	;read ISPEC string from tree
	zstr=read_atomic_name(z)
	mdsopen,tree,shot
	istr=mdsvalue(node+'.'+strupcase(zstr[0])+':ISPEC',/quiet,status=status)
	mdsclose,tree,shot
	IF NOT status THEN RETURN,-1

	ispec=read_ispec_file(istr,/str)

	;write to file /tmp/tmp_z.ispec
	;path='/tmp/tmp_'+num2str(z,1)+'.ispec'
	;openw,lun,path,/get_lun
	;lines=strsplit(istr,string(10B),/extract)
	;FOR i=0,n(lines) DO printf,lun,lines[i]
	;close,lun
	;free_lun,lun	

	;send path to read_ispec_file
	;ispec=read_ispec_file(path)	

	RETURN,ispec
END

;+
;NAME:
;	IS_SPEC
;
;PURPOSE:
;	This function is used to check if a given shot has been setup
;	for a given Z for use with IMPSPEC.
;
;CALLING SEQUENCE:
;	result=IS_SPEC(shot,z)
;	
;MODIFICATION HISTORY:
;	M.L. Reinke	10/20/15 - modified to use environment variables for tree/nodes
;-

FUNCTION is_ispec,shot,z
	tree=getenv('IMPSPEC_MDS_TREE')
	node=getenv('IMPSPEC_MDS_PATH')
	mdsopen,tree,shot
	nchk=n(z)+1
	chk=intarr(nchk)
	FOR i=0,nchk-1 DO BEGIN
		zstr=read_atomic_name(z[i])
		dummy=mdsvalue(node+'.'+strupcase(zstr[0])+':ISPEC',/quiet,status=status)
		IF status THEN chk[i]=1
        ENDFOR
	mdsclose,tree,shot
	RETURN,chk

END

PRO copy_ispec,tshot,fshot,z=z
	tree=getenv('IMPSPEC_MDS_TREE')
	node=getenv('IMPSPEC_MDS_PATH')
	FOR i=0,n(z) DO BEGIN
		zstr=read_atomic_name(z[i])
		mdsopen,tree,fshot
		istr=mdsvalue(node+'.'+strupcase(zstr[0])+':ISPEC')
		mdsclose,tree,fshot
		mdsopen,tree,tshot
		mdsput,node+'.'+strupcase(zstr[0])+':ISPEC','$',istr
		mdsclose,tree,tshot	
	ENDFOR
END

;+
;NAME:
;	ADD_IMPSPEC
;
;PURPOSE:
;	This procedure adds the nodes to a given shot that allow
;	IMPSPEC to be run
;
;CALLING SEQUENCE:
;	ADD_IMPSPEC,shot
;
;INPUTS:
;	shot	LONG	shot number
;	
;OPTIONAL INPUTS:
;	z	INTARR	of the the atomic numbers to add
;
;KEYWORD PARAMETERS:
;	force	/force will delete a given element subnode before writing
;	rm	/rm will delete the \SPECTROSCOPY::TOP.IMPSPEC node 
;
;PROCEDURE:
;	The default behavior is to simply add the 
;	\SPECTROSCOPY::TOP.IMPSPEC node if it is not already there for a
;	given shot.  If the z optional input is specified then both
;	WRITE_ISPEC2NODE and WRITE_ISPEC2TREE are run for that shot
;
;MODIFICATION HISTORY:
;	Written by:	M.L. Reinke - 3/12
;	7/3/12		M.L. Reinke - added the PROF structure and nodes for use with GENTRAN_WRITE2TREE
;	10/20/15	M.L. Reinke - modified to use environment variables for tree/nodes
;
;-

PRO add_impspec,shot,z=z,force=force,rm=rm
	tree=getenv('IMPSPEC_MDS_TREE')
	node=getenv('IMPSPEC_MDS_PATH')
	mdsopen,tree,shot
	dummy=mdsvalue(node+'.N:ISPEC',/quiet,status=status)
	mdsclose,tree,shot
	IF status AND keyword_set(rm) THEN BEGIN						;removes the IMPSPEC directory
		mdstcl, "set verify"
		mdstcl, 'edit spectroscopy /shot='+num2str(shot,1)
		mdstcl, 'DELETE NODE '+node
		mdstcl, 'write'
		mdstcl, 'close' 	
	ENDIF
	IF NOT status THEN BEGIN
		mdstcl, "set verify"
		mdstcl, 'edit '+tree+' /shot='+num2str(shot,1)
		mdstcl, 'ADD NODE '+node
		mdstcl, 'ADD NODE '+node+'.PROF'
		mdstcl, 'ADD NODE/USAGE=SIGNAL '+node+'.PROF:DENS'
		mdstcl, 'ADD NODE/USAGE=SIGNAL '+node+'.PROF:TEMP'
		mdstcl, 'ADD NODE/USAGE=SIGNAL '+node+'.PROF:NEUT'
		mdstcl, 'ADD NODE/USAGE=TEXT '+node+'.PROF:VERSION'
		mdstcl, 'write'
		mdstcl, 'close' 
        ENDIF ELSE print, node+' already exists, use /rm to delete'
	IF keyword_set(z) AND NOT keyword_set(rm) THEN BEGIN
		FOR i=0,n(z) DO BEGIN
			zstr=read_atomic_name(z[i])
			IF is_ispec(shot,z[i]) THEN BEGIN
				IF keyword_set(force) THEN BEGIN
					mdstcl, "set verify"
					mdstcl, 'edit spectroscopy /shot='+num2str(shot,1)
					mdstcl, 'DELETE NODE '+node+'.'+strupcase(zstr)
					mdstcl, 'write'
					mdstcl, 'close' 	
                                ENDIF ELSE print, node+'.'+strupcase(zstr)+' already exists, use /force to overwrite'
                        ENDIF ELSE write_ispec2node,shot,path,z=z[i]
			write_ispec2tree,shot,path,z=z[i]
                ENDFOR
	ENDIF
END


;+
;NAME:
;	APPEND_IMPSEPC_MODELING
;
;PURPOSE:
;	This procedure will add the IMSPEC.PROF nodes for use with GENTRAN_TENE_WRITE2TREE
;	in order to modify pre-existing IMPSPEC trees
;
;-

PRO append_impspec_modeling,shot,rm=rm
	tree=getenv('IMPSPEC_MDS_TREE')
	node=getenv('IMPSPEC_MDS_PATH')
	mdsopen,tree,shot
	dummy=mdsvalue(node+'.PROF.VERSION',/quiet,status=status)
	mdsclose,tree,shot

	IF status AND keyword_set(rm) THEN BEGIN						;removes the IMPSPEC directory
		mdstcl, "set verify"
		mdstcl, 'edit '+tree+' /shot='+num2str(shot,1)
		mdstcl, 'DELETE NODE '+node+'.PROF'
		mdstcl, 'write'
		mdstcl, 'close' 	
	ENDIF
	IF NOT status THEN BEGIN
		mdstcl, "set verify"
		mdstcl, 'edit '+tree+' /shot='+num2str(shot,1)
		mdstcl, 'ADD NODE '+node+'.PROF'
		mdstcl, 'ADD NODE/USAGE=SIGNAL '+node+'.PROF:DENS'
		mdstcl, 'ADD NODE/USAGE=SIGNAL '+node+'.PROF:TEMP'
		mdstcl, 'ADD NODE/USAGE=SIGNAL '+node+'.PROF:NEUT'
		mdstcl, 'ADD NODE/USAGE=TEXT '+node+'.PROF:VERSION'
		mdstcl, 'write'
		mdstcl, 'close' 
        ENDIF ELSE print, node+'.PROF already exists, use /rm to delete'
END

;+
;NAME:
;	WRITE_ISPEC2NODE
;
;PURPOSE:
;	This procedure takes an ISPEC file and creates a subnode under
;	\SPECTROSCOPY::TOP.IMPSPEC for that element, allowing
;	RUN_IMPSPEC to be used on that shot.
;
;CALLING SEQUENCE:
;	WRITE_ISPEC2NODE
;	
;INPUTS:
;	shot	LONG	shot number
;	path	STRING	of the path to the ISPEC file (sent to READ_ISPEC_FILE)
;		LONG	of the shot number copy from (sent to READ_ISPEC_TREE)
;
;OPTIONAL INPUTS:
;	z	INT	of the atomic number of interest (sent READ_ISPEC_PATH)		
;
;KEYWORD PARAMETERS:
;	rm	/rm will remove the existing IMPSPEC.Z node and overwrite with specified file
;	base 	/base will write only the base IMPSPEC nodes, not those for GENTRAN modeling (as if that will ever occur?)
;
;OUTPUTS:
;	Nodes are added to the tree based on the data in the ISPEC
;	file.  A subnode \SPECTROSCOPY::TOP.IMPSPEC.LINE# is added for each
;	line specified where the brightness, coefficients and information for
;	the IMPSPEC scope are located.
;
;MODIFICATION HISTORY:
;	Written by:	M.L. Reinke 3/12
;	7/3/12		M.L. Reinke - added the DIFF/CONV/CSDEN to ELEM structure and BR_MOD, NZ 
;					and ETA to LINE# nodes for use with GENTRAN_WRITE2TREE
;	8/1/15		M.L. Reinke - added the rm and base keywords
;  	10/20/15	M.L. Reinke - modified to use environment variables for tree/nodes                     
;
;-

PRO write_ispec2node,shot,path,z=z,rm=rm,base=base,force=force
	tree=getenv('IMPSPEC_MDS_TREE')
	node=getenv('IMPSPEC_MDS_PATH')
	IF size(path,/type) EQ 3 THEN ispec=read_ispec_tree(path,z) ELSE ispec=read_ispec_file(path,z=z)
	mdsopen,tree,shot
	dummy=mdsvalue(node+'.'+ispec.elem+':ISPEC',/quiet,status=status)
	mdsclose,tree,shot
	IF status AND keyword_set(rm) OR keyword_set(force) THEN BEGIN						;removes the IMPSPEC directory
		mdstcl, "set verify"
		mdstcl, 'edit '+tree+' /shot='+num2str(shot,1)
		mdstcl, 'DELETE NODE '+node+'.'+ispec.elem+' /confirm'
		mdstcl, 'write'
		mdstcl, 'close' 	
	ENDIF
	mdstcl, 'set verify'
	mdstcl, 'edit '+tree+' /shot='+num2str(shot,1)
	mdstcl, 'ADD NODE '+node+'.'+ispec.elem
	mdstcl, 'ADD NODE/USAGE=NUMERIC '+node+'.'+ispec.elem+':NLINES'
	mdstcl, 'ADD NODE/USAGE=TEXT '+node+'.'+ispec.elem+':ISPEC'
	mdstcl, 'ADD NODE/USAGE=ACTION '+node+'.'+ispec.elem+':RUN_ISPEC'
	IF NOT keyword_set(base) THEN BEGIN
		mdstcl, 'ADD NODE/USAGE=SIGNAL '+node+'.'+ispec.elem+':DIFF'
		mdstcl, 'ADD NODE/USAGE=SIGNAL '+node+'.'+ispec.elem+':CONV'
		mdstcl, 'ADD NODE/USAGE=SIGNAL '+node+'.'+ispec.elem+':CSDEN'
	ENDIF
	FOR i=0,ispec.nlines-1 DO BEGIN
		lpath=ispec.elem+'.LINE'+num2str(i,1)
		mdstcl, 'ADD NODE '+node+'.'+lpath
		mdstcl, 'ADD NODE/USAGE=SIGNAL '+node+'.'+lpath+':BR'
		IF NOT keyword_set(base) THEN mdstcl, 'ADD NODE/USAGE=SIGNAL '+node+'.'+lpath+':BR_MOD'
		IF NOT keyword_set(base) THEN mdstcl, 'ADD NODE/USAGE=SIGNAL '+node+'.'+lpath+':NZ'
		IF NOT keyword_set(base) THEN mdstcl, 'ADD NODE/USAGE=NUMERIC '+node+'.'+lpath+':ETA'
		mdstcl, 'ADD NODE/USAGE=SIGNAL '+node+'.'+lpath+':COEFS'
		mdstcl, 'ADD NODE/USAGE=NUMERIC '+node+'.'+lpath+':LAM'
		mdstcl, 'ADD NODE/USAGE=TEXT '+node+'.'+lpath+':LABEL'
		mdstcl, 'ADD NODE/USAGE=NUMERIC '+node+'.'+lpath+':DLAM'
        ENDFOR	
	mdstcl, 'write'
	mdstcl, 'close' 
END

;+
;NAME:
;	APPEND_ISPEC_MODELING
;
;PURPOSE:
;	This procedure will add the IMSPEC.XX.DIFF/CONV/CSDEN nodes for use with GENTRAN_WRITE2TREE
;	and the LINE#.BR_MOD, NZ and ETA.  This is used to modify pre-existing IMPSPEC trees
;
;-

PRO append_ispec_modeling,shot,path,z=z
	tree=getenv('IMPSPEC_MDS_TREE')
	node=getenv('IMPSPEC_MDS_PATH')
	IF size(path,/type) EQ 3 THEN ispec=read_ispec_tree(path,z) ELSE ispec=read_ispec_file(path,z=z)
	mdstcl, "set verify"
	mdstcl, 'edit '+tree+' /shot='+num2str(shot,1)
	mdstcl, 'ADD NODE/USAGE=SIGNAL '+node+'.'+ispec.elem+':DIFF'
	mdstcl, 'ADD NODE/USAGE=SIGNAL '+node+'.'+ispec.elem+':CONV'
	mdstcl, 'ADD NODE/USAGE=SIGNAL '+node+'.'+ispec.elem+':CSDEN'
	FOR i=0,ispec.nlines-1 DO BEGIN
		lpath=ispec.elem+'.LINE'+num2str(i,1)
		mdstcl, 'ADD NODE/USAGE=SIGNAL '+node+'.'+lpath+':BR_MOD'
		mdstcl, 'ADD NODE/USAGE=SIGNAL '+node+'.'+lpath+':NZ'
		mdstcl, 'ADD NODE/USAGE=NUMERIC '+node+'.'+lpath+':ETA'
        ENDFOR
	mdstcl, 'write'
	mdstcl, 'close' 
END


;+
;NAME:
;	WRITE_ISPEC2TREE
;
;PURPOSE:
;	This procedure writes and ISPEC file for a specified atomic number
;	to the tree for given shot.
;
;CALLING SEQUENCE:
;	WRITE_ISPEC2TREE,shot,path
;
;INPUTS:
;	shot	LONG	shot number
;	path	STRING	of the path to the ISPEC file (sent to READ_ISPEC_FILE)
;		LONG	of the shot number copy from (sent to READ_ISPEC_TREE)
;
;OPTIONAL INPUTS:
;	z	INT	of the atomic number of interest (sent READ_ISPEC_PATH)		
;
;MODIFICATION HISTORY:
;	Written by	M.L. Reinke - 3/12
;       fixed path overwrite when passing a config file and no z ANJ 20120615
;	7/14/12 	M.L. Reinke - fixed bug where read_ispec_tree was using optional input of z instead of input
;  	10/20/15	M.L. Reinke - modified to use environment variables for tree/nodes 
;
;-

PRO write_ispec2tree,shot,path,z=z
	tree=getenv('IMPSPEC_MDS_TREE')
	node=getenv('IMPSPEC_MDS_PATH')
	;IF size(path,/type) EQ 3 THEN ispec=read_ispec_tree(path,z) ELSE ispec=read_ispec_file(path,z=z)
	IF size(path,/type) NE 7 THEN BEGIN ; type 7 is a string
	 ispec=read_ispec_tree(path,z) 
	 istr=ispec2string(path,z=z)
        ENDIF ELSE BEGIN
	 ispec=read_ispec_file(path,z=z)
	 istr=ispec2string(path) ; don't pass z or path is overwritten...
	ENDELSE

	mdsopen,tree,shot
	mdsput,node+'.'+ispec.elem+':ISPEC','$',istr
	mdsclose,tree,shot
END

;+
;NAME:
;	LORENTZIAN_INTEG FITS
;
;PURPOSE:
;	This function calculates the lorentzian line profile for the sum of 
;	an arbitrary number of lines plus an optional DC offset.  The format
;	is made to be used with MPFIT. Each point is calculated as the integral 
;       of a lorentzian from midway to neighboring points.
;
;CALLING_SEQUENCE:
;	result=LORENTZIAN_INTEG_FITS(x,p)
;
;INPUTS:
;	x	FLTARR	[n_x] of the points for which to calculate the spectra
;	p	FLTARR 	[n_gauss*3+nb] where nb will determine the order of the baseline fit
;
;OUTPUTS:
;	result:	FLTARR  [n_x] of the sum of all the lines at each point
;
;PROCEDURE:
;	Lorentzians are specified as a/pi*c/2/((x-b)^2+(c/2)^2), for the amplitude a (height of the line, not area), center b, and FWHM c.
;	The integral is then a difference of arctangents:
;         int_x1^x2 lor dx = a/pi*( atan((x1-b)*2/c) - atan((x2-b)*2/c) )
;       If the parameter values of p are a double precision array, then the result is returned as a double.
;	
;	IF n_elements(p) MOD 3 = 0 then the last three are assume to be a quadraditc baseline
;	IF n_elements(p) MOD 3 = 1 then the last three are assume to be a constant baseline
;	IF n_elements(p) MOD 3 = 2 then the last three are assume to be a linear baseline
;
;MOFICATION HISTORY:
;	Written by: 	A. N. James 20120712
;                       copied from hirexsr_fit_spectra.pro in THACO
;	3/12		modified the baseline reference to use the
;			first element of the x scaling vector
;
;-

FUNCTION lorentzian_integ_fits, x, p,base=base
	CASE n_elements(p) MOD 3 OF
		0 : BEGIN
			IF n_elements(p) EQ 3 THEN BEGIN ;no baseline
				n_line=1
				basen=3
			ENDIF ELSE BEGIN		;quadradic baseline
				n_line=n_elements(p)/3-1
				basen=2
			ENDELSE
		END
		1 : BEGIN		;constant baseline
			n_line=n_elements(p)/3
			basen=0
		END

		2 : BEGIN		;linear baseline
			n_line=n_elements(p)/3
			basen=1	
		END
	ENDCASE
			
	type=size(p,/type)
	IF type EQ 5 THEN L = dblarr(3,n_line) ELSE L=fltarr(3,n_line)
	FOR i = 0, n_line-1 DO FOR j = 0,2 DO L[j,i] = p[3*i+j]
	
	nx=n(x)+1
	y=fltarr(nx)
	dx=x(1:n(x)-1)-x(0:n(x)-1)
	dx=[dx[1],dx]

	;FOR i = 0, n_line-1 DO y += L[0,i]*exp(-(x-L[1,i])^2/(2.*L[2,i]^2))
	;FOR i = 0, n_line-1 DO y += L[0,i]/!PI*L[2,i]/2/( (x-L[1,i])^2 + (L[2,i]/2)^2 )
	FOR i = 0, n_line-1 DO y += L[0,i]*L[2,i]/dx*( atan((x+dx/2-L[1,i])/L[2,i]) - atan((x-dx/2-L[1,i])/L[2,i]) )
;         int_x1^x2 lor dx = a/pi*( atan((x1-b)*2/c) - atan((x2-b)*2/c) )
	
	x0=x[0]	
	CASE basen OF 
		3 : base = 0
		2 : base = p[n(p)] + p[n(p)-1]*(x-x0) + p[n(p)-2]*(x-x0)^2	;
		1 : base = p[n(p)] + p[n(p)-1]*(x-x0)				; 
		0 : base = p[n(p)]
		ELSE :	
	ENDCASE
	y+=base

	RETURN, y
END

;+
;NAME:
;	GAUSSIAN_FITS
;
;PURPOSE:
;	This function calculates the gaussian line profile for the sum of 
;	an arbitrary number of gaussians plus an optional DC offset.  The format
;	is made to be used with MPFIT
;
;CALLING_SEQUENCE:
;	result=GAUSSIAN_FITS(x,p)
;
;INPUTS:
;	x	FLTARR	[n_x] of the points for which to calculate the spectra
;	p	FLTARR 	[n_gauss*3+nb] where nb will determine the order of the baseline fit
;
;OUTPUTS:
;	result:	FLTARR  [n_x] of the sum of all the gaussians at each point
;
;PROCEDURE:
;	Gaussians are specified as a*exp(-(x-b)^2/(2*c^2)).  If the parameter values of p
;	are a double precision array, then the result is returned as a double.
;	
;	IF n_elements(p) MOD 3 = 0 then the last three are assume to be a quadraditc baseline
;	IF n_elements(p) MOD 3 = 1 then the last three are assume to be a constant baseline
;	IF n_elements(p) MOD 3 = 2 then the last three are assume to be a linear baseline
;
;MOFICATION HISTORY:
;	Written by: 	copied from hirexsr_fit_spectra.pro in THACO
;	3/12		modified the baseline reference to use the
;			first element of the x scaling vector
;
;-

FUNCTION gaussian_fits, x, p,base=base
	CASE n_elements(p) MOD 3 OF
		0 : BEGIN
			IF n_elements(p) EQ 3 THEN BEGIN
				n_line=1
				basen=3
			ENDIF ELSE BEGIN		;quadradic baseline
				n_line=n_elements(p)/3-1
				basen=2
			ENDELSE
		END
		1 : BEGIN		;constant baseline
			n_line=n_elements(p)/3
			basen=0
		END

		2 : BEGIN		;linear baseline
			n_line=n_elements(p)/3
			basen=1	
		END
	ENDCASE
			
	type=size(p,/type)
	IF type EQ 5 THEN L = dblarr(3,n_line) ELSE L=fltarr(3,n_line)

	FOR i = 0, n_line-1 DO FOR j = 0,2 DO L[j,i] = p[3*i+j]
	nx=n(x)+1
	x0=x[0]	
	y=fltarr(nx)
	FOR i = 0, n_line-1 DO y = y + L[0,i]*exp(-(x-L[1,i])^2/(2.*L[2,i]^2))
	CASE basen OF 
		3 : base=0
		2 : base=p[n(p)]+p[n(p)-1]*(x-x0)+p[n(p)-2]*(x-x0)^2	;
		1 : base=p[n(p)]+p[n(p)-1]*(x-x0)				; 
		0 : base=p[n(p)]
		ELSE :	
	ENDCASE
	y+=base

	RETURN, y
END

;+
;NAME:
;	ISPEC_FIT_LINE
;
;PURPOSE:
;	This function is a lower level code that actually does the
;	fitting of a spectrum based on the LINE information from the ISPEC file
;
;CALLING SEQUENCE:
;	result=IMPSPEC_FIT_LINE(spec,lam,sig,line)
;
;INPUTS:
;	spec	FLTARR	[nlam] of the spectral brightness vs. wavelength
;	lam	FLTARR	[nlam] of the wavelength values
;	sig	FLTARR	[nlam] of the photon statistic uncertainties in spec
;	line	STRUC	substructure of the ISPEC file (see READ_ISPEC_FILE)
;			(*.ilines specifies the # of lines in the fit)
;
;OPTIONAL INPUTS:
;	fitz	INTARR	[nz] of the atomic numbers to include in fit DEFAULT: all listed
;
;KEYWORD PARAMETERS:
;	plot	/plot will display a graphical output of the data, seed and fit
;
;OUTPUTS:
;	result	FLTARR	[ilines*3+2) of the fit coefficients along
;			linear background
;
;OPTIONAL OUTPUTS:
;	error	FLTARR	[ilines*3+2] of the uncertainty in fit coefficients
;
;PROCEDURE:
;	The line.dlam array defines the subset of spectrum that is
;	used in the fit.  Non-linear least-squares fitting is done using
;	MPFIT using the following constraints.
;		line intensitity and DC offset > 0
;		line width > 0 and < line.inst
;		all line widths are the same
;		all lines are tied to the main peak, shifting together
;		constrains shift to < 0.5 angstroms
;		
;	Those impurities which are not specified with the
;	fitz optional input have their peaks fixed at 0.0, thus removing them
;	from the fit.
;	
;MODIFICATION HISTORY:
;	Written by:	M.L. Reinke - 3/12
;	M.L. Reinke	9/6/2014 - modified initial guess for line to remove baseline to help fit low signal lines
;
;-

FUNCTION impspec_fit_line,spec,lam,sig,line,error=error,plot=plot,status=status,fitz=fitz

	tmp=where(lam GE line.dlam[0] AND lam LE line.dlam[1])
	y=spec[tmp]
	x=lam[tmp]
	ysig=sig[tmp]
	chk=where(ysig LE 0)
	IF chk[0] NE -1 THEN ysig[chk]=abs(mean(y))

	shift=0.5			;allowable shift (put in ISPEC file?)
	n_lines=line.ilines
	L=fltarr(3,n_lines)
	L[2,*]=line.inst
	L[1,*]=line.lam
	L[0]=y(ipt(x,L[1]))-min(y)
	base_line=min(y)
	estimate=fltarr(n_lines*3+2)
	FOR i = 0,n_lines-1 DO FOR j = 0,2 DO estimate[3*i+j] = L[j,i]
	estimate[n_lines*3+1]=base_line
	
	parinfo = replicate({value:0., fixed:0, limited:[0,0],limits:[0.,0], tied:''}, n_elements(estimate))
	parinfo[*].value=estimate

	;sets the lower intensity bound to be 0.0
	FOR i=0,n_lines-1 DO parinfo[3*i].limited[0]=1		;states that all intensity coefs will have a lower bound
	FOR i=0,n_lines-1 DO parinfo[3*i].limits[0]=0.0      	;sets the lower bound to be 0.0

	;sets the lower width bound to be 0.0
	FOR i=0,n_lines-1 DO parinfo[3*i+2].limited[0]=1	;all width coefs will have a lower bound
	FOR i=0,n_lines-1 DO parinfo[3*i+2].limits[0]=0.0      	;sets the lower bound to be 0.0

	;sets the upper width bound to be specified from the line file
	FOR i=0,n_lines-1 DO parinfo[3*i+2].limited[1]=1		;all width coefs will have an upper bound
	FOR i=0,n_lines-1 DO parinfo[3*i+2].limits[1]=line.inst      	;sets the upper bound

	;set the +/- shift to 0.5 Angstroms to prevent gross wander
	parinfo[1].limited[0]=1
	parinfo[1].limits[0]=line.lam[0]-shift
	parinfo[1].limited[1]=1
	parinfo[1].limits[1]=line.lam[0]+shift

	IF keyword_set(fitz) THEN BEGIN
		FOR i=1,n_lines-1 DO BEGIN
			chk=where(fitz EQ line.z[i])
			IF chk[0] EQ -1 THEN parinfo[3*i].fixed=1	;if not in the fitted-z list then fix background line to 0.0
		ENDFOR
	ENDIF

	;fixes DC offset to positive definite
	parinfo[n_lines*3+1].limited[0]=1
	parinfo[n_lines*3+1].limits[0]=0.0

	;tie all wavelength shifts together
	FOR i=1,n_lines-1 DO BEGIN
		IF line.lam[0]-line.lam[i] LT 0 THEN parinfo(3*i+1).tied = 'P(1)+'+strtrim(abs(line.lam[0]-line.lam[i]),1) ELSE $
			parinfo(3*i+1).tied = 'P(1)-'+strtrim(abs(line.lam[0]-line.lam[i]),1)
	ENDFOR

	;tie all widths  together
	FOR i=1,n_lines-1 DO parinfo(3*i+2).tied = 'P(2)'

	coefs = mpfitfun('gaussian_fits', x,y,ysig, estimate,perror=error,parinfo=parinfo,status=status,niter=niter,/quiet)
	IF keyword_set(plot) AND status NE 0 THEN BEGIN
		plot,x,y,/xst,/ysty,xtit='Wavelength [Ang]',ytit='Spectral Brightness [AU]',yr=[0,max(y)*1.1],tit=line.label[0],/nodata
		oploterror,x,y,ysig,psym=8
		xplot=make(x[0],last(x),1000)
		oplot,xplot,gaussian_fits(xplot,estimate),color=100,linestyle=1
		oplot,xplot,gaussian_fits(xplot,coefs),color=200
		FOR i=0,n_lines-1 DO oplot,xplot,gaussian_fits(xplot,coefs[i*3:(i+1)*3-1]),color=30,linestyle=2
		oplot,xplot,coefs[n_lines*3+1]+(xplot-xplot[0])*coefs[n_lines*3],color=30,linestyle=2
        ENDIF
	RETURN,coefs
END

;+
;NAME:
;	VUV_LOAD_SPEC
;
;PURPOSE:
;	This procedure loads the VUV spectrometer data from the tree
;	for a given shot
;
;CALLING SEQUENCE:
;	VUV_LOAD_SPEC,shot,spec,specbr,lam,time
;
;INPUTS:
;	shot	LONG 	shot number
;	spec	STRING	identifier of the spectrometer to be loaded
;
;OUTPUTS:
;	specbr	FLTARR	[nlam,ntime] of the spectral brightnes [ph/s/m^2/Ang] 
;	lam	FLTARR	[nlam] of the wavelengths [Angstroms]	
;	time	FLTARR	[ntime] of the time points [seconds]
;
;OPTIONAL OUTPUTS:
;	sigbr	FLTARR	[nlam,ntime] of the photon statistics error in specbr
;
;MODIFICATION HISTORY:
;	Written by:	M.L. Reinke (taken from VUV_TREE_ANALYSIS.PRO)
;	M.L. Reinke	10/20/15 - made to look for machine specific loading routines
;
;-

PRO vuv_load_spec,shot,spec,specbr,lam,time,sigbr=sigbr
	mach=getenv('MACHINE')
	CASE mach OF
		'cmod' : BEGIN
			CASE spec OF 
				'xeus' : path='\SPECTROSCOPY::TOP.XEUS:SPEC'
				'loweus' : path='\SPECTROSCOPY::TOP.LOWEUS:SPEC'
			ENDCASE
			mdsopen,'spectroscopy',shot
			specbr=mdsvalue('_sig='+path)
			lam=mdsvalue('dim_of(_sig,0)')
			time=mdsvalue('dim_of(_sig,1)')
			sigbr=mdsvalue('dim_of(_sig,2)')
			mdsclose,'spectroscopy',shot
                END

		'nstxu' : BEGIN
			CASE spec OF 
				'xeus' : path='\PSPEC_PC::TOP.XEUS:IMAGE'
				'loweus' : path='\PSPEC_PC::TOP.LOWEUS:IMAGE'
			ENDCASE
			mdsopen,'pspec_pc',shot
			ispecbr=float(mdsvalue('_sig='+path+':IMAGE',/quiet,status=status))
			lam=mdsvalue('dim_of(_sig,0)',/quiet)
			time=float(mdsvalue('dim_of(_sig,2)',/quiet))/1.0e3
			mdsclose,'pspec_pc',shot
			x=size(ispecbr)
			specbr=fltarr(x[1],x[3])
			FOR i=0,x[1]-1 DO FOR j=0,x[3]-1 DO specbr[i,j]=median(ispecbr[i,*,j])
			sigbr=specbr*0.0	;hardcode for now
                END

	ENDCASE

END

;+
;NAME:
;	CHROMEX_LOAD_SPEC
;
;-

PRO chromex_load_spec,shot,chr,raw=raw,lam_auto=lam_auto
	if not keyword_set(raw) then raw=0
	if not keyword_set(lam_auto) then lam_auto=0

	mdsopen,'spectroscopy',shot
	 IF NOT keyword_set(raw) THEN BEGIN
	  int=mdsvalue('\SPECTROSCOPY::TOP.CHROMEX.ANALYSIS:ANAL_SPECTRA')
	  lam=mdsvalue('dim_of(\SPECTROSCOPY::TOP.CHROMEX.ANALYSIS:ANAL_SPECTRA,0)')
	  ch=mdsvalue('dim_of(\SPECTROSCOPY::TOP.CHROMEX.ANALYSIS:ANAL_SPECTRA,1)')
	  t=mdsvalue('dim_of(\SPECTROSCOPY::TOP.CHROMEX.ANALYSIS:ANAL_SPECTRA,2)')
	 ENDIF ELSE BEGIN
	  int=mdsvalue('\SPECTROSCOPY::TOP.CHROMEX.ANALYSIS:RAW_SPECTRA')
	  lam=mdsvalue('\SPECTROSCOPY::TOP.CHROMEX.CALIBRATION:LAMBDA_TRUE')
	  dt = mdsvalue('\SPECTROSCOPY::TOP.CHROMEX.CAMAC:DELTA_TIME')
	  t_start=mdsvalue('\SPECTROSCOPY::TOP.CHROMEX.CAMAC:TIME_START')
	  t = t_start + dt * findgen(n_elements(int[*,4,4])) +dt/2
	  ch=indgen(n(int[0,0,*])+1)
	 ENDELSE
	 cwl=mdsvalue('\SPECTROSCOPY::TOP.CHROMEX.SETTINGS:LAMBDA_CTR')
	 cfg_per=mdsvalue('\SPECTROSCOPY::TOP.CHROMEX.SETTINGS:CFG_PERISCOP')
	 cfg_fib=mdsvalue('\SPECTROSCOPY::TOP.CHROMEX.SETTINGS:CFG_PER_FIBR')
	mdsclose,'spectroscopy',shot

	lab=strarr(16)
	FOR i=0,n(cfg_per) DO lab[i]=cfg_per[i]+'_'+num2str(cfg_fib[i],1)

	; automatic wavelength calibration to bypass occasional shifts
	shift=0
	if keyword_set(lam_auto) then begin 
		case cwl of ; use omnipresent Mo triplets to fit the spectra
			410 : begin
				lines=[379.8,386.4,390.2]
			end
			550 : begin
				lines=[550.649,553.305,557.045]
			end
			else : print, 'cannot do automatic wl cal for cwl '+num2str(cwl,1)
		endcase

		x=lam
		y=total(total( chr.int ,1),2) ; the spectrum averaged over all time
		tmp=where(t lt 0.02)
		ysig=total(total( chr.int(tmp,*,*) ,1),2)

		nl=n(lines)+1
		L=fltarr(3,nl)
		wid=.25 ; [nm] min width
		L[2,*]=wid
		L[1,*]=lines ; center
		L[0,*]=max(y); amplitude
		estimate=fltarr(nl*3+2)
		FOR i = 0,nl-1 DO $
		 FOR j = 0,2 DO $ 
		  estimate[3*i+j] = L[j,i]

		base_line=min(y)
		estimate[nl*3+1]=base_line
	
		parinfo = replicate({value:0., fixed:0, limited:[0,0],limits:[0.,0], tied:''}, n_elements(estimate))
		parinfo[*].value=estimate

		FOR i=0,nl-1 DO BEGIN
			parinfo[3*i].limited=[1,	1] ; states that all intensity coefs will have upper and lower bounds
			parinfo[3*i].limits =[0.0,	max(y)*1.5]
		
			mxshift=2 ; allowable shift in nm
			parinfo[3*i+1].limited=[1,1] ; set the +/- shift to prevent gross wander
			parinfo[3*i+1].limits =lines[0] + mxshift*[-1,1]
			IF lines[0]-lines[i] LT 0 $ ; tie all wavelength shifts together
			 THEN parinfo(3*i+1).tied = 'P(1)+'+strtrim(abs(lines[0]-lines[i]),1) $
			 ELSE parinfo(3*i+1).tied = 'P(1)-'+strtrim(abs(lines[0]-lines[i]),1)
		
			parinfo[3*i+2].limited=[1,	1]	; all width coefs will have upper and lower bounds
			parinfo[3*i+2].limits =[wid,	5*wid]
			parinfo(3*i+2).tied   ='P(2)' ; tie all widths together
		ENDFOR

		;fixes DC offset to positive definite
		parinfo[nl*3+1].limited=[1,	1]
		parinfo[nl*3+1].limits= [0.0,	mean(y)]

		coefs = mpfitfun('lorentzian_integ_fits', x,y,ysig, estimate,perror=error,parinfo=parinfo,status=status,niter=niter,/quiet)
		shift=lines[0]-coefs[1]
		print,'auto shifted chromex wavelength by '+num2str(shift,1),coefs

		print,status
		plot,x,y
		oplot,x,lorentzian_integ_fits(x,estimate),color=96, linestyle=1
		FOR i=0,nl-1 DO $
                 oplot,x,lorentzian_integ_fits(x,coefs[i*3:(i+1)*3-1]),color=164,linestyle=2
	endif

	chr={shot:shot,int:int,lam:lam,ch:ch,t:t,lab:lab,shift:shift}
	stop
END

;+
;NAME:
;	IMPSPEC
;
;PURPOSE:
;	This is a high level procedure which takes in a shot number
;	and an ISPEC file and computes fit and calculates line brightness
;
;CALLING SEQUENCE:
;	IMPSPEC,shot,ispec,br,coefs
;
;INPUTS:
;	shot		LONG	shot number
;	ispec		STRUC	formatted as per READ_ISPEC_FILE
;	
;OPTIONAL INPUTS:
;	kline		INT	select one line from ISPEC to run (use -1 to select 0th line)
;	fitz		INTARR	of the Z's to include in fit (sent to IMPSPEC_FIT_LINE)
;
;KEYWORD PARAMETERS:
;	plot	/plot sent to IMPSOEC_FIT_LINE to plot results of fit
;	verb	/verb not used yet
;	debug	/debug stops the code at the end of each line fitting
;
;OUTPUTS:
;	br	STRUC	containing the line brightness
;		*.LINE#		FLTARR	[ntime,3] of the brightnes [*,0],
;					time [*,1] and uncertainty [*,2] 
;		*.ELEM		STRING	elemental symbol
;		*.Z		INT	atomic number
;		*.NLINES	INT	number of lines specified in the file		
;	coefs	STRUC	containing the fit coefficients
;		*.LINE#		FLTARR	[line.iline*3+2,3] of the
;					spectral fit coefficients
;		*.ELEM		STRING	elemental symbol
;		*.Z		INT	atomic number
;		*.NLINES	INT	number of lines specified in the file
;
;OPTIONAL OUTPUTS:
;	data	STRUC of the data and spectrometer labels for use w/ repeated calls if necessary
;	
;PROCEDURE:
;	
;
;MODIFICATION HISTORY:
;	Written by:	M.L. Reinke - 3/12
;	10/20/15	M.L. Reinke - modified to use a non-hardcoded named spectrometer data loading system
;                                    
;-

PRO impspec,shot,ispec,br,coefs,data=data,plot=plot,debug=debug,kline=kline,fitz=fitz,verb=verb
	
	nlines=ispec.nlines
	instr=strarr(nlines)
	FOR i=0,nlines-1 DO BEGIN
		j=nlines-1-i
		instr[i]=ispec.(j).spec
	ENDFOR
	
	;load data for 0th line
	IF NOT keyword_set(data) THEN BEGIN
		vuv_load_spec,shot,instr[0],specbr,lam,time,sigbr=sigbr
		idat={specbr:specbr,sig:sigbr,lam:lam,time:time,nt:n(time)+1}
		spec=instr[0]
		data=create_struct('d0',idat)
		FOR i=0,nlines-1 DO BEGIN
			tmp=where(spec EQ instr[i])
			IF tmp[0] EQ -1 THEN BEGIN
				vuv_load_spec,shot,instr[i],specbr,lam,time,sigbr=sigbr
				idat={specbr:specbr,sig:sigbr,lam:lam,time:time,nt:n(time)+1}
				spec=[spec,instr[i]]
				data=create_struct(data,'d'+num2str(n(spec),1),idat)
			ENDIF
        	ENDFOR
		data=create_struct(data,'spec',spec)
	ENDIF

	br={elem:ispec.elem,z:ispec.z,nlines:ispec.nlines}
	coefs={elem:ispec.elem,z:ispec.z,nlines:ispec.nlines}
	IF NOT keyword_set(kline) THEN BEGIN
		start=0
		stop=nlines-1
        ENDIF ELSE BEGIN
		IF iline EQ -1 THEN iline=0
		start=nlines-1-kline
		stop=nlines-1-kline
	ENDELSE
	FOR i=start,stop DO BEGIN
		j=nlines-1-i
		index=where(data.spec EQ instr[i])
		jcoefs=fltarr(3*ispec.(j).ilines+2,data.(index).nt,2)
		jbr=fltarr(data.(index).nt,3)
		jstatus=intarr(data.(index).nt)
		FOR k=0,data.(index).nt-1 DO BEGIN
			kcoefs=impspec_fit_line(data.(index).specbr[*,k],data.(index).lam,data.(index).sig[*,k],ispec.(j),plot=plot,status=status,error=kerror,fitz=fitz)
			IF status EQ 6 OR status EQ 7 THEN BEGIN
				jcoefs[*,k,0]=kcoefs
				jcoefs[*,k,1]=kerror
				jbr[k,0]=sqrt(2*!pi)*kcoefs[0]*kcoefs[2]
				jbr[k,2]=jbr[k,0]*sqrt((kerror[0]/kcoefs[0])^2+(kerror[2]/kcoefs[2])^2)
			ENDIF
			jstatus[k]=status
		ENDFOR
		jbr[*,1]=data.(index).time
		addstr="'line"+num2str(j,1)+"',jcoefs"
		result=execute("coefs=create_struct("+addstr+",coefs)")	
		addstr="'line"+num2str(j,1)+"',jbr"
		result=execute("br=create_struct("+addstr+",br)")	
		IF keyword_set(debug) THEN stop
	ENDFOR
END

;+
;NAME:
;	RUN_IMPSPEC
;	
;PURPOSE:
;	This procedure runs IMPSPEC and stores the data in the tree
;
;CALLING SEQUENCE:
;	RUN_IMPSPEC,shot,z
;
;INPUTS:
;	shot	LONG 	shot number
;	z	INT	atomic number of impurity to look at
;
;OPTIONAL INPUTS:
;	fitz	INTARR	of the other impurities to include in the fit (sent to IMPSPEC)
;
;OUTPUTS:
;	all data output is stored into the \SPECTROSCOPY::TOP.IMPSPEC tree
;
;PROCEDURE:
;	The ISPEC file is loaded from the tree using READ_ISPEC_TREE
;
;MODIFICAION HISTORY:
;	Written by:	M.L. Reinke - 3/12
;
;-

PRO run_impspec,shot,z,fitz=fitz
	tree=getenv('IMPSPEC_MDS_TREE')
	node=getenv('IMPSPEC_MDS_PATH')
	IF NOT keyword_set(fitz) THEN fitz=0

	;load ISPEC configuration file
	ispec=read_ispec_tree(shot,z)
	
	;run IMPSPEC
	impspec,shot,ispec,br,coefs,fitz=fitz

	mdsopen,tree,shot
	mdsput,node+'.'+ispec.elem+':NLINES','build_with_units($,"")',ispec.nlines
	FOR i=0,ispec.nlines-1 DO BEGIN
		path=node+'.'+ispec.elem+'.LINE'+num2str(i,1)
		mdsput,path+':BR','build_signal(build_with_units($1,"ph/s/m^2"),*,build_with_units($2,"seconds"),build_with_units($3,"ph/s/m^2"))',$
			br.(i)[*,0],br.(i)[*,1],br.(i)[*,2]
		mdsput,path+':COEFS','build_signal(build_with_units($1," "),*,build_with_units($2,"seconds"),build_with_units($3," "))',$
			coefs.(i)[*,*,0],br.(i)[*,1],coefs.(i)[*,*,1]
		mdsput,path+':DLAM','build_with_units($,"Angstroms")',ispec.(i).dlam
		mdsput,path+':LAM','build_with_units($,"Angstroms")',ispec.(i).lam
		mdsput,path+':LABEL','build_with_units($," ")',ispec.(i).label
        ENDFOR

	mdsclose,tree,shot
END
