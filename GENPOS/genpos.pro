;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;FUNCTIONS FOR GENERAL POSTION VECTOR;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

;+
;NAME:
;	_GENPOS
;
;PURPOSE:
;	The GENPOS_* suite of procedures and functions are the
;	application of formulas and methods in 'General Equations for
;	Radiometry in Tokamak Plasmas' for pinhole cameras.  Due to
;	the specifics of geometry contained in these functions it is
;	highly recommened to obtain a copy of this paper while using
;	these functions since its FIGURES will help dramatically
;
;	Version 1.0 - original 
;	Version 1.1 - Updated speed on GENPOS_VOL_COEFS
;	Version 1.2 - Added GENPOS_PLANAR tools
;	Version 1.3 - Added the matrix inversion tools
;
;MODIFICATION HISTORY:
;	Version 1.3
;	ML Reinke 
;	July 2006
;-

;+
;NAME:
;	_TUTORIAL  - CREATING A VOLUME WEIGHTING VECTOR
;
;PURPOSE:
;
;	The volume weightings are referred to as GPV (short for GENPOS Volumes) throughout
;	this program.  These are generated for a given grid defined by using GENPOS_GRID.  The
;	center points of this grid are found using the /center keyword in GENPOS_GRID.  For
;	C-Mod, the GRID_VES function can be used to generate grides that cover most of the
;	tokamak cross-section.  If the emissivity (Power/Volume - NO STERRADIANS!) is known at
;	each of these grid centers, then the total power desposited onto the dector for which
;	the GPV was calculated is the sum of the GPV times the emissivity at every grid point.
;	Note that when making such sums, it is quicker to use WHERE(GPV NE 0) to isolate
;	elements from which emission is actually collected.
;
;	TEST_SCRIPT:					;assume upos and du are given
;	nx=100
;	ny=100
;	ves_grid=GRID_VES(nx=nx, ny=ny) 		; grid for input to GENPOS_VOL_COEFS
;	ves_cent=GRID_VES(nx=nx, ny=ny, /center)	; 100x100 (R,Z) points of the grid centers for general use
;	gpv_raw=GENPOS_VOL_COEFS(ves_grid,upos,du)	; generate a raw GPV given a upos array [4,n] and du [n] vector
;	gpv=SUM_ARRAY(gpv_raw,/i)			; creates the GPV array of intereste 
;	GENPOS_GPV2CONTOUR,gpv,ves_cent			; make a contour plot of the gpv
;
;
;-

;+
;NAME:
;	GENPOS_PLANAR_INFO
;
;PURPOSE:
;	This function creates an detector information structure that can be used with other
;	functions in GENIE to create POS vectors and GPV data.  This function can only handle
;	a planar type of information file.  Data is entered in an ASCII format
;	
;CALLING SEQUENCE:
;	result=GENPOS_PLANAR_INFO(path)
;	
;INPUTS:
;	path	STR path to info file
;
;KEYWORD_PARAMETERS:
;	debug:	/debug stops the code before returning
;
;OUTPUTS:
;	resut:	STRUC of information describing the detector view.
;		*.name		STR name of the view
;		*.ap		STRUC aperture substructure
;		*.*.vec		FLTARR [ro,0.0,zo] of aperture [meters]
;		*.*.rot		FLTARR [alpha,beta,gamma] of aperture normal [radians]
;		*.*.size	FLTARR [y,z] of aperture dimensions [meters]
;		*.det		STRUC detector substructure
;		*.*.x0
;		*.*.x1		x3 FLTARR that locate the detector plane in aperture coordinates [meters]
;		*.*.x2
;		*.*.xi		FLTARR of xi locations of pixels on the detector [meters]
;		*.*.zeta	FLTARR of zeta locations of pixels on the detector [meters]
;		*.*.size	FLTARR [xi,zeta] of detector dimensions [meters]
;		*.type		STR of type (must be planar)
;		*.author	STR name of the author of the info file
;
;	If path points to an info file that doesn't have (on the second line) 'type = planar' then the
;	function will RETURN, -1
;
;
;RESTRICTIONS:
;	There are formatting issues with the INFO file.  See /home/mlreinke/idl/genie/data/info/planar_template.info 
;	to get a template to enter the data.  Contact the author with questions.
;
;MODIFICATION HISTORY:
;	Written by:	ML Reinke: 7-28-06
;	8-9-06:		ML Reinke - adjusted the reading of lengths from millimeters to centimeters
;       10-16-15:	ML Reinke - modified cntr variable to be 0L to allow for large pixel # detectors
;
;-

FUNCTION genpos_planar_info,path,debug=debug
	openr,lun,path, /get_lun
	line=strarr(1)

	;HEADER
	;---------------------------------------------------------------------------------------
	readf,lun,line
	tmp=strsplit(line, '=',/extract)
	name=tmp[1]
	readf,lun,line
	tmp=strsplit(line, '=',/extract)
	type=tmp[1]
	IF strmatch (type, '*planar*') EQ 0 THEN RETURN,-1	
	readf,lun,line
	tmp=strsplit(line, '=',/extract)
	author=tmp[1]

	readf,lun,line
	WHILE eof(lun) NE 1 AND strmatch(line,'*DATA START*',/fold_case) EQ 0 DO readf, lun, line


	;APERTURE
	;---------------------------------------------------------------------------------------
	;create ap_vector
	point_lun,-lun,data_start	
	readf,lun,line
	WHILE eof(lun) NE 1 AND strmatch(line,'ap_ro *',/fold_case) EQ 0 DO readf, lun, line
	ro=float(str_dbl_ext(line,'=',';'))
	point_lun,lun,data_start
	readf,lun,line
	WHILE eof(lun) NE 1 AND strmatch(line,'ap_zo *',/fold_case) EQ 0 DO readf, lun, line
	zo=float(str_dbl_ext(line,'=',';'))
	ap_vec=[ro,0.0,zo]

	;create ap_rot
	point_lun,lun,data_start
	readf,lun,line
	WHILE eof(lun) NE 1 AND strmatch(line,'ap_alpha *',/fold_case) EQ 0 DO readf, lun, line
	alpha=float(str_dbl_ext(line,'=',';'))
	point_lun,lun,data_start
	readf,lun,line
	WHILE eof(lun) NE 1 AND strmatch(line,'ap_beta *',/fold_case) EQ 0 DO readf, lun, line
	beta=float(str_dbl_ext(line,'=',';'))
	point_lun,lun,data_start
	readf,lun,line
	WHILE eof(lun) NE 1 AND strmatch(line,'*ap_gamma *',/fold_case) EQ 0 DO readf, lun, line
	gamma=float(str_dbl_ext(line,'=',';'))
	ap_rot=[alpha,beta,gamma]

	;create ap_size
	point_lun,lun,data_start
	readf,lun,line
	WHILE eof(lun) NE 1 AND strmatch(line,'*ap_y *',/fold_case) EQ 0 DO readf, lun, line
	y=float(str_dbl_ext(line,'=',';'))
	point_lun,lun,data_start
	readf,lun,line
	WHILE eof(lun) NE 1 AND strmatch(line,'*ap_z *',/fold_case) EQ 0 DO readf, lun, line
	z=float(str_dbl_ext(line,'=',';'))
	ap_size=[y,z]*1.0e-2

	;form apeture structure
	ap_struc=create_struct('vec',ap_vec,'rot',ap_rot,'size',ap_size)

	;DETECTOR
	;---------------------------------------------------------------------------------------
	;create det_size
	point_lun,lun,data_start
	readf,lun,line
	WHILE eof(lun) NE 1 AND strmatch(line,'det_xi *',/fold_case) EQ 0 DO readf, lun, line
	xi=float(str_dbl_ext(line,'=',';'))
	point_lun,lun,data_start
	readf,lun,line
	WHILE eof(lun) NE 1 AND strmatch(line,'det_zeta *',/fold_case) EQ 0 DO readf, lun, line
	zeta=float(str_dbl_ext(line,'=',';'))
	det_size=[xi,zeta]*1.0e-2

	;create detector positioning
	point_lun,lun,data_start
	readf,lun,line
	WHILE eof(lun) NE 1 AND strmatch(line,'*x0 *',/fold_case) EQ 0 DO readf, lun, line
	tmp=str_dbl_ext(line,'=',';')
	tmp=strsplit(tmp, ',', /extract)
	x0=float(tmp)*1.0e-2	
	point_lun,lun,data_start
	readf,lun,line
	WHILE eof(lun) NE 1 AND strmatch(line,'*x1 *',/fold_case) EQ 0 DO readf, lun, line
	tmp=str_dbl_ext(line,'=',';')
	tmp=strsplit(tmp, ',', /extract)
	x1=float(tmp)*1.0e-2
	point_lun,lun,data_start
	readf,lun,line
	WHILE eof(lun) NE 1 AND strmatch(line,'*x2 *',/fold_case) EQ 0 DO readf, lun, line
	tmp=str_dbl_ext(line,'=',';')
	tmp=strsplit(tmp, ',', /extract)
	x2=float(tmp)*1.0e-2

	;create detector locations
	point_lun,lun,data_start
	readf,lun,line
	WHILE eof(lun) NE 1 AND strmatch(line,'xi *',/fold_case) EQ 0 DO readf, lun, line
	xi=str_dbl_ext(line,'=',';')
	point_lun,lun,data_start
	readf,lun,line
	WHILE eof(lun) NE 1 AND strmatch(line,'zeta *',/fold_case) EQ 0 DO readf, lun, line
	zeta=str_dbl_ext(line,'=',';')
	IF strmatch(xi, '*-1*') AND strmatch(zeta, '*-1*') THEN BEGIN
		point_lun,lun,data_start
		readf,lun,line
		WHILE eof(lun) NE 1 AND strmatch(line,'*xi_ch *',/fold_case) EQ 0 DO readf, lun, line
		xi_ch=float(str_dbl_ext(line,'=',';'))
		point_lun,lun,data_start
		readf,lun,line
		WHILE eof(lun) NE 1 AND strmatch(line,'*n_xi *',/fold_case) EQ 0 DO readf, lun, line
		n_xi=float(str_dbl_ext(line,'=',';'))	
		point_lun,lun,data_start
		readf,lun,line
		WHILE eof(lun) NE 1 AND strmatch(line,'*xi_o *',/fold_case) EQ 0 DO readf, lun, line
		xi_o=float(str_dbl_ext(line,'=',';'))
		point_lun,lun,data_start
		readf,lun,line
		WHILE eof(lun) NE 1 AND strmatch(line,'*zeta_ch *',/fold_case) EQ 0 DO readf, lun, line
		zeta_ch=float(str_dbl_ext(line,'=',';'))
		point_lun,lun,data_start
		readf,lun,line
		WHILE eof(lun) NE 1 AND strmatch(line,'*n_zeta *',/fold_case) EQ 0 DO readf, lun, line
		n_zeta=float(str_dbl_ext(line,'=',';'))	
		point_lun,lun,data_start
		readf,lun,line
		WHILE eof(lun) NE 1 AND strmatch(line,'*zeta_o *',/fold_case) EQ 0 DO readf, lun, line
		zeta_o=float(str_dbl_ext(line,'=',';'))
		
		cntr=0L
		num=n_xi*n_zeta
		IF n_xi EQ 0 THEN num=n_zeta
		IF n_zeta EQ 0 THEN num=n_xi
		xi=fltarr(num)
		zeta=fltarr(num)
		IF n_zeta EQ 0 THEN n_zeta=1
		IF n_xi EQ 0 THEN n_xi = 1
		FOR i=0,n_zeta-1 DO BEGIN
			FOR j=0,n_xi-1 DO BEGIN
				xi[cntr]=xi_ch*(j-xi_o)*1.0e-2
				zeta[cntr]=zeta_ch*(i-zeta_o)*1.0e-2
				cntr+=1
			ENDFOR
		ENDFOR
		
	ENDIF ELSE BEGIN
		xi=float(strsplit(xi,',', /extract))*1.0e-2
		zeta=float(strsplit(zeta,',',/extract))*1.0e-2
	ENDELSE
	det_struc=create_struct('x0',x0,'x1',x1,'x2',x2,'xi',xi, 'zeta',zeta,'size',det_size)	

	close,lun
	free_lun,lun
	
	output={name:name, ap:ap_struc, det:det_struc, type:type, author:author}
	IF keyword_set(debug) THEN stop
	RETURN, output
END


;+
;NAME:
;	GENPOS_INVERT_POS
;	
;PURPOSE:
;	This function takes in POS vectors who's psi values are outside of the -!pi/2 < psi < !pi/2 range
;	where the functionality of GENPOS is defined and converts them to an equivlent POS vector looking
;	throught the plasma in the opposite direction.
;
;CALLING SEQUENCE:
;	result=GENPOS_INVERT_POS(pos)
;	
;INPUTS:
;	pos:	FLTARR of 4 x m where m is the number of pos vectors.
;
;OPTIONAL INPUTS:
;	rzbnd:	FLTARR of [r_min,r_max,abs(z)] of the plasma boundary.  See PROCEDURE for use.  DEFAULT: [0.44,1.0,0.5]
;
;KEYWORD PARAMETERS:
;	debug:	/debug stops the code before the RETURN statement
;
;OUTPUTS:
;	result:	FLTARR of 4 x m of coverted POS vectors.  If pos[3,i] was not out of bounds, the input POS is output without
;		any manipulation
;
;PROCEDURE:
;	First off, I'm pissed that I didn't realized I'd need this function until now.  Although the POS vector still uniquely
;	determines the hyperbola, the plasma view can vary from being l > 0 to l < 0 depending on the points used.  
;	Although an extra bit (or encoded on R_TANG since it's always > 0) to tell the direciton in l towards the plasma could solve the
;	problem it wouldn't work with existing codes (GENPOS_VOL_COEFS).  Instead, a POS vector can be found that is along the 
;	same hyperbola but has l > 0 describing the path through plasma.  To do this, a boundary in R and Z is needed that will
;	outside the plasma.  The inverted pos vector will have an (R1,Z1) on this boundary, the same R_TANG and a new PSI that
;	is inbounds.  The LINE_PATH_PLOTS of a set of POS vectors where some have been inverted will look awkward but their utility
;
;MODIFICATION HISTORY
;	Written by:	ML Reinke: 8-28-06
;
;-	

FUNCTION genpos_invert_pos,pos,rzbnd=rzbnd,debug=debug
	IF NOT keyword_set(rzbnd) THEN rzbnd=[0.44,1.0,0.6]
	in_pos=pos
	new_pos=fltarr(4,n(pos[0,*])+1)

	FOR i=0,n(pos[0,*]) DO BEGIN
		pos=reform(in_pos[*,i])
		IF pos[3] GT !pi/2.0 OR pos[3] LT -!pi/2.0 THEN BEGIN
			l_r=1.0-sqrt(1.0-(pos[0]^2-rzbnd[1]^2)/(pos[0]^2-pos[2]^2))
			IF pos[3] GT 0 THEN l_z=(pos[1]+rzbnd[2])/(sqrt(pos[0]^2-pos[2]^2)*tan(pos[3]))
			IF pos[3] LT 0 THEN l_z=(pos[1]-rzbnd[2])/(sqrt(pos[0]^2-pos[2]^2)*tan(pos[3]))		
			l=max([l_r,l_z])
			r1=sqrt(pos[0]^2+l*(l-2.0)*(pos[0]^2-pos[2]^2))
 
			z1=pos[1]-l*sqrt(pos[0]^2-pos[2]^2)*tan(pos[3])

			r2=pos[0]
			z2=pos[1]
			rt=pos[2]	
			cosphi=-((sqrt(r1^2-rt^2)-sqrt(r2^2-rt^2))^2-r1^2-r2^2)/(2.0*r1*r2)
			psi=atan((z1-z2)/sqrt(r1^2+r2^2-2.0*r1*r2*cosphi))
			new_pos[*,i]=[r1,z1,rt,psi]
		ENDIF ELSE new_pos[*,i]=pos
	ENDFOR
	
	pos=in_pos
	output=new_pos
	IF keyword_set(debug) THEN stop
	RETURN,output
END

;+
;NAME:
;	GENPOS_ROT_MATRIX
;
;PURPOSE:
;	This function generates the 3 x 3 rotation matrix used to
;	transform coordinates from the (xyz) aperature coordinate system to
;	the (x''' y''' z''') coordinate system whos axes are parallel
;	to the XYZ coordinates in which the tokamak center is located.
;
;CALLING SEQUENCE
;	result=GENPOS_ROT_MATRIX(a,b,g)
;
;INPUTS:
;	a:	FLT angle alpha (in radians) in FIGURE 7
;	b:	FLT angle beta (in radians) in FIGURE 7
;	c:	FLT angle gamma (in radians) in FIGURE 7
;
;OUTPUTS:
;	result:	FLTARR 3 x 3 matrix for which (xyz)''' can be determined
;		using result##[x,y,z].
;
;PROCEDURE:
;	The rotation angles alpha, beta and gamma locate the view of
;	the aperature normal w/r/t the tokamak center.  More
;	description is given FIGURE 7 of the General Equations paper.
;
;MODIFICATION HISTORY:
;	Written by:  	ML Reinke  - June 2006
;	ML Reinke	6-09-06 fixed an embarasses matrix multiplication error
;-

FUNCTION genpos_rot_matrix,a,b,g
	Ra=[[1,0,0],$
	    [0,cos(a),sin(a)],$
	    [0,-sin(a),cos(a)]]
	Rb=[[cos(b),0,sin(b)],$
	    [0,1,0],$
	    [-sin(b),0,cos(b)]]
	Rg=[[cos(g), sin(g),0],$
	    [-sin(g),cos(g),0],$
	    [0,0,1]]
	R=Rg##Rb##Ra
	RETURN,R
END


;+
;NAME:
;	GENPOS_CYL2TOK
;
;PURPOSE:
;	This function generates the tangency radius (R_t) and
;	declination angle (psi) given two RZ positions and cosine of
;	the angle between them
;
;CALLING SEQUENCE:
;	result=GENPOS_CYL2TOK(rz1,rz2,cosphi)
;
;INPUTS:
;	rz1:		FLTARR [r1,z1]
;	rz2:		FLTARR [r2,z2]
;	cosphi:		FLT cos(phi)
;
;OUTPUTS:
;	result:		FLTARR [R_t, psi]
;
;PROCEDURE:
;	This function uses equations (13) and (15) in Section 3.1 -
;	Tokamak Coordinates.  Note that if R_2 > R_1 the psi returned
;	will not be accurate.
;
;MODIFICATION HISTORY:
;	Written by:  	ML Reinke - June 2006
;   	8-30-06:	ML Reinke - added an if statement for when cosphi approaches one, the
;                                   singularity is avoided in computing r_tang
;                                                      
;-

FUNCTION genpos_cyl2tok,rz1,rz2,cosphi
	
	r1=double(rz1[0])
	z1=double(rz1[1])
	r2=double(rz2[0])
	z2=double(rz2[1])

	rt=r2*r1*sqrt(1.0-double(cosphi)^2)/sqrt(r1^2+r2^2-2.0*r1*r2*double(cosphi))
        IF 1.0-cosphi LT 1.0e-4 AND (r1-r2) LT 1.0e-4 THEN rt=0.0
	psi=atan((z1-z2)/sqrt(r1^2+r2^2-2.0*r1*r2*cosphi))
        IF r1 LT rt OR r2 LT rt THEN print, 'error in computing r_tang' 
	output=[rt,psi]
	RETURN, output

END


;+
;NAME:
;	GENPOS_DBLXYZ2TOK
;
;PURPOSE:
;	This function generates the tangency radius (R_t) and
;	declination angle (psi) given two positions in the xyz
;	(aperate) coordinate system and information about the aperture
;	position and viewing direction
;
;CALLING SEQUENCE:
;	result=GENPOS_DBLXYZ2TOK(ap_vec,xyz1,xyz2,rot)
;
;INPUTS:	
;	ap_vec:		FLTARR [ro, 0.0, zo] where ro and zo is the aperture position
;	xyz1:		FLTARR [x,y,z] location of point 1 in aperature coordinates
;	xyz2:		FLTARR [x,y,z] location of point 2 in aperature coordinates
;	rot: 		FLTARR rotation matrix generated using GENPOS_ROT_MATRIX
;
;OUTPUTS:
;	result:		FLTARR [R_t, psi]
;
;OPTIONAL OUTPUTS:
;	bigXYZ_vec:	FLTARR [2,3] where [i,*] is the XYZ vector of each point
;	cyl_vec:	FLTARR [2,3] where [i,*] is the cylindrical vector of each point 
;
;PROCEDURE:
;	This function generates the XYZ positions of each point then
;	forms the RZ positions and the cosphi and uses GENPOS_CYL2TOK
;	to find R_t and psi.
;
;MODFICATION HISTORY:
;	Written by:	ML Reinke - June 2006
;
;-

FUNCTION genpos_dblxyz2tok,ap_vec,xyz1,xyz2,rot,bigXYZ_vec=bigXYZ_vec,cyl_vec=cyl_vec
	;determine XYZ positions
	bigXYZ1=[0.0, ap_vec[0],ap_vec[2]]+rot##xyz1
        bigXYZ2=[0.0, ap_vec[0],ap_vec[2]]+rot##xyz2

        ;convert to cyl variables
        rz1=[sqrt(bigXYZ1[0]^2+bigXYZ1[1]^2),bigXYZ1[2]]
        phi1=asin(bigXYZ1[0]/sqrt(bigXYZ1[0]^2+bigXYZ1[1]^2))
        rz2=[sqrt(bigXYZ2[0]^2+bigXYZ2[1]^2),bigXYZ2[2]]
        phi2=asin(bigXYZ2[0]/sqrt(bigXYZ2[0]^2+bigXYZ2[1]^2))
        cosphi=cos(phi1-phi2)

        bigXYZ_vec=[[bigXYZ1],[bigXYZ2]]
        cyl_vec=[[rz1[0],phi1,rz1[1]],[rz2[0],phi2,rz2[1]]]
        tok_vec=genpos_cyl2tok(rz1,rz2,cosphi) ;turn into tokamak variables
        RETURN, tok_vec
END


;+
;NAME:
;	GENPOS_XYZ2TOK
;
;CALLING SEQUENCE:
;	result=GENPOS_XYZ2TOK(ap_vec,xyz_vec,rot)
;
;PURPOSE:
;	This function generates the tangency radius (R_t) and
;	declination angle (psi) given a single point in aperature
;	coordinates and assumes the aperate is the 2nd point.
;
;INPUTS:
;	ap_vec:		FLTARR [ro, 0.0, zo] where ro and zo is the aperture position
;	xyz_vec:	FLTARR [x,y,z] location of a point in aperature coordinates
;	rot: 		FLTARR rotation matrix generated using GENPOS_ROT_MATRIX
;
;OUTPUTS:
;	result:		FLTARR [R_t, psi]
;
;OPTIONAL OUTPUTS:
;	bigXYZ_vec:	FLTARR of the XYZ vector of the input xyz_vec
;	cyl_vec:	FLTARR of the cylindrical vector of the input xyz_vec
;
;PROCEDURE:
;	This function generates the XYZ positions of xyz_vec and forms
;	its RZ position.  The cosphi is generated assume the aperture
;	is at theta = 0 and has RZ of the aperature.  Then GENPOS_CYL2TOK
;	to find R_t and psi.
;
;MODFICATION HISTORY:
;	Written by:	ML Reinke - June 2006
;	8-25-06:	ML Reinke - changed it so that call to GENPOS_CYL2TOK uses detector as
;                                   rz1 and aperture as rz2
;
;-

FUNCTION genpos_xyz2tok,ap_vec,xyz_vec,rot,bigXYZ_vec=bigXYZ_vec,cyl_vec=cyl_vec
	bigXYZ_vec=[0.0, ap_vec[0],ap_vec[2]]+rot##xyz_vec

	z2=bigXYZ_vec[2]
	r2=sqrt(bigXYZ_vec[0]^2+bigXYZ_vec[1]^2)
	cosphi=bigXYZ_vec[1]/r2

	cyl_vec=[r2,acos(cosphi),z2]

	tok_vec=genpos_cyl2tok([r2,z2],[ap_vec[0],ap_vec[2]],cosphi)
	RETURN, tok_vec
END

;+
;NAME:
;	GENPOS_DET2XYZ
;
;PURPOSE:
;	This function transforms a position located on the detector
;	plane (xi, zeta) to a position in the aperature coordinate
;	system (xyz)
;
;CALLING SEQUENCE:
;	result=GENPOS_DET2XYZ(x0,x1,x2,det_pos)
;
;INPUTS:
;	x0:
;	x1:		x3 FLTARR locations of the detector's fiducial points in aperture coordinates.
;	x2:
;	det_pos:	FLTARR [2,n] where [*,i] are [xi,zeta] positions on the detector plane
;
;KEYWORD PARAMETERS:
;	debug:	/debug stops the program before RETURN
;
;OUTPUTS:
;	result:		FLTARR [3,n] where [*,i] are [x,y,z] positions in the aperature coordinate system
;
;PROCEDURE
;	The detector coordinates are located by x0,x1,x2 where x0->x1
;	is in the zeta direction and x0->x2 is in the xi direction.
;	The detector plane is define by xi and zeta and xi cross zeta
;	is the detector normal which points towards the aperature (ie
;	n_hat dot x_hat is > 0).  See section 3.2 'Aperture Coordinates'
;
;MODIFICATION HISTORY:
;	Written by:	ML Reinke, June 2006
;
;-

FUNCTION genpos_det2xyz,x0,x1,x2,det_pos,debug=debug	

	xi_uv=(x2-x0)/sqrt(total((x2-x0)^2))
	zeta_uv=(x1-x0)/sqrt(total((x1-x0)^2))

        n=n(det_pos[0,*])+1
        xyz_vec=fltarr(3,n)
        FOR i=0L,n-1 DO xyz_vec[*,i]=det_pos[0,i]*xi_uv+det_pos[1,i]*zeta_uv+x0
        IF keyword_set(debug) THEN stop

	RETURN,xyz_vec
END


;+
;NAME:
;	GENPOS_DET2U
;
;PURPOSE:
;	This function calculates the etendue, U, given a detectors
;	position and det/ap areas.
;
;CALLING SEQUENCE:
;	result=GENPOS_DET2U(x0,x1,x2,det_pos)
;
;INPUTS:
;	x0:
;	x1:		x3 FLTARR locations of the detector's fiducial points in aperture coordinates.
;	x2:
;	det_pos:	FLTARR [2,n] where [*,i] are [xi,zeta] positions on the detector plane	
;
;OPTIONAL INPUTS:
;	a_det:		FLT of the detector area
;	a_ap:		FLT of the aperture area
;
;KEYWORD PARAMETERS:
;	debug:	/debug stops the program before RETURN
;
;OUTPUTS:
;	result:		FLTARR of length n, with each element containing the exact etendue, including det/ap
;			cosine terms.  If a_det and/or a_ap are not included then output is per unit detector
;			and/or aperature area.
;
;OPTIONAL OUTPUTS:
;	det_cos:	FLTARR [n] of the n_hat dotted into xyz_hat. This is the cosine of the angle of incidence of the ray onto
;				the detector.  Useful for angular dependent filter transmissions at detector
;	det_cos:	FLTARR [n] of the x_hat dotted into xyz_hat. This is the cosine of the angle of incidence of the ray onto
;				through the aperture.  Useful for angular dependent filter transmissions at aperature
;
;PROCEDURE:
;	The exact relation for etendue is calculated (see Section 3.2 'Aperture Coordinates') by finding the projection of the
;	detector and aperture normals onto the line connecting the center of the aperature and detector.
;
;MODIFICATION HISTORY:
;	Written by:	ML Reinke - June 2006
;	12/2/08:	ML Reinke - added det_cos and ap_cos optional output
;
;-

FUNCTION genpos_det2U,x0,x1,x2,det_pos,a_det=a_det,a_ap=a_ap,debug=debug,det_cos=det_cos,ap_cos=ap_cos
	IF NOT keyword_set(a_det) THEN a_det=1.0
        IF NOT keyword_set(a_ap) THEN a_ap=1.0

	xi_uv=(x2-x0)/sqrt(total((x2-x0)^2))
	zeta_uv=(x1-x0)/sqrt(total((x1-x0)^2))

        n=n(det_pos[0,*])+1
        U=fltarr(n)
        det_cos=fltarr(n)
        ap_cos=fltarr(n)

        FOR i=0L,n-1 DO BEGIN
        	xyz=det_pos[0,i]*xi_uv+det_pos[1,i]*zeta_uv+x0
                xyz_mag=sqrt(total(xyz*xyz))
                norm=crossp((x2-x0),(x1-x0))
                n_hat=norm/sqrt(total(norm*norm))
;                IF i EQ 0 THEN BEGIN
;                	print, xyz
;                        print, n_hat
;                        print, (total(n_hat*xyz)*a_det/xyz_mag)
;                        print, (total([1.0,0.0,0.0]*xyz)*a_ap/xyz_mag)
;                        print, xyz_mag^2
;                ENDIF
                U[i]=(total(n_hat*xyz)*a_det/xyz_mag)*(total([1.0,0.0,0.0]*xyz)*a_ap/xyz_mag)/xyz_mag^2
                det_cos[i]=abs(total(n_hat*xyz)/xyz_mag)
                ap_cos[i]=abs(total([1.0,0.0,0.0]*xyz)/xyz_mag)
        ENDFOR
        
        IF keyword_set(debug) THEN stop

	RETURN,U
END




;+
;NAME:
;	GENPOS_GRID
;
;PURPOSE:
;	This function creates a regular a grid of points for use with
;	GENPOS_UPOS and GENPOS_VOL_COEFS.
;
;CALLING SEQUENCE:
;	result=GENPOS_GRID(delx,dely,n)
;
;INPUTS:
;	delx:	FLT extent in the horizontal direction
;	dely:	FLT extend in the vertical direction
;	n:	FLT number of rectangular grid locations (square spacing)
;		FLTARR [nx,ny] points in the [hor,ver] directions
;
;OPTIONAL INPUTS:
;	xo:	FLT x-position of the grid center
;	yo:	FLT y=position of the grid center
;
;KEYWORD_PARAMETERS:
;	center:	/gives points that are at the center of the grid points instead of their boundaries.
;	circle: /circle makes r=delx/2.0 and spacing is square with nx number of points in both directions.
;	debug:	/debug stops the program before RETURN
;
;OUTPUTS:
;	result: STRUC containing the grid points and the areas
;		*.pnts FLTARR [2,(nx+1)*(ny+1)] of points that bound a grid location
;		*.area FLT of the areas of grid elements
;		*.n    FLTARR [nx,ny] of input
;
;	(note that *.pnts is [2, nx*ny] if /center is used)
;
;RESTRICTIONS:
;	To use /circle you must have GENIE_LINE compiled since
;	LINE_INLCFS is use.  Also /circle isn't perfect yet.
;
;	This functions requires MLR_FUNCTIONS as well.
;
;EXAMPLE:
;	Using /center, GENPOS_GRID can be used to make grids that are pixel boundaries and
;	grids of points that are centers of those pixels.
;
;	out=GENPOS_GRID(1.0,2.0,[5,10])
;	out_c=GENPOS_GRID(1.0,2.0,[5,10],/center)
;	plot,out.pnts[0,*],out.pnts[1,*],psym=5,xr=[-0.75,0.75],/xsty,yr=[-1.25,1.25],/ysty
;	oplot,out_c.pnts[0,*],out_c.pnts[1,*],psym=4
;
;	This will plot triangles that are the locations of the grid
;	element boundaries and diamonds that are the centers of these
;	pixels.  Both have *.area that are the same.
;
;PROCEDURE:
;
;	The orientation of the grid points is
;
;	         DEFAULT				 /CENTER
;
;      4 --------------------------              --------------------------              
;	 |    |    |    |    |    |           3  |    |    |    |    |    |              
;      3 --------------------------              --------------------------
;        |    |    |    |    |    |           2  |    |    |    |    |    |              
;  y   2 --------------------------              --------------------------
;        |    |    |    |    |    |           1  |    |    |    |    |    |              
;      1 --------------------------              --------------------------
;        |    |    |    |    |    |           0  |    |    |    |    |    |              
;      0 --------------------------              --------------------------
;        0    1    2    3    4    5                 0    1    2    3    4    
;	             x
;      (4,4)-----------------------              --------------------------
;	 |    |    |    |    |    |              |3,3 |7,7 |    |    |    |
;      (3,3)-----------------------              --------------------------
;        |    |    |    |    |    |              |2,2 |6,6 |    |    |    |
;  y   (2,2)(7,7)------------------              --------------------------
;        |    |    |    |    |    |              |1,1 |5,5 |    |    |    |
;      (1,1)(6,6)------------------              --------------------------
;        |    |    |    |    |    |              |0,0 |4,4 |    |    |    |
;      (0,0)(5,5)------------------              --------------------------
;       	             
;	So to get the ygr and xgr vectors use:
;		ygr=grid.pnts[1,0:grid.n[1]]
;		xgr=grid.pnts[0,indgen(out.n[0]+1)*(out.n[1]+1)]
;
;	If using a /center grid then
;		ygr=grid.pnts[1,0:grid.n[1]-1]
;		xgr=grid.pnts[0,indgen(out.n[0])*(out.n[1])]
;
;MODIFICATION HISTORY:
;	Written by: ML Reinke June 2006
;	7-27-06:	ML Reinke - added the input n [#x, #y] to output structure for reasons of portability
;
;-

FUNCTION genpos_grid,delx,dely,n,xo=xo,yo=yo,circle=circle,debug=debug,center=center

	IF n(n) NE 0 THEN BEGIN
            nx=long(n[0])
            ny=long(n[1])
        ENDIF ELSE BEGIN
            nx=long(n)
            ny=long(n)
        ENDELSE
        n_new=[nx,ny]
        IF keyword_set(circle) THEN BEGIN 
        	ny=nx
                dely=delx
        ENDIF
        IF NOT keyword_set(xo) THEN xo=0.0
        IF NOT keyword_set(yo) THEN yo=0.0
        
        IF keyword_set(center) THEN BEGIN
            	area=delx/nx*dely/ny
        	xpts=make(xo+delx/2.0*(1.0/nx-1.0), xo+delx/2.0*(-1.0/nx+1.0),nx)
	        ypts=make(yo+dely/2.0*(1.0/ny-1.0), yo+dely/2.0*(-1.0/ny+1.0),ny)
        ENDIF ELSE BEGIN
            	area=delx/nx*dely/ny
                nx+=1
                ny+=1
		xpts=make(xo-delx/2.0,xo+delx/2.0,nx)
                ypts=make(yo-dely/2.0,yo+dely/2.0,ny)
 
        ENDELSE 

        cntr=0L
        pnts=fltarr(2,nx*ny)
        FOR i=0,nx-1 DO BEGIN
        	FOR j=0,ny-1 DO BEGIN
                	pnts[*,cntr]=[xpts[i],ypts[j]]
                        cntr+=1
                ENDFOR
        ENDFOR
        

        IF keyword_set(circle) THEN BEGIN
        	in_circ=fltarr(nx*ny)
                th=make(0.0,2.0*!pi,50)
                th=th[0:n(th)-1]
                xbnd=delx/2.0*cos(th)
                ybnd=delx/2.0*sin(th)
                FOR i=0,nx*ny-1 DO in_circ[i]=line_inlcfs(xbnd,ybnd,[xo,yo],pnts[*,i])
                good=where(in_circ EQ 1.0)
                pnts=pnts[*,good]
        ENDIF
               

        output={pnts:pnts, area:area, n:n_new}             	
        IF keyword_set(debug) THEN stop  


        RETURN, output
END

;+	
;NAME:
;	GENPOS_GRID2CENT
;	
;PURPOSE:
;	This function converts a GENPOS_GRID structure that is used in GENPOS_VOL_COEFS into a GENPOS_GRID structure that
;	was generated using the /center keyword.  A "centers" structure is what is used in the rest of GENPOS and this function
;	can be used in GENPOS_VOL_COEFS and GENPOS_PLANAR2GPV to convert the ves_grid into ves_cent to store with the GPV and
;	lhat data.
;
;CALLING SEQUENCE:
;	result=GENPOS_GRID2CENT(grid)
;
;INPUTS:
;	grid: 	STRUC containing the grid edge locations (see OUTPUT of GENPOS_GRID)
;
;KEYWORD PARAMETERS:
;	debug:	/debug stops the function before the RETURN statement
;
;OUTPUTS:
;	result:	STRUC containing the grid center locations (see OUTPUT of GENPOS_GRID using /center keyword)
;
;MODFICIATION HISTORY:
;	Written by:	ML Reinke - 8/8/07
;
;-	

FUNCTION genpos_grid2cent,grid,debug=debug
        nx=grid.n[0]
        ny=grid.n[1]
        area=grid.area
	n=[nx,ny]

        ;reduce from full set of grid points
	rgr=reform(grid.pnts[0,*])
	rgr=rgr[indgen(nx+1)*(ny+1)]
	zgr=reform(grid.pnts[1,*])
	zgr=zgr[0:ny]

	new_pnts=fltarr(2,nx*ny)
        
        cntr=0L
        FOR i=0,nx-1 DO BEGIN
        	FOR j=0,ny-1 DO BEGIN
                	new_pnts[0,cntr]=0.5*(rgr[i]+rgr[i+1])
                        new_pnts[1,cntr]=0.5*(zgr[j]+zgr[j+1])
                        cntr+=1
                ENDFOR
        ENDFOR
        
        cent={pnts:new_pnts, area:area[0], n:n}
        output=cent
        IF keyword_set(debug) THEN stop
        RETURN, output
END

;+
;NAME:
;	GENPOS_GRID2GR
;
;PURPOSE:
;	This procedure takes in a grid structure generated from
;	GENPOS_GRID and then returns the 1D vectors of the xgrid and
;	ygrid values.
;
;CALLING SEQUENCE:
;	GENPOS_GRID2GR,grid,xgr,ygr
;
;INPUTS:
;	grid:	STRUC from GENPOS_GRID (can be generated as default or using /center)
;
;OUTPUTS:
;	xgr:	FLTARR of the x-locations of the grids
;		length = grid.n[0] if /center, grid.n[0]+1 if default
;	ygr:	FLTARR of the y-location of the grids
;		length = grid.n[1] if /center, grid.n[1]+1 if default
;
;MODIFICATION HISTORY:
;	Written by:	M.L. Reinke 10-21-06
;
;-

PRO genpos_grid2gr,grid,xgr,ygr
	IF long(grid.n[0])*grid.n[1] LT n(grid.pnts[0,*])+1 THEN del=1 ELSE del=0
	ygr=grid.pnts[1,0:grid.n[1]-1+del]
	xgr=grid.pnts[0,indgen(grid.n[0]+del)*(grid.n[1]+del)]
END



;+
;NAME:
;	GENPOS_UPOS
;
;PURPOSE:
;	This function takes a grided detector and grided aperture and
;	determines the positions vectors and etendues for all the
;	possible connections.
;
;CALLING SEQUENCE:
;	result=GENPOS_UPOS(ap_vec,ap_rot,x0,x1,x2,ap_grid,det_grid)
;
;INPUTS:
;	ap_vec:		FLTARR [ro, 0.0, zo] where ro and zo is the aperture position
;	ap_rot: 	FLTARR [alpha,beta,gamma] set of orientation angles of the aperture [in radians]
;	x0:
;	x1:		x3 FLTARR locations of the detector's fiducial points in aperture coordinates.
;	x2:
;	ap_grid:	STRUC grid of aperature subarray centers (use GENPOS_GRID with /center)
;	det_grid:	STRUC grid of detector subarray centers (use GENPOS_GRID with /center)
;
;OPTIONAL INPUTS
;	rot_matrix:	The output of GENPOS_ROT_MATRIX for ap_rot can
;			be input to save computation time if processing multiple detectors
;
;KEYWORD_PARAMETERS:
;	debug:		/debug to stop before RETURN
;
;OUTPUTS:
;	result:		STRUC containing the position vectors and etendues of the n_ap*n_det different combinations 
;			*.upos FLTARR [4,n_ap*n_det] position vectors r1,z1,rt,th] with point 1 being the detector
;			*.du   FLTARR [n_ap*n_det] values of the etendue linking the detector and aperature subarray
;
;RESTRICTIONS:
;	The aperature and detector grids must be generated using GENPOS_GRID so the ap_grid points will be in the y-z plane
;	with the center being the origin.  The det_grid points will be in the xi-zeta plane located by x0,x1 and x3.  The center of
;	the grid (input to GENPOS_GRID) should be the center of the detector prior to being subdivided.
;
;PROCEDURE
;	See the 'Beyone the Line Integral Approximation' section for info on why this is a useful thing to do.
;
;MODIFICATION HISTORY:
;	Written by:	M.L. Reinke June 2006
;	6-27-06:	M.L. Reinke: corrected calc of dU to include cosine terms
;	8-25-06:	ML Reinke - adjusted call to GENPOS_DBLXYZ2TOK to use detector points as XYZ1 and
;                                   aperture as XYZ2
;	8-28-06:	ML Reinke - added the ability to handle normally out of bounds psi values
;                                   using GENPOS_INVERT_POS
;
;-

FUNCTION genpos_upos,ap_vec,ap_rot,x0,x1,x2,ap_grid,det_grid,rot_matrix=rot_matrix,debug=debug,t=t,ap=ap,clear=clear

	IF NOT keyword_set(rot_matrix) THEN rot_matrix=genpos_rot_matrix(ap_rot[0],ap_rot[1],ap_rot[2])

        apnts=ap_grid.pnts
        dpnts=det_grid.pnts
        n_ap=n(apnts[0,*])+1
        n_det=n(dpnts[0,*])+1
        upos=fltarr(4,n_ap*n_det)
        dU=fltarr(n_ap*n_det)
       	clear=fltarr(n_ap*n_det)+1.0
                

        ap_xyz=[transpose(fltarr(n_ap)),apnts]
        det_xyz=genpos_det2xyz(x0,x1,x2,dpnts)
        norm=crossp((x2-x0),(x1-x0))
	n_hat=norm/sqrt(total(norm*norm))

        cntr=0L
        FOR i=0,n_ap-1 DO BEGIN
        	FOR j=0,n_det-1 DO BEGIN
                	out=genpos_dblxyz2tok(ap_vec,det_xyz[*,j],ap_xyz[*,i],rot_matrix,cyl_vec=cyl_vec)
                        upos[*,cntr]=[cyl_vec[0,0],cyl_vec[2,0],out[0],out[1]]
                        IF cyl_vec[0,0] LT cyl_vec[0,1] THEN BEGIN
                        	IF cyl_vec[2,0] GT cyl_vec[2,1] THEN upos[3,cntr]=!pi-out[1]
                        	IF cyl_vec[2,0] LT cyl_vec[2,1] THEN upos[3,cntr]=-!pi+out[1]
                        	upos[*,cntr]=genpos_invert_pos(reform(upos[*,cntr]))
                        ENDIF
			;find dU
	        	xyz=det_xyz[*,j]-ap_xyz[*,i]
        	        xyz_mag=sqrt(total(xyz*xyz))
        	        dU[cntr]=(total(n_hat*xyz)*det_grid.area/xyz_mag)*(total([1.0,0.0,0.0]*xyz)*ap_grid.area/xyz_mag)/xyz_mag^2.0
			;dU[cntr]=det_grid.area*ap_grid.area/xyz_mag^2.0
			;print, 1-total(n_hat*xyz)*total([1.0,0.0,0.0]*xyz)/xyz_mag^2.0
                        IF keyword_set(t) AND keyword_set(ap) THEN BEGIN
                        	l_int=[(ap.size[1]/2.0-ap_xyz[2,i])/(det_xyz[2,j]-ap_xyz[2,i]), (-ap.size[1]/2.0-ap_xyz[2,i])/(det_xyz[2,j]-ap_xyz[2,i]),$
                                	(ap.size[0]/2.0-ap_xyz[1,i])/(det_xyz[1,j]-ap_xyz[1,i]), (-ap.size[0]/2.0-ap_xyz[1,i])/(det_xyz[1,j]-ap_xyz[1,i])]
                                tmp=where(l_int GE 0 AND l_int LE -t/det_xyz[0,j])
                                IF tmp[0] NE -1 THEN BEGIN
                                    	;stop
                                	clear[cntr]=0.0
                                ENDIF
                                
                        ENDIF                       	
                        cntr+=1
                ENDFOR
        ENDFOR

        output={upos:upos, dU:dU}
        IF keyword_set(debug) THEN stop
        RETURN, output
END

;+
;NAME:
;	GENPOS
;
;PURPOSE:
;	This function calls several lower level functions to calculate,from user input, 
;	a set of position vectors ([ro,zo,rt,psi]) that can be used in other GENIE programs.
;
;CALLING SEQUENCE
;	result=GENPOS(ap_vec,ap_rot,x0,x1,x2,det_pos)
;	
;INPUTS:
;	ap_vec: 	FLTARR [ro, 0.0, zo] where ro and zo is the aperture position
;	ap_rot:		FLTARR [alpha,beta,gamma] set of orientation angles of the aperture [in radians]
;	x0:
;	x1:		x3 FLTARR locations of the detector's fiducial points in aperture coordinates.
;	x2:
;	det_pos:	FLTARR [xi, zeta] x m where m is the number of points in the deteector plane
;				for which GENPOS will cacluate position vectors
;
;OPTIONAL INPUTS:
;	a_det:		FLT area of detector [DEFAULT = 1.0]
;	a_det:		FLTARR [xi zeta] (width,height) of the detector
;	a_ap:		FLT area of aperature [DEFAULT = 1.0]
;	a_ap:		FLT [y,z] (width, height) of the aperature.  Circle is defined by [rad, -1.0]
;	
;KEYWORD PARAMETERS:
;	debug:		/debug stops before the RETURN statements
;
;OUTPUTS:
;	result:		FLTARR [ro,zo,rt,psi] x m position vectors
;
;OPTIONAL OUTPUTS:
;	bigXYZ_vec: 	FLTARR [X,Y,Z] vector of points in the XYZ system (tokamak)
;	cyl_vec:	FLTARR [R,phi,Z] points in R-PSI-Z (tokamak cylindrical)
;	xyz_vec:	FLTARR [x,y,z] vector of points in the xyz system (aperture)
;	etendue:	FLTARR [n] of etendue values for the det_pos and x0,x1,x2.  If optional inputs a_det
;				and a_ap are not input, this is per unit aperature and detector area.
;	det_cos:	FLTARR [n] of the n_hat dotted into xyz_hat. This is the cosine of the angle of incidence of the ray onto
;				the detector.  Useful for angular dependent filter transmissions at detector
;	det_cos:	FLTARR [n] of the x_hat dotted into xyz_hat. This is the cosine of the angle of incidence of the ray onto
;				through the aperture.  Useful for angular dependent filter transmissions at aperature
;
;RESTRICTIONS:
;	This function requires MLR_FUNCTIONS and other functions in GENPOS.PRO
;
;PROCEDURE:
;	See ML Reinke's "General Equations for Radiometry in Tokamak Plasmas" for an in-depth
;	descriptions of the analysis that is performed in this and other GENPOS functions.
;
;	GENPOS was formely used to find field of view bounds, but this functionality has been removed since 
;	GENPOS_VOL_COEFS has been developed to do the job more accurately.
;
;	The POS vector is found using GENPOS_XYZ2TOK and both points are checked to see if psi is in bounds.  If not, it is
;	corrected and the pos vector is sent to GENPOS_INVERT_POS to form the equivilent line of sight.
;
;MODIFICATION HISTORY:
;	Written by: 	ML Reinke 4-17-06
;	5-18-06:	ML Reinke - Overhaulled to included the field of view position vectors
;	8-28-06:	ML Reinke - Overhaulded removing full,fovpos type crap since it's covered by
;                                   GENPOS_VOL_COEFS.  Also fixed the -!pi/2 < psi < !pi/2 error using GENPOS_INVERT_POS
;	12-2-08:	ML Reinke - Added links to the ap_cos and det_cos optional outputs from GENPOS_DET2U
;
;-
 
FUNCTION genpos,ap_vec,ap_rot,x0,x1,x2,det_pos,debug=debug,a_ap=a_ap,a_det=a_det,bigXYZ_vec=bigXYZ_vec,cyl_vec=cyl_vec,$
                xyz_vec=xyz_vec,etendue=etendue,ap_cos=ap_cos,det_cos=det_cos

	;setup sizing demanded by det_pos
        xvec=size(det_pos)
	IF xvec[1] NE 2 THEN RETURN, -1
	IF xvec[0] EQ 1 THEN npos=1 ELSE npos=xvec[2]
	pos=fltarr(4,npos)

        ;generate rotation matrix and center detector positions
	rot_matrix=genpos_rot_matrix(ap_rot[0],ap_rot[1],ap_rot[2])
        xyz_vec=genpos_det2xyz(x0,x1,x2,det_pos)


        ;generate tokamak position vector for each detector
	FOR i=0L,npos-1 DO BEGIN
		tok_vec=genpos_xyz2tok(ap_vec,xyz_vec[*,i],rot_matrix,bigXYZ_vec=bigXYZ_vec,cyl_vec=cyl_vec)
                ;create the pos vector defined with detector as point 1 (R1,Z1) see section 3.3
                pos[0,i]=cyl_vec[0]  ;R of detector point
		pos[1,i]=cyl_vec[2]  ;Z of detector point
		pos[2,i]=tok_vec[0]  ;tangency radius of view
                pos[3,i]=tok_vec[1]  ;psi
                IF cyl_vec[0] LT ap_vec[0] THEN BEGIN
                	IF cyl_vec[2] GT ap_vec[2] THEN pos[3,i]=!pi-tok_vec[1]
                        IF cyl_vec[2] LT ap_vec[2] THEN pos[3,i]=-!pi+tok_vec[1]
                        pos[*,i]=genpos_invert_pos(reform(pos[*,i]))
                ENDIF
        ENDFOR

        IF keyword_set(a_ap) THEN IF n(a_ap) NE 0 THEN ap_area=a_ap[0]*a_ap[1] ELSE ap_area=a_ap 
        IF keyword_set(a_det) THEN IF n(a_det) NE 0 THEN det_area=a_det[0]*a_det[1] ELSE det_area=a_det
        etendue=genpos_det2U(x0,x1,x2,det_pos,a_ap=ap_area,a_det=det_area,ap_cos=ap_cos,det_cos=det_cos)  ;etendue of each channel

        IF keyword_set(debug) THEN stop

	RETURN,pos
END


;+
;NAME:
;	GENPOS_VOL_COEFS
;
;PURPOSE:
;	This function takes an array of tokamak positions vectors, and their etendues and determines the volume 
;	element w/r/t to a user supplied grid.  It is assumed that the inputs are sub-arrays of a single detector channel.
;
;CALLING SEQUENCE:
;	result=GENPOS_VOL_COEFS(grid,upos,du)
;
;INPUTS:
;	grid:		STRUC boundaries of rectangular pixels.  It's necessary to use the output of GENPOS_GRID (DO NOT use /center)
;	upos:		FLTARR of size 4 x m of pos [ro,zo,rt,psi] For planar detectors, use the output of GENPOS_UPOS *.upos
;	du:		FLTARR of length m of the differential etendues [m^2 str] for each pos vector in upos.  For planar detectors,
;				use the output of GENPOS_UPOS *.du
;
;OPTIONAL INPUTS:
;	vol2d:		If vol2d is set to a non-zero variable then it will be filled with the volume 
;			information availale in 2D (see optional outputs)
;	lhat:		Set lhat to a non-zero named variable and it will be filled with the optional output lhat structure	
;	
;KEYWORD PARAMETERS:
;	phineg:		/phineg changes the sign on the toroidal direction unit vector
;	debug:		/debug stops the code in various places
;	kdebug:		/kdebug invokes loop_debug and plots for upos[*,k]
;	loop_debug:	/loop_debug stops the code during the UPOS loop
;	plots:		/plots outputs plots of the grid, intersection points and line of sight.
;	contour:	/contour generates a contour plot of the volume elements on the (R,Z) grid input.
;	ps:		/ps suppresses calls to windows and allows PS files to be generated
;	 		NOTE: There will be UPOS number of plots if this is invoked, so do with care.
;
;OUTPUTS:
;	result:		FLTARR of length n x m, where n is the number of grid points and m is the number of upos elements.  The
;				values are volumes associated with elements described by grid.  To get the total, assuming upos 
;				and du describe one detector, use vol=SUM_ARRAY(result, /i) (SUM_ARRAY is a
;				function in MLR_FUNCTIONS which is necessary to run GENPOS)
;
;OPTIONAL OUTPUTS:
;	vol2d:		STRUC of the volume information in 2D
;			*.vol 	FLTARR [nx,ny] of the volume elements
;			*.r   	FLTARR of length nx of the radial grid points (centers)
;			*.z   	FLTARR of length ny of the vertical grid points (centers)
;	ctime:		FLT ctime will return the total computation time of GENPOS_VOL_COEFS
;	lhat:		STRUC of the weighted unit vectors in (r,phi,z).  For each UPOS, dV*unit vector is added at each
;			voxel intersection.  See PROCEDURE for info on creating the useful input to GENSPEC functions
;			*.r		STRUC of the weighted radial unit vectors and weightings
;			*.*.in		FLTARR [n, m] of the weighted radial unit vector going in
;			*.*.dvin	FLTARR [n, m] of the volume weightings going in
;			*.*.out		FLTARR [n, m] of the weighted radial unit vector coming out
;			*.*.dvout	FLTARR [n, m] of the volume weightings coming out
;			*.phi		FLTARR [n, m] of the weighted toroidal unit vector
;			*.z		FLTARR [n, m] of the weighted vertical unit vector
;
;	IF /ps is used then postscript outputs will be stored in /home/username/idl/plots
;		/contour - gpv_genplt_*.ps
;		/plots - gpv_grids.ps
;
;PROCEDURE:
;	GENPOS_VOL_COEFS must be run for each detector element.  The output of this function need only be calculated once for
;	installed detectors and then can be stored for subsequent use in an inversion algorithm.  To model detectors from a
;	known/modeled emissivity, the power deposited on a detector can be calculated from these volume estimates as total(vol*emiss)
;
;	The extents of GRID are used to determine bounds in the parameter-l that will intersect the grid (l_min,l_max).  
;	These are turned into boundaries in R and Z which are used to create a subset of the R and Z locations of the grid.  
;	The intersection points on the grid are found by solving for the parameter-l where the hyperbola definded by POS
;	intersects either an R or Z grid location.  The intersection points are then sorted in parameter-l and the volume element
;	is added to the corresponding pixel.  dV+=dU/(4pi)*delta_l*sqrt[(R1^2-RT^2)*(1+(tan(psi))^2)].  See 'Beyond the Line
;	Integral Approximation' for why this value is useful
;
;	NOTE:   it is up to the end user to show convergance of the volume elements as the aperture and detector meshing
;		is increased.
;
;	Like GPV, to form the lhat array that should be stored you must use sum_array to add over all the UPOS vectors.  Then
;	divided out the total GPV, limiting to non-zero values, to get the average unit vector at each spatial point.
;
;
;MODIFICATION HISTORY:
;	Written by:	ML Reinke - June 2006
;	6-30-06:	ML Reinke - large overhaul in the loop for each POS vector to decrease calculation time
;				    intersections are done in parameter-l and (R,Z) points are only calculated
;				    when plotting.
;	7-31-06:	ML Reinke - adjusted the finding of zgr intersections so that duplicates (line intersects a grid point) are
;                                   not counted.
;	8-09-06:	ML Reinke - added a flag to contour plot to check to see if there were any volume elements added.  This
;                                   allows arrays that have channels off the grid to be contour plotted in a simple loop.
;	8-30-06:	ML Reinke - discovered that GPV eats it when psi->!pi/2.0
;	8-31-06:	ML Reinke - fixed a bug when i_int[tmp] and j_int[tmp] are computed that made values = -2 when slightly over/under
;				    the  max/min of the grid size. Values are now forced to edge if over the edge.
;	8-08-07:	ML Reinke - added the lhat optional output to calculate weighted unit vectors for GENSPEC
;	8-14-07:	ML Reinke - changed the formatting of lhat.r to be a sub structure for inward and outward views.
;-

FUNCTION genpos_vol_coefs,grid,upos,du,debug=debug,plots=plots,ps=ps,contour=contour,loop_debug=loop_debug,vol2d=vol2d,ctime=ctime,$
                          kdebug=kdebug,vessel=vessel,lhat=lhat,phineg=phineg
	
	start_time=systime(/seconds)
        du/=grid.area*1.0e-6		;for small gridings numbers becomes small, removed out after summation

	;determine grid sizing from grid variable (create using GENPOS_GRID)
	tmp=where(grid.pnts[0,*] EQ grid.pnts[0,0])
	ny=tmp[n(tmp)]+1
	nx=int((n(grid.pnts[0,*])+1)/ny)
        npos=n(upos[0,*])+1		;number of views
	vol=fltarr((nx-1)*(ny-1),npos)	;volume elements for each pixel
					;using GENPOS_GRID pixel=0 is min(r) and min(z)
					;and count increase in r-direction first
        ;setup lhat arrays if called for
	IF keyword_set(lhat) THEN BEGIN
        	dV_in=fltarr((nx-1)*(ny-1),npos)
                dV_out=fltarr((nx-1)*(ny-1),npos)
           	lhat_r_in=fltarr((nx-1)*(ny-1),npos)
                lhat_r_out=fltarr((nx-1)*(ny-1),npos)
        	lhat_phi=fltarr((nx-1)*(ny-1),npos)
                lhat_z=fltarr((nx-1)*(ny-1),npos)
                IF keyword_set(phineg) THEN sign=-1.0 ELSE sign=1.0
        ENDIF

	;reduce from full set of grid points
	rgr=reform(grid.pnts[0,*])
	rgr=rgr[indgen(nx)*ny];rgr[int(make(0,nx-1,nx))*ny]
	zgr=reform(grid.pnts[1,*])
	zgr=zgr[0:ny-1]

        minr=min(rgr)
        maxr=max(rgr)
        minz=min(zgr)
        maxz=max(zgr)

	;load limiter data from TREE
	mdsopen, 'analysis',-1
	node_year='2002'
	r_lim=mdsvalue('.LIMITERS.WALL.WALL_'+node_year+':RLIM')
	z_lim=mdsvalue('.LIMITERS.WALL.WALL_'+node_year+':ZLIM')
	r_rf=mdsvalue('.LIMITERS.RF_LIMITER:R')
	z_rf=mdsvalue('.LIMITERS.RF_LIMITER:Z')
	mdsclose

	zoff=0.0001*(zgr[1]-zgr[0])
        roff=0.0001*(rgr[1]-rgr[0])
        IF keyword_set(kdebug) THEN kpt=kdebug ELSE kpt=-1
        IF keyword_set(loop_debug) THEN ldebug=1 ELSE ldebug=0
        IF keyword_set(plots) THEN pl=1 ELSE pl=0

	;start loop to fill volume dataset
	FOR k=0L,npos-1 DO BEGIN
		IF k EQ kpt THEN loop_debug=1 ELSE loop_debug=ldebug
                IF k EQ kpt THEN plots=1 ELSE plots=pl
                pos=reform(upos[*,k])	;form current pos vector
                ;if pos describes a horizontal line on a grid boundary move it up/down and change up/down direction of next occurance.
		IF pos[3] EQ 0 AND total(where(pos[1] EQ zgr)) NE -1 THEN BEGIN
			pos[1]+=zoff
			zoff*=-1.0
                ENDIF
                ;if pos describes a vertical line on a grid boundary move it up/down and change up/down direction of next occurance.
                IF pos[3] EQ !pi/2.0 and total(where(pos[0] EQ rgr)) NE -1 THEN BEGIN
                	pos[0]+=roff
                        roff*=-1.0
                ENDIF
		IF keyword_set(plots) THEN BEGIN
			;open window and plot vessal
			IF NOT keyword_set(ps) THEN BEGIN
				device, window_state=var
				IF var[19] EQ 0 THEN window,19,xsize=500,ysize=690,xpos=0,ypos=670,title='vessel cx,19' ELSE wset,19
                        ENDIF ELSE BEGIN
                        	d_old=!d
                                device, xsize=7.0, ysize=7.0*675.0/500.0, /inches
                        ENDELSE

			IF NOT keyword_set(vessel) THEN plot,r_lim,z_lim,title='C-MOD: ',chars=1.3 ELSE $
                        	vessel_plot,ps=ps,title='C-MOD'

			;plot grid and line
			makesym,10
			oplot, grid.pnts[0,*],grid.pnts[1,*],psym=8
                        FOR i=0,n(rgr) DO oplot, [rgr[i],rgr[i]],[zgr[0],zgr[n(zgr)]]
                        FOR i=0,n(zgr) DO oplot, [rgr[0],rgr[n(rgr)]],[zgr[i],zgr[i]]
			l=make(0.0,2.0,100)
			r_hyp=sqrt(pos[0]^2+l*(l-2.0)*(pos[0]^2-pos[2]^2))
			z_hyp=pos[1]-l*sqrt(pos[0]^2-pos[2]^2)*tan(pos[3])
			oplot, r_hyp,z_hyp,color=100
		ENDIF
	
		;find the interval in parameter-l where the line intersects the edges of the grid
		;this will bound calculations of intersection points, speeding up the calc
		;for large plasma, ap and det grids.
		a=sqrt(pos[0]^2-pos[2]^2)*tan(pos[3])
		aa=(pos[0]^2-pos[2]^2)
		IF pos[2] GT 0 THEN BEGIN		
			l_top=(pos[1]-maxz)/a	
			l_bot=(pos[1]-minz)/a
		ENDIF ELSE BEGIN
			l_top=(pos[1]-maxz)/(pos[0]*tan(pos[3]))
			l_bot=(pos[1]-minz)/(pos[0]*tan(pos[3]))
		ENDELSE
		b=(pos[0]^2-minr^2)/aa
		c=(pos[0]^2-maxr^2)/aa
		l_right=1-sqrt(1-c)
		IF pos[2] LE minr THEN l_left=1-sqrt(1-b) ELSE l_left=1+sqrt(1-c)

		lrange=[l_top,l_bot,l_right,l_left]
		order=sort(lrange)
		lmin=lrange[order[1]]*0.999999
		lmax=lrange[order[2]]*1.000001

		;calculate the l values of the intersection
		l_int=[0]	;insert intial zero so array is not empty if either tmp[0]=-1
		i_int=[0]
		j_int=[0]
		rbounds=[sqrt(pos[0]^2+lmin*(lmin-2.0)*aa),sqrt(pos[0]^2+lmax*(lmax-2.0)*aa)]
                zo=pos[1]-a
		IF pos[2] GT minr AND (zo LT maxz AND zo GT minz) THEN rbounds[1]=pos[2]
		max=max(rbounds,min=min)
		tmp=where(rgr LE max AND rgr GE min)
		IF tmp[0] NE -1 THEN BEGIN
			c=(pos[0]^2-rgr[tmp]^2)/(pos[0]^2-pos[2]^2)	;extra calculations since must first calc +/- R's
			ipts=[tmp,tmp]
			l_tmp=[1+sqrt(1-c),1-sqrt(1-c)]
			inside=where(l_tmp LE lmax AND l_tmp GE lmin)
			l_int=[l_int,l_tmp[inside]]
			i_int=[i_int,ipts[inside]]
			j_int=[j_int,intarr(n(inside)+1)-1]
		ENDIF
		zbounds=[pos[1]-lmin*a,pos[1]-lmax*a]
		max=max(zbounds,min=min)
		tmp=where(zgr LE max AND zgr GE min)
		IF tmp[0] NE -1 THEN BEGIN
                        l_int_z=(pos[1]-zgr[tmp])/a
                        check=setintersection(l_int_z,l_int)
                        IF check[0] NE -1 THEN BEGIN
                        	FOR ii=0,n(tmp) DO BEGIN
                            		IF where(check EQ ii) EQ -1 THEN BEGIN
                                        	l_int=[l_int,l_int_z[ii]]
                                                i_int=[i_int, -1]
                                                j_int=[j_int,tmp[ii]]
                                            ENDIF
                                        ENDFOR
                        ENDIF ELSE BEGIN
				l_int=[l_int,l_int_z]
                                num=n(tmp)+1
                                i_int=[i_int, intarr(num)-1]
                                j_int=[j_int,tmp]
                        ENDELSE
		ENDIF
		num=n(l_int)
                IF num NE 0 THEN BEGIN
                	l_int=l_int[1:num] ;remove initial zero
			i_int=i_int[1:num]
			j_int=j_int[1:num]
			
			order=sort(l_int)
			l=l_int[order]
			i_int=i_int[order]
			j_int=j_int[order]
            
			;if plotting plot the intersection points
			IF keyword_set(plots) THEN BEGIN
				int=fltarr(2,n(l)+1)
				int[0,*]=sqrt(pos[0]^2+l*(l-2.0)*aa)
				int[1,*]=pos[1]-l*a
				makesym,10
				oplot,int[0,*],int[1,*],color=200,psym=8
			ENDIF

			IF keyword_set(loop_debug) THEN stop		

			;fill the i_int, j_int values that are unknown with the lesser
			;of the two bounding indices
                        tmp=where(i_int EQ -1) 
                        IF keyword_set(loop_debug) THEN i_int_tmp=float(i_int)		
 			IF tmp[0] NE -1 THEN BEGIN
 				rpts=sqrt(pos[0]^2+l[tmp]*(l[tmp]-2.0)*aa)
                                IF keyword_set(loop_debug) THEN i_int_tmp[tmp]=interp_vec_reform(rgr,rpts)
 				FOR m=0,n(tmp) DO BEGIN
 					loc=where(rgr-rpts[m] GE 0)
                                        IF rpts[m] LT rgr[0] THEN loc=-1
                                        IF loc[0] EQ -1 THEN IF rpts[m] GT maxr THEN i_int[tmp[m]]=rgr[nx-1] ELSE i_int[tmp[m]]=rgr[0] ELSE $
                                        	i_int[tmp[m]]=loc[0]-1
 				ENDFOR
 			ENDIF
 			tmp=where(j_int EQ -1)
                        IF keyword_set(loop_debug) THEN j_int_tmp=float(j_int)
 			IF tmp[0] NE -1 THEN BEGIN
 				zpts=pos[1]-l[tmp]*a
                                IF keyword_set(loop_debug) THEN j_int_tmp[tmp]=interp_vec_reform(zgr,zpts)
 				FOR m=0,n(tmp) DO BEGIN
 					loc=where(zgr-zpts[m] GE 0)
                                        IF zpts[m] LT zgr[0] THEN loc=-1
                                        IF loc[0] EQ -1 THEN IF zpts[m] GT maxz THEN j_int[tmp[m]]=zgr[ny-1] ELSE j_int[tmp[m]]=zgr[0] ELSE $
                                        	j_int[tmp[m]]=loc[0]-1
 				ENDFOR
 			ENDIF

 	                IF keyword_set(loop_debug) THEN stop
                        
                        
         	        ;start loop to calc volume elements for each point
 			FOR i=0,n(l)-1 DO BEGIN
 				ipt=i_int[i] < i_int[i+1]
 				jpt=j_int[i] < j_int[i+1]
 				IF (i_int[i] EQ i_int[i+1] AND j_int[i] EQ j_int[i+1]) AND (l[i] LE 1.0 AND l[i+1] GE 1.0)  THEN ipt-=1
                                IF keyword_set(loop_debug) THEN BEGIN
                                	print, '(i,j)   i',i_int_tmp[i],j_int_tmp[i]
                                        print, '(i,j) i+1',i_int_tmp[i+1],j_int_tmp[i+1]
                                        print, '(ipt,jpt)',ipt,jpt
                                ENDIF
				IF keyword_set(plots) THEN BEGIN
 					makesym,13
 					oplot, [(rgr[ipt]+rgr[ipt+1])*0.5],[(zgr[jpt]+zgr[jpt+1])*0.5],psym=8,symsize=2.0
                                ENDIF
                                dV=(l[i+1]-l[i])*du[k]/(4.0*!pi)*sqrt((pos[0]^2-pos[2]^2)*(1.0+(tan(pos[3]))^2))
                              	vol[jpt+ipt*(ny-1),k]+=dV 						;eq in section 4
                                IF keyword_set(lhat) THEN BEGIN
                                        l_ave=0.5*(l[i]+l[i+1])						;find average l-value
                                	lhat_th=atan((1.0-l_ave)/pos[2]*sqrt(pos[0]^2-pos[2]^2))	;find relative angle
					IF l_ave LT 1.0 THEN BEGIN
                                        	lhat_r_in[jpt+ipt*(ny-1),k]+=dV*mean(-sin(lhat_th)*cos(pos[3]))		;l_r going in
                                                dV_in[jpt+ipt*(ny-1),k]+=dV
                                        ENDIF ELSE BEGIN
                                        	lhat_r_out[jpt+ipt*(ny-1),k]+=dV*mean(-sin(lhat_th)*cos(pos[3]))	;l_r coming out
                                                dV_out[jpt+ipt*(ny-1),k]+=dV
                                        ENDELSE
                                        
                                        lhat_phi[jpt+ipt*(ny-1),k]+=dV*mean(cos(lhat_th)*cos(pos[3]))*sign		;l_phi
					lhat_z[jpt+ipt*(ny-1),k]+=dV*mean(fltarr(2)-sin(pos[3]))			;l_z
                                ENDIF
        	        ENDFOR
                  
                        IF keyword_set(loop_debug) THEN stop
		ENDIF
                	
        ENDFOR
        ;rescale data
        non_zero=where(vol NE 0)
        IF non_zero[0] NE -1 THEN BEGIN
        	vol[non_zero]*=grid.area*1.0e-6
                IF keyword_set(lhat) THEN BEGIN
                	lhat_r_in*=grid.area*1.0e-6
                        lhat_r_out*=grid.area*1.0e-6
                        lhat_phi*=grid.area*1.0e-6
                        lhat_z*=grid.area*1.0e-6
                        dV_in*=grid.area*1.0e-6
                        dV_out*=grid.area*1.0e-6
                ENDIF
        ENDIF
        du*=grid.area*1.0e-6

        IF keyword_set(lhat) THEN lhat={r:{in:lhat_r_in, dVin:dV_in, out:lhat_r_out, dVout:dV_out}, phi:lhat_phi,z:lhat_z}
      
        IF keyword_set(ps) THEN BEGIN
           	device, xsize=float(d_old.x_size)/d_old.x_px_cm,ysize=float(d_old.y_size)/d_old.y_px_cm
        	device, /close
                spawn, 'cp idl.ps /home/'+logname()+'/idl/plots/gpv_grids.ps'
        ENDIF

	IF keyword_set(debug) THEN stop

	IF keyword_set(vol2d) OR keyword_set(contour) THEN BEGIN
		r=fltarr(nx-1)
		FOR i=0,nx-2 DO r[i]=0.5*(rgr[i]+rgr[i+1])
		z=fltarr(ny-1)
		FOR i=0,ny-2 DO z[i]=0.5*(zgr[i]+zgr[i+1])
		vol2d=fltarr(nx-1,ny-1)
		FOR k=0L,npos-1 DO BEGIN
			cntr=0L
			FOR i=0,nx-2 DO BEGIN
				FOR j=0,ny-2 DO BEGIN
					vol2d[i,j]+=vol[cntr,k]
					cntr+=1
				ENDFOR
			ENDFOR
		ENDFOR
		vol2d={vol:vol2d, r:r,z:z}
	ENDIF


	IF keyword_set(contour) THEN IF mean(vol2d.vol) GT 0 THEN BEGIN
                labels={ilab:'Major Radius [m]',jlab:'Z [m]', klab:'Volume [m!u3!n]',ctit:'n_r: '+num2str(n(rgr),1)+' n_z: '$
                	+num2str(n(zgr),1)+' n_view: '+num2str(npos,1),itit:'',jtit:''}
		genplt,vol2d.vol,r,z,ps=ps,labels=labels,prefix='gpv_'
	ENDIF

	IF keyword_set(debug) THEN stop
	ctime=systime(/seconds)-start_time

	RETURN, vol	
END

FUNCTION genpos_planar2umod,info,t
	xyzpos=genpos_det2xyz(info.det.x0,info.det.x1,info.det.x2,transpose([[info.det.xi],[info.det.zeta]]))
        umod=fltarr(n(xyzpos[0,*])+1)+1.0
        FOR i=0,n(xyzpos[0,*]) DO BEGIN
            ;IF abs(xyzpos[1,i]) GT info.ap.size[0]/2.0 THEN umod[i]=1.0-t/info.ap.size[0]*(abs(xyzpos[1,i])-info.ap.size[0]/2.0)/(abs(xyzpos[0,i])-t)
            IF xyzpos[1,i] EQ 0 THEN dist=xyzpos[2,i] ELSE dist=xyzpos[1,i]	;adjust for poloidally viewing arrays
            umod[i]=1.0-t*info.ap.size[1]/(info.ap.size[1]*info.ap.size[0])*abs(dist)/abs(xyzpos[0,i])
            
        ENDFOR
        ;stop
        output=umod
        RETURN, output
END


;+
;NAME:
;	GENPOS_PLANAR2POS
;	
;PURPOSE:
;	This function takes a planar information structure and returns the corresponding POS vectors.
;
;CALLING SEQUENCE:
;	result=GENPOS_PLANAR2POS(info)
;
;INPUTS:
;	info:		STRUC of detector information (see GENPOS_PLANAR_INFO)
;	t:		FLT of the aperture thickness [meters]
;
;KEYWORD PARAMETERS:
;	debug:		/debug stops before the RETURN statement
;
;OUTPUTS:
;	result:		FLTARR 4 x m where m is the number of detector positions in info.det.xi and info.det.zeta
;
;OPTIONAL OUTPUTS:
;	etendue:	FLTARR of length m of the etendue for each detector
;	umod:		FLTARR of length m of the etendue mod for each detector due to finite thickness (u=etendue*umod)
;
;RESTRICTIONS:
;	This function calls GENPOS to generate the POS vectors as well as calculate the etendues.
;
;MODIFICATION HISTORY:
;	Written by:	ML Reinke 7-27-06
;	11-2-07:	ML Reinke - added the t optional input and umod optional output to compute
;                                   finite aperture thickness etendues
;
;-

FUNCTION genpos_planar2pos,info,etendue=etendue,ap_cos=ap_cos,det_cos=det_cos,debug=debug,t=t,umod=umod

	IF strmatch(info.type ,'*planar*') EQ 0 THEN RETURN, -1
	num=n(info.det.xi)+1
	det_pos=fltarr(2,num)
	det_pos[0,*]=info.det.xi
	det_pos[1,*]=info.det.zeta
	pos=genpos(info.ap.vec,info.ap.rot,info.det.x0,info.det.x1,info.det.x2,det_pos,etendue=etendue,$
		a_det=info.det.size[0]*info.det.size[1],a_ap=info.ap.size[0]*info.ap.size[1],ap_cos=ap_cos,det_cos=det_cos)
        IF keyword_set(t) THEN umod=genpos_planar2umod(info,t) ELSE umod=fltarr(num)+1.0
	IF keyword_set(debug) THEN stop
	RETURN,pos
END

FUNCTION genpos_planar2xyz,info
	xyzpos=genpos_det2xyz(info.det.x0,info.det.x1,info.det.x2,transpose([[info.det.xi],[info.det.zeta]]))
	RETURN,xyzpos
END

FUNCTION genpos_planar2th1,info,th_ap=th_ap
	IF NOT keyword_set(th_ap) THEN th_ap=0.0
	ro=info.ap.vec[0]
	xyzpos=genpos_det2xyz(info.det.x0,info.det.x1,info.det.x2,transpose([[info.det.xi],[info.det.zeta]]))
        th=fltarr(n(xyzpos[0,*])+1)
        FOR i=0,n(xyzpos[0,*]) DO BEGIN
		a=sqrt(xyzpos[0,i]^2+xyzpos[1,i]^2+xyzpos[2,i]^2)
                th[i]=acos((a^2-pos[0,i]^2-ro^2)/(-2.0*ro*pos[0,i]))
        ENDFOR
        th-=th_ap
        output=-th
        RETURN, output
END



;+
;NAME:
;	GENPOS_PLANAR2GPV
;
;PURPOSE:
;	This function takes in a planar INFO structure and a grid and uses GENPOS_VOL_COEFS to determine
;	the volume weights on a poloidal grid.
;	
;CALLING SEQUENCE:
;	result=GENPOS_PLANAR2GPV(info,ves_grid)
;
;INPUTS:
;	info:		STRUC of detector information (see GENPOS_PLANAR_INFO)
;	ves_grid:	STRUC of pixel grid boundaries (see GENPOS_GRID) that will be passed to GENPOS_VOL_COEFS
;				which means /center should NOT be used.
;
;OPTIONAL INPUTS:
;	n_ap:		INT/INTARR of aperture griding.  If n_ap is INT then there will be n_ap x n_ap grid (square spacing)
;				If n_ap is INTARR [nx,ny], then there will be n_ap[0]*n_ap[1] grids.  DEFAULT = 3
;	n_det:		INT/INTARR of detector griding (same as n_ap) DEFAULT = 3
;	path:		STR of path where the output will be saved DEFAULT = '/home/username/genie/data/gpv/gpv.dat'
;	lpath:		STR of the path where the lhat optional output, if selected, will be stored.  
;				DEFAULT = '/home/username/genie/data/lhat/lhat.dat'
;	detector:	INT of detector out of the array to for which GENPOS_VOL_COEFS will be run.  This is primarily for debugging purposes.
;	lhat:		Set lhat to a non-zero named variable and it will be filled with the optional output lhat structure	
;
;KEYWORD PARAMETERS:
;	raw:		/raw returns the non-summed output of GENPOS_VOL_COEFS (see optional outputs below)
;	gpv_debug:	/gpv_debug sends a /debug to GENPOS_VOL_COEFS
;	gpv_plots:	/gpv_plots sends a /plots to GENPOS_VOL_COEFS
;	gpv_contour:	/gpv_contour sends a /contour to GENPOS_VOL_COEFS
;	gpv_ps:		/gpv_ps sends a /ps to GENPOS_VOL_COEFS
;	quiet:		/quiet doesn't display status information about current channel being processed
;	debug:		/debug stops before the RETURN statement
;
;OUTPUTS:
;	result:		FLTARR of m x n volumes [m^3] where m is the number of detectors described in the info file and n is the
;				number of grid points (ves_grid.n[0]*ves_grid.n[1]).
;
;OPTIONAL OUTPUTS:
;	ves_cent:	STRUC the ves_cent structure calculated using GENPOS_GRID2CENT that is stored with the GPV and LHAT data.
;	lhat:		STRUC of the average unit vectors in (r,phi,z) (w/o /raw invoked)
;			*.lr	STRUC sub structure for the radial unit vector data
;				*.in	FLTARR [n, m] of the average radial unit vector going in
;				*.dvin	FLTARR [n, m] of the volume weightings going in
;				*.out	FLTARR [n, m] of the average radial unit vector coming out
;				*.dvout	FLTARR [n, m] of the volume weightings coming out
;			*.lphi	FLTARR [n, m] of the average toroidal unit vector
;			*.lz	FLTARR [n, m] of the average vertical unit vector
;
;	===================================================R A W================================================================
;
;	If a /raw is used then the output array will be STRUC that holds the unsummed GPV data.  If the array was stored as a 
;	3D data set the size would be unreasonable, even for a modest number of channels.  Honestly, I crashed IDL cause I ran out
;	of memory. Yeah what's up now punk!
;
;	result: 	STRUC containing gpv for each UPOS for each detector only for non-zero elements
;			*.d0	gpv and spatial information for detector 0
;				*.tmp	LONARR 	[#non_zero] of array elements
;				*.gpv	FLTARR 	[#non_zero,l] of gpv values [1/m^3]
;			*.d1	gpv and spatial information for detector 0
;			  .			 .
;			  .			
;			*.dn	gpv and spatial information for detector n	
;
;	The *.dx.tmp array correspond to array locations (R,Z)=(ves_cent.pnts[0,tmp] and ves_cent.pnts[1,tmp]) while
;	the *.dx.gpv array are the volume coefficients at those points.  Similarly, the lhat structure will be stored in a condensed
;	form when /raw is invoked.
;
;	lhat: 		STRUC containing lhat structures for each UPOS for each detector only for non-zero elements
;			*.d0	lhat and spatial information for detector 0	
;				*.tmp	LONARR 	[#non_zero] of array elements
;				*.lhat	STRUC	of lhat data
;					*.lr	STRUC of unit vectors in the radial (major radius) direction (in and out)
;						*.in	FLTARR [#non_zero] of the average radial unit vector going in
;						*.dvin	FLTARR [#non_zero] of the volume weightings going in
;						*.out	FLTARR [#non_zero] of the average radial unit vector coming out
;						*.dvout	FLTARR [#non_zero] of the volume weightings coming out
;					*.lphi 	FLTARR [#non_zero] average unit vector in the polar angle (toroidal) direction
;					*.lz	FLTARR [#non_zero] average unit vector in the vertical (vertical) direction
;       ========================================================================================================================
;
;RESTRICTIONS:
;	This codes uses GENPOS_GRID to generate the detector and aperture grids, GENPOS_UPOS to process grids and
;	then GENPOS_VOL_COEFS to find the volume coefficents.  Not all optional inputs are carried through to
;	GENPOS_VOL_COEFS.
;
;MODIFICATION HISTORY:
;	Written by:	ML Reinke 7-27-06
;	8-12-06:	ML Reinke - added the RAW output option
;	8-31-06:	ML Reinke - added the detector optional input
;	8-08-07:	ML Reinke - added the ability to get lhat as an optional output and select a
;                                   save path for lhat data.
;	8-08-07:	ML Reinke - now calcualted ves_cent using GENPOS_GRID2CENT and saves with gpv and lhat	
;	8-08-07:	ML Reinke - overhauld how /raw data is stored to save space, now only non_zero values are stored
;-

FUNCTION genpos_planar2gpv,info,ves_grid,n_ap=n_ap,n_det=n_det,path=path,gpv_debug=gpv_debug,gpv_plots=gpv_plots,kdebug=kdebug,$
		gpv_contour=gpv_contour,gpv_ps=gpv_ps,quiet=quiet,debug=debug,raw=raw,detector=detector,lhat=lhat_tmp,lpath=lpath,$
		ves_cent=ves_cent

	IF NOT keyword_set(n_ap) THEN n_ap=3
	IF NOT keyword_set(n_det) THEN n_det=3
	IF NOT keyword_set(path) THEN path='/home/'+logname()+'/idl/genie/data/gpv/gpv.dat'
        IF NOT keyword_set(lpath) THEN lpath='/home/'+logname()+'/idl/genie/data/lhat/lhat.dat' 
	IF strmatch(info.type ,'*planar*') EQ 0 THEN RETURN, -1
	num_det=n(info.det.xi)+1
	num_grid=ves_grid.n[0]*ves_grid.n[1]
        IF NOT keyword_set(detector) THEN BEGIN
        	start=0
                stop=num_det-1
        ENDIF ELSE BEGIN
            	start=detector
                stop=detector
        ENDELSE
        ves_cent=genpos_grid2cent(ves_grid)

        n_ap_bak=n_ap
        n_det_bak=n_det
        x_ap=size(n_ap)
        x_det=size(n_det)
        IF x_ap[0] EQ 0 THEN BEGIN
        	ap_pts=n_ap^2
                det_pts=n_det^2
                tmp=intarr(2,num_det)
                FOR i=0,num_det-1 DO tmp[*,i] = [n_ap,n_ap]
                n_ap=tmp
                tmp=intarr(2,num_det)
                FOR i=0,num_det-1 DO tmp[*,i] = [n_det,n_det]
                n_det=tmp
                redo_apgrid=0
        ENDIF

        IF x_ap[0] EQ 1 THEN BEGIN
        	IF x_ap[1] EQ 2 AND x_det[1] EQ 2 THEN BEGIN
                	ap_pts=n_ap[0]*n_ap[1]
                        det_pts=n_det[0]*n_det[1]
                        tmp=intarr(2,num_det)
                        FOR i=0,num_det-1 DO tmp[*,i] = n_ap
                        n_ap=tmp
                        tmp=intarr(2,num_det)
                        FOR i=0,num_det-1 DO tmp[*,i] = n_det
                        n_det=tmp
                        redo_apgrid=0
                ENDIF 
                IF x_ap[1] EQ num_det AND x_det[1] EQ num_det THEN BEGIN
                	IF keyword_set(raw) THEN RETURN,-2
                        n_ap=[transpose(n_ap),transpose(n_ap)]
                        n_det=[transpose(n_det),transpose(n_det)]
                        redo_apgrid=1
                ENDIF 
        ENDIF

        IF x_ap[0] EQ 2 THEN BEGIN
        	IF keyword_set(raw) THEN RETURN,-2
                redo_apgrid=1
        ENDIF
	IF keyword_set(lhat_tmp) THEN lhatstr=' LHAT ' ELSE lhatstr=''
        IF keyword_set(raw) THEN rawstr=' RAW ' ELSE rawstr=''
        IF keyword_set(raw) THEN BEGIN
         	;not used at this time
                
        ENDIF ELSE BEGIN 
        	gpv=fltarr(num_det,num_grid)
                IF keyword_set(lhat_tmp) THEN BEGIN
                	lhat_r_in=fltarr(num_det,num_grid)
                        v_in=fltarr(num_det,num_grid)
                        lhat_r_out=fltarr(num_det,num_grid)
                        v_out=fltarr(num_det,num_grid)
                        lhat_phi=fltarr(num_det,num_grid)
                        lhat_z=fltarr(num_det,num_grid)
                ENDIF
        ENDELSE

	IF NOT redo_apgrid THEN ap_grid=genpos_grid(info.ap.size[0],info.ap.size[1],n_ap[*,start],/center)
	rot_matrix=genpos_rot_matrix(info.ap.rot[0],info.ap.rot[1],info.ap.rot[2])

        FOR i=start,stop DO BEGIN
		IF NOT keyword_set(quiet) THEN print, 'GPV: Channel - '+num2str(i,1)+' of '+num2str(num_det-1,1)+rawstr+lhatstr
                IF redo_apgrid THEN ap_grid=genpos_grid(info.ap.size[0],info.ap.size[1],n_ap[*,i],/center)
		det_grid=genpos_grid(info.det.size[0],info.det.size[1],n_det[*,i],/center,xo=info.det.xi[i],yo=info.det.zeta[i])
		upos=genpos_upos(info.ap.vec,info.ap.rot,info.det.x0,info.det.x1,info.det.x2,ap_grid,det_grid,rot_matrix=rot_matrix)
                IF keyword_set(gpv_contour) THEN vol2d=1.0 ELSE vol2d=0
                tmp=genpos_vol_coefs(ves_grid,upos.upos,upos.du,debug=gpv_debug,plots=gpv_plots,ps=gpv_ps,contour=gpv_contour,vol2d=vol2d,$
                	kdebug=kdebug,ctime=ctime,lhat=lhat_tmp)
                IF NOT keyword_set(quiet) THEN print, '     CTIME = '+num2str(ctime,dp=2)+' sec'
		IF keyword_set(raw) THEN BEGIN
                    	det_str='d'+num2str(i,1)
                        non_zero=where(sum_array(tmp[*,*],/i) NE 0)		;find values where there is non-zero GPV overall UPOS
                        det_struc={tmp:non_zero, gpv:tmp[non_zero,*]}
                        IF i EQ start THEN gpv={d0:det_struc} ELSE gpv=create_struct(gpv, det_str,det_struc)
                        IF keyword_set(lhat_tmp) THEN BEGIN
                                lhat_r_in=fltarr(num_grid,ap_pts*det_pts)
                                lhat_r_out=fltarr(num_grid,ap_pts*det_pts)
                                lhat_phi=fltarr(num_grid,ap_pts*det_pts)
                                lhat_z=fltarr(num_grid,ap_pts*det_pts)
                            	FOR j=0,ap_pts*det_pts-1 DO BEGIN		;form average for each lhat_XXX vector at each UPOS
                                	non_in=where(lhat_tmp.r.dvin[*,j] NE 0)
                                        IF non_in[0] NE -1 THEN lhat_r_in[non_in,j]=lhat_tmp.r.in[non_in,j]/lhat_tmp.r.dvin[non_in,j]
                                        non_out=where(lhat_tmp.r.dvout[*,j] NE 0)
                                        IF non_out[0] NE -1 THEN lhat_r_out[non_out,j]=lhat_tmp.r.out[non_out,j]/lhat_tmp.r.dvout[non_out,j]
                                        non_all=where(tmp[*,j] NE 0)
                                        IF non_all[0] NE -1 THEN BEGIN
                                        	lhat_phi[non_all,j]=lhat_tmp.phi[non_all,j]/tmp[non_all,j]
                                                lhat_z[non_all,j]=lhat_tmp.z[non_all,j]/tmp[non_all,j]
                                        ENDIF
                                ENDFOR
                                det_struc={tmp:non_zero, lhat:{lr:{in:lhat_r_in[non_zero,*], vin:lhat_tmp.r.dvin[non_zero,*], $
                                	out:lhat_r_out[non_zero,*],vout:lhat_tmp.r.dvout[non_zero,*]}, lphi:lhat_phi[non_zero,*], lz:lhat_z[non_zero,*]}}
                                IF i EQ start THEN lhat={d0:det_struc} ELSE lhat=create_struct(lhat, det_str,det_struc)
                        ENDIF
                ENDIF ELSE BEGIN 
                	gpv[i,*]=sum_array(tmp,/i)
                        IF keyword_set(lhat_tmp) THEN BEGIN
                            	;form lhat_r_in from raw data
                        	lhat_r_in[i,*]=reform(sum_array(lhat_tmp.r.in,/i))
                                v_in[i,*]=reform(sum_array(lhat_tmp.r.dvin,/i))
                                non_zero=where(v_in[i,*] NE 0)
                                lhat_r_in[i,non_zero]/=v_in[i,non_zero]
                                
                                ;form lhat_r_out from raw data
                                lhat_r_out[i,*]=reform(sum_array(lhat_tmp.r.out,/i))
                                v_out[i,*]=reform(sum_array(lhat_tmp.r.dvout,/i))
                                non_zero=where(v_out[i,*] NE 0)
                                IF non_zero[0] NE -1 THEN lhat_r_out[i,non_zero]/=v_out[i,non_zero]
                                
                                ;form lhat_phi and lhat_z for channel i
                                lhat_phi[i,*]=sum_array(lhat_tmp.phi,/i)
                                lhat_z[i,*]=sum_array(lhat_tmp.z,/i)
                                non_zero=where(gpv[i,*] NE 0)
                                lhat_phi[i,non_zero]/=gpv[i,non_zero]
                                lhat_z[i,non_zero]/=gpv[i,non_zero]
                        ENDIF                      
                ENDELSE
                IF keyword_set(debug) THEN stop
        ENDFOR    
        n_det=n_det_bak
        n_ap=n_ap_bak

	output=gpv
	save, gpv,ves_cent,filename=path
        IF keyword_set(lhat_tmp) THEN BEGIN
        	IF NOT keyword_set(raw) THEN lhat={lr:{in:lhat_r_in, vin:v_in, out:lhat_r_out, vout:v_out}, lphi:lhat_phi, lz:lhat_z}
                lhat_tmp=lhat
                save, lhat, ves_cent,filename=lpath
        ENDIF
	IF keyword_set(debug) THEN stop
	RETURN,output
		
 END

;+
;NAME:
;	GENPOS_PLANAR2UPOS
;
;PURPOSE:
;	This function takes in a planar INFO structure and a grid and uses GENPOS_VOL_COEFS to determine
;	the volume weights on a poloidal grid.
;	
;CALLING SEQUENCE:
;	result=GENPOS_PLANAR2UPOS(info)
;
;INPUTS:
;	info:		STRUC of detector information (see GENPOS_PLANAR_INFO)
;
;OPTIONAL INPUTS:
;	n_ap:		INT/INTARR of aperture griding.  If n_ap is INT then there will be n_ap x n_ap grid (even spacing)
;				If n_ap is INTARR [nx,ny], then there will be n_ap[0]*n_ap[1] grids.  DEFAULT = 3
;	n_det:		INT/INTARR of detector griding (same as n_ap) DEFAULT = 3
;	path:		STR of path where the output will be saved DEFAULT = '/home/username/genie/data/gpv/upos.dat'
;
;KEYWORD PARAMETERS:
;	quiet:		/quiet doesn't display status information about current channel being processed
;	debug:		/debug stops before the RETURN statement
;
;OUTPUTS:
;	result:		FLTARR of 4 x n x m position vectors [m^3] where m is the number of detectors described in the info file and n is the
;				number of connecting points (n_ap[0]*n_ap[1]*n_det[0]*n_det[1])
;
;RESTRICTIONS:
;	This codes uses GENPOS_GRID to generate the detector and aperture grids, GENPOS_UPOS to process grids.
;
;MODIFICATION HISTORY:
;	Written by:	ML Reinke 8-18-06 (adapted from GENPOS_PLANAR2GPV)
;	
;-

FUNCTION genpos_planar2upos,info,n_ap=n_ap,n_det=n_det,t=t,clear=clear,path=path,quiet=quiet,debug=debug

	IF NOT keyword_set(n_ap) THEN n_ap=3
	IF NOT keyword_set(n_det) THEN n_det=3
	IF NOT keyword_set(path) THEn path='/home/'+logname()+'/idl/genie/data/gpv/upos.dat'
	
	IF strmatch(info.type ,'*planar*') EQ 0 THEN RETURN, -1
	num_det=n(info.det.xi)+1

        IF n(n_ap) EQ 0 THEN ap_pts=n_ap^2 ELSE ap_pts=n_ap[0]*n_ap[1]
        IF n(n_det) EQ 0 THEN det_pts=n_det^2 ELSE det_pts=n_det[0]*n_det[1]

	ap_grid=genpos_grid(info.ap.size[0],info.ap.size[1],n_ap,/center)
	rot_matrix=genpos_rot_matrix(info.ap.rot[0],info.ap.rot[1],info.ap.rot[2])
        upos=fltarr(4,ap_pts*det_pts,num_det)
        du=fltarr(ap_pts*det_pts,num_det)
        clear=fltarr(ap_pts*det_pts,num_det)+1.0
	FOR i=0,num_det-1 DO BEGIN
		IF NOT keyword_set(quiet) THEN print, 'UPOS: Channel - '+num2str(i,1)+' of '+num2str(num_det,1)
		det_grid=genpos_grid(info.det.size[0],info.det.size[1],n_det,/center,xo=info.det.xi[i],yo=info.det.zeta[i])
		tmp=genpos_upos(info.ap.vec,info.ap.rot,info.det.x0,info.det.x1,info.det.x2,ap_grid,det_grid,rot_matrix=rot_matrix,t=t, ap=info.ap, clear=clear_i)
                upos[*,*,i]=tmp.upos
                du[*,i]=tmp.du
                clear[*,i]=clear_i
                IF keyword_set(debug) THEN stop
	ENDFOR

	upos={upos:upos,du:du}
        output=upos
	save,upos ,filename=path
	IF keyword_set(debug) THEN stop
	RETURN,output
		
 END

;+
;NAME:
;	GENPOS_GRID2RMID
;
;PURPOSE:
;	This function takes a vessel grid and, using EFIT flux surfaces, creates a vector of points
;	on the outboard midplane. 
;
;CALLING SEQUENCE:
;	result = GENPOS_GRID2RMID(grid,shot)
;
;INPUTS:
;	grid:		STRUC 	of grid points generated from GENPOS_GRID using /center
;	shot:		LONINT 	shot number
;
;OPTIONAL INPUTS:
;	tpts:		FLTARR 	of time points that the output is generated for DEFAULT: all EFIT time points
;	efit_times:	FLTARR 	of the the EFIT time points can be input to save on load time.  Otherwise LINE_GETTIMES 
;				is uses to load the EFIT data.
;
;KEYWORD PARAMETERS:
;	sol:		/sol allows values outside the LCFS to be returned as their EFIT_RZ2RMID values and not = -1
;	rho:		/rho returns the points as normalized radius, rho.
;	debug:		/debug stops the code before the RETURN
;	
;OUTPUTS:
;	result:		FLTARR 	n x t where n is the number of grid points (grid.n[0]*grid.n[1]) and t is
;				the number of time points.  Each value of of the array is the return of
;				EFIT_RZ2RMID for that (R,Z) point. Values outside of the LCFS = -1 (unless /sol is invoked)
;
;RESTRICTIONS:
;	This function uses EFIT_RZ2RMID.  In fact this is little more than a script that abuses the privledge that
;	is...EFIT_RZ2RMID.  In fact, most of GENIE would not be possible without that wonderful script.
;
;PROCEDURE:
;	EFIT_RZ2RMID is used to find the locations of all the pixels. This channel map can be created once for each shot,
;	and any number of views which use this gridding can use the same channel map.  Because only one EFIT_RZ2RMID call 
;	is used, then even large grids (100 x 100) can be processed for an entire shot in less than a minute 
;	on one of the Alcator workstations.
;
;
;MODIFICATION HISTORY:
;	Written by:	ML Reinke 7-27-06
;	8-01-06:	ML Reinke - Changed output format to just send back the radius values (EFIT_RZ2RMID output)
;	8-08-07:	ML Reinke - Added the /sol keyword
;
;-

FUNCTION genpos_grid2rmid,grid,shot,tpts=tpts,debug=debug,efit_times=efit_times,rho=rho,sol=sol
	IF NOT keyword_set(efit_times) THEN efit_times=line_gettimes(shot)
	IF NOT keyword_set(tpts) THEN tpts=efit_times
	rmid=line_getrmid(shot)
	num_rmid=n(rmid[0,*])+1
	num_time=n(tpts)+1
	rmidpts=efit_rz2rmid(grid.pnts[0,*],grid.pnts[1,*],tpts,shot=shot,rho=rho)
        IF NOT keyword_set(sol) THEN FOR j=0,num_time-1 DO BEGIN
		IF NOT keyword_set(rho) THEN BEGIN
			rout=rmid[ipt(efit_times,tpts[j]),num_rmid-1]
			tmp=where(rmidpts[*,j] GT rout)
			IF tmp[0] NE -1 THEN rmidpts[tmp,j]=-1
		ENDIF ELSE BEGIN
			tmp=where(rmidpts[*,j] GT 1.0)	
			IF tmp[0] NE -1 THEN rmidpts[tmp,j]=-1
		ENDELSE
	ENDFOR        
	IF keyword_set(debug) THEN stop
	output=rmidpts
	RETURN,output
END



;+
;NAME:
;	GENPOS_GRID_INVESSEL
;	
;PURPOSE:
;	This function takes a vessel grid and outputs an INTARR describing whether each point on the grid
;	is insides (1) or outside (0) of the vessel.  This way, grids can be truncated prior to inversion to remove
;	pixels that would be outside the vessel.
;
;CALLING SEQUENCE:
;	result=GENPOS_GRID_INVESSEL(ves_grid)
;
;INPUTS:
;	ves_cent:	STRUC output of GENPOS_GRID that defines the vessel (use /center)
;	
;OPTIONAL INPUTS:
;	del:		FLT if del is specified then LINE_VESSEL is run with /new at the given del.  Otherwise
;				LINE_VESSEL is run with /load.
;
;KEYWORD PARAMETERS:
;	plot:		/plot will produce a plot of the vessel with filled in circles for the grid points that
;				are inside the vessel and hollow circles for those outside.  This can be used
;				to verify the grid has been properly truncated (suggested).  Grids larger
;				than 3600 points switch to just plotting dots inside the vessel for clarity
;	debug:		/debug stops the code before the RETURN command
;
;OUTPUTS:
;	result:		INTARR of length ves_grid.n[0]*ves_grid.n[1] that contains 1 if the pixel is inside the vessel
;				and 0 if it is outside.  The order is that of ves_grid.pnts.
;
;PROCEDURE:
;	This function uses a tightly meshed contour of the vessel wall (NOT LIMITER) as the boundary points for
;	LINE_INLCFS.  Since LINE_INLCFS was designed for boundaries that are convex to interior points, some hacks
;	have been implimented to make it work with the divertor geometry.  It works pretty well, but since it's
;	not gaurenteed, a /plot should be used when generating the output for a new grid meshing.
;		
;MODIFICATION HISTORY:
;	Written by:	ML Reinke: 8-4-06
;	10-26-06:	ML Reinke - updated so that it works for old/new divertor geometry and in the divertor itself.
;	1-20-06:	ML Reinke - updated so that the /plot feature calls VESSEL_PLOT
;-
	

FUNCTION genpos_grid_invessel,ves_cent,del=del,plot=plot,debug=debug,shot=shot
	IF keyword_set(del) OR keyword_set(shot) THEN ves=line_vessel(/new, del=del,shot=shot) ELSE ves=line_vessel(/load)
	num_pts=long(ves_cent.n[0])*ves_cent.n[1]
	good_pts=intarr(num_pts)

;	edit 10-26-06
;	FOR i=0L,num_pts-1 DO BEGIN
;		IF ves_cent.pnts[1,i] LT -0.2 AND ves_cent.pnts[0,i] LT 0.7 THEN axis=[0.55,-0.35] ELSE axis=[0.68,0.2]
;		IF line_inlcfs(ves.wall.r,ves.wall.z,axis,[ves_cent.pnts[0,i],ves_cent.pnts[1,i]]) EQ 1 THEN IF $
;			(ves_cent.pnts[1,i] LT -0.4 AND ves_cent.pnts[0,i] GT 0.65) EQ 0 THEN good_pts[i] = 1
;                
;	ENDFOR

	FOR i=0L,num_pts-1 DO BEGIN
		IF ves_cent.pnts[1,i] LT -0.35 AND ves_cent.pnts[0,i] LT 0.75 THEN BEGIN
                    	r=ves_cent.pnts[0,i] 
                        z=ves_cent.pnts[1,i]
                        CASE 1 OF
                            	r GT 0.525 AND r LT 0.6 AND z GT -0.48 AND z LE -0.35 : good_pts[i] = 1
                                r GT 0.575 AND r LT 0.615 AND z GE -0.575 AND z LE -0.48 : good_pts[i] = 1
                        	r LE 0.525 AND z GE -0.425 : BEGIN
                                	axis=[0.45,-0.4]
                                        IF line_inlcfs(ves.wall.r,ves.wall.z,axis,[ves_cent.pnts[0,i],ves_cent.pnts[1,i]]) EQ 0 $
                                        	THEN good_pts[i] = 1
                                END
                        	r LE 0.525 AND z LT -0.425 AND z GT -0.48: BEGIN
                                	axis=[0.4,-0.5]
                                        IF line_inlcfs(ves.wall.r,ves.wall.z,axis,[ves_cent.pnts[0,i],ves_cent.pnts[1,i]]) EQ 0 $
                                        	THEN good_pts[i] = 1
                                END
                                r GE 0.6 AND z LE -0.35 AND z GT -0.575 : BEGIN
                                	axis=[0.7,-0.45]
                                        IF line_inlcfs(ves.wall.r,ves.wall.z,axis,[ves_cent.pnts[0,i],ves_cent.pnts[1,i]]) EQ 0 $
                                        	THEN good_pts[i] = 1
                                END
                                r LT 0.6 AND z LE -0.48 AND z GT -0.575 : BEGIN
                                	axis=[0.525,-0.525]
                                        IF line_inlcfs(ves.wall.r,ves.wall.z,axis,[ves_cent.pnts[0,i],ves_cent.pnts[1,i]]) EQ 0 $
                                       		THEN good_pts[i] = 1
                                END
                                ELSE:
                        ENDCASE
                ENDIF ELSE BEGIN 
                	axis=[0.68,0.2]
			IF line_inlcfs(ves.wall.r,ves.wall.z,axis,[ves_cent.pnts[0,i],ves_cent.pnts[1,i]]) EQ 1 THEN IF $
                          (ves_cent.pnts[1,i] LT -0.4 AND ves_cent.pnts[0,i] GT 0.65) EQ 0 THEN good_pts[i] = 1
                ENDELSE
                
        ENDFOR
        tmp=where(ves_cent.pnts[0,*] GT 0.7 AND ves_cent.pnts[1,*] LT -0.5)
        IF tmp[0] NE -1 THEN good_pts[tmp] = 0.0
        tmp=where(ves_cent.pnts[0,*] LT 0.455 AND ves_cent.pnts[1,*] LT -0.45)
        IF tmp[0] NE -1 THEN good_pts[tmp] = 0.0      
	
	IF keyword_set(plot) THEN BEGIN
;		device,window_stat=var
;		IF var[19] EQ 0 THEN wreset,19 ELSE wset,19	
;		plot, ves.wall.r,ves.wall.z
            	vessel_plot,title='GENPOS_GRID_INVESSEL'
		makesym,9
		IF num_pts LE 3600 THEN oplot, ves_cent.pnts[0,*],ves_cent.pnts[1,*],color=200,psym=8
		makesym,10
		IF num_pts GT 3600 THEN use_sym=3 ELSE use_sym=8
		FOR i=0L,num_pts-1 DO IF good_pts[i] EQ 1 THEN oplot, [ves_cent.pnts[0,i]],[ves_cent.pnts[1,i]],color=200,psym=use_sym
		oplot, ves.wall.r,ves.wall.z
	ENDIF

	output=good_pts
	IF keyword_set(debug) THEN stop
	RETURN,output
END

FUNCTION genpos_gpv2vol2d,ves_cent,gpv
	x=size(gpv)
	IF x[0] EQ 1 THEN num_det = 1 ELSE num_det=x[1]	

	genpos_grid2gr,ves_cent,r,z
        nr=n(r)+1
        nz=n(z)+1
	vol2d=fltarr(nr,nz,num_det)
        FOR k=0,num_det-1 DO BEGIN
		FOR i=0,nr-1 DO BEGIN
			FOR j=0,nz-1 DO BEGIN
				vol2d[i,j]+=gpv[k,i*nz+j]
			ENDFOR
		ENDFOR
	ENDFOR
	vol2d={vol:vol2d, r:r,z:z}
	output=vol2d
        RETURN,output
END

;+
;NAME:
;	GENPOS_GPV2COEF_MATRIX
;
;PURPOSE:
;	This function is the rosetta stone!  It takes in a set of volume coefs (output of GENPOS_VOL_COEFS) and
;	a set of rho values corresponding to these (R,Z) points (output of GENPOS_GRID2RMID) and generates
;	the coeffcient matrix to be used in a linear least squares fit to a radial emissivity profile.  Since this process assumes
;	flux surface symmetry, the matrix is generated to force the slope at rho=0 to be zero.
;
;CALLING SEQUENCE:
;	result=GENPOS_GPV2COEF_MATRIX(gpv,rhopts,order)
;
;INPUTS:
;	gpv:	FLTARR 	n x m where n is the number of detectors and m is the number of pixels.  These are the volume coefficents
;				generated using GENPOS_VOL_COEFS (or GENPOS_PLANAR2POS if a planar array was used) [meters^3]
;	rhopts:	FLTARR	m x t where t is the number of time points.  These are normalized radius values that correspond to the m
;				pixels and should be generated using GENPOS_GRID2RMID.
;	order:	INT	order of the polynomial (DEFAULT) or bessel function expansion
;
;KEYWORD PARAMETERS:
;	bessel:		/bessel uses a bessel function of the first kind (BESELJ) as the basis function instead of a polynomial
;	debug:		/debug stops the code before the return statement
;
;OUTPUTS:
;	result:	FLTARR n x order x t where [*,*,i] is the matrix, A, to be used in the least squares inversion of (A*x=b) where
;			b is the vector of length n and x is the vector of length order and are the coeffcients of the
;			See http://en.wikipedia.org/wiki/Linear_least_squares for a brief overview and Numerical Recipes for
;			something more indepth.
;
;PROCEDUE:
;	There will be a note in 'General Equations for Radiometry in Tokamak Plasmas' regarding how the expansion works but
;	for now just talk to ML Reinke directly if you've got questions on the method.
;
;	The expansion for poly is f=a_0+a_1*x^2+x_2*x^3 so that when using POLY on the coefs found, you need to use
;	coefs=[coefs[0],0.0,coefs[1:n(coefs)]].  This means that the expansion has x^order terms.  For the
;	Bessel expansion, only even Bessels are used to fix the slope at rho=0 to be zero.  Also, the Bessel expansion is done in
;	orthogonal functions so f=a_0*J_0(xi_00*rho)+a_1*J_2(xi_20*rho)+a_2*J_4*(xi_40*rho).... in order to force emiss to zero at rho=1
;
;MODIFICATION HISTORY:
;	Written by:	ML Reinke 8-1-06
;	8-4-06:		ML Reinke - adjusted bessel/poly expansion so that there slope at rho=0 is zero.
;	6-7-07:		ML Reinke - made it so the WHERE command selects points where rhopts > 0 and GPV > 0
;
;-


FUNCTION genpos_gpv2coef_matrix,gpv,rhopts,order,debug=debug,bessel=bessel
	
	num_det=n(gpv[*,0])+1
	num_time=n(rhopts[0,*])+1

	coef_matrix=dblarr(num_det,order,num_time)
        bessel_zeros=bessel_zeros()

	FOR i=0,num_det-1 DO BEGIN
		FOR j=0,order-1 DO BEGIN
			FOR k=0,num_time-1 DO BEGIN
				tmp=where(rhopts[*,k] GE 0 AND gpv[i,*] GT 0)
				IF NOT keyword_set(bessel) THEN BEGIN 
                                    	IF j GT 0 THEN jplus=1 ELSE jplus=0
                                	coef_matrix[i,j,k]=total(gpv[i,tmp]*rhopts[tmp,k]^(j+jplus))
                                ENDIF ELSE BEGIN                                      
					coef_matrix[i,j,k]=total(gpv[i,tmp]*beselj(rhopts[tmp,k]*bessel_zeros[2*j],2*j))
                                ENDELSE
			ENDFOR
		ENDFOR
	ENDFOR
	

	output=coef_matrix
	IF keyword_set(debug) THEN stop
	RETURN,output
END	


;+
;NAME:
;	GENPOS_COEFS2PROFILE
;
;PURPOSE:
;	This function converts a set of coefficients determined by programs like GENSPEC_INVERT_MOMENTS into
;	profiles in normalized major radius, rho.
;
;CALLING SEQUENCE:
;	result=GENPOS_COEFS2PROFILE(coefs)
;	
;INPUTS:
;	coefs:		FLTARR of size order x n_time of the coefficents to a polynomial [DEFAULT] or orthogonal
;				bessel function expansion
;
;OPTIONAL INPUTS:
;	nr:		FLT of number of points on which to create profile DEFAULT: 50
;
;KEYWORD PARAMETERS:
;	bessel:		/bessel uses an expansion of even orthogonal bessel functions (see PROCEDURE)
;
;OUTPUTS:
;	result:		FLTARR of size nr x n_time of the profile.
;
;PROCEDURE:
;	Let a = coefs,x=rho m=order-1 and result=y to simplify the explanation
;	For the polynomial 		y=a[0]+0.0*x+a[1]*x^2+a[2]*x^3....+a[m]*x^(m+1)	
;	For the bessel fucntion		y=a[0]*J_0(X_00*x)+a[1]*J_2(X_02*x)....+a[m]*J_(2*m)(X_0(2*m)*x)
;	(X_ij are the zeros)
;
;	These expansions both force the profiles to have zero derivitave on axis which must be true when
;	flux surface symmetry is assumed.  The bessel function expansion forces the profile to be zero
;	at rho=1 which will be valid for peaked or slightly hollow profiles.  For edge dominated profiles
;	an alternative expansion will be necessary
;
;MODIFICATION HISTORY:
;	Written by:	ML Reinke: 8-21-06
;
;-

FUNCTION genpos_coefs2profile,coefs,bessel=bessel,nr=nr
	IF NOT keyword_set(nr) THEN nr=50
	num_time=n(coefs[0,*])+1
	order=n(coefs[*,0])+1
	
	rho=make(0.0,1.0,nr)
	profile=fltarr(nr,num_time)	
	bessel_zeros=bessel_zeros()

        FOR k=0,num_time-1 DO BEGIN
		FOR j=0,order-1 DO BEGIN
			IF NOT keyword_set(bessel) THEN BEGIN 
                        	IF j GT 0 THEN jplus=1 ELSE jplus=0
                               	profile[*,k]+=coefs[j,k]*rho^(j+jplus)
                        ENDIF ELSE BEGIN                                      
				profile[*,k]+=coefs[j,k]*beselj(rho*bessel_zeros[2*j],2*j)
                        ENDELSE
		ENDFOR
	ENDFOR
        
	output=profile
	IF keyword_set(debug) THEN stop
	RETURN,output
END


;+
;NAME:
;	GENPOS_GRID_BFIELD
;
;PURPOSE:
;	This function takes in an output of GENPOS_GRID or GRID_VES and finds the magnitudes of the
;	magnetic field in cyldrical coordinates.  This is used in GENSPEC programs that need these
;	values when invoking incompressibility for the velocity.
;
;CALLING SEQUENCE:
;	result=GENPOS_GRID_BFIELD(gridpts,shot)
;	
;INPUTS:
;	gridpts:	STRUC that is the output of GENPOS_GRID or GRID_VES using /center
;	shot:		LONG shot number
;
;OPTIONAL INPUTS:
;	t_pts:		FLTARR of time points to perform operations DEFAULT: EFIT time points
;
;OUTPUTS:
;	result:		STRUC of data relating to the magnetic field components at the grid points
;			*.br		FLTARR [n_grid, n_time] radial component of the magnetic field [T]
;			*.bphi 		FLTARR [n_grid, n_time] toroidal component of the magnetic field [T]
;			*.bz		FLTARR [n_grid, n_time] veritical component of the magnetic field [T]
;			*.psi_norm	FLTARR [n_grid, n_time] normalized psi values
;			*.t_pts		FLTARR [n_time] time points [s]
;			*.shot		LONG shot number
;
;PROCEDURE:
;	The radial and veritical field componenets are computed by EFIT_RZ2PSI.  The values of psi returned by that function
;	are used, along with PSI0 and PSIBDY to computer psi_norm.  The toroidal field is calculated using the experimental
;	toroidal field and assume the spatial variation is Bphi=Bt0*R0/R.  This ignores the diamagnetic term, but a comparison
;	with Bt=F/R shows this is negligbable.  See BTOR_COMPARE for evidence.
;
;MODIFICATION HISTORY:
;	Written by: 	ML Reinke - 6/6/07
;
;-

FUNCTION genpos_grid_bfield,gridpts,shot,t_pts=t_pts
	
	n_grid=gridpts.n[0]*gridpts.n[1]
	
	mdsopen, 'magnetics', shot
        Bt=abs(mdsvalue("\magnetics::BTOR",/quiet,status=status1))
  	Bt_time=mdsvalue("dim_of(\magnetics::BTOR)",/quiet,status=status2)
	mdsclose,'magnetics',shot
	
	mdsopen,'analysis',shot
	raxis=mdsvalue('\ANALYSIS::TOP:EFIT.RESULTS.G_EQDSK:RMAXIS')
	zaxis=mdsvalue('\ANALYSIS::TOP:EFIT.RESULTS.G_EQDSK:RMAXIS')
	efit_time=mdsvalue('dim_of(\ANALYSIS::TOP:EFIT.RESULTS.G_EQDSK:RMAXIS)')
	psi0=mdsvalue('\ANALYSIS::TOP:EFIT.RESULTS.A_EQDSK:PSI0')
	psibndry=mdsvalue('\ANALYSIS::TOP:EFIT.RESULTS.A_EQDSK:PSIBDY')
	mdsclose,'analysis',shot

	IF NOT keyword_set(t_pts) THEN t_pts = efit_time
	n_tpts=n(t_pts)+1
	Bt0=interpol(bt,bt_time,t_pts)
	R0=interpol(raxis,efit_time,t_pts)
	Z0=interpol(zaxis,efit_time,t_pts)
	psi0=interpol(psi0,efit_time,t_pts)
	psibndry=interpol(psibndry,efit_time,t_pts)
	psi=efit_rz2psi(reform(gridpts.pnts[0,*]),reform(gridpts.pnts[1,*]),t_pts,bz,br,shot=shot)
	
	bphi=fltarr(n_grid,n_tpts)
	psi_norm=fltarr(n_grid,n_tpts)
	FOR i=0, n(t_pts) DO BEGIN
		psi_norm[*,i]=(psi[*,i]-psi0[i])/(psibndry[i]-psi0[i])
		bphi[*,i]=Bt0[i]*R0[i]/gridpts.pnts[0,*]
	ENDFOR

	output={br:br, bphi:bphi, bz:bz, psi_norm:psi_norm, t_pts:t_pts, shot:shot}

	RETURN, output
END


;+
;NAME:
;	GENPOS_LHAT
;
;PURPOSE:
;	This function computes the unit vector components for a single line of sight at
;	certain points on a GENPOS_GRID.
;
;CALLING SEQUENCE:
;	result=GENPOS_LHAT(pos,gridpts)
;
;INPUTS:
;	pos		FLTARR [Ro,Zo,Rt,psi] pos vector paramterizing line of sight
;	gridpts		STRUC of plasma grid.  Use output of GENPOS_GRID or GRID_VES with /center
;
;OPTIONAL INPUTS:
;	good:		INTARR of the points (R,Z) = (gridpts.pnts[0,good], gridpts.pnts[1,good]) that you
;			want the calculation done fore.  DEFAULT is all of gridpts but this shouldn't be
;			used.  Typically good=where(gpv[i,*] GT 0) should input.  See procedure for why.
;
;KEYWORD PARAMETERS:
;	debug:		/debug will stop the code before the output is returned
;	rev:		/rev assumes that the POS vector has been reversed (due to position problems) and *= lz and lphi and lr by -1.0
;	diff:		/diff calculates lr from sqrt(1.0-lz^2-lphi^2).  Good for poloidal arrays
;
;OUTPUTS:
;	result:		STRUC of unit vectors in cylindrical coordinates (R,phi,Z)
;			*.lr	unit vector in the radial (major radius) direction ***see restrictions***
;			*.lphi 	unit vector in the polar angle (toroidal) direction
;			*.lz	unit vector in the vertical (vertical) direction
;
;RESTRICTIONS:
;	The values of result.lr will be valid when pos[3] is well away from 0.  At psi=0 the unit vector is degenerate since
;	on the way in it points inward and on the way out it points outward.  For gridded plasmas there will be a range above and below
;	psi=0 where the view will enter and exit along the same horizontal slice of pixels.
;
;	****NOTE 8/08/07: LHAT VALUES CAN NOW BE CALCULTED USING GENPOS_VOL_COEFS****
;
;PROCEDURE:
;	This function isn't complete but it's getting the job done right now.  If you take the total of the unit vector it will
;	be slightly greater or less than one due to approximations.  This error should decrease as smaller voxel sizes are used.  
;
;
;MODIFICATION HISTORY:
;	Written by: 	ML Reinke - 6/04/07
;	8-08-07:	ML Reinke - fixed a bug that was calculating lr with the wrong sign.
;	2-12-08:	ML Reinke - further patched this crumbling shitty function (rev and diff)
;-

FUNCTION genpos_lhat,pos,gridpts,good=good,debug=debug,rev=rev,diff=diff
	n_grid=gridpts.n[0]*gridpts.n[1]
	IF NOT keyword_set(good) THEN good=indgen(n_grid)
	tmp=good
	IF tmp[0] EQ -1 THEN RETURN,-1
	n_good=n(tmp)+1

	lr=fltarr(n_good)
	rz=fltarr(3,2)
	dist=fltarr(3)
	lrz=fltarr(3)
	FOR i=0, n(tmp) DO BEGIN
		l1=line_s(pos,r=gridpts.pnts[0,tmp[i]])	;will be [outer, inner]
		IF l1[0] EQ -1 THEN lrz[0:1]=[l1,l1] ELSE lrz[0:1]=l1
		rz[0,*]=[gridpts.pnts[0,tmp[i]],line_z(lrz[0],pos)]
		rz[1,*]=[gridpts.pnts[0,tmp[i]],line_z(lrz[1],pos)]
		lrz[2]=line_s(pos,z=gridpts.pnts[1,tmp[i]])
		rz[2,*]=[line_r(lrz[2],pos),gridpts.pnts[1,tmp[i]]]
		dist=sqrt((rz[*,0]-gridpts.pnts[0,tmp[i]])^2+(rz[*,1]-gridpts.pnts[1,tmp[i]])^2)
		l=lrz[minloc(dist)]
		lr[i]=-(1.0-l)/gridpts.pnts[0,tmp[i]]*sqrt(pos[0]^2-pos[2]^2)
	ENDFOR
        lphi=reform(pos[2]/gridpts.pnts[0,tmp]*cos(pos[3]))
        lz=fltarr(n_good)-sin(pos[3])
        
        ;setup shitty patches
        IF keyword_set(rev) THEN sign=-1.0 ELSE sign=1.0
        IF keyword_set(diff) THEN lr=-1.0*sqrt(1.0-lphi^2-lz^2)*sign
        IF keyword_set(rev) THEN BEGIN
        	lz*=sign
                lphi*=sign
        ENDIF

	output={lr:lr,lphi:lphi,lz:lz}
	IF keyword_set(debug) THEN stop
	RETURN,output
END

;+
;NAME:
;	GENPOS_GPV2VOXEL_MATRIX
;
;PURPOSE
;	This function will take a set of GPVs and the rho values of their locations and
;	turn it into a voxel weighting matrix for emissivitiy profile inversion (think L matrix for Abel inversion)
;
;CALLING SEQUENCE
;	result=GENPOS_GPV2VOXEL_MATRIX(gpv, rhopts)
;
;INPUTS:
;	gpv		FLTARR [n_ch, n_pts] of the volume weightings generated using GENPOS_VOL_COEFS
;	rhophts		FLTARR [n_pts,n_time] of the rho locations of the (R,Z) points that locate each GPV at each time point
;	
;OPTIONAL INPUTS:
;	n_rho:		INT number of points from 0.0 -> 1.0 (inclusive) for the result
;	rho_min:	FLT instead of 0.0 the rho values are taken from some input, minimum rho
;	rho_vec:	FLTARR [n_rho] of the rho points that each column of the result
;
;OUTPUTS
;	result:		FLTARR [n_ch, n_rho, n_time] of the voxel weightings for each channel at each rho to be used in profile inversion
;
;OPTIONAL OUTPUTS:
;	rho_vec:	FLTARR if set to an unnamed variable it will be filled with the rho.
;
;
;MODIFIATION HISTORY:
;	Written by:	ML Reinke 6/11/07
;	7-12-07:	ML Reinke - fixed a WHERE bug that would crash if no weighting in the plasma was found
;       7-17-07:	ML Reinke - allowed 2D rhopts so that multiple time slices can be calculated
;	8-02-07:	ML Reinke - added the rho_min optional input
;	1-03-08:	ML Reinke - modified the voxel filling to interpolate between two nearest bins instead of placing all the
;                                   gpv in the nearest bin
;	1-04-08:	ML Reinke - fixed a bug in the previous mod that would cause a bin to goto NAN if mapped rho value EQ a bin rho
;
;-

FUNCTION genpos_gpv2voxel_matrix,gpv,rhopts,rho_vec=rho_vec,n_rho=n_rho,rho_min=rho_min,exp_rho=exp_rho

	IF NOT keyword_set(n_rho) THEN IF keyword_set(rho_vec) THEN n_rho=n(rho_vec) ELSE n_rho=20
        IF NOT keyword_set(rho_min) THEN rho_min=0.0
	IF NOT keyword_set(rho_vec) THEN rho_vec=make(rho_min,1.0,n_rho)
        IF NOT keyword_set(exp_rho) THEN exp_rho=42.5*n_rho/23.0
	n_ch=n(gpv[*,0])+1
	n_grid=n(rhopts[*,0])+1
        n_time=n(rhopts[0,*])+1

        ;stop
	voxel=fltarr(n_ch,n_rho,n_time)
        FOR k=0,n_time-1 DO BEGIN
		FOR i=0,n_ch-1 DO BEGIN
			tmp=where(gpv[i,*] GT 0 AND rhopts[*,k] GT 0)
			IF tmp[0] NE -1 THEN FOR j=0,n(tmp) DO BEGIN
;                            	ipt=ipt(rho_vec,rhopts[tmp[j],k])
;                               ibnd=ipt
                            	rho_j=rhopts[tmp[j],k]
                                ibnd=ibound(rho_vec,rho_j)
                        	IF ibnd[0] NE -1 THEN BEGIN
                                 	del_rho=rho_vec[ibnd[1]]-rho_vec[ibnd[0]]
					del_high=rho_vec[ibnd[1]]-rho_j
                                        del_low=rho_j-rho_vec[ibnd[0]]   	
;                                	voxel[i,ipt(rho_vec,rhopts[tmp[j],k]),k]+=gpv[i,tmp[j]]
                                        IF ibnd[0] NE ibnd[1] THEN BEGIN
                                        	voxel[i,ibnd[0],k]+=gpv[i,tmp[j]]*(1.0-del_low/del_rho)
	                                        voxel[i,ibnd[1],k]+=gpv[i,tmp[j]]*(1.0-del_high/del_rho)
                                        ENDIF ELSE voxel[i,ibnd[0],k]+=gpv[i,tmp[j]]
                                ENDIF ELSE IF rho_j GE rho_vec[n_rho-1] THEN  voxel[i,n_rho-1,k]+=gpv[i,tmp[j]]*exp(-(rho_j-rho_vec[n_rho-1])*exp_rho)
                        ENDFOR
		ENDFOR
        ENDFOR
        ;stop    
	output=voxel
	RETURN,output
END

;+
;NAME:
;	GENPOS_POS2VOXEL_MATRIX
;
;PURPOSE:
;	This function creates a voxel matrix using a pos vector and an etendue vector using an EFIT reconstruction
;
;CALLING SEQUENCE:
;	result=genpos_pos2voxel_matrix(pos,u,shot)
;
;INPUTS:
;	pos:		FLTARR	[4,n_ch] of the [Ro,Zo,Rt,Psi] for each channel 
;	u		FLTARR	[n_ch] of the etendue for each channel [m^2-str]
;	shot:		LONG	of the shot number for EFIT reference
;
;OPTIONAL INPUTS:
;	tpts:		FLTARR	of the time points to calculate the voxel matrix DEFAULT = all EFIT times
;	n_s:		INT	number of points along line of sight to divide view into DEFAULT: 300
;	rzbnd		FLTARR 	[r_min,r_max,|z|] of the rectangular bounding box for the lines DEFAULT: [0.44,1.0,0.45]
;	rho_min		FLT	minimum value of rho to be used DEFAULT: 0.0
;	n_rho		INT	number of rho points to use DEFAULT: 20
;	rho_vec		FLTARR	of the rho points for which the voxel matrix is created: DEFAULT: make(rho_min,1.0,n_rho)
;	rz_ap		FLTARR  [R_ap, Z_ap] of the aperature to start the line of sight at instead of detector [Ro,Zo]
;	r_ap		FLT	of the R value to translate the pos vector to to start the tracing.
;	exp_rho:	FLT	of the scaling past the last rho point to include in the last bin DEFAULT: 1.85*n_rho
;	kdebug:		INT	channel number to stop the voxel filling code at for debugging DEFAULT: OFF
;	ro:		FLT	value of the major radius of the magnetic axis [m] DEFAULT: 0.68
;	a_rho:		FLTARR	[n_rho] of the exponential asymmetry term.  Used in iteration 1 and beyond (see PROCEDURE)
;	m:		INT	of the poloidal m number of the sine/cosine weighting matrix (if selected) DEFAULT: 1
;	tree:		STRING	of the EFIT tree to use for calculating weighting matrix [DEFAULT: 'analysis']
;
;KEYWORD PARAMETERS:
;	debug:		/debug stops the function before the RETURN statement
;	inout:		/inout multiplies each voxel weighting by an asymmetry term (see PROCEDURE)
;	sine:		/sine multiples each voxel weighting by sine (using (R,Z) pt and magnetic axis at time point)
;	cosine:		/cosine multiples each voxel weighting by cosine (using (R,Z) pt and magnetic axis at time point)

;OUTPUTS
;	result:		FLTARR [n_ch, n_rho, n_time] of the voxel weightings for each channel at each rho to be used in profile inversion
;				at each time point.
;
;OPTIONAL OUTPUTS:
;	rho_vec:	FLTARR if set to an unnamed variable it will be filled with the rho vector used if it is not input through rho_vec.
;	rhopts:		FLTARR [n_ch*n_s,n_time] of the rho values for all the channels along the lines of sight for each time
;				slice.  This output is tied to a specific pos,tpts and n_sinput and is useful for iterating on rho_vec.
;
;PROCEDURE:
;	When /inout is invoked without a_rho, the asymmetry term R^2/Ro^2-1.0) is multiplied into the volume term, which is
;	used to derive a_rho.  If a_rho is given and /inout set, then the asymmetry term exp(a_interp*(r_ch[j]^2/ro^2-1.0)) 
;	is used which will allow a better first order emissivity profile to be calculated.
;
;MODFICATION HISTORY:
;	Written by:	ML Reinke 12/21/07
;	1-03-08:	ML Reinke - updated the voxel filling to interpolate the volume element to the two nearest bins.
;					added the exp_rho optional input to handle the edge bin
;	1-04-08:	ML Reinke - fixed a bug in the previous mod that would cause a bin to goto NAN if mapped rho value EQ a bin rho
;	1-09-08:	ML Reinke - added the rhopts optional output so iterations can be done on rho_vec and save time
;	1-31-08:	ML Reinke - added the ability to fit a rotation induced in/out asymmetry (see inout, ro, a_rho)
;	7-15-08:	ML Reinke - added the ability to generate sine or cosine asymmetry weighted voxel matrices
;	8-12-10:	ML Reinke - added the r_ap optional input which allows the n_s to start from a defined major radius rather than the aperature.
;   	4-21-11:	ML Reinke - updated the verb output to specify the type of matrix
;	6-13-11:	ML Reinke - added the tree optional input for use with EFIT_RZ2XXX codes.	
;
;-

FUNCTION genpos_pos2voxel_matrix,pos,u,shot,tpts=tpts,n_s=n_s,rzbnd=rzbnd,rho_vec=rho_vec,n_rho=n_rho,rho_min=rho_min,rz_ap=rz_ap,r_ap=r_ap,exp_rho=exp_rho,$
		kdebug=kdebug,debug=debug,verb=verb,rhopts=rhopts,inout=inout,ro=ro,a_rho=a_rho,sine=sine,cosine=cosine,m=m,psinorm=psinorm,tree=tree

	IF NOT keyword_set(n_rho) THEN IF keyword_set(rho_vec) THEN n_rho=n(rho_vec)+1 ELSE n_rho=20
        IF NOT keyword_set(rho_min) THEN rho_min=0.0
	IF NOT keyword_set(rho_vec) THEN rho_vec=make(rho_min,1.0,n_rho)
	IF NOT keyword_set(rzbnd) THEN rzbnd=[0.44,1.0,0.45]	;set path's [r_min,r_max,abs(z_max)] limit
	IF NOT keyword_set(n_s) THEN n_s=300
	IF NOT keyword_set(tpts) THEN tpts=line_gettimes(shot)
	IF NOT keyword_set(exp_rho) THEN exp_rho=42.5*n_rho/23.0
        IF NOT keyword_set(ro) THEN ro=0.68
        IF NOT keyword_set(m) THEN m=1.0
	n_ch=n(pos[0,*])+1
	n_time=n(tpts)+1
	n_s=long(n_s)
        
        ;find the magnetic axis if sine or cosine is set so theta can be found
        IF keyword_set(cosine) OR keyword_set(sine) THEN BEGIN
        	efit_axis=line_getaxis(shot)
                efit_time=line_gettimes(shot)
                IF total(efit_time-tpts) EQ 0 THEN BEGIN	;skip interp if using all times
                	axis_ro=reform(efit_axis[*,0])
                	axis_zo=reform(efit_axis[*,1])
                ENDIF ELSE BEGIN
                	axis_ro=interpol(reform(efit_axis[*,0]),efit_time,tpts)
                        axis_zo=interpol(reform(efit_axis[*,1]),efit_time,tpts)
                ENDELSE
        ENDIF

	;intiate (r,z) and s points as vectors for fast EFIT_RZ2RMID and determine r_z points
	pos_r=fltarr(n_ch*n_s)
	pos_z=fltarr(n_ch*n_s)
	pos_s=fltarr(n_ch*n_s)
	FOR i=0,n_ch-1 DO BEGIN
		;determine where the line of sight terminates in s due to radial boundries
		IF pos[2,i] LT rzbnd[0] THEN BEGIN	;if r_tang < r_min
			svec=line_s(pos[*,i],r=rzbnd[0])
			smin_r=svec[1]	 ;choose the smallest s where r=rmin
		ENDIF ELSE BEGIN			;if r_tang > r_min
			svec=line_s(pos[*,i],r=rzbnd[1])
			smin_r=svec[0]	 ;choose the largest (and non-negative) s where r=rmax
 		ENDELSE

		;determine where the line of sight terminates in s due to z boundries
		IF pos[3,i] NE 0 THEN BEGIN
			IF line_s(pos[*,i],z=rzbnd[2]) LT 0 THEN $
				smin_z=line_s(pos[*,i],z=-1.0*rzbnd[2]) ELSE $
				smin_z=line_s(pos[*,i],z=rzbnd[2])
		ENDIF ELSE smin_z=2.0*smin_r	;if no inclination, smin_z set > smin_r	
		IF smin_r GT smin_z THEN smin = smin_z ELSE smin=smin_r		;set minimum s
                s_ap=0
		IF keyword_set(rz_ap) THEN s_ap=line_s(pos[*,i],z=rz_ap[1]) > 0
                IF keyword_set(r_ap) THEN s_ap=min(line_s(pos[*,i],r=r_ap))> 0
		s=make(s_ap,smin,n_s)	;n_s points going from the aperture to the boundary intersecion
		pos_s[i*n_s:(i+1)*n_s-1]=s
		pos_r[i*n_s:(i+1)*n_s-1]=line_r(s,pos[*,i])
		pos_z[i*n_s:(i+1)*n_s-1]=line_z(s,pos[*,i])
        ENDFOR
        IF keyword_set(psinorm) THEN inv_str='RZ2RHO' ELSE inv_str='RZ2RMID'
        IF keyword_set(verb) THEN print, '(R,Z) points found, calling '+inv_str
	IF NOT keyword_set(rhopts) THEN BEGIN
        	IF keyword_set(psinorm) THEN rhopts=efit_rz2rho(pos_r,pos_z,tpts,shot=shot,/psinorm,tree=tree) ELSE rhopts=efit_rz2rmid(pos_r,pos_z,tpts,shot=shot,/rho,tree=tree)
        ENDIF
        IF keyword_set(verb) THEN print, inv_str+' done, filling voxel matrix'
        IF keyword_set(cosine) THEN print, 'filling cosine m='+num2str(m,1)
        IF keyword_set(sine) THEN print, 'filling sine m='+num2str(m,1)
	voxel=fltarr(n_ch,n_rho,n_time)
        FOR k=0,n_time-1 DO BEGIN
		FOR i=0,n_ch-1 DO BEGIN
			rho_ch=rhopts[i*n_s:(i+1)*n_s-1,k]
			s_ch=pos_s[i*n_s:(i+1)*n_s-1]
			r_ch=pos_r[i*n_s:(i+1)*n_s-1]
			z_ch=pos_z[i*n_s:(i+1)*n_s-1]
			IF keyword_set(kdebug) THEN if i EQ kdebug THEN stop
			FOR j=1,n_s-2 DO BEGIN			;using +/- points to determine delta_s 
                                ;ipt=ipt(rho_vec,rho_ch[j])
                                ibnd=ibound(rho_vec,rho_ch[j])
                                del_s=0.5*(s_ch[j+1]-s_ch[j-1])	
				vol=u[i]/(4.0*!pi)*sqrt(((pos[0,i])^2-(pos[2,i])^2)*(1.0+(tan(pos[3,i]))^2))*del_s
                                asym_term=1.0
                        	IF ibnd[0] NE -1 THEN BEGIN
                                	IF keyword_set(inout) THEN BEGIN
                                		IF keyword_set(a_rho) THEN BEGIN 
                                                	m_a=(a_rho[ibnd[1]]-a_rho[ibnd[0]])/rho_vec[ibnd[1]]-rho_vec[ibnd[0]]
                                                	a_interp=a_rho[ibnd[0]]+m_a*(rho_ch[j]-rho_vec[ibnd[0]])
                                                	asym_term=exp(a_interp*(r_ch[j]^2/ro^2-1.0)) ;use to get better 1st order emissivity
	                                        ENDIF ELSE asym_term=(r_ch[j]^2/ro^2-1.0) ;use to get a_rho
                                        ENDIF
                                        IF keyword_set(cosine) OR keyword_set(sine) THEN BEGIN
                                        	deltaR=r_ch[j]-axis_ro[k]
                                                deltaZ=z_ch[j]-axis_zo[k]
                                                IF deltaR GE 0 THEN BEGIN
                                                	IF deltaZ GE 0 THEN sign=1.0 ELSE sign=-1.0
                                                        IF deltaZ GE 0 THEN dth=0 ELSE dth=2.0*!pi
                                                ENDIF
                                                IF deltaR LT 0 THEN BEGIN
                                                	IF deltaZ GE 0 THEN sign=-1.0 ELSE sign=1.0
                                                        IF deltaR LT 0 THEN dth=!pi
                                                ENDIF
                                                th=sign*atan(abs(deltaZ/deltaR))+dth
                                                ;length=sqrt(deltaR^2+deltaZ^2)
                                                ;IF keyword_set(sine) THEN asym_term=(deltaZ/length)
                                                ;IF keyword_set(cosine) THEN asym_term=(deltaR/length)
                                            	IF keyword_set(sine) THEN asym_term=sin(m*th)
                                                IF keyword_set(cosine) THEN asym_term=cos(m*th)
                                        ENDIF
					del_rho=rho_vec[ibnd[1]]-rho_vec[ibnd[0]]
					del_high=rho_vec[ibnd[1]]-rho_ch[j]
					del_low=rho_ch[j]-rho_vec[ibnd[0]]
					;voxel[i,ipt[0],k]+=vol
                                        IF ibnd[0] NE ibnd[1] THEN BEGIN
						voxel[i,ibnd[0],k]+=vol*(1.0-del_low/del_rho)*asym_term
						voxel[i,ibnd[1],k]+=vol*(1.0-del_high/del_rho)*asym_term
                                        ENDIF ELSE voxel[i,ibnd[0],k]+=vol*asym_term
                                ENDIF ELSE IF rho_ch[j] GE rho_vec[n_rho-1] THEN voxel[i,n_rho-1,k]+=vol*$
                                	exp(-(rho_ch[j]-rho_vec[n_rho-1])*exp_rho)*asym_term ;this will have a forced asym_term=1.0
                        ENDFOR
                        IF keyword_set(debug) THEN stop
                        IF i EQ 0 AND k EQ 0 AND keyword_set(verb) THEN print, 'voxel matrix filled for ch '+num2str(i,1)+' of '+num2str(n_ch-1,1)
                        IF i EQ 0 AND k EQ 1 AND keyword_set(verb) THEN print, 'voxel matrix filled for time '+num2str(k,1)+' of '+num2str(n_time-1,1)
		ENDFOR
        ENDFOR
	output=voxel
	IF keyword_set(debug) THEN stop
	RETURN,output
END


;+
;NAME:
;	GENPOS_PROFILE_INVERT:
;	
;PURPOSE:
;	This function inverts a 1D moment profile given the moment, the voxel weighting matrix and
;	a smoothing parameter to get the local moment/volume.  Basically the ill-conditioned b = V*a is inverted
;	to calculate a given measurements of b using linear-regularization.
;
;CALLING SEQUENCE:
;	result=GENPOS_PROFILE_INVERT(mom,voxel,eps)
;
;INPUTS:
;	mom:		FLTARR of length m of "brightnesses" to be inverted
;	voxel:		FLTARR of size m x n of the spatial weighting coefficients [inside-->outside]
;	eps:		FLT weighting of the regularization matrix
;
;OPTIONAL INPUTS:
;	eta:		FLT weighting coefficient for the edge zero DEFAULT: 0.0 
;	err:		FLTARR of length m of the statistical error in the input mom
;	n_err_inv:	INT number of inversions to complete when caclulation error propigation DEFAULT: 1000 [NO LONGER DOES ANYTHING 2-16-09]
;	nprof:		INT the number of spatial profiles that voxel is describing (radial+sine+cosine for example).  
;				This will adjust regularization matrices accordingly.
;	inv_matrix:	DBLARR of the inverted least squares/smoothing matrix
;
;KEYWORD PARAMETERS:
;	svdc:		/svdc inverts using singular value decomposition (SSVDC) where the default is LA_INVERT
;	nofirst:	/nofirst removes the weighting of the weighting of the first derivative to be used when profile does not goto rho=0.
;
;OUTPUTS:
;	result:		FLTARR of length n of the "emissivity" from the inversion
;
;OPTIONAL OUTPUTS:
;	svdout:		STRUC if /svdc is called then this struc is filled with the outputs of SSVDC
;	brchk:		STRUC of the brightness check
;				*.mom FLTARR of length m of voxel#result which is a check on the quality of the inversion
;				*.inverr FLTARR [n] of the emissivity error using propigation of error on the matrix equation.
;				  DEFAULT will be zeros if optional input err is not set
;	inv_matrix:	DBLARR of the inverted least squares/smoothing matrix
;
;PROCEDURER:
;	This is routine method for inverting noisy profile that have suffcient spatial sampling of a profile.  It is assumed
;	that the output emissivity profile will be in rho or psi so that the derivative on axis should be zero.  Use
;	GENPOS_POS2VOXEL_MATRIX, GENPOS_GPV2VOXEL_MATRIX or GENPOS_GPV2VOXEL_VEL_MATRIX to general voxel.  The correct ordering is for the
;	lower values of mom to be from the innermost choords.
;
;	Inversion errors can be calculated if the error of the mom "brightness" profile is included.  Many (n_err_inv) inversions
;	are done with the same length/weighting matrix varying the mom profile by random noise who magnitude is perscribed by the err optional input.
;
;MODIFICATION HISTORY:
;	Written by: 	ML Reinke - June 2007
;	7-17-07:	ML Reinke - added the inv_matrix optional output for use in GENSPEC if EPS_EM and EPS_TI are the same.
;	8-15-07:	ML Reinke - added the edge zero controlled by weighting factor eta
;	10-30-08:	ML Reinke - added the err and n_err_inv optional inputs and the inverr
;                                   output in the brchk structure. These allow statistical error to be calculated.
;	11-04-08:	ML Reinke - added the nprof optional input which adjust the regularization
;                                   matrix to reflect multiple spatial profiles are being fit at the same time.                                  
;	2-16-08:	ML Reinke - switched the inversion error calculation over to propigation error from doing multiple inversions
;	4-21-11:	ML Reinke - added the ability to use the inv_matrix as an optional input so user can speed up inversions
;
;-

FUNCTION genpos_profile_invert,mom,voxel,eps,svdout=svdout,brchk=brchk,verb=verb,svdc=svdc,nofirst=nofirst,inv_matrix=inv_matrix,eta=eta,err=err,n_err_inv=n_err_inv,nprof=nprof

        IF keyword_set(verb) THEN print, ' Performing inversion with smoothing'
        IF NOT keyword_set(n_err_inv) THEN n_err_inv=1000
        IF NOT keyword_set(nprof) THEN nprof=1.0
	x=size(voxel)
        npts=x[2]/nprof

	;2nd derivative matrix over whole profile
	vec=[-1.0,2.0,-1.0]
	d=fltarr(npts-2,npts)
	FOR i=0,npts-3 DO d[i,i:i+2]=vec

	;1st derivative matrix at inner two points
	IF keyword_set(nofirst) THEN vec = [0.0,0.0] ELSE vec=[-1.0,1.0]
	f=fltarr(npts-1,npts)
	FOR i=0,1 DO f[i,i:i+1]=vec

        ;make edge zero
        IF NOT keyword_set(eta) THEN eta=0.0
        z=fltarr(npts,npts)
        z[npts-1,npts-1]=1.0
        u_o=fltarr(npts)

        ;copy the weighting matrices to cover the whole profile
	IF nprof GT 1.0 THEN BEGIN
                d_tot=fltarr((npts-2)*nprof,npts*nprof)
                FOR i=0,nprof-1 DO d_tot[i*(npts-2):(i+1)*(npts-2)-1,i*npts:(i+1)*npts-1]=d
                f_tot=fltarr((npts-1)*nprof,npts*nprof)
                FOR i=0,nprof-1 DO f_tot[i*(npts-1):(i+1)*(npts-1)-1,i*npts:(i+1)*npts-1]=f
                z_tot=fltarr(npts*nprof,npts*nprof)
                FOR i=0,nprof-1 DO z_tot[i*npts:(i+1)*npts-1,i*npts:(i+1)*npts-1]=z
                d=d_tot
                f=f_tot
                z=z_tot
                u_o=fltarr(npts*nprof)
        ENDIF
  
	;create the least-squares + second derivative matrix + first derivative
	mat = (transpose(voxel)#voxel + eps*transpose(d)#d+eps*transpose(f)#f+eta*transpose(z)#z)

	;invert directly or using SVDC
	IF NOT keyword_set(svdc) THEN BEGIN
        	IF NOT keyword_set(inv_matrix) THEN inv_matrix=la_invert(mat,/double)
                emiss=inv_matrix#(transpose(voxel)#mom+eta*z#u_o)
                IF keyword_set(err) THEN BEGIN
                    	IF keyword_set(verb) THEN print, 'Calculating Inversion Error Using '+num2str(n_err_inv,1)+' Iterations'
                	;emiss_err=fltarr(n(emiss)+1,n_err_inv)
                        ;inverr=fltarr(n(emiss)+1)
                        ;FOR i=0,n_err_inv-1 DO BEGIN
                        	;mom_err=mom+randomn(seed,n(mom)+1)*err
                                ;emiss_err[*,i]=inv_matrix#(transpose(voxel)#mom_err+eta*z#u_o)
                        ;ENDFOR
                        ;FOR i=0,n(emiss) DO inverr[i]=stdev(emiss_err[i,*])
                    	errmatrix=inv_matrix#transpose(voxel)
                        inverr=fltarr(n(emiss)+1)
                        FOR i=0,n(emiss) DO inverr[i]=sqrt(total((errmatrix[i,*]*err)^2))
                ENDIF ELSE inverr=fltarr(n(emiss)+1)
        ENDIF ELSE BEGIN

		;do the singular value decomposition
		ssvdc,mat,w,u,v

		;obtain the inverse matrix
		w_recip=fltarr(npts,npts)
		FOR i=0,npts-1 DO IF w[i] NE 0 THEN w_recip[i,i]=(1.0/w[i])
		IF abs(min(w)/max(w)) LT 1.0e-5 THEN print, 'Inversion is ill-conditioned, |w_min/w_max| = '$
			+num2str(abs(min(w)/max(w)))
		emiss= v#w_recip#transpose(u)#(transpose(voxel)#mom+eta*z#u_o)
		svdout={u:u,v:v,w:w}
	ENDELSE
        brchk={mom:voxel#emiss,inverr:inverr}

	RETURN,emiss
END


;+
;NAME:
;	GENPOS_EMISS_INVERT
;
;PURPOSE:
;	This function is an upper level function to calculate the emissivity profile from an array of power deposited
;	and and GPV values.
;
;CALLING SEQUENCE:
;	result=GENPOS_EMISS_INVERT(power,gpv,shot,time,ves_cent)
;
;INPUTS:
;	power:		FLTARR 	[n_ch, n_time] of power deposited (ph/s or Watts) NOTE: this is NOT brightness (W/m^2)
;	gpv:		FLTARR 	[n_ch, n_pnts] of volume weightings for each channel [m^3] (see GENPOS_VOL_COEFS)
;	shot:		LONG	shot number
;	time:		FLTARR	[n_time] of time scale of powers [seconds]
;	ves_cent	STRUC	of information for the voxel gridding for which gpv was calculated.
;	
;OPTIONAL INPUTS:
;	good:		INTARR	[n_ch] with 1's (use) or 0's (do not use) indicating which channels to use during the inversion
;	rhopts:		FLTARR	[n_grid, n_time] of the rho values corresponding to the (R,Z) grid points.  This will be calculated
;				by GENSPEC_MATRIX_INVERT but is an optional output that can be reused on different spectral lines to
;				save computation time.
;	n_rho:		INT	number of points from minimum < rho < 1.0 (inclusive) to be used in inversion DEFAULT: 25
;	rho_vec:	FLTARR 	the actually rho points to be used DEFAULT: not used35
;	eta:		FLT	edge zero weighting factor DEFAULT: 0.0
;	eps:		FLT	 profile smoothing factor. DEFAULT: 1.0
;
;KEYWORD PARAMETERS:
;	nofirst:	/nofirst prevents the matrix inversion from weighting the profile so that it has zero derivative at rho=0
;	debug:		/debug stops the code in various places and just before the RETURN statement
;	quiet:		/quiet suppresses terminal messages displaying computation times of various procedures
;
;OUTPUTS:
;	result:		STRUC	containing the kinetic profiles and their moment profile checks
;			*.rho		FLTARR 	[n_rho] of the rho values at which each profile has been calculated
;			*.time		FLTARR 	[n_time] of the time points (copy of time input for convience)
;			*.emiss 	FLTARR 	[n_rho, n_time] of the emissivity "power"/m^3 [ph/s/m^3 or W/m^3]
;			*.ch		INTARR 	[n_ch] of only the GOOD channels used in the inversion
;			*.brchk		FLTARR 	[3, n_ch, n_time] of the moments calculated from the calculated profiles.  This can be
;						be used to check the quality of the output versus the input data.
;
;PROCEDURE:
;	This program calculates the spatial weighting matrices using GENPOS_GPV2VOXEL_MATRIX and then inverts them using 
;	GENPOS_PROFILE_INVERT.  The eps value at each time point are weighted to the max value of
;	the voxel matrix.  For example, the GENPOS_PROFILE_INVERT eps at each time point i is: 
;		
;		double(eps_em*(max(voxel[*,*,i]))^2)
;
;RESTRICTIONS:
;	Use GENIE_INI.bat to compile the GENPOS functions in the correct order.  You should know what you're doing before
;	using this function. 
;
;MODFICATION HISTORY:
;	Written by:	ML Reinke - 8/21/07
;
;-

FUNCTION genpos_emiss_invert,power,gpv,shot,time,ves_cent,good=good,rhopts=rhopts,n_rho=n_rho,rho_vec=rho_vec,eta=eta,$
		eps=eps,nofirst=nofirst,debug=debug,quiet=quiet,sol=sol,err=err

	n_ch=n(gpv[*,0])+1
	n_time=n(time)+1
	ch=indgen(n_ch)

	IF NOT keyword_set(eps) THEN eps=1.0
        IF NOT keyword_set(eta) THEN eta=0.0

	IF NOT keyword_set(n_rho) THEN IF NOT keyword_set(rho_vec) THEN n_rho=22 ELSE n_rho=n(rho_vec)+1
	IF NOT keyword_set(good) THEN good=intarr(n_ch)+1
	start_time=systime(/seconds)
	IF NOT keyword_set(rhopts) THEN rhopts=genpos_grid2rmid(ves_cent,shot,tpts=time,/rho,sol=sol)
	ctime=systime(/seconds)-start_time
	IF NOT keyword_set(quiet) THEN print, 'VES_CENT converted to RHOPTS: '+num2str(ctime,dp=2)
	IF keyword_set(debug) THEN stop
	tmp=where(good EQ 1)
        IF tmp[0] EQ -1 THEN BEGIN
        	print, 'No good channels selected - FAILURE'
                RETURN,-1
        ENDIF
	n_good=n(tmp)+1
	good_ch=ch[tmp]
        rho_min=1.0
        FOR i=0,n_good-1 DO BEGIN
            tmp2=where(gpv[tmp[i],*] GT 0 AND rhopts GT 0)
            IF tmp2[0] NE -1 THEN rho_min = rho_min < min(rhopts[tmp2])
        ENDFOR
        IF NOT keyword_set(quiet) THEN print, 'RHO values from '+num2str(rho_min,dp=3)+' < rho < '+num2str(max(rhopts),dp=3)+' selected'
	start_time=systime(/seconds)
	voxel=genpos_gpv2voxel_matrix(gpv[tmp,*],rhopts,rho_vec=rho_vec,n_rho=n_rho,rho_min=rho_min)
	ctime=systime(/seconds)-start_time
	IF NOT keyword_set(quiet) THEN print, 'VOXEL arrays calculated: '+num2str(ctime,dp=2)
	IF keyword_set(debug) THEN stop

	emiss=fltarr(n_rho,n_time)
	brchk=fltarr(n_good,n_time)
        error=fltarr(n_rho,n_time)
        

	start_time=systime(/seconds)
	FOR i=0,n_time-1 DO BEGIN
		a=max(voxel[*,*,i])
		emiss[*,i]=genpos_profile_invert(power[tmp,i],voxel[*,*,i],double(eps*a^2),brchk=bright,nofirst=nofirst,eta=double(eta*a^2),err=err)
		brchk[*,i]=bright.mom
                error[*,i]=bright.inverr
	ENDFOR
	ctime=systime(/seconds)-start_time
	IF NOT keyword_set(quiet) THEN print, 'inversions completed for '+num2str(n_time)+' time points: '+num2str(ctime,dp=2)
	output={rho:rho_vec,time:time,emiss:emiss,ch:good_ch,brchk:brchk,err:error}
	IF keyword_set(debug) THEN stop
	RETURN,output	
END


;+
;NAME:
;	GENPOS_LINE_BR
;
;PURPOSE:
;	This function calculates the line-integrated brightness for an emissivity profile and given set of POS vectors.  This
;	duplicates, in part the LINE_BR functionality.
;
;CALLING SEQUENCE:
;	result=GENPOS_LINE_BR(pos,emiss,r,t,shot,tau)
;
;INPUTS:
;	pos	FLTARR	[4,nch] of the POS vectors
;	emiss	FLTARR	[nr,nt] of the emissivity [X/m^3]
;	r	FLTARR	[nr] of the radial scale vector [m]
;	t	FLTARR	[nt] of the time scale [sec]
;	shot	LONG	shot number
;	tau	FLTARR	[ntau] of the time points at which to calculate the brightness [sec]
;
;OPTIONAL INPUTS:
;	emerr	FLTARR	[nr,nt] of the uncertainty in emissivity [X/m^3]
;	n_s	INT	number of points to discretize along line of sight DEFAULT=75
;	rmap	FLTARR	of the mapping data to reinsert for different emissivities
;	efit	STRUC	of EFIT data to resinsert for different emissivities
;	
;KEYWORD PARAMETERS:
;	rho	/rho assums that the r input is actually r/a
;	verb	/verb displays (a lot) of status messages
;	debug	/debug stops the code at various points for inspection
;	
;OUTPUTS:
;	result	FLTARR	[nch,nt] of the line-integrated brightness [X/m^2]
;	
;OPTIONAL OUTPUTS:
;	brerr	FLTARR	[nch,nt] of the uncertainty in result [X/m^2]
;	efit	STRUC	of EFIT data to be used as optional input
;	rmap	FLTARR	of mapping data to be used as optional input
;
;MODIFICATION HISTORY:
;	Written by:	M.L. Reinke - 4/18/12 (adapted from LINE_BR)
;
;-

FUNCTION genpos_line_br,pos,emiss,r,t,shot,tau,emerr=emerr,n_s=n_s,verb=verb,debug=debug,brerr=brerr,rmap=rmap,efit=efit,rho=rho

	IF NOT keyword_set(n_s) THEN n_s=75	;		set number of segments to define line of sight
	IF NOT keyword_set(rzbnd) THEN rzbnd=[0.44,1.0,0.5]	;set path's [r_min,r_max,abs(z_max)] limit
	IF NOT keyword_set(emerr) THEN emerr=emiss*0.0
	
	npos=n(pos[0,*])+1	;determine the number of views
	ntau=n(tau)+1		;determine the number of time points
	br=fltarr(npos, ntau)	;initialize output brightness
	brerr=fltarr(npos,ntau)	;initialize output uncertainty [filled if emerr specified]

	;load EFIT data or decomile structure
	IF NOT keyword_set(efit) THEN BEGIN
		efit_time=line_gettimes(shot)
		efit_lcfs=line_getlcfs(shot)
		efit_axis=line_getaxis(shot)
		efit_rmid=line_getrmid(shot)
		efit={time:efit_time,lcfs:efit_lcfs,axis:efit_axis,rmid:efit_rmid}
	ENDIF ELSE BEGIN
		efit_rmid=efit.rmid
		efit_time=efit.time
		efit_lcfs=efit.lcfs
		efit_axis=efit.axis
	ENDELSE
	
	;debug information
	IF keyword_set(verb) THEN print, ' npos = '+num2str(npos,1)+'  ntau = '+num2str(ntau,1)
	IF keyword_set(debug) THEN stop

	;start filling the brightness array cycling through each view
	FOR i=0,npos-1 DO BEGIN
		IF keyword_set(verb) THEN print, ' i = '+num2str(i,1)

		;determine where the line of sight terminates in s due to radial boundries
		IF pos[2,i] LT rzbnd[0] THEN BEGIN	;if r_tang < r_min
			svec=line_s(pos[*,i],r=rzbnd[0])
			smin_r=svec[1]	 ;choose the smallest s where r=rmin
		ENDIF ELSE BEGIN			;if r_tang > r_min
			svec=line_s(pos[*,i],r=rzbnd[1])
			smin_r=svec[0]	 ;choose the largest (and non-negative) s where r=rmax
 		ENDELSE

		;determine where the line of sight terminates in s due to z boundries
		IF pos[3,i] NE 0 THEN BEGIN
			IF line_s(pos[*,i],z=rzbnd[2]) LT 0 THEN $
				smin_z=line_s(pos[*,i],z=-1.0*rzbnd[2]) ELSE $
				smin_z=line_s(pos[*,i],z=rzbnd[2])
		ENDIF ELSE smin_z=2.0*smin_r	;if no inclination, smin_z set > smin_r	
		IF smin_r GT smin_z THEN smin = smin_z ELSE smin=smin_r		;set minimum s
		s=make(0,smin,n_s)
		IF keyword_set(verb) THEN print, s		
		;there is now array of points that go from the (pos[0],pos[1]) to the closest limiting surface and is divided into n_s parts.  

		;cycle through all the time points
		FOR j=0,ntau-1 DO BEGIN
			IF keyword_set(verb) THEN print, ' j = '+num2str(j,1)
			emiss_s=fltarr(n_s)	;initialize emissivity as a function of path length array
			emerr_s=fltarr(n_s)
			efit_i=ipt(efit_time,tau[j])
			IF efit_i NE -1 THEN BEGIN
				;load relevant LCFS boundary and prepare for format needed for LINE_INCLFS
				rbdry=efit_lcfs[0,where(efit_lcfs[0,*,efit_i] NE 0 AND efit_lcfs[1,*,efit_i] NE 0),efit_i]
				zbdry=efit_lcfs[1,where(efit_lcfs[0,*,efit_i] NE 0 AND efit_lcfs[1,*,efit_i] NE 0),efit_i]
				rbdry=rbdry[0:n(rbdry)-1]
				zbdry=zbdry[0:n(zbdry)-1]
				axis=[efit_axis[efit_i,0],efit_axis[efit_i,1]]
				
				;define an array of 1's (inside) and 0's (outside) where s is in/out LCFS
				s_inlcfs=intarr(n_s)
				FOR k=0,n_s-1 DO s_inlcfs[k]=line_inlcfs(rbdry,zbdry,axis,$
					[line_r(s[k],pos[*,i]), line_z(s[k],pos[*,i])])
				IF keyword_set(sol) THEN BEGIN
					goodpts=where(line_r(s,pos[*,i]) GT 0.44 AND line_r(s,pos[*,i]) LT 0.92)
				ENDIF ELSE goodpts=where(s_inlcfs NE 0) 
				IF keyword_set(verb) THEN print, goodpts
				
				;if any points are in the LCFS fill emiss_s and integrate
				IF goodpts[0] NE -1 THEN BEGIN
					;print debug info
					IF keyword_set(verb) THEN BEGIN
						print, 'Good R'
						print, line_r(s[goodpts], pos[*,i])
						print, 'Good Z'
						print, line_z(s[goodpts], pos[*,i])
						print, 'Time Point'
						print, efit_time[efit_i]
					ENDIF
				
					;determine the set of (R,Z) points for the view that are inside LCFS
					good_r=line_r(s[goodpts], pos[*,i])
					good_z=line_z(s[goodpts], pos[*,i])
					ngood=n(good_r)+1
                                        
                                	;if specified in RMAJ, convert to r/a
					IF NOT keyword_set(rho) THEN BEGIN
						irmid=reform(efit_rmid[efit_i,*])
						irho=(irmid-irmid[0])/(last(irmid)-irmid[0])
						rhovec=interpol(irho,irmid,r)
                                        ENDIF ELSE rhovec=r

					;find their equivilent midplane radii
					IF NOT keyword_set(rmap) THEN rhomid=efit_rz2rmid(good_r,good_z,efit_time[efit_i],shot=shot,/rho) ELSE rhomid=rmap
					IF npos EQ 1 AND ntau EQ 1 THEN rmap=rhomid
					IF keyword_set(verb) THEN BEGIN
						print, 'Good RMID'
						print, rhomid
					ENDIF

					;if any, set values of r_on_mid that are < min(r) to min(r)
					tmp=where(rhomid LT min(rhovec))
					IF tmp[0] NE -1 THEN rhomid[tmp]=min(rhovec)
		
					;use IDL function INTERPOLATE to finde emiss_s values
					;special case is isolated of emiss is a 1D array
					rho_interp=interp_vec_reform(rhovec, rhomid)
					tau_interp=interp_vec_reform(t,fltarr(ngood)+tau[j])
					IF n(t) EQ 0 THEN t_interp=fltarr(ngood)	;if only 1 time slice, default to it
					IF keyword_set(verb) THEN BEGIN
						print, 'R_INTERP'
						print, rho_interp
						print, 'T_INTERP'
						print, tau_interp
					ENDIF
					IF n(t) NE 0 THEN BEGIN
						emiss_good=interpolate(emiss,rho_interp,tau_interp, missing=1.0e-4) 
						emerr_good=interpolate(emerr,rho_interp,tau_interp, missing=1.0e-4)
                                        ENDIF  ELSE BEGIN
						emiss_good=interpolate(emiss,rho_interp,missing=1.0e-4)
						emerr_good=interpolate(emerr,rho_interp,missing=1.0e-4)
					ENDELSE
					emiss_s[goodpts]=emiss_good
					emerr_s[goodpts]=emerr_good
					IF keyword_set(verb) THEN BEGIN
						print, 'Good EMISS_S'
						print, emiss_good
					ENDIF

					;the scaling of l is arbitrary with no physical signifigance except that
					;l=0 is pt1 and l=1 is pt2 if the l values are multiplied by the distance
					;between 1 and 2 then the qunatitative units are applied to the scaling.
					lscale=sqrt((pos[0,i]^2-pos[2,i]^2)*(1.0+(tan(pos[3,i]))^2))
					br[i,j]=int_tabulated(s, emiss_s*lscale)
					FOR m=0,n_s-2 DO brerr[i,j]+=(0.5*(s[m+1]-s[m]))^2*(emerr_s[m]*lscale)^2+(0.5*(s[m+1]-s[m]))^2*(emerr_s[m+1]*lscale)^2
					brerr[i,j]=sqrt(brerr[i,j])
			
				ENDIF ELSE br[i,j]=0.0	;chord outside of LCFS
			ENDIF ELSE br[i,j]=-1	;value out of efit time range	
		ENDFOR
		
	ENDFOR	
	output=br
	IF keyword_set(debug) THEN stop
	RETURN, output
END

FUNCTION genpos_line_brtau,pos,emiss,r,t,shot,tau,emerr=emerr,n_s=n_s,verb=verb,debug=debug,brerr=brerr,rmap=rmap,efit=efit,rho=rho

	IF NOT keyword_set(n_s) THEN n_s=75	;		set number of segments to define line of sight
	IF NOT keyword_set(rzbnd) THEN rzbnd=[0.44,1.0,0.5]	;set path's [r_min,r_max,abs(z_max)] limit
	IF NOT keyword_set(emerr) THEN emerr=emiss*0.0
	
	npos=n(pos[0,*])+1	;determine the number of views
	ntau=n(tau)+1		;determine the number of time points
	br=fltarr(npos, ntau)	;initialize output brightness
	brerr=fltarr(npos,ntau)	;initialize output uncertainty [filled if emerr specified]

	;load EFIT data or decomile structure
	IF NOT keyword_set(efit) THEN BEGIN
		efit_time=line_gettimes(shot)
		efit_lcfs=line_getlcfs(shot)
		efit_axis=line_getaxis(shot)
		efit_rmid=line_getrmid(shot)
		efit={time:efit_time,lcfs:efit_lcfs,axis:efit_axis,rmid:efit_rmid}
	ENDIF ELSE BEGIN
		efit_rmid=efit.rmid
		efit_time=efit.time
		efit_lcfs=efit.lcfs
		efit_axis=efit.axis
	ENDELSE
	
	;debug information
	IF keyword_set(verb) THEN print, ' npos = '+num2str(npos,1)+'  ntau = '+num2str(ntau,1)
	IF keyword_set(debug) THEN stop

	;start filling the brightness array cycling through each view
	FOR i=0,npos-1 DO BEGIN
		IF keyword_set(verb) THEN print, ' i = '+num2str(i,1)

		;determine where the line of sight terminates in s due to radial boundries
		IF pos[2,i] LT rzbnd[0] THEN BEGIN	;if r_tang < r_min
			svec=line_s(pos[*,i],r=rzbnd[0])
			smin_r=svec[1]	 ;choose the smallest s where r=rmin
		ENDIF ELSE BEGIN			;if r_tang > r_min
			svec=line_s(pos[*,i],r=rzbnd[1])
			smin_r=svec[0]	 ;choose the largest (and non-negative) s where r=rmax
 		ENDELSE

		;determine where the line of sight terminates in s due to z boundries
		IF pos[3,i] NE 0 THEN BEGIN
			IF line_s(pos[*,i],z=rzbnd[2]) LT 0 THEN $
				smin_z=line_s(pos[*,i],z=-1.0*rzbnd[2]) ELSE $
				smin_z=line_s(pos[*,i],z=rzbnd[2])
		ENDIF ELSE smin_z=2.0*smin_r	;if no inclination, smin_z set > smin_r	
		IF smin_r GT smin_z THEN smin = smin_z ELSE smin=smin_r		;set minimum s
		s=make(0,smin,n_s)
		IF keyword_set(verb) THEN print, s		
		;there is now array of points that go from the (pos[0],pos[1]) to the closest limiting surface and is divided into n_s parts.  
	
		rpts=fltarr(n_s*ntau)
		zpts=fltarr(n_s*ntau)
		rzgood=intarr(n_s*ntau)
		FOR j=0,ntau-1 DO BEGIN
			efit_i=ipt(efit_time,tau[j])
			IF efit_i NE -1 THEN BEGIN
				;load relevant LCFS boundary and prepare for format needed for LINE_INCLFS
				rbdry=efit_lcfs[0,where(efit_lcfs[0,*,efit_i] NE 0 AND efit_lcfs[1,*,efit_i] NE 0),efit_i]
				zbdry=efit_lcfs[1,where(efit_lcfs[0,*,efit_i] NE 0 AND efit_lcfs[1,*,efit_i] NE 0),efit_i]
				rbdry=rbdry[0:n(rbdry)-1]
				zbdry=zbdry[0:n(zbdry)-1]
				axis=[efit_axis[efit_i,0],efit_axis[efit_i,1]]
				
				;define an array of 1's (inside) and 0's (outside) where s is in/out LCFS
				s_inlcfs=intarr(n_s)
				FOR k=0,n_s-1 DO s_inlcfs[k]=line_inlcfs(rbdry,zbdry,axis,[line_r(s[k],pos[*,i]), line_z(s[k],pos[*,i])])
				goodpts=where(s_inlcfs NE 0) 
				IF goodpts[0] NE -1 THEN BEGIN
					rpts[n_s*j+goodpts]=line_r(s[goodpts], pos[*,i])
					zpts[n_s*j+goodpts]=line_z(s[goodpts], pos[*,i])
					rzgood[n_s*j+goodpts]=1
				ENDIF
                	ENDIF
		ENDFOR	
		IF NOT keyword_set(rmap) THEN rmap=efit_rz2rmid(rpts,zpts,tau,shot=shot,/rho)
 
	
		;cycle through all the time points
		FOR j=0,ntau-1 DO BEGIN
			IF keyword_set(verb) THEN print, ' j = '+num2str(j,1)
			emiss_s=fltarr(n_s)	;initialize emissivity as a function of path length array
			emerr_s=fltarr(n_s)
			efit_i=ipt(efit_time,tau[j])
			IF efit_i NE -1 THEN BEGIN

				goodpts=where(rzgood[j*n_s:(j+1)*n_s-1] EQ 1)
				ngood=total(rzgood[j*n_s:(j+1)*n_s-1])
				;if any points are in the LCFS fill emiss_s and integrate
				IF goodpts[0] NE -1 THEN BEGIN
					;print debug info
					IF keyword_set(verb) THEN BEGIN
						print, 'Good R'
						print, line_r(s[goodpts], pos[*,i])
						print, 'Good Z'
						print, line_z(s[goodpts], pos[*,i])
						print, 'Time Point'
						print, efit_time[efit_i]
					ENDIF
				
	                             	;if specified in RMAJ, convert to r/a
					IF NOT keyword_set(rho) THEN BEGIN
						irmid=reform(efit_rmid[efit_i,*])
						irho=(irmid-irmid[0])/(last(irmid)-irmid[0])
						rhovec=interpol(irho,irmid,r)
                                        ENDIF ELSE rhovec=r
					rhomid=rmap[n_s*j+goodpts,j]
					
					IF keyword_set(verb) THEN BEGIN
						print, 'Good RMID'
						print, rhomid
					ENDIF

					;if any, set values of r_on_mid that are < min(r) to min(r)
					tmp=where(rhomid LT min(rhovec))
					IF tmp[0] NE -1 THEN rhomid[tmp]=min(rhovec)
		
					;use IDL function INTERPOLATE to finde emiss_s values
					;special case is isolated of emiss is a 1D array
					rho_interp=interp_vec_reform(rhovec, rhomid)
					tau_interp=interp_vec_reform(t,fltarr(ngood)+tau[j])
					IF n(t) EQ 0 THEN t_interp=fltarr(ngood)	;if only 1 time slice, default to it
					IF keyword_set(verb) THEN BEGIN
						print, 'R_INTERP'
						print, r_interp
						print, 'T_INTERP'
						print, t_interp
					ENDIF
					IF n(t) NE 0 THEN BEGIN
						emiss_good=interpolate(emiss,rho_interp,tau_interp, missing=1.0e-4) 
						emerr_good=interpolate(emerr,rho_interp,tau_interp, missing=1.0e-4)
                                        ENDIF  ELSE BEGIN
						emiss_good=interpolate(emiss,rho_interp,missing=1.0e-4)
						emerr_good=interpolate(emerr,rho_interp,missing=1.0e-4)
					ENDELSE
					emiss_s[goodpts]=emiss_good
					emerr_s[goodpts]=emerr_good
					IF keyword_set(verb) THEN BEGIN
						print, 'Good EMISS_S'
						print, emiss_good
					ENDIF

					;the scaling of l is arbitrary with no physical signifigance except that
					;l=0 is pt1 and l=1 is pt2 if the l values are multiplied by the distance
					;between 1 and 2 then the qunatitative units are applied to the scaling.
					lscale=sqrt((pos[0,i]^2-pos[2,i]^2)*(1.0+(tan(pos[3,i]))^2))
					br[i,j]=int_tabulated(s, emiss_s*lscale)
					FOR m=0,n_s-2 DO brerr[i,j]+=(0.5*(s[m+1]-s[m]))^2*(emerr_s[m]*lscale)^2+(0.5*(s[m+1]-s[m]))^2*(emerr_s[m+1]*lscale)^2
					brerr[i,j]=sqrt(brerr[i,j])
			
				ENDIF ELSE br[i,j]=0.0	;chord outside of LCFS
			ENDIF ELSE br[i,j]=-1	;value out of efit time range	
		ENDFOR
		
	ENDFOR	
	output=br
	IF keyword_set(debug) THEN stop
	RETURN, output
END
