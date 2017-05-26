;+
;NAME:
;	GENPOS_SPHERICAL_INFO
;
;PURPOSE:
;	This function creates an detector information structure that can be used with other
;	functions in GENIE to create POS vectors and GPV data.  This function can only handle
;	a spherical bragg reflector type of information file.  Data is entered in an ASCII format
;	
;CALLING SEQUENCE:
;	result=GENPOS_SPERICAL_INFO(path)
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
;		*.m		STRUC mirror substructure
;		*.*.vec		FLTARR [ro,0.0,zo] of aperture [meters]
;		*.*.rot		FLTARR [alpha,beta,gamma] of aperture normal [radians]
;		*.*.size	FLTARR [y,z] of aperture dimensions [meters]
;		*.*.bragg	STRUC of the bragg information
;		*.*.*.twod	FLOAT of the 2d spacing [Ang]
;		*.*.*.iref	FLOAT of the integrated reflectivity
;		*.*.*.rwid	FLOAT of the rocking curve width [mRad]
;		*.*.rad		FLOAT of the mirror radius [m]
;		*.det		STRUC detector substructure
;		*.*.x0
;		*.*.x1		x3 FLTARR that locate the detector plane in aperture coordinates [meters]
;		*.*.x2
;		*.*.xi		FLTARR of xi locations of pixels on the detector [meters]
;		*.*.zeta	FLTARR of zeta locations of pixels on the detector [meters]
;		*.*.size	FLTARR [xi,zeta] of detector dimensions [meters]
;		*.*.n_xi	INT number of xi pts
;		*.*.n_zeta	INT number of zeta pts
;		*.type		STR of type (must be spherical)
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
;	Written by:	ML Reinke: 1/10
;
;-

FUNCTION genpos_spherical_info,path,debug=debug
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
	IF strmatch (type, '*spherical*') EQ 0 THEN RETURN,-1	
	readf,lun,line
	tmp=strsplit(line, '=',/extract)
	author=tmp[1]

	readf,lun,line
	WHILE eof(lun) NE 1 AND strmatch(line,'*DATA START*',/fold_case) EQ 0 DO readf, lun, line


	;MIRROR
	;---------------------------------------------------------------------------------------
	;create m_vector
	point_lun,-lun,data_start	
	readf,lun,line
	WHILE eof(lun) NE 1 AND strmatch(line,'m_ro *',/fold_case) EQ 0 DO readf, lun, line
	ro=float(str_dbl_ext(line,'=',';'))
	point_lun,lun,data_start
	readf,lun,line
	WHILE eof(lun) NE 1 AND strmatch(line,'m_zo *',/fold_case) EQ 0 DO readf, lun, line
	zo=float(str_dbl_ext(line,'=',';'))
	m_vec=[ro,0.0,zo]

	;create m_rot
	point_lun,lun,data_start
	readf,lun,line
	WHILE eof(lun) NE 1 AND strmatch(line,'m_alpha *',/fold_case) EQ 0 DO readf, lun, line
	alpha=float(str_dbl_ext(line,'=',';'))
	point_lun,lun,data_start
	readf,lun,line
	WHILE eof(lun) NE 1 AND strmatch(line,'m_beta *',/fold_case) EQ 0 DO readf, lun, line
	beta=float(str_dbl_ext(line,'=',';'))
	point_lun,lun,data_start
	readf,lun,line
	WHILE eof(lun) NE 1 AND strmatch(line,'*m_gamma *',/fold_case) EQ 0 DO readf, lun, line
	gamma=float(str_dbl_ext(line,'=',';'))
	m_rot=[alpha,beta,gamma]

	;create m_size
	point_lun,lun,data_start
	readf,lun,line
	WHILE eof(lun) NE 1 AND strmatch(line,'*m_y *',/fold_case) EQ 0 DO readf, lun, line
	y=float(str_dbl_ext(line,'=',';'))
	point_lun,lun,data_start
	readf,lun,line
	WHILE eof(lun) NE 1 AND strmatch(line,'*m_z *',/fold_case) EQ 0 DO readf, lun, line
	z=float(str_dbl_ext(line,'=',';'))
	m_size=[y,z]*1.0e-2

	;read Bragg properties
	point_lun,lun,data_start
	readf,lun,line
	WHILE eof(lun) NE 1 AND strmatch(line,'*m_twod *',/fold_case) EQ 0 DO readf, lun, line
	m_twod=float(str_dbl_ext(line,'=',';'))
	point_lun,lun,data_start
	readf,lun,line
	WHILE eof(lun) NE 1 AND strmatch(line,'*m_iref *',/fold_case) EQ 0 DO readf, lun, line
	m_iref=float(str_dbl_ext(line,'=',';'))
	point_lun,lun,data_start
	readf,lun,line
	WHILE eof(lun) NE 1 AND strmatch(line,'*m_rwid *',/fold_case) EQ 0 DO readf, lun, line
	m_rwid=float(str_dbl_ext(line,'=',';'))
	bragg=create_struct('twod',m_twod,'iref',m_iref,'rwid',m_rwid)

	;read radius of curvature
	point_lun,lun,data_start
	readf,lun,line
	WHILE eof(lun) NE 1 AND strmatch(line,'*m_rad *',/fold_case) EQ 0 DO readf, lun, line
	m_rad=float(str_dbl_ext(line,'=',';'))

	;form mirror structure
	m_struc=create_struct('vec',m_vec,'rot',m_rot,'size',m_size,'bragg',bragg,'rad',m_rad)

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
	det_struc=create_struct('x0',x0,'x1',x1,'x2',x2,'xi',xi, 'zeta',zeta,'size',det_size,'n_xi',int(n_xi),'n_zeta',int(n_zeta))	

	close,lun
	free_lun,lun
	
	output={name:name, m:m_struc, det:det_struc, type:type, author:author}
	IF keyword_set(debug) THEN stop
	RETURN, output
END


;+
;NAME:
;	GENPOS_SPHERICAL2QUADCURVE
;
;PURPOSE:
;	This function generate the coefficients for the quadradic curve on the detector plane
;	in pixel space (xi,zeta) given an info file and a wavelength.
;
;CALLING SEQUENCE:
;	result=GENPOS_SPHERICAL2QUADCURVE(info,lam_o)
;
;INPUTS:
;	info	STRUC	output of GENPOS_SPHERICAL_INFO
;	lam_o	FLOAT	wavelength [Ang]
;
;OUTPUTS:
;	result:	FLTARR [6] of the quaddradic curve coefficients [a,b,c,d,f,g]
;		such that xi and zeta to be specified in pixel #'s not physical length
;
;PROCEDURE:
;	The coffiencts specify the curve in the equation:
;		a*xi^2+2*b*xi*zeta+c*zeta^2+2*d*xi+2*f*zeta+g=0
;	and the derivation of the dependance of the coefficients on the INFO file parameters can
;	be found in "General Equations for X-Ray Crystal Imaging Spectroscopy in Tokamaks" Section 3.
;
;	The outputs can be used in functions like ELLIPSE_XPT which given zeta's calculate xi's of the
;	(usually) ellipse defined by the detector/mirror alignment in the INFO file.
;
;MODFICATION HISTORY:
;	Written by: 	M.L. Reinke 1/10
;	8/11		M.L. Reinke - added the phi input to move the origin in x-y plane
;-

FUNCTION genpos_spherical2quadcurve,info,lam_o,phi=phi

	x0=info.det.x0
	x1=info.det.x1
	x2=info.det.x2
	dxi=info.det.size[0]
	dzeta=info.det.size[1]

	th_bragg=asin(lam_o/info.m.bragg.twod)
	cbsqr=1.0/(tan(th_bragg))^2

	;rotate coordinate system to use different parts of the crystal
	IF keyword_set(phi) THEN BEGIN
		R=info.m.rad

		;rotate x0
		x=x0[0]
		y=x0[1]
		x0[0]=(x-R*(1-cos(phi)))*cos(phi)-(y-R*sin(phi))*sin(phi)
		x0[1]=(x-R*(1-cos(phi)))*sin(phi)+(y-R*sin(phi))*cos(phi)

		;rotate x1
		x=x1[0]
		y=x1[1]
		x1[0]=(x-R*(1-cos(phi)))*cos(phi)-(y-R*sin(phi))*sin(phi)
		x1[1]=(x-R*(1-cos(phi)))*sin(phi)+(y-R*sin(phi))*cos(phi)

		;rotate x2
		x=x2[0]
		y=x2[1]
		x2[0]=(x-R*(1-cos(phi)))*cos(phi)-(y-R*sin(phi))*sin(phi)
		x2[1]=(x-R*(1-cos(phi)))*sin(phi)+(y-R*sin(phi))*cos(phi)
	ENDIF
 	
	xo=x0[0]
	yo=x0[1]
	zo=x0[2]

	xi=(x2-x0)
	xi_hat=xi/sqrt(total(xi*xi))
	xi_x=xi_hat[0]
	xi_y=xi_hat[1]
	xi_z=xi_hat[2]

	zeta=(x1-x0)
	zeta_hat=zeta/sqrt(total(zeta*zeta))
	zeta_x=zeta_hat[0]
	zeta_y=zeta_hat[1]
	zeta_z=zeta_hat[2]
	
	a=(1.0-xi_x^2*(1.0+cbsqr))*dxi^2
	b=-2.0*xi_x*zeta_x*(1.0+cbsqr)*dxi*dzeta
	c=(1.0-zeta_x^2*(1.0+cbsqr))*dzeta^2
	d=(2.0*(xo*xi_x+yo*xi_y+zo*xi_z-xo*xi_x*(1.0+cbsqr)))*dxi
	f=(2.0*(xo*zeta_x+yo*zeta_y+zo*zeta_z-xo*zeta_x*(1.0+cbsqr)))*dzeta
	g=yo^2+zo^2-xo^2*cbsqr
	
	out=[a,b,c,d,f,g]
	RETURN,out
END

;+
;NAME:
;	GENPOS_SPHERICAL_REFL
;
;-

FUNCTION genpos_spherical_refl,x1,xs,n
	x2=2.0*xs-x1-2.0*total((xs-x1)*n)*n
	RETURN,x2
END

;+
;NAME:
;	GENPOS_SPHERE
;
;PURPOSE:
;	This function calculates the POS vector and etendue for a set of detector points
;	given detector position and mirror position information usually specified by an INFO file.
;	This traces rays from the detector through the center of the mirror.
;
;CALLING SEQUENCE:
;	result=GENPOS_SPHERE(m_vec,m_rot,m_rad,x0,x1,x2,det_pos)
;
;INPUTS:
;	m_vec: 		FLTARR 	[ro, 0.0, zo] where ro and zo is the mirror position [meter]
;	m_rot:		FLTARR 	[alpha,beta,gamma] set of orientation angles of the mirror [in radians]
;	m_rad:		FLOAT 	of the mirror radius of curvature [meter]
;	x0:
;	x1:		x3 FLTARR locations of the detector's fiducial points in the mirror coordinate system
;	x2:
;	det_pos:	FLTARR 	[2 x m] where m is the number of points in the deteector plane
;				for which GENPOS_SPHERE will cacluate POS vectors
;
;OPTIONAL INPUTS:
;	a_det:		FLOAT area of detector DEFAULT = 1.0
;	a_m:		FLOAT area of mirror DEFAULT = 1.0
;	
;KEYWORD PARAMETERS:
;	debug:		/debug stops before the RETURN statements
;
;OUTPUTS:
;	result:		FLTARR [ro,zo,rt,psi] x m position vectors
;
;OPTIONAL OUTPUTS:
;	etendue:	FLTARR [m] of etendue values
;
;MODIFICATION HISTORY:
;	Written by:	M.L. Reinke 1/10
;	8/5/10:		M.L. Reinke - Changed variable labeling that was screwing up the etendue calculation
;-

FUNCTION genpos_sphere,m_vec,m_rot,m_rad,x0,x1,x2,det_pos,etendue=etendue,thb=thb,a_det=a_det,a_m=a_m,debug=debug
	IF NOT keyword_set(a_det) THEN a_det=1.0
	IF NOT keyword_set(a_m) THEN a_m=1.0		

	size=size(det_pos)
	num=size[2]
	rot=genpos_rot_matrix(m_rot[0],m_rot[1],m_rot[2])
	xyz_vec=genpos_det2xyz(x0,x1,x2,det_pos)

	pos=fltarr(4,num)
	etendue=fltarr(num)
	thb=fltarr(num)
	FOR i=0L,num-1 DO BEGIN
		xa=xyz_vec[*,i]
		xs=[0.0,0.0,0.0]	;if gridding mirror this is non-zero
		n_vec=[m_rad,0,0]-xs
		n_mag=sqrt(total(n_vec*n_vec))
		n=n_vec/n_mag
		xb=genpos_spherical_refl(xa,xs,n)
		
		xyz1=xs		;position on mirror
		xyz2=xb		;reflected position
		rt_psi=genpos_dblxyz2tok(m_vec,xyz1,xyz2,rot)
		pos[0,i]=m_vec[0]
		pos[1,i]=m_vec[2]
		pos[2,i]=rt_psi[0]
		pos[3,i]=rt_psi[1]
		;stop

		;adapted to spherical mirror from GENPOS_DET2U
              	norm=crossp((x2-x0),(x1-x0))
                n_hat_det=norm/sqrt(total(norm*norm))
		n_hat_m=n
		xyz=xa-xs
		xyz_mag=sqrt(total(xyz*xyz))
		thb[i]=asin(total(xyz_mag*n_hat_m))	
		etendue[i]=-1.0*total(n_hat_det*xyz)*a_det/xyz_mag*total(n_hat_m*xyz)*a_m/xyz_mag/xyz_mag^2
	ENDFOR

	IF keyword_set(debug) THEN stop
	output=pos
	RETURN,output		

END

;+
;NAME:
;	GENPOS_SPHERE_UPOS
;
;PURPOSE:
;	This function calculates POS vectors for a finite sized mirror.
;
;CALLING SEQUENCE:
;	result=GENPOS_SPHERE_UPOS(m_vec,m_rot,m_rad,x0,x1,x2,m_grid,det_grid)
;
;INPUTS
;	m_vec: 		FLTARR 	[ro, 0.0, zo] where ro and zo is the mirror position [meters]
;	m_rot:		FLTARR 	[alpha,beta,gamma] set of orientation angles of the mirror [in radians]
;	m_rad:		FLOAT 	of the mirror radius of curvature [meter]
;	x0:
;	x1:		x3 FLTARR locations of the detector's fiducial points in the mirror coordinate system
;	x2:
;	m_grid:		FLTARR 	[3,num_m] of the <x,y,z> points on the mirror
;	det_grid:	FLTARR 	[2,num_det] of the <xi,zeta> points on the detector plane
;
;KEYWORD PARAMETERS:
;	debug		/debug will stop the code before the RETURN 
;	
;OUTPUTS:
;	result:		FLTARR	[4,num_m,num_det] of the POS vector for each det/mirror grid point
;
;OPTIONAL OUTPUTS:
;	du:		FLTARR	[num_m,num_det]	of the etendue for each det/mirror grid point
;	thb:		FLTARR	[num_m,num_det] of the bragg angle for each det/mirror grid point [radians]
;	bigXYZ:		FLTARR	[3,num_m,num_det,3] of the XYZ (tokamak) coordiates of:
;				[*,*,*,0] - the detector point
;				[*,*,*,1] - the mirror point
;				[*,*,*,2] - the reflected point
;	
;MODIFICATION HISTORY:
;	Written by:	ML Reinke 1/10
;	2/27/10		ML Reinke - changed the output array formatting to seperate out each mirror/det combination
;	8/20/10		ML Reinke - updated so that the du and thb optional outputs are calculated
;
;- 

FUNCTION genpos_sphere_upos,m_vec,m_rot,m_rad,x0,x1,x2,m_grid,det_grid,a_det=a_det,a_m=a_m,du=du,thb=thb,debug=debug,bigXYZ=bigXYZ,xyz_vec=xyz_vec
	IF NOT keyword_set(a_det) THEN a_det=1.0
	IF NOT keyword_set(a_m) THEN a_m=1.0		


	;determine size of arrays
	size=size(det_grid)
	IF size[0] EQ 1 THEN num_det=size[1] ELSE num_det=size[2]
	IF n(a_det)+1 NE num_det THEN dA_det=fltarr(num_det)+a_det[0] ELSE dA_det=a_det
	size=size(m_grid)
	IF size[0] EQ 1 THEN num_m=1 ELSE num_m=size[2]
	IF n(a_m)+1 NE num_m THEN dA_m=fltarr(num_m)+a_m[0] ELSE dA_m=a_m

	rot=genpos_rot_matrix(m_rot[0],m_rot[1],m_rot[2])	;generate rotation matrix
	xyz_vec=genpos_det2xyz(x0,x1,x2,det_grid)		;tranlate det coords to xyz

	upos=fltarr(4,num_m,num_det)
	bigXYZ=fltarr(3,num_m,num_det,3)
	du=fltarr(num_m,num_det)
	thb=fltarr(num_m,num_det)

	cntr=0L
	FOR i=0L,num_m-1 DO BEGIN
		FOR j=0L,num_det-1 DO BEGIN
			xa=xyz_vec[*,j]
			xs=m_grid[*,i]
			n_vec=[m_rad,0,0]-xs
			n_mag=sqrt(total(n_vec*n_vec))
			n=n_vec/n_mag
			xb=genpos_spherical_refl(xa,xs,n)
			xyz1=xs									;position on mirror
			xyz2=xb									;reflected position
			rt_psi=genpos_dblxyz2tok(m_vec,xyz1,xyz2,rot,cyl_vec=cyl_vec)
			upos[0,i,j]=cyl_vec[0,0]
			upos[1,i,j]=cyl_vec[2,0]
			upos[2,i,j]=rt_psi[0]
			upos[3,i,j]=rt_psi[1]

			;generate XYZ coordinates of points
			bigXYZ[*,i,j,0]=[0.0, m_vec[0],m_vec[2]]+rot##x1	;det
			bigXYZ[*,i,j,1]=[0.0, m_vec[0],m_vec[2]]+rot##xs	;mirror
			bigXYZ[*,i,j,1]=[0.0, m_vec[0],m_vec[2]]+rot##x2	;reflection

			;adapted to spherical mirror from GENPOS_DET2U
		    	norm=crossp((x2-x0),(x1-x0))
                	n_hat_det=norm/sqrt(total(norm*norm))
			n_hat_m=n
			xyz=xa-xs						;line connecting mirror to detector
			xyz_mag=sqrt(total(xyz*xyz))
			thb[i,j]=asin(total(xyz/xyz_mag*n_hat_m))			;calculate the bragg angle for this mirror/detector point
			du[i,j]=-1.0*total(n_hat_det*xyz)*dA_det[j]/xyz_mag*total(n_hat_m*xyz)*dA_m[i]/xyz_mag/xyz_mag^2
		ENDFOR
	ENDFOR
		
	IF keyword_set(debug) THEN stop
	output=upos
	RETURN,output		
		

END

;+
;NAME:
;	GENPOS_GRID_SPHERE
;
;PURPOSE:
;	This function is used to grid (centers) a spherical mirror in the xyz coordinate system.
;
;CALLING SEQUENCE:
;	result=GENPOS_GRID_SPHERE(info)
;
;INPUTS:
;	info	STRUC 	of info file
;
;OPTIONAL INPUTS:
;	ny	INT	of the number of grid points in the y-hat direction DEFAULT: 5	
;	nz	INT	of the number of grid points in the z-hat direction DEFAULT: 5
;
;KEYWORD PARAMETERS:
;	debug	/debug will stop the program before the RETURN statement
;
;OUTPUTS:
;	result:	FLTARR [3,nz*ny] of the <x,y,z> points on the mirror [m]
;
;OPTIONAL OUTPUS:
;	dA	FLTARR	[nz*ny] of the area [m^2] of each grided element
;
;PROCEDURE:
;	This assumes that the center of the mirror is at the origin of the xyz coordinate system with the x-hat
;	being along the mirror normal.
;
;MODIFICATION HISTORY:
;	Written by: 	M. L. Reinke 1/10
;	8/20/10		M.L. Reinke - modified the rectangular gridder to use centers and calculate the area
;-

FUNCTION genpos_grid_sphere,info,ny=ny,nz=nz,circle=circle,dA=dA,debug=debug
	IF NOT keyword_set(nz) THEN nz=5
	IF NOT keyword_set(ny) THEN ny=5
	mdy=info.m.size[0]	;dy of mirror [m]
	mdz=info.m.size[1]	;dz of mirror [m]
	IF ny EQ 1 AND nz EQ 1 THEN BEGIN
		dA=mdy*mdz
		RETURN,fltarr(3)	;return immediately if the center is chosen
	ENDIF
	R=info.m.rad
	nz=float(nz)
	ny=float(ny)
	IF NOT keyword_set(circle) THEN BEGIN
		m_grid=fltarr(3,nz*float(ny))
		dA=fltarR(nz*float(ny))
		psi=make(-atan(mdz/2.0/R),atan(mdz/2.0/R),nz+1)	;psi boundaries
		dpsi=psi[1]-psi[0]
		phi=make(-atan(mdy/2.0/R),atan(mdy/2.0/R),ny+1)	;psi boundaries
		dphi=phi[1]-phi[0]
		FOR i=0,ny-1 DO BEGIN
			iphi=phi[i]+(phi[i+1]-phi[i])/2.0
			FOR j=0,nz-1 DO BEGIN
				jpsi=psi[j]+(psi[j+1]-psi[j])/2.0
				dA[i*nz+j]=R^2*cos(jpsi)*dphi*dpsi
				m_grid[*,i*nz+j]=[R*(1.0-cos(iphi))*cos(jpsi),R*cos(jpsi)*sin(iphi),R*sin(jpsi)]
			ENDFOR
		ENDFOR
	ENDIF ELSE BEGIN
;		m_grid=fltarr(3,nz*ny+1)
;		thvec=make(0,2.0*!pi,ny+1)
;		thvec=thvec[1:*]
;		rvec=make(0,(mdz/2.0)^2,nz)
;		rvec=sqrt(rvec)
;		cntr=1.0	;sets m_grid[*,0]=<0,0,0>
;		FOR i=0,ny-1 DO BEGIN	;variation in theta
;			FOR j=0,nz-1 DO BEGIN	;variation in radius
;				yvec=rvec[j]*cos(thvec[i])
;				zvec=rvec[j]*sin(thvec[i])
;				xvec=R-sqrt(R^2-yvec^2-zvec^2)				
;				m_grid[0,cntr]=xvec
;				m_grid[1,cntr]=yvec
;				m_grid[2,cntr]=zvec
;				cntr+=1
;			ENDFOR
;		ENDFOR
		m_grid=fltarr(3,nz*float(ny))
		yvec=make(-1.0*mdy/2.0,mdy/2.0,ny) 
		zvec=make(-1.0*mdz/2.0,mdz/2.0,nz)
		IF ny EQ 1 THEN yvec=[0.0]
		IF nz EQ 1 THEN zvec=[0.0]
		cntr=0
		FOR i=0L,ny-1 DO BEGIN
			xvec=R-sqrt(R^2-yvec[i]^2-zvec^2)
			m_grid[0,cntr*nz:(cntr+1)*nz-1]=xvec
			m_grid[1,cntr*nz:(cntr+1)*nz-1]=fltarr(nz)+yvec[i]
			m_grid[2,cntr*nz:(cntr+1)*nz-1]=zvec
			cntr+=1
		ENDFOR
		tmp=where(sqrt(m_grid[1,*]^2+m_grid[2,*]^2) LE mdy/2.0)
		m_grid=m_grid[*,tmp]
	ENDELSE

	IF keyword_set(debug) THEN stop
	RETURN,m_grid
END

;+
;NAME:
;	GENPOS_SPHERICAL2POS
;	
;PURPOSE:
;	This function takes a spherical information structure and returns the corresponding POS vectors.
;
;CALLING SEQUENCE:
;	result=GENPOS_SPERICAL2POS(info)
;
;INPUTS:
;	info:		STRUC of detector information (see GENPOS_SPHERICAL_INFO)
;
;KEYWORD PARAMETERS:
;	debug:		/debug stops before the RETURN statement
;
;OUTPUTS:
;	result:		FLTARR 4 x m where m is the number of detector positions in info.det.xi and info.det.zeta
;
;OPTIONAL OUTPUTS:
;	etendue:	FLTARR of length m of the etendue for each detector
;
;RESTRICTIONS:
;	This function calls GENPOS to generate the POS vectors as well as calculate the etendues.
;
;MODIFICATION HISTORY:
;	Written by:	ML Reinke 1/10
;
;-

FUNCTION genpos_spherical2pos,info,etendue=etendue,debug=debug,xi=xi

	IF strmatch(info.type ,'*spherical*') EQ 0 THEN RETURN, -1
	IF keyword_set(xi) THEN BEGIN
		xi_pts=info.det.xi[0:info.det.n_xi-1]
		ipt=ipt(xi_pts,xi*info.det.size[0])
		index=indgen(info.det.n_zeta)*info.det.n_xi+ipt
		det_pos=fltarr(2,info.det.n_zeta)
		det_pos[0,*]=info.det.xi[index]
		det_pos[1,*]=info.det.zeta[index]
	ENDIF ELSE BEGIN
		num=n(info.det.xi)+1
		det_pos=fltarr(2,num)
		det_pos[0,*]=info.det.xi
		det_pos[1,*]=info.det.zeta
	ENDELSE
	a_det=info.det.size[0]*info.det.size[1]
	a_m=info.m.size[0]*info.m.size[1]
	pos=genpos_sphere(info.m.vec,info.m.rot,info.m.rad,info.det.x0,info.det.x1,info.det.x2,det_pos,etendue=etendue,a_det=a_det,a_m=a_m)
	IF keyword_set(debug) THEN stop
	RETURN,pos
END



;+
;NAME:
;	GENPOS_SPHERICAL2UPOS
;
;-

FUNCTION genpos_spherical2upos,info,debug=debug,corners=corners,ny=ny,nz=nz,bigXYZ=bigXYZ,du=du,xyz_m=xyz_m,xyz_det=xyz_det,thb=thb
	IF strmatch(info.type ,'*spherical*') EQ 0 THEN RETURN, -1

	xyz_m=genpos_grid_sphere(info,ny=ny,nz=nz,dA=dA)
	
	IF keyword_set(corners) THEN BEGIN
		n_xi=info.det.n_xi
		n_zeta=info.det.n_zeta
		det_pos=fltarr(2,4)
		det_pos[*,0]=[info.det.xi[0],info.det.zeta[0]]
		det_pos[*,1]=[info.det.xi[n_xi-1],info.det.zeta[n_xi-1]]
		det_pos[*,2]=[info.det.xi[n_xi*(n_zeta-1)],info.det.zeta[(n_zeta-1)*n_xi]]
		det_pos[*,3]=[info.det.xi[n_zeta*n_xi-1],info.det.zeta[n_zeta*n_xi-1]]
		
	ENDIF ELSE BEGIN
		num=n(info.det.xi)+1
		det_pos=fltarr(2,num)
		det_pos[0,*]=info.det.xi
		det_pos[1,*]=info.det.zeta
	ENDELSE
	a_det=info.det.size[0]*info.det.size[1]		;pixel area in [m^2]

	upos=genpos_sphere_upos(info.m.vec,info.m.rot,info.m.rad,info.det.x0,info.det.x1,info.det.x2,xyz_m,det_pos,a_det=a_det,a_m=dA,du=du,bigXYZ=bigXYZ,xyz_vec=xyz_det,thb=thb)
	IF keyword_set(debug) THEN stop
	RETURN,upos
END

;+
;NAME:
;	GENPOS_ALIGN_ROT
;
;-

FUNCTION genpos_align_rot,x6,a,b,c,Rabc=Rabc
	Ra=[[cos(a),sin(a),0],$
	    [-1.0*sin(a),cos(a),0],$
	    [0,0,1]]
	Rb=[[1,0,0],$
	    [0,cos(b),sin(b)],$
	    [0,-sin(b),cos(b)]]
	Rc=[[cos(c),0,-sin(c)],$
	    [0,1,0],$
	    [sin(c),0,cos(c)]]
	Rabc=(Ra##Rb##Rc)
	x3=reform(Rabc##x6)
	RETURN,x3
END

FUNCTION genpos_align2info,align,info,plotpts=plotpts
	pix=info.det.size[0]
	det=[info.det.n_xi,info.det.n_zeta]*pix
	genpos_align2xyz,align[0],align[1],align[2],align[3]*pix,align[4]*pix,align[5],align[6],align[7],x0,x1,x2,pix=pix,det=det,plotpts=plotpts
	newinfo=info
	newinfo.det.x0=x0
	newinfo.det.x1=x1
	newinfo.det.x2=x2
	RETURN,newinfo
END

;+
;NAME:
;	GENPOS_ALIGN2XYZ
;
;-

PRO genpos_align2xyz,l,h,q,x,z,a,b,c,x0,x1,x2,pix=pix,det=det,plotpts=plotpts
	;set optional input of pix size and detector size
	IF NOT keyword_set(pix) THEN pix=172.0e-6		;[m]
	IF NOT keyword_set(det) THEN det=[195.0,487.0]*pix	;[m]
	
	Rq=[[cos(q),-sin(q),0],$
	    [sin(q),cos(q),0],$
	    [0, 0, 1]]

	l3=[0,-l,h]
	d06=[x,0,-1.0*z]
	d16=[x,0,-1.0*z+det[1]]
	d26=[x-det[0],0,-1.0*z]
	
	x03=genpos_align_rot(d06,a,b,c)
	x13=genpos_align_rot(d16,a,b,c)
	x23=genpos_align_rot(d26,a,b,c)
	orig=reform(Rq##l3)
	plotpts=[[orig],[x03],[x13],[x23]]

	x0=reform(Rq##(l3+x03))
	x1=reform(Rq##(l3+x13))
	x2=reform(Rq##(l3+x23))
END

;+
;NAME:
;	GENPOS_XYZ2ALIGN
;
;-

PRO genpos_xyz2align,x0,x1,x2,l,h,q,x,z,a,b,c,pix=pix,det=det
	;set optional input of pix size and detector size
	IF NOT keyword_set(pix) THEN pix=172.0e-6	;[m]
	IF NOT keyword_set(det) THEN det=[195.0,487.0]*pix ;in [m]

	;find theta
	tq_num=(-1.0*x0[0]*det[1]+(x0[0]-x1[0])*z)*det[0]+(x0[0]-x2[0])*x*det[1]
	tq_dem=(-1.0*det[1]*x0[1]+(x0[1]-x1[1])*z)*det[0]+(x0[1]-x2[1])*x*det[1]
	q=atan(-1.0*tq_num/tq_dem)
	sq=sin(q)
	cq=cos(q)
	tq=tq_num/tq_dem
	
	;find l
	l=1.0/cq*(-1.0*x0[1]+(x0[1]-x1[1])*z/det[1]+(x0[1]-x2[1])*x/det[0])
	
	;find h
	h=x0[2]+(x1[2]-x0[2])*z/det[1]+(x0[2]-x2[2])*x/det[0]


	;find c
	tc=-1.0*(x0[2]-x2[2])/(x0[2]-x1[2])*det[1]/det[0]
	c=atan(tc)
	sc=sin(c)
	cc=cos(c)

	;find bigC
	bigC_num=det[1]*(L-sq*x2[0]+cq*x2[1])+x*(sq*(x2[0]-x1[0])+cq*(x1[1]-x2[1]))
	bigC_dem=det[0]*(det[1]-z)-det[0]*x
	bigC=-1.0*bigC_num/bigC_dem
	
	;find bigD
	bigD_num=det[0]*(l+cq*x1[1]-sq*x1[0])+x*(sq*(x1[0]-x2[0])+cq*(x2[1]-x1[1]))
	bigD_dem=bigC_dem
	bigD=bigD_num/bigD_dem

	;find a
	sa=bigD*sc-bigC*cc
	a=asin(sa)
	ca=cos(a)

	;find b
	sb=(sc*bigC+cc*bigD)/ca
	b=asin(sb)

END

;+
;NAME:
;	GENPOS_INFO2ALIGN
;
;-

FUNCTION genpos_info2align,info,x=x,z=z
	IF NOT keyword_set(x) THEN x=0.0
	IF NOT keyword_set(z) THEN z=0.0
	x0=info.det.x0
	x1=info.det.x1
	x2=info.det.x2
	pix=info.det.size[0]
	det=[info.det.n_xi,info.det.n_zeta]*pix
	genpos_xyz2align,x0,x1,x2,l,h,q,x,z,a,b,c,pix=pix,det=det
	align=[l,h,q,x,z,a,b,c]
	RETURN,align
END

;+
;NAME:
;	LINE_PLANE_INT
;
;-

FUNCTION line_plane_int,l,p
	;see wikipedia entry on line-plane intersection
	b=(l[*,0]-p[*,0])
	A=[[l[*,0]-l[*,1]],[p[*,1]-p[*,0]],[p[*,2]-p[*,0]]]
	tuv=la_invert(A)#b
	RETURN,tuv
END
