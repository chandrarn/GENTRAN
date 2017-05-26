;+
FUNCTION is_xtomo_emiss,shot
	mdsopen,'xtomo',shot
	dummy=mdsvalue('\XTOMO::TOP.RESULTS.CORE.CONFIG:EPS',/quiet,status=status)
	mdsclose,'xtomo',shot
	IF status THEN output=1 ELSE output=0
	RETURN,output
END

FUNCTION is_xtomo_asym,shot
	mdsopen,'xtomo',shot
	dummy=mdsvalue('dim_of(\XTOMO::TOP.RESULTS.CORE:EMISS,5)',/quiet,status=status)
	IF status THEN RETURN,1 ELSE RETURN,0
END

;NAME:
;	XRAY_TOMO_INFO
;
;PURPOSE:
;	This function loads geometry data from the tree to create an planar INFO file for use with GENPOS.
;
;CALLING SEQUENCE:
;	result=XRAY_TOMO_INFO(array)
;
;INPUTS:
;	array:		INT of the array number (see procedure)
;
;OPTIONAL INPUTS:
;	shot:		LONG of the shot number to load from DEFAULT: -1
;
;KEYWORD PARAMETERS:
;	debug:		/debug stops the code before the output is returned
;
;OUTPUTS:
;	result:		STRUC in the INFO file format.  (See GENPOS help)
;	
;PROCEDURE:
;	Data is loaded from the \MHD::TOP.XTOMO.GEOMETRY node.  The array value uses the tree
;	notation.
;			1 - core top
;			2 - edge top
;			3 - core outboard
;			4 - edge outboard	
;
;MODIFICATION HISTORY:
;	Written by:	ML Reinke - Spring 2006
;
;-


FUNCTION xray_tomo_info,array_num,debug=debug,shot=shot
	a=array_num-1
	IF NOT keyword_set(shot) THEN shot=-1
	mdsopen, 'xtomo',shot
	angle_aper=mdsvalue('\XTOMO::TOP.GEOMETRY:ANGLE_APER')
	angle_array=mdsvalue('\XTOMO::TOP.GEOMETRY:ANGLE_ARRAY')
	length_aper=mdsvalue('\XTOMO::TOP.GEOMETRY:LENGTH_APER')
	length_det=mdsvalue('\XTOMO::TOP.GEOMETRY:LENGTH_DET')
	num_in_array=mdsvalue('\XTOMO::TOP.GEOMETRY:NUM_IN_ARRAY')
	r_aperture=mdsvalue('\XTOMO::TOP.GEOMETRY:R_APERTURE')
	r_array=mdsvalue('\XTOMO::TOP.GEOMETRY:R_ARRAY')
	spacing_det=mdsvalue('\XTOMO::TOP.GEOMETRY:SPACING_DET')
	width_aper=mdsvalue('\XTOMO::TOP.GEOMETRY:WIDTH_APER')
	width_det=mdsvalue('\XTOMO::TOP.GEOMETRY:WIDTH_DET')
	z_aperture=mdsvalue('\XTOMO::TOP.GEOMETRY:Z_APERTURE')
	z_array=mdsvalue('\XTOMO::TOP.GEOMETRY:Z_ARRAY')
	mdsclose,'xtomo',shot

	name='xray_tomo_'+num2str(a,1)
	ap_vec=[r_aperture[a],0.0,z_aperture[a]]
	ap_rot=[0.0,(angle_aper[a]-90.0)*!pi/180.0,!pi/2.0]
	ap_size=[length_aper[a],width_aper[a]]
	ap={vec:ap_vec,rot:ap_rot,size:ap_size}

	th=angle_aper[a]*!pi/180.0
	phi=angle_array[a]*!pi/180.0
	del_r=r_array[a]-r_aperture[a]
	del_z=z_array[a]-z_aperture[a]

	IF th GT !pi/2.0 THEN x0=[-del_r*sin(!pi-th)-del_z*cos(!pi-th),0.0,-del_r*cos(!pi-th)+del_z*sin(!pi-th)] ELSE $
		 x0=[-del_r*sin(th)+del_z*cos(th),0.0,del_r*cos(th)+del_z*sin(th)]
	x1=x0+[-sin(th-phi),0.0,cos(th-phi)]*0.01
	x2=x0+[0.0,0.01,0.0]
	
	xi=fltarr(num_in_array[a])
	zeta=fltarr(num_in_array[a])

	FOR i=0,num_in_array[a]-1 DO zeta[i]=((num_in_array[a]-1)/2.0-i)*spacing_det[a]
	det_size=[length_det[a],width_det[a]]

	det={x0:x0,x1:x1,x2:x2,xi:xi,zeta:zeta,size:det_size}
	type='planar'
	author='mlreinke/granetz'

	info={name:name,ap:ap,det:det,type:type,author:author}

	output=info
	IF keyword_set(debug) THEN stop
	RETURN,output
END

;+
;NAME:
;	XRAY_TOMO_POS
;	
;PURPOSE:
;	This function uses XRAY_TOMO_INFO AND GENPOS_PLANAR2POS to find the pos vector and etendue for one or more arrays.
;
;CALLING SEQUENCE:
;	result=XRAY_TOMO_POS()
;
;KEYWORD PARAMETERS:
;	reform:		/reform will run GENPOS_POS_REFORM, necessary for useful POS data from ARRAY 1.  RZBND=[0.44,1.0,0.4,-0.4]
;
;OPTIONAL INPUTS:
;	shot:		LONG shot number DEFAULT: -1
;	array:		INT  array number [1,2,3,4] DEFAULT: 1
;
;OUTPUTS:
;	result:		FLTARR [4,n_ch] of the POS vector for the array(s)
;
;OPTIONAL OUTPUTS:
;	etendue:	FLTARR [n_ch] of the etendue for each channel [m^2-str]
;	rz_ap:		FLTARR [Ro,Zo] of the aperture location
;
;MODIFICATION HISTORY:
;	Written by:	ML Reinke - Spring 2006
;	1-31-08:	ML Reinke - added the reform keyword [use with ARRAY 1]
;
;-

FUNCTION xray_tomo_pos,etendue=etendue,array=array,shot=shot,rz_ap=rz_ap,reform=reform,zshift=zshift,rshift=rshift,dbeta=dbeta,t=t
	IF NOT keyword_set(array) THEN array=1
	IF NOT keyword_set(zshift) THEN zshift=0.0
	IF NOT keyword_set(rshift) THEN rshift=0.0
	rz_ap=fltarr(2,n(array)+1)
	FOR i=0,n(array) DO BEGIN
		info=xray_tomo_info(array[i],shot=shot)
		IF keyword_set(dbeta) THEN info.ap.rot[1]+=dbeta*!pi/180.0
		rz_ap[*,i]=[info.ap.vec[0],info.ap.vec[2]]
		pos=genpos_planar2pos(info,etendue=U,t=t,umod=umod)
		pos[0,*]+=rshift
		pos[1,*]+=zshift
		IF i EQ 0 THEN out_pos=pos ELSE out_pos=[[out_pos],[pos]]
		IF i EQ 0 THEN etendue=U*umod ELSE etendue=[etendue,U*umod]
	ENDFOR
	IF keyword_set(reform) THEN BEGIN
		rzbnd=[0.44,1.0,0.4,-0.4]
		genpos_pos_reform,out_pos,rzbnd
	ENDIF
	RETURN,out_pos
END

FUNCTION xray_tomo_gpv,ves_grid,n_ap=n_ap,n_det=n_det,gpv_contour=gpv_contour,load=load,detector=detector,kdebug=kdebug,$
		array=array,shot=shot,ves_cent=ves_cent
	IF NOT keyword_set(array) THEN array=0
	out_path='/home/mlreinke/idl/genie/data/gpv/xtomo/'
	file_path='xtomo_'
	FOR i=0,n(array) DO BEGIN
		info=xray_tomo_info(array[i],shot=shot)
		path=out_path+file_path+num2str(array[i],1)+'.dat'
		IF keyword_set(load) THEN restore, path, /verb ELSE $
			gpv=genpos_planar2gpv(info,ves_grid,n_ap=n_ap,n_det=n_det,path=path,gpv_contour=gpv_contour,$
				detector=detector,kdebug=kdebug)
		IF i EQ 0 THEN output=gpv ELSE output=[output,gpv]
	ENDFOR
	
	RETURN,output
END


;+
;NAME:
;	XTOMO_GET_GOOD
;
;PURPOSE:
;	This function loads the GOOD vector for use with the brightness profiles.
;
;CALLING SEQUENCE:
;	result=XTOMO_GET_GOOD(array)
;
;INPUTS:
;	array:	INT of the array [1,2,3,4]
;
;OPTIONAL INPUTS:
;	shot:	LONG of the shot number DEFAULT: -1
;
;OUTPUTS:
;	result:	INTARR [n_ch] of 1's and 0's showing which channels are good
;
;PROCEDURE:
;	Currently this is hardcoded...staytuned.
;
;MODIFICATION HISTORY:
;	Written by:	ML Reinke - 12/07
;	1-30-07:	ML Reinke - updated good for array 3
;	7-31-07:	ML Reinke - updated shot evolving channel maps for FY05 and FY06
;
;-

FUNCTION xtomo_get_good,array,shot=shot,core=core
	good=intarr(38)+1
	IF NOT keyword_set(shot) THEN shot = -1
	CASE array OF 
		1 : 	BEGIN
				good=good
			END	
		2 : 	BEGIN
				good[12]=0
				IF shot GT 1060309000 AND shot LT 1060729000 THEN good[[12,make(26,37,12)]]=0

			END
		3 : 	BEGIN	
				IF shot GT 1040309000 AND shot LT 1040311000 THEN good[[14,25]]=0
				IF shot GT 1050301000 AND shot LT 1060309000 THEN good[[14, make(20,37,18)]]=0
				IF shot GT 1060309000 AND shot LT 1060729000 THEN good[[0,1,2,3,4,14,25]]=0
				good[[25]]=0
			END
		4 : 	BEGIN
				good=good
				IF shot GT 1101022000 THEN good[7]=0
			END
		5 :	BEGIN
				good=good
			END
	ENDCASE

	IF keyword_set(core) THEN BEGIN
		mdsopen,'xtomo',shot
		good=mdsvalue('\XTOMO::TOP.RESULTS.CORE.CONFIG:GOOD')
		mdsclose,'xtomo',shot
	ENDIF

	RETURN,good
END

PRO xtomo_put_good,shot,good,array=array
	mdsopen,'xtomo',shot
	IF keyword_set(array) THEN BEGIN
		array_str=num2str(array,1)
		mdsput,'\XTOMO::TOP.RESULTS.ARRAY_'+array_str+'.CONFIG:GOOD','build_with_units($,"")',float(good)
	ENDIF ELSE mdsput,'\XTOMO::TOP.RESULTS.CORE.CONFIG:GOOD','build_with_units($,"")',float(good)
	mdsclose,'xtomo',shot
END


;+
;NAME:
;	CALC_XTOMOBR_DATA
;
;PURPOSE:
;	This procedure loads brightness data from the XTOMO tree, time averages and forms it into a 2D array.
;
;CALLING SEQUENCE
;
;	CALC_XTOMOBR_DATA,shot,array,br,ch,t
;
;INPUTS:
;	shot:	LONG shot number
;	array:	INT array number [1,2,3,4]
;	
;OPTIONAL INPUTS:
;	del_i	INT number of time points to average over DEFAULT = 5 (50 kHz-> 10 kHz)
;	good:	INTARR [n_ch] of 1's and 0's to toggle channels DEFAULT: good from XTOMO_GET_GOOD
;
;KEYWORD PARAMETERS:
;	rt:	/rt fills ch with rtang values instead of channel numbers
;
;OUTPUTS:
;	br:	FLTARR [n_ch, n_time] of the brightness [W/m^2] 
;	ch:	INTARR [n_ch] of the channel number (not really useful) (see optional output)
;	t:	FLTARR [n_time] of the time points
;
;OPTIONAL OUTPUTS:
;	ch:	FLTARR [n_ch,n_time] of the midplane "tangency radius" found by tracing the line of
;			sight through the flux surfaces.  Uses GENPOS_POS2RMIDTANG.
;
;PROCEDURE:
;	Data is loaded from the \XTOMO::TOP.SIGNALS nodes and the etendue is calculated from XRAY_TOMO_POS.
;	Background subtraction is done for each channel seperately using an average from the first point to
;	the point where t=-20 ms.  The data is averaged over DEL_I time points and the new time point is
;	the center.  DEL_I is forced to DEL_I-1 if it is an even number.
;
;MODIFICATION HISTORY:
;	Written by:	ML Reinke - 12/07
;	1-17-08		ML Reinke - added the rt keyword which fills CH with RTANG values
;	7-23-08:	ML Reinke - included the shot optional input in the call to XRAY_TOMO_POS
;	7-31-08:	ML Reinke - included the shot optional input in the call to XRAY_GET_GOOD
;				    if timing data is not consistent across array then interpolate 
;				    to finer time grid (bad for fast events)
;
;-

PRO calc_xtomobr_data,shot,array,br,ch,t,del_i=del_i,good=good,rt=rt,brerr=brerr,dbeta=dbeta,nopos=nopos
	IF NOT keyword_set(del_i) THEN IF shot LT 1071201000 THEN del_i=3 ELSE del_i=5
	pos=xray_tomo_pos(array=array,etendue=u,shot=shot,dbeta=dbeta,t=0.00005)
	mdsopen,'xtomo',shot
	det=mdsvalue('\XTOMO::TOP.GEOMETRY:NUM_IN_ARRAY',/quiet,status=status)
	IF NOT status THEN num_ch=38 ELSE num_ch=det[array-1]
	IF NOT keyword_set(good) THEN good=xtomo_get_good(array,shot=shot)
	data=mdsvalue('\XTOMO::TOP.SIGNALS.ARRAY_'+num2str(array,1)+':CHORD_01')
	t=mdsvalue('dim_of(\XTOMO::TOP.SIGNALS.ARRAY_'+num2str(array,1)+':CHORD_01)')
	IF min(t) GT -0.02 THEN min_tpt=t[0]/2.0 ELSE min_tpt=-0.02
	i_bl=ipt(t,min_tpt)
	br=fltarr(num_ch,n(t)+1)
	brerr=fltarr(num_ch,n(t)+1)
	FOR i=1,num_ch DO BEGIN
		array_str=num2str(array,1)
		IF i LT 10 THEN ch_str='0'+num2str(i,1) ELSE ch_str=num2str(i,1)
		br_ch=mdsvalue('\XTOMO::TOP.SIGNALS.ARRAY_'+array_str+':CHORD_'+ch_str,status=status,/quiet)
		IF NOT status THEN br_ch=fltarr(n(t)+1)
		br_ch*=good[i-1]/u[i-1]*4.0*!pi
		IF n(br_ch) EQ n(t) THEN br[i-1,*]=br_ch ELSE BEGIN
			ch_t=mdsvalue('dim_of(\XTOMO::TOP.SIGNALS.ARRAY_'+array_str+':CHORD_'+ch_str+')')
			br[i-1,*]=interpol(br_ch,ch_t,t)
		ENDELSE
		br[i-1,*]-=mean(br[i-1,0:i_bl])
	ENDFOR
	;rad=mdsvalue('\XTOMO::TOP.BRIGHTNESSES.ARRAY_'+num2str(array,1)+':CHORD_RADII')
	;th=mdsvalue('\XTOMO::TOP.BRIGHTNESSES.ARRAY_'+num2str(array,1)+':CHORD_ANGLES')
	mdsclose,'xtomo',shot
		
	;define output data
	IF del_i NE 1 THEN BEGIN
		i = long((del_i-1)/2)
		cntr=0
		time_new = fltarr(floor((n(t)+1)/del_i))
		diodebr_new = fltarr(num_ch, n(time_new)+1)
		tmp = fltarr(num_ch)
		FOR cntr=0L,n(time_new) DO BEGIN
			time_new[cntr]=t[i]
			FOR j=i-(del_i-1)/2,i+(del_i-1)/2 DO tmp+=br[*,j]
			diodebr_new[*,cntr]=tmp/del_i
			i+=del_i
			tmp*=0.0
		ENDFOR
		t = time_new						;replace arrays				
		br=diodebr_new
	ENDIF 

	IF keyword_set(rt) THEN ch=genpos_pos2rmidtang(pos,shot,t,/efit) ELSE ch=indgen(num_ch)
END

;+
;NAME:
;	CALC_XTOMO_EMISS_DATA
;	
;PURPOSE:
;	This procedure calculates the emissivity profiles for the edge xray tomography system assuming flux
;	surface symmetry.  GENPOS_POS2VOXEL_MATRIX is used to generate the spatial weightings and
;	GENPOS_PROFILE_INVERT is used to invert the data.
;
;INPUTS:
;	shot:		LONG 	shot number
;	array:		INT	edge array number [2,4]
;
;OPTIONAL INPUTS:
;	bright:		FLTARR 	[n_ch,n_t] of the brightness [W/m^2] DEFAULT calls CALC_XTOMOBR_DATA
;	t:		FLTARR	[n_t] of the time points [sec] DEFAULT calls CALC_XTOMOBR_DATA
;	del_i:		INT	used in call to CALC_XTOMOBR_DATA
;	good:		INTARR	of 1's and 0's for operational channels DEFAULT calls XTOMO_GET_GOOD
;	eps:		FLT	of the smoothing weighting factor DEFAULT = 1.0
;	eta:		FLT 	of the edge zero weighting factor DEFAULT = 1.0
;	n_rho:		INT	number of rhopts (final data will be n_rho-2) DEFAULT = 30
;	rho_vec:	FLTARR	of rho values to use for inversion (see PROCEDURE for array specific DEFAULTS)
;	n_s:		INT	number of points to divide line of sight into for GENPOS_POS2VOXEL_MATRIX
;	
;OUTPUTS:
;	emiss:		FLTARR	[n_rho-2,n_time] of the x-ray emissivity [W/m^3]
;	rho:		FLTARR	[n_rho-2] of the rho values for the emissivity
;	time:		FLTARR	[n_time] of the time points [sec]
;	brchk:		FLTARR	[n_good,n_time] of the brightnesses calculated from the emissivity [W/m^2]
;
;OPTIONAL OUTPUTS:
;	rmaj:		FLTARR	[n_rho-2,n_time] of the outboard major radii of the rho values [m]
;	
;	If unassigned variables are used for bright,ch and t then they become optional outputs.
;
;PROCEDURE:
;	If the brightness data is not input via the optional inputs than CALC_XTOMOBR_DATA is used to load the data.
;	GENPOS_POS2VOXEL_MATRIX is used to calculated the weighting matrix using rho_vec specific to each array:
;		2 : rho_vec=make(0.88,1.05,n_rho)
;		4 : rho_vec=make(0.86,1.05,n_rho)
;	The voxel matrix is run on the EFIT time base and then interpolated for use with the much faster x-ray data.
;	The bad channels identified by GOOD are removed from bright and voxel and the brightness is truncated to the EFIT
;	time window.  GENPOS_PROFILE_INVERT is used to do the inversion.
;
;MODFICATION HISTORY:
;	Written by:	ML Reinke - 1/22/08
;
;-		

PRO calc_xtomo_emiss_data,shot,array,emiss,rho,time,brchk,rmaj=rmaj,bright=bright,ch=ch,t=t,del_i=del_i,good=good,$
		debug=debug,eps=eps,eta=eta,n_rho=n_rho,rho_vec=rho_vec,n_s=n_s,tr=tr

	;load brightness data if not provided
	IF NOT keyword_set(bright) AND NOT keyword_set(t) THEN calc_xtomobr_data,shot,array,bright,ch,t,del_i=del_i,good=good,rt=rt
	print, 'brightness data loaded'
	IF NOT keyword_set(good) THEN good=xtomo_get_good(array)		

	;setup inversion data
	IF NOT keyword_set(n_rho) THEN IF keyword_set(rho_vec) THEN n_rho=n(rho_vec)+1 ELSE n_rho=30
	IF NOT keyword_set(eta) THEN eta=1.0
	IF NOT keyword_set(eps) THEN eps=1.0
	IF NOT keyword_set(rho_vec) THEN BEGIN
		CASE array OF 
			2 : rho_vec=make(0.88,1.05,n_rho)
			4 : rho_vec=make(0.86,1.05,n_rho)
		     else : RETURN
		ENDCASE
	ENDIF
	goodpts=where(good EQ 1)
	bright=bright[goodpts,*]
	n_good=total(good)	

	;generate voxel matrix on subset of EFIT time grid points
	t_vox=line_gettimes(shot)
	IF keyword_set(tr) THEN BEGIN
		min_t=min(t_vox) > tr[0]
		max_t=max(t_vox) < tr[1]
	ENDIF ELSE BEGIN
		min_t=min(t_vox)
		max_t=max(t_vox)
	ENDELSE
	i_low=ipt(t_vox,min_t)
	i_high=ipt(t_vox,max_t)
	t_vox=t_vox[i_low:i_high]
	print, 'doing inversion for '+num2str(min_t,dp=2)+' < t < '+num2str(max_t,dp=2)	

	;generate voxel data
	pos=xray_tomo_pos(array=array,etendue=u)
	info=xray_tomo_info(array)
	voxel=genpos_pos2voxel_matrix(pos,u,shot,tpts=t_vox,rho_vec=rho_vec,rz_ap=[info.ap.vec[0],info.ap.vec[2]],/verb,rhopts=rhopts,n_s=n_s)
	voxel=voxel[goodpts,*,*]	;remove bad channels
	vox_inner_view=reform(voxel[*,0,*])	
	voxel=voxel[*,1:n_rho-2,*]		;remove outer and inner voxel
	voxel[*,0,*]+=vox_inner_view		;add old inner to new inner
	rho=rho_vec[1:n_rho-2]		;remove outer and inner rho
	n_rho=n(rho)+1

	IF keyword_set(debug) THEN stop

	;truncate brightness and time to be within EFIT time
	i_low=ipt(t,min(t_vox))+1
	i_high=ipt(t,max(t_vox))-1
	t=t[i_low:i_high]
	bright=bright[*,i_low:i_high]
	n_time=n(t)+1

	;generate output data
	rmaj=line_getrmaj(shot,rho,t)
	time=t
	emiss=fltarr(n_rho,n_time)
	brchk=fltarr(n_good,n_time)
	print, 'inverting profiles'
	FOR i=0L,n_time-1 DO BEGIN
		;interpolate voxel matrix that was generated on EFIT time points
		ibnd=ibound(t_vox,t[i])
		vox_i=voxel[*,*,ibnd[0]]+(voxel[*,*,ibnd[1]]-voxel[*,*,ibnd[0]])/(t_vox[ibnd[1]]-t_vox[ibnd[0]])*(t[i]-t_vox[ibnd[0]])
		a=max(vox_i)
		pow=bright[*,i]*u[goodpts]/(4.0*!pi)
		IF total(vox_i) NE 0 THEN BEGIN
			emiss[*,i]=genpos_profile_invert(pow,vox_i,double(eps*a^2),brchk=brchk_i,/nofirst,eta=double(eta*a^2))	
			brchk[*,i]=brchk_i.mom*4.0*!pi/u[goodpts]
		ENDIF
	ENDFOR
	print, 'inversions done'
	IF keyword_set(debug) THEN stop
END


;+
;NAME:
;	CALC_XTOMOCORE_EMISS_DATA
;	
;PURPOSE:
;	This procedure calculates the emissivity profile for the core xray tomography system assuming flux
;	surface symmetry.  Deviations seen on the BRCHK are useful for observing in/out and up/down asymmetries
;	GENPOS_POS2VOXEL_MATRIX is used to generate the spatial weightings and GENPOS_PROFILE_INVERT is used to 
;	invert the data.
;
;INPUTS:
;	shot:		LONG 	shot number
;
;OPTIONAL INPUTS:
;	del_i:		INT	used in call to CALC_XTOMOBR_DATA
;	good:		INTARR	[76] of 1's and 0's for operational channels DEFAULT calls XTOMO_GET_GOOD for array [1,3]
;	eps:		FLT	of the smoothing weighting factor DEFAULT = 0.75
;	eta:		FLT 	of the edge zero weighting factor DEFAULT = 1.0
;	n_rho:		INT	number of rhopts DEFAULT = 30
;	rho_vec:	FLTARR	of rho values to use for inversion (see PROCEDURE for array specific DEFAULTS)
;	n_s:		INT	number of points to divide line of sight into for GENPOS_POS2VOXEL_MATRIX
;	tr:		FLTARR 	[t_low,t_high] of the EFIT subset data to do the inversion
;	
;KEYWORD PARAMETERS:
;	array2		/array2 will use array2 in calculated the emissivity profile, brightness vector 1,3,2
;	array5		/array5 will use array5 in calculated the emissivity profile, brightness vector 1,3,5
;
;OUTPUTS:
;	emiss:		FLTARR	[n_rho,n_time] of the x-ray emissivity [W/m^3]
;	rho:		FLTARR	[n_rho] of the rho values for the emissivity
;	time:		FLTARR	[n_time] of the time points [sec]
;	brchk:		FLTARR	[n_good,n_time] of the brightnesses calculated from the emissivity [W/m^2]
;
;OPTIONAL OUTPUTS:
;	rmaj:		FLTARR	[n_rho,n_time] of the outboard major radii of the rho values [m]
;	bright:		FLTARR 	[n_good,n_time] of the brightness [W/m^2] DEFAULT calls CALC_XTOMOBR_DATA
;	ch:		INTARR	[n_good] of channel numbers [0-76] of the good channels
;	etree:		STRING	of the EFIT tree to use for inversion [DEFAULT: ANALYSIS]
;
;PROCEDURE:
;	CALC_XTOMOBR_DATA is used to load the data and a 76 channel [array_1,array_3] super array is formed (114 if /array2 or /array5)
;	GENPOS_POS2VOXEL_MATRIX is used to calculated the weighting matrix using its default rho_vec with n_rho as defined above.
;	The voxel matrix is run on the EFIT time base (or subset if tr is used) and then interpolated for use with the much faster x-ray data.
;	GENPOS_PROFILE_INVERT is used to do the inversion.
;
;MODFICATION HISTORY:
;	Written by:	ML Reinke - 1/31/08 (adapated from CALC_XTOMO_EMISS_DATA)
;	2/28/08:	ML Reinke - Fixed a bug that would mix up etendues if good was not the default.
;	1/27/10:	ML Reinke - Fixed a bug so that if GOOD is sent to it, it will break it up
;				    and send it to CALC_XTOMOBR_DATA
;	6/18/11		M.L. Reinke - added the ETREE optional input to allow for alternate EFIT trees
;	1/10/12		ML Reinke - added the ability to use array5 keyword
;	5/30/12		M.L. Reinke - fixed a bug that caused inv_matrix to be held constant over inversions that was put in during the /fast update
;-		

PRO calc_xtomocore_emiss_data,shot,emiss,rho,time,brchk,rmaj=rmaj,bright=bright,ch=ch,del_i=del_i,good=good,emerr=emerr,$
		debug=debug,eps=eps,eta=eta,n_rho=n_rho,rho_vec=rho_vec,n_s=n_s,tr=tr,array2=array2,array5=array5,psinorm=psinorm,fast=fast,dbeta=dbeta,etree=etree
	acase=0
	IF keyword_set(array2) AND NOT keyword_set(array5) THEN acase=1			;3 arrays 1,2,3
	IF NOT keyword_set(array2) AND keyword_set(array5) THEN acase=2			;3 arrays 1,3,5
	IF keyword_set(array2) AND keyword_set(array5) THEN acase=3			;4 arrays 1,2,3,5	;expand if/when need x4 system
	CASE acase OF 
		0 : BEGIN
			acstr='0, ARRAYS: 1,3'
			good1=good[0:38*1-1]
			good3=good[38:38*2-1]
		END
		1 : BEGIN
			acstr='1, ARRAYS: 1,2,3'
			good1=good[0:38*1-1]
			good3=good[38:38*2-1]
			good2=good[38*2:38*3-1]	
		END
		2 : BEGIN
			acstr='2, ARRAYS: 1,3,5'
			good1=good[0:38*1-1]
			good3=good[38:38*2-1]
			good5=good[38*2:38*3-1]	
		END
		ELSE : RETURN
	ENDCASE

	;load brightness data if not provided
	calc_xtomobr_data,shot,1,bright1,ch1,t_bright,del_i=del_i,good=good1,brerr=brerr1,dbeta=dbeta
	calc_xtomobr_data,shot,3,bright3,ch3,t,del_i=del_i,good=good3,brerr=brerr3
	IF keyword_set(array2) THEN calc_xtomobr_data,shot,2,bright2,ch2,t,del_i=del_i,good=good2,brerr=brerr2
	IF keyword_set(array5) THEN calc_xtomobr_data,shot,5,bright5,ch5,t,del_i=del_i,good=good5,brerr=brerr5
	print, 'ACASE='+acstr+' brightness data loaded - SHOT: '+num2str(shot,1)

	;load line of sight and etendue data
	pos1=xray_tomo_pos(array=1,etendue=u1,/reform,shot=shot)
	pos3=xray_tomo_pos(array=3,etendue=u3,shot=shot)
	IF keyword_set(array2) THEN pos2=xray_tomo_pos(array=2,etendue=u2,shot=shot)
	IF keyword_set(array5) THEN pos5=xray_tomo_pos(array=5,etendue=u5,shot=shot)

	CASE acase OF 
		0 : BEGIN
			bright=[bright1,bright3]
			err=[brerr1,brerr3]
			ch=[ch1,ch3+max(ch1)+1]
			good=[good1,good3] 
			pos=[[pos1],[pos3]]
			u=[u1,u3]
		END
		1 : BEGIN
			bright=[bright1,bright3,bright2]
			err=[brerr1,brerr3,brerr2]
			ch=[ch1,ch3+max(ch1)+1,ch2+max(ch1)+max(ch3)+2]	
			good=[good1,good3,good2]
			pos=[[pos1],[pos3],[pos2]]
			u=[u1,u3,u2]	
		END
		2 : BEGIN
			bright=[bright1,bright3,bright5]
			err=[brerr1,brerr3,brerr5]
			ch=[ch1,ch3+max(ch1)+1,ch5+max(ch1)+max(ch5)+2]	
			good=[good1,good3,good5]
			pos=[[pos1],[pos3],[pos5]]
			u=[u1,u3,u5]
		END
	END

	;setup inversion data
	IF NOT keyword_set(n_rho) THEN IF keyword_set(rho_vec) THEN n_rho=n(rho_vec)+1 ELSE n_rho=30
	IF keyword_set(rho_vec) THEN rho=rho_vec
	IF NOT keyword_set(eta) THEN eta=0.75
	IF NOT keyword_set(eps) THEN eps=1.0

	IF keyword_set(debug) THEN stop	

	;remove bad channels
	goodpts=where(good EQ 1)
	bright=bright[goodpts,*]
	err=err[goodpts,*]
	u=u[goodpts]
	pos=pos[*,goodpts]
	n_good=total(good)
	ch=ch[goodpts]

	;generate voxel matrix on subset of valid EFIT time grid points
	chk=efit_check(t_vox,ntvox,shot=shot,tree=etree)
	IF keyword_set(tr) THEN BEGIN
		min_t=min(t_vox) > tr[0]
		max_t=max(t_vox) < tr[1]
	ENDIF ELSE BEGIN
		min_t=min(t_vox)
		max_t=max(t_vox)
	ENDELSE
	i_low=ipt(t_vox,min_t)
	i_high=ipt(t_vox,max_t)
	t_vox=t_vox[i_low:i_high]
	print, 'doing inversion for '+num2str(min_t,dp=2)+' < t < '+num2str(max_t,dp=2)	
	voxel=genpos_pos2voxel_matrix(pos,u,shot,tpts=t_vox,rho_vec=rho,/verb,rhopts=rhopts,n_s=n_s,n_rho=n_rho,psinorm=psinorm,tree=etree)

	IF keyword_set(debug) THEN stop

	;truncate brightness and time to be within EFIT time
	i_low=ipt(t,min(t_vox))+1
	i_high=ipt(t,max(t_vox))-1
	t=t[i_low:i_high]
	bright=bright[*,i_low:i_high]
	n_time=n(t)+1

	;generate output data
	rmaj=line_getrmaj(shot,rho,t,psinorm=psinorm)
	time=t
	emiss=fltarr(n_rho,n_time)
	emerr=fltarr(n_rho,n_time)
	brchk=fltarr(n_good,n_time)
	start_time=systime(/seconds)
	IF keyword_set(fast) THEN print, 'inverting profiles - FAST' ELSE print,'inverting profiles'
	inv_matrix=0
	FOR i=0L,n_time-1 DO BEGIN
		ibnd=ibound(t_vox,t[i])
		IF keyword_set(fast) THEN BEGIN
			IF t[i]-t_vox[ibnd[0]] GE t_vox[ibnd[1]]-t[i] THEN BEGIN
				index=ibnd[1]
				IF t[i-1]-t_vox[ibnd[0]] LT t_vox[ibnd[1]]-t[i-1] THEN inv_matrix=0	;reset the inv_matrix if for the first time selecting new index
			ENDIF ELSE index=ibnd[0]
			vox_i=voxel[*,*,index]
		ENDIF ELSE BEGIN 	;interpolate voxel matrix that was generated on EFIT time points
			vox_i=voxel[*,*,ibnd[0]]+(voxel[*,*,ibnd[1]]-voxel[*,*,ibnd[0]])/(t_vox[ibnd[1]]-t_vox[ibnd[0]])*(t[i]-t_vox[ibnd[0]])
			inv_matrix=0
		ENDELSE
		a=max(vox_i)
		pow=bright[*,i]*u/(4.0*!pi)
		perr=err[*,i]*u/(4.0*!pi)
		IF total(vox_i) NE 0 THEN BEGIN
			emiss[*,i]=genpos_profile_invert(pow,vox_i,double(eps*a^2),brchk=brchk_i,eta=double(eta*a^2),err=perr,inv_matrix=inv_matrix)	
			brchk[*,i]=brchk_i.mom*4.0*!pi/u
			emerr[*,i]=brchk_i.inverr
		ENDIF
	ENDFOR
	ctime=systime(/seconds)-start_time
	print, num2str(n_time,1)+' inversions done in '+num2str(ctime,dp=2)
	IF keyword_set(debug) THEN stop
END

;+
;NAME:
;	CALC_XTOMOCORE_ASYMEMISS_DATA
;
;MODIFICATION HISTORY:
;	4/6/11		M.L. Reinke - allowed the inclusion of the m=2 cosine term
;	6/18/11		M.L. Reinke - added the ETREE optional input to allow for alternate EFIT trees
;-

PRO calc_xtomocore_asymemiss_data,shot,emiss,rho,time,emr,emc1,emc2,ems1,brchk,brchkasym,rmaj=rmaj,bright=bright,ch=ch,del_i=del_i,good=good,$
		debug=debug,eps=eps,eta=eta,n_rho=n_rho,rho_vec=rho_vec,n_s=n_s,tr=tr,psinorm=psinorm,fast=fast,dbeta1=dbeta1,dbeta3=dbeta3,$
		no_br1=no_br1, no_br3=no_br3,no_sin=no_sin,no_cos=no_cos,no_m2=no_m2,etree=etree
	
	IF NOT keyword_set(etree) THEN etree='ANALYSIS'
	IF keyword_set(no_br1) THEN no_cos=1		;if not using top view, cannot resolve cosine term
	IF keyword_set(no_br3) THEN no_sin=1		;if not using outboard view, cannot resolve sine term

	;load brightness data if not provided
	IF keyword_set(good) THEN BEGIN
		good1=good[0:37]
		good3=good[38:75]
	ENDIF
	calc_xtomobr_data,shot,1,bright1,ch1,t_bright,del_i=del_i,good=good1
	calc_xtomobr_data,shot,3,bright3,ch3,t,del_i=del_i,good=good3
	print, 'brightness data loaded - SHOT: '+num2str(shot,1)

	;load line of sight data
	pos1=xray_tomo_pos(array=1,etendue=u1,/reform,shot=shot,dbeta=dbeta1)
	pos3=xray_tomo_pos(array=3,etendue=u3,shot=shot,dbeta=dbeta3)
	IF keyword_set(no_br1) OR keyword_set(no_br3) THEN BEGIN
		IF keyword_set(no_br1) THEN BEGIN
			pos=pos3
			u=u3
			ch=ch3
			good=good3
			bright=bright3
		ENDIF
		IF keyword_set(no_br3) THEN BEGIN
			pos=pos1
			u=u1
			ch=ch1
			good=good1
			bright=bright1
		ENDIF
	ENDIF ELSE BEGIN
		pos=[[pos1],[pos3]]
		u=[u1,u3]
		bright=[bright1,bright3]
		ch=[ch1,ch3+max(ch1)+1]
		IF NOT keyword_set(good) THEN good=[good1,good3]

	ENDELSE

	;setup inversion configuration data
	IF NOT keyword_set(n_rho) THEN IF keyword_set(rho_vec) THEN n_rho=n(rho_vec)+1 ELSE n_rho=30
	IF keyword_set(rho_vec) THEN rho=rho_vec
	IF NOT keyword_set(eta) THEN eta=0.75
	IF NOT keyword_set(eps) THEN eps=1.0
	IF keyword_set(debug) THEN stop	

	;remove bad channels
	goodpts=where(good EQ 1)
	bright=bright[goodpts,*]
	u=u[goodpts]
	pos=pos[*,goodpts]
	n_good=total(good)
	ch=ch[goodpts]

	;generate voxel matrix on subset of EFIT time grid points
	IF shot EQ 1110201023 AND etree EQ 'EFIT09' THEN BEGIN		;temp hardcore because this shot has a "good" EFIT that isn't converged (apparently)
		t_vox=0.9
		ntvox=1
	ENDIF ELSE chk=efit_check(t_vox,ntvox,shot=shot,tree=etree)
		

	IF keyword_set(tr) THEN BEGIN
		min_t=min(t_vox) > tr[0]
		max_t=max(t_vox) < tr[1]
	ENDIF ELSE BEGIN
		min_t=min(t_vox)
		max_t=max(t_vox)
	ENDELSE
	i_low=ipt(t_vox,min_t)
	i_high=ipt(t_vox,max_t)
	t_vox=t_vox[i_low:i_high]
	print, 'doing inversion for '+num2str(min_t,dp=2)+' < t < '+num2str(max_t,dp=2)	
	voxel_r=genpos_pos2voxel_matrix(pos,u,shot,tpts=t_vox,rho_vec=rho,/verb,rhopts=rhopts,n_s=n_s,n_rho=n_rho,psinorm=psinorm,tree=etree)
	IF NOT keyword_set(no_cos) THEN BEGIN
		voxel_c1=genpos_pos2voxel_matrix(pos,u,shot,tpts=t_vox,rho_vec=rho,/verb,rhopts=rhopts,n_s=n_s,n_rho=n_rho,psinorm=psinorm,/cosine,m=1,tree=etree)
		;voxel_c1=genpos_pos2voxel_matrix(pos,u,shot,tpts=t_vox,rho_vec=rho,/verb,rhopts=rhopts,n_s=n_s,n_rho=n_rho,psinorm=psinorm,/inout)
		IF NOT keyword_set(no_m2) THEN voxel_c2=genpos_pos2voxel_matrix(pos,u,shot,tpts=t_vox,rho_vec=rho,/verb,rhopts=rhopts,n_s=n_s,n_rho=n_rho,$
							psinorm=psinorm,/cosine,m=2,tree=etree)
	ENDIF
	IF NOT keyword_set(no_sin) THEN voxel_s1=genpos_pos2voxel_matrix(pos,u,shot,tpts=t_vox,rho_vec=rho,/verb,rhopts=rhopts,n_s=n_s,n_rho=n_rho,psinorm=psinorm,/sine,m=1,tree=etree)
	
	IF keyword_set(debug) THEN stop

	;truncate brightness and time to be within EFIT time
	IF ntvox NE 1 THEN BEGIN
		i_low=ipt(t,min(t_vox))+1
		i_high=ipt(t,max(t_vox))-1
		t=t[i_low:i_high]
		bright=bright[*,i_low:i_high]
	ENDIF ELSE BEGIN	;assume that single EFIT is to be used for all points in trange
		tmp=where(t GE tr[0] AND t LE tr[1])
		t=t[tmp]
		bright=bright[*,tmp]
	ENDELSE
	n_time=n(t)+1

	;generate output data
	rmaj=line_getrmaj(shot,rho,t,psinorm=psinorm)
	time=t
	emiss=fltarr(n_rho,n_time)
	emr=fltarr(n_rho,n_time)
	emc1=fltarr(n_rho,n_time)
	ems1=fltarr(n_rho,n_time)
	emc2=fltarr(n_rho,n_time)
	brchk=fltarr(n_good,n_time)
	brchkasym=fltarr(n_good,n_time)
	start_time=systime(/seconds)
	IF keyword_set(fast) THEN print, 'inverting profiles - FAST' ELSE print,'inverting profiles'

	invcase=0									;include all terms (default)
	IF keyword_set(no_m2) AND NOT keyword_set(no_sin) THEN invcase=1		;include only m=1 terms
	IF keyword_set(no_cos) THEN invcase=2						;include only odd terms
	IF keyword_set(no_sin) THEN invcase=3						;include only even terms
	IF keyword_set(no_m2) AND keyword_set(no_sin) THEN invcase=4			;include only m=1 cosine 

	inv_matrix_r=0
	inv_matrix_a=0
	FOR i=0L,n_time-1 DO BEGIN
		;interpolate voxel matrix that was generated on EFIT time points
		ibnd=ibound(t_vox,t[i])
		IF keyword_set(fast) THEN BEGIN
			IF t[i]-t_vox[ibnd[0]] GE t_vox[ibnd[1]]-t[i] THEN BEGIN
				index=ibnd[1]
				IF t[i-1]-t_vox[ibnd[0]] LT t_vox[ibnd[1]]-t[i-1] THEN BEGIN
					inv_matrix_r=0							;reset the inv_matrix if for the first time selecting new index
					inv_matrix_a=0
				ENDIF
			ENDIF ELSE index=ibnd[0]
			vox_r=voxel_r[*,*,index]
			IF invcase NE 2 THEN vox_c1=voxel_c1[*,*,index]
			IF invcase NE 1 AND invcase NE 2 AND invcase NE 4 THEN vox_c2=voxel_c2[*,*,index]
			IF invcase NE 3 AND invcase NE 4 THEN vox_s1=voxel_s1[*,*,index]
		ENDIF ELSE BEGIN
			IF ntvox NE 1 THEN BEGIN
				vox_r=voxel_r[*,*,ibnd[0]]+(voxel_r[*,*,ibnd[1]]-voxel_r[*,*,ibnd[0]])/(t_vox[ibnd[1]]-t_vox[ibnd[0]])*(t[i]-t_vox[ibnd[0]])
				IF invcase NE 2 THEN vox_c1=voxel_c1[*,*,ibnd[0]]+(voxel_c1[*,*,ibnd[1]]-voxel_c1[*,*,ibnd[0]])/(t_vox[ibnd[1]]-t_vox[ibnd[0]])*(t[i]-t_vox[ibnd[0]])
				IF invcase NE 1 AND invcase NE 2 AND invcase NE 4 THEN vox_c2=voxel_c2[*,*,ibnd[0]]+(voxel_c2[*,*,ibnd[1]]-voxel_c1[*,*,ibnd[0]])/$
					(t_vox[ibnd[1]]-t_vox[ibnd[0]])*(t[i]-t_vox[ibnd[0]])
				IF invcase NE 3 AND invcase NE 4 THEN vox_s1=voxel_s1[*,*,ibnd[0]]+(voxel_s1[*,*,ibnd[1]]-voxel_s1[*,*,ibnd[0]])/(t_vox[ibnd[1]]-t_vox[ibnd[0]])*(t[i]-t_vox[ibnd[0]])
				inv_matrix_r=0
				inv_matrix_a=0
			ENDIF ELSE BEGIN
				vox_r=voxel_r
				IF invcase NE 2 THEN vox_c1=voxel_c1
				IF invcase NE 1 AND invcase NE 2 AND invcase NE 4 THEN vox_c2=voxel_c2
				IF invcase NE 3 AND invcase NE 4 THEN vox_s1=voxel_s1
				inv_matrix_r=0
				inv_matrix_a=0
			ENDELSE
		ENDELSE

		;invert for just the radial profile
		vox_i=vox_r
		a=max(vox_i)
		pow=bright[*,i]*u/(4.0*!pi)
		IF total(vox_i) NE 0 THEN BEGIN
			emiss[*,i]=genpos_profile_invert(pow,vox_i,double(eps*a^2),brchk=brchk_i,eta=double(eta*a^2),inv_matrix=inv_matrix_r)	
			brchk[*,i]=brchk_i.mom*4.0*!pi/u
		ENDIF

		;invert assuming selected asymmetry parameters
		CASE invcase OF 
			0 : BEGIN
				vox_i=[[vox_r],[vox_c1],[vox_c2],[vox_s1]]			;use all terms
				nprof=4
			END
			1 : BEGIN
				vox_i=[[vox_r],[vox_c1],[vox_s1]]				;use only m=1	
				nprof=3	
			END
			2 : BEGIN
				vox_i=[[vox_r],[vox_s1]]					;use only sine
				nprof=2
			END
			3 : BEGIN
				vox_i=[[vox_r],[vox_c1],[vox_c2]]				;use only m=1,m=2 cosine
				nprof=3
			END
			4 : BEGIN
				vox_i=[[vox_r],[vox_c1]]					;use only cosine
				nprof=2
			END
		ENDCASE
		IF keyword_set(debug) THEN stop
		a=max(vox_i)
		pow=bright[*,i]*u/(4.0*!pi)
		IF total(vox_i) NE 0 THEN BEGIN
			em_tot=genpos_profile_invert(pow,vox_i,double(eps*a^2),brchk=brchk_i,eta=double(eta*a^2),nprof=nprof,inv_matrix=inv_matrix_a)	
			emr[*,i]=em_tot[0:n_rho-1]
			CASE invcase OF 
				0 : BEGIN
					emc1[*,i]=em_tot[n_rho:2.0*n_rho-1]
					emc2[*,i]=em_tot[2.0*n_rho:3.0*n_rho-1]
					ems1[*,i]=em_tot[3.0*n_rho:4.0*n_rho-1]
				END
				1 : BEGIN
					emc1[*,i]=em_tot[n_rho:2.0*n_rho-1]
					ems1[*,i]=em_tot[2.0*n_rho:3.0*n_rho-1]
				END
				2 : BEGIN
					ems1[*,i]=em_tot[n_rho:2.0*n_rho-1]
				END
				3 : BEGIN
					emc1[*,i]=em_tot[n_rho:2.0*n_rho-1]
					emc2[*,i]=em_tot[2*n_rho:3.0*n_rho-1]
				END
				4 : BEGIN
					emc1[*,i]=em_tot[n_rho:2.0*n_rho-1]
				END
			ENDCASE
			brchkasym[*,i]=brchk_i.mom*4.0*!pi/u
		ENDIF
	ENDFOR
	ctime=systime(/seconds)-start_time
	print, num2str(n_time,1)+' inversions done in '+num2str(ctime,dp=2)
	IF keyword_set(debug) THEN stop
END



;old script for adding the RESULTS node to the tree
PRO xtomo_add_results2tree,shot,force=force,debug=debug
	mdsopen,'xtomo',shot
	dummy=mdsvalue('\XTOMO::TOP.RESULTS:COMMENT',/quiet,status=status)
	mdsclose,'xtomo',shot
	IF keyword_set(debug) THEN stop
	
	;if force then delete the node before adding (if present)
	IF keyword_set(force) AND status THEN BEGIN
		mdstcl, "set verify"
		mdstcl, 'edit xtomo /shot='+num2str(shot,1)
		mdstcl, 'delete node \XTOMO::TOP.RESULTS /noconfirm'
		mdstcl, 'write'
		mdstcl, 'close'
	ENDIF
	
	;add the nodes using the TCL script
	IF NOT status OR keyword_set(force) THEN BEGIN
		mdstcl, "set verify"
		mdstcl, 'edit xtomo /shot='+num2str(shot,1)
		mdstcl, 'add node \XTOMO::TOP.RESULTS'
		mdstcl, '@/home/mlreinke/xtomo/XTOMO_ADD_RESULTS.TCL'
		mdstcl, 'close'
	ENDIF ELSE IF NOT keyword_set(quiet) THEN print, 'SHOT: '+num2str(shot,1)+' has RESULT nodes'
END

;old script for adding the RESULTS.CORE node to the tree
PRO xtomo_add_core_results2tree,shot,force=force,debug=debug,new=new
	mdsopen,'xtomo',shot
	dummy=mdsvalue('\XTOMO::TOP.RESULTS.CORE.CONFIG:DEL_I',/quiet,status=status)
	mdsclose,'xtomo',shot
	IF keyword_set(debug) THEN stop
	
	;if force then delete the node before adding (if present)
	IF keyword_set(force) AND status THEN BEGIN
		mdstcl, "set verify"
		mdstcl, 'edit xtomo /shot='+num2str(shot,1)
		mdstcl, 'delete node \XTOMO::TOP.RESULTS.CORE /noconfirm'
		mdstcl, 'write'
		mdstcl, 'close'
	ENDIF
	
	;add the nodes using the TCL script
	IF NOT status OR keyword_set(force) THEN BEGIN
		mdstcl, "set verify"
		mdstcl, 'edit xtomo /shot='+num2str(shot,1)
		IF keyword_set(new) THEN mdstcl, 'add node \XTOMO::TOP.RESULTS'
		mdstcl, 'add node \XTOMO::TOP.RESULTS.CORE'
		mdstcl, '@/home/mlreinke/xtomo/XTOMO_ADD_CORE_RESULTS.TCL'
		mdstcl, 'close'
	ENDIF ELSE IF NOT keyword_set(quiet) THEN print, 'SHOT: '+num2str(shot,1)+' has core RESULT nodes'
END



;+
;NAME:
;	XTOMO_WRITE2TREE
;
;PURPOSE:
;	This procedure calculates emissivity and brightness data and inserts it into the tree
;
;CALLING SEQUENCE
;	XTOMO_WRITE2TREE,shot,array
;
;INPUTS:
;	shot:		LONG	shot number
;	array:		INT	array number [2,4]
;
;OPTIONAL INPUTS:
;	del_i:		INT	number of points to average brightness over before inversion DEFAULT: 5
;	good:		INTARR	of 1's and 0's describing good channels DEFAULT: see CALC_XTOMO_EMISS_DATA
;	eps:		FLT	smoothing weighting factor DEFAULT: 1.0
;	eta:		FLT	edge zero weighting factor DEFAULT: 1.0	
;	n_rho		INT	number of rho points DEFAULT: see CALC_XTOMO_EMISS_DATA
;	rho_vec:	FLTARR	values of the rho points to do the inversion on DEFAULT: see CALC_XTOMO_EMISS_DATA
;	n_s:		INT	number of points to divide line of sight into for calculating VOXEL matrix
;
;KEYWORD PARAMETERS:
;	no_emiss:	/no_emiss does not calculate emissivity but just stores BRIGHT, DEL_I and GOOD
;	debug:		/debug stops the code before the end for debugging purposes
;
;OUTPUTS:
;	All data is stored to the tree in array specific nodes in \XTOMO::TOP.RESULTS.ARRAY_X
;		\XTOMO::TOP.RESULTS.ARRAY_X:BRCHK - brightness checks [W/m^2]
;		\XTOMO::TOP.RESULTS.ARRAY_X:BRIGHT - brightness [W/m^2]
;		\XTOMO::TOP.RESULTS.ARRAY_X:EMISS - emissivity [W/m^3]
;		\XTOMO::TOP.RESULTS.ARRAY_X:RMAJ - major radius values for the inversion (time evolving)
;		\XTOMO::TOP.RESULTS.ARRAY_X:RTMID - outboard midplane values for the brightness and brightness check
;		\XTOMO::TOP.RESULTS.ARRAY_X.CONFIG:DEL_I
;		\XTOMO::TOP.RESULTS.ARRAY_X.CONFIG:EPS
;		\XTOMO::TOP.RESULTS.ARRAY_X.CONFIG:ETA
;		\XTOMO::TOP.RESULTS.ARRAY_X.CONFIG:GOOD
;		\XTOMO::TOP.RESULTS.ARRAY_X.CONFIG:N_S
;		\XTOMO::TOP.RESULTS.ARRAY_X.CONFIG:RHO_VEC
;
;PROCEDURE:
;	If /no_emiss is used, then the brightness data will be stored where the GOOD=0 channels will be forced to zero.  When the emissivity
;	is stored as well, then the brightness data has the bad channels removed.  Also, since EFIT is required for the inversion, the brightness
;	and emissivity data is only run between the first and last EFIT times and the datasets are truncated.
;
;MODIFICATION HISTORY:
;	Written by:	ML Reinke - 1/22/08
;
;-

PRO xtomo_write2tree,shot,array,del_i=del_i,good=good,debug=debug,no_emiss=no_emiss,eps=eps,eta=eta,n_rho=n_rho,rho_vec=rho_vec,n_s=n_s,tr=tr
	IF NOT keyword_set(del_i) THEN BEGIN
		CASE array OF
			2: 		del_i=5
			4: 		del_i=5
			ELSE :		del_i=1
		ENDCASE
	ENDIF

	IF NOT keyword_set(eps) THEN BEGIN
		CASE array OF
			2: 		eps=1.0
			4: 		eps=1.0
			ELSE :		eps=1.0
		ENDCASE
	ENDIF

	IF NOT keyword_set(eta) THEN BEGIN
		CASE array OF
			2: 		eta=1.0
			4: 		eta=1.0
			ELSE :		eta=1.0
		ENDCASE
	ENDIF

	array_str=num2str(array,1)

	pos=xray_tomo_pos(array=array)
	IF keyword_set(no_emiss) THEN BEGIN
		calc_xtomobr_data,shot,array,br,ch,time,del_i=del_i,good=good 
		rtmid=genpos_pos2rmidtang(pos,shot,time,/efit)	
	ENDIF ELSE BEGIN	
		calc_xtomo_emiss_data,shot,array,emiss,rho,time,brchk,rmaj=rmaj,bright=br,ch=ch,del_i=del_i,good=good,$
			eps=eps,eta=eta,n_rho=n_rho,rho_vec=rho_vec,n_s=n_s,tr=tr
		rtmid=genpos_pos2rmidtang(pos[*,where(good EQ 1)],shot,time,/efit)
	ENDELSE
	
	mdsopen,'xtomo',shot
	mdsput,'\XTOMO::TOP.RESULTS.ARRAY_'+array_str+'.CONFIG:DEL_I','build_with_units($,"")',del_i
	mdsput,'\XTOMO::TOP.RESULTS.ARRAY_'+array_str+'.CONFIG:GOOD','build_with_units($,"")',float(good)
	mdsput,'\XTOMO::TOP.RESULTS.ARRAY_'+array_str+':BRIGHT',$
        	'build_signal(build_with_units($1, "W/m^2"),*,'+$
                'build_with_units($2,"CH#"),'+$
          	'build_with_units($3,"seconds"))',br,ch,time
	mdsput,'\XTOMO::TOP.RESULTS.ARRAY_'+array_str+':RTMID',$
        	'build_signal(build_with_units($1, "m"),*,'+$
                'build_with_units($2,"CH#"),'+$
          	'build_with_units($3,"seconds"))',rtmid,ch,time
	IF NOT keyword_set(no_emiss) THEN BEGIN
		mdsput,'\XTOMO::TOP.RESULTS.ARRAY_'+array_str+':EMISS',$
        		'build_signal(build_with_units($1, "W/m^3"),*,'+$
	                'build_with_units($2,""),'+$
        	  	'build_with_units($3,"seconds"))',emiss,rho,time
		mdsput,'\XTOMO::TOP.RESULTS.ARRAY_'+array_str+':RMAJ',$
        		'build_signal(build_with_units($1, "m"),*,'+$
	                'build_with_units($2,""),'+$
        	  	'build_with_units($3,"seconds"))',rmaj,rho,time
		mdsput,'\XTOMO::TOP.RESULTS.ARRAY_'+array_str+':BRCHK',$
        		'build_signal(build_with_units($1, "W/m^2"),*,'+$
	                'build_with_units($2,"CH#"),'+$
        	  	'build_with_units($3,"seconds"))',brchk,ch,time
		mdsput,'\XTOMO::TOP.RESULTS.ARRAY_'+array_str+'.CONFIG:EPS','build_with_units($,"")',eps
		mdsput,'\XTOMO::TOP.RESULTS.ARRAY_'+array_str+'.CONFIG:ETA','build_with_units($,"")',eta
		mdsput,'\XTOMO::TOP.RESULTS.ARRAY_'+array_str+'.CONFIG:N_S','build_with_units($,"")',n_s
		mdsput,'\XTOMO::TOP.RESULTS.ARRAY_'+array_str+'.CONFIG:RHO_VEC','build_with_units($,"")',rho_vec
	ENDIF	
	mdsclose,'xtomo',shot
	
	IF keyword_set(debug) THEN stop
END

;+
;NAME:
;	XTOMO_CORE_WRITE2TREE
;
;PURPOSE:
;	This procedure calculates emissivity and brightness data and inserts it into the tree for the core xray arrays.
;	Flux surface symmetry is assumed (erroneously past 1st order) and the resulting brightness checks can be used to see
;	what type of asymmetry is present (in/out or up/down) so further analysis can be done.
;
;CALLING SEQUENCE
;	XTOMO_CORE_WRITE2TREE,shot
;
;INPUTS:
;	shot:		LONG	shot number
;
;OPTIONAL INPUTS:
;	del_i:		INT	number of points to average brightness over before inversion DEFAULT: from TREE
;	good:		INTARR	of 1's and 0's describing good channels DEFAULT: from TREE
;	eps:		FLT	smoothing weighting factor DEFAULT: see CALC_XTOMOCORE_EMISS_DATA
;	eta:		FLT	edge zero weighting factor DEFAULT: see CALC_XTOMOCORE_EMISS_DATA
;	n_rho		INT	number of rho points DEFAULT: see CALC_XTOMOCORE_EMISS_DATA
;	rho_vec:	FLTARR	values of the rho points to do the inversion on DEFAULT: see CALC_XTOMOCORE_EMISS_DATA
;	n_s:		INT	number of points to divide line of sight into for calculating VOXEL matrix
;	tr:		FLTARR 	[t_low,t_high] of the time interval to do inversions [sec] DEFAULT: from TREE
;	etree:		STRING	of the EFIT tree to use for inversion [DEFAULT: ANALYSIS]
;
;KEYWORD PARAMETERS:
;	no_emiss:	/no_emiss does not calculate emissivity but just stores BRIGHT, DEL_I and GOOD
;	debug:		/debug stops the code before the end for debugging purposes
;
;OUTPUTS:
;	All data is stored to the tree in array specific nodes in \XTOMO::TOP.RESULTS.ARRAY_X
;		\XTOMO::TOP.RESULTS.CORE:BRCHK - brightness checks [W/m^2]
;		\XTOMO::TOP.RESULTS.CORE:BRIGHT - brightness [W/m^2]
;		\XTOMO::TOP.RESULTS.CORE:EMISS - emissivity [W/m^3]
;		\XTOMO::TOP.RESULTS.CORE:RMAJ - major radius values for the inversion (time evolving) [m]
;		\XTOMO::TOP.RESULTS.CORE.CONFIG:DEL_I
;		\XTOMO::TOP.RESULTS.CORE.CONFIG:EPS
;		\XTOMO::TOP.RESULTS.CORE.CONFIG:ETA
;		\XTOMO::TOP.RESULTS.CORE.CONFIG:GOOD
;		\XTOMO::TOP.RESULTS.CORE.CONFIG:N_S
;		\XTOMO::TOP.RESULTS.CORE.CONFIG:TR
;
;PROCEDURE:
;	If /no_emiss is used, then the brightness data will be stored where the GOOD=0 channels will be forced to zero.  When the emissivity
;	is stored as well, then the brightness data has the bad channels removed.  Also, since EFIT is required for the inversion, the brightness
;	and emissivity data is only run between the first and last EFIT times and the datasets are truncated.
;
;MODIFICATION HISTORY:
;	Written by:	ML Reinke - 1/31/08
;	2/11		ML Reinke - updated to calc the m=1 sine and cosine terms
;	4/6/11		ML Reinke - updated to calc and store the m=2 cosine term
;	1/10/12		ML Reinke - added array5 for shots > 1111201001
;-

PRO xtomo_core_write2tree,shot,del_i=del_i,good=good,eps=eps,eta=eta,n_rho=n_rho,rho_vec=rho_vec,n_s=n_s,tr=tr,debug=debug,no_emiss=no_emiss,$
		asym=asym,no_cos=no_cos,no_sin=no_sin,no_m2=no_m2,no_br1=no_br1,no_br3=no_br3,no_br5=no_br5,psinorm=psinorm,dbeta1=dbeta1,dbeta3=dbeta3,etree=etree

	IF shot GT 1111201001 THEN isbr5=1 ELSE isbr5=0
	IF keyword_set(no_br5) THEN isbr5=0
	
	mdsopen,'xtomo',shot
	IF NOT keyword_set(del_i) THEN del_i=mdsvalue('\XTOMO::TOP.RESULTS.CORE.CONFIG:DEL_I')
	IF NOT keyword_set(good) THEN good=mdsvalue('\XTOMO::TOP.RESULTS.CORE.CONFIG:GOOD')
	IF NOT keyword_set(tr) THEN tr=mdsvalue('\XTOMO::TOP.RESULTS.CORE.CONFIG:TR')
	mdsclose,'xtomo',shot
	;eps,eta,n_s are taken from the inversion program
	emlogic=0
	IF NOT keyword_set(no_emiss) THEN emlogic=1
	IF keyword_set(asym) THEN emlogic=2

	IF keyword_set(no_emiss) THEN BEGIN
		calc_xtomobr_data,shot,1,bright1,ch1,time,del_i=del_i,good=good[0:37]
		calc_xtomobr_data,shot,3,bright3,ch3,time,del_i=del_i,good=good[38:75]
		IF isbr5 THEN calc_xtomobr_data,shot,5,bright5,ch5,time,del_i=del_i,good=good[76:113]
		print, 'brightness data loaded'
		IF isbr5 THEN BEGIN
			bright=[bright1,bright3,bright5]
			ch=[ch1,ch3+ch1[n(ch1)],ch5+last(ch1)+last(ch3)]
		ENDIF ELSE BEGIN
			bright=[bright1,bright3]
			ch=[ch1,ch3+ch1[n(ch1)]]	
		ENDELSE
	ENDIF ELSE BEGIN
		IF NOT keyword_set(asym) THEN BEGIN
			calc_xtomocore_emiss_data,shot,emiss,rho,time,brchk,rmaj=rmaj,bright=bright,ch=ch,del_i=del_i,good=good,array5=isbr5,$
				eps=eps,eta=eta,n_rho=n_rho,rho_vec=rho_vec,n_s=n_s,tr=tr,psinorm=psinorm,dbeta=dbeta,etree=etree
		ENDIF ELSE BEGIN
			calc_xtomocore_asymemiss_data,shot,emiss,rho,time,emr,emc1,emc2,ems1,brchk,brchkasym,rmaj=rmaj,bright=bright,ch=ch,del_i=del_i,good=good,$
				eps=eps,eta=eta,n_rho=n_rho,rho_vec=rho_vec,n_s=n_s,tr=tr,psinorm=psinorm,dbeta1=dbeta1,dbeta3=dbeta3,$
				no_cos=no_cos,no_sin=no_sin,no_br1=no_br1,no_br3=no_br3,no_m2=no_m2,etree=etree
		ENDELSE
	ENDELSE
	
	mdsopen,'xtomo',shot
	mdsput,'\XTOMO::TOP.RESULTS.CORE.CONFIG:DEL_I','build_with_units($,"")',int(del_i)
	mdsput,'\XTOMO::TOP.RESULTS.CORE.CONFIG:GOOD','build_with_units($,"")',int(good)
	mdsput,'\XTOMO::TOP.RESULTS.CORE:BRIGHT',$
        	'build_signal(build_with_units($1, "W/m^2"),*,'+$
                'build_with_units($2,"CH#"),'+$
          	'build_with_units($3,"seconds"))',bright,ch,time
	CASE emlogic OF
		0 : BEGIN ;do nothing
			
		END
		1 : BEGIN ;write 1-D profiles to the tree
			mdsput,'\XTOMO::TOP.RESULTS.CORE:EMISS',$
		       		'build_signal(build_with_units($1, "W/m^3"),*,'+$
	        	        'build_with_units($2,""),'+$
	        	        'build_with_units($3,"seconds"),'+$
        	  		'build_with_units($4,""))',emiss,rho,time,etree
			mdsput,'\XTOMO::TOP.RESULTS.CORE:BRCHK',$
        			'build_signal(build_with_units($1, "W/m^2"),*,'+$
		                'build_with_units($2,"CH#"),'+$
	        	  	'build_with_units($3,"seconds"))',brchk,ch,time
		END
		2 : BEGIN ;writed 2-D profiles to the tree
			mdsput,'\XTOMO::TOP.RESULTS.CORE:EMISS',$
		       		'build_signal(build_with_units($1, "W/m^3"),*,'+$
	        	        'build_with_units($2,""),'+$
        	  		'build_with_units($3,"seconds"),'+$
    				'build_with_units($4,"W/m^3"),'+$
				'build_with_units($5,"W/m^3"),'+$
				'build_with_units($6,"W/m^3"),'+$
				'build_with_units($7,"W/m^3"),'+$
  				'build_with_units($8,""))',emiss,rho,time,emr,emc1,emc2,ems1,etree
			mdsput,'\XTOMO::TOP.RESULTS.CORE:BRCHK',$
        			'build_signal(build_with_units($1, "W/m^2"),*,'+$
		                'build_with_units($2,"CH#"),'+$
	        	  	'build_with_units($3,"seconds"),'+$
				'build_with_units($4,"W/m^2"))',brchk,ch,time,brchkasym
		END
		ELSE : print, 'emlogic error'
	ENDCASE
	IF emlogic NE 0 THEN BEGIN	;write data independent of emlogic 
		mdsput,'\XTOMO::TOP.RESULTS.CORE:RMAJ',$
	        	'build_signal(build_with_units($1, "m"),*,'+$
		        'build_with_units($2,""),'+$
        	  	'build_with_units($3,"seconds"))',rmaj,rho,time
		mdsput,'\XTOMO::TOP.RESULTS.CORE.CONFIG:EPS','build_with_units($,"")',float(eps)
		mdsput,'\XTOMO::TOP.RESULTS.CORE.CONFIG:ETA','build_with_units($,"")',float(eta)
		mdsput,'\XTOMO::TOP.RESULTS.CORE.CONFIG:N_S','build_with_units($,"")',int(n_s)
		mdsput,'\XTOMO::TOP.RESULTS.CORE.CONFIG:TR','build_with_units($,"sec")',float(tr)	
	ENDIF
	mdsclose,'xtomo',shot
	
	IF keyword_set(debug) THEN stop
END






	