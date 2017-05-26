;+
;NAME:
;	GRID_VES
;
;PURPOSE:
;	This function returns a grid of arbitrary density that effectively
;	describes the Alcator C-Mod geometry.
;
;CALLING SEQUENCE:
;	result=GRID_VES()
;
;OPTIONAL INPUTS:
;	nx:	INT number of grids in the radial direction 	DEFAULT: 10
;	ny:	INT number of grids in the vertical direction 	DEFAULT: 10
;	rrange:	FLTARR of the radial range DEFAULT [0.425, 0.925]
;	zrange: FLTARR of the z range DEFAULT [-0.45,0.45]
;
;KEYWORD PARAMETERS:
;	center:		/center is sent to GENPOS_GRID
;	div:		/div sets the grid to parameters that describe the divertor
;			The rrange and zrange optional inputs do not work with divertor
;
;OUTPUTS:
;	Output formatted the same as GENPOS_GRID.  See it's helpfile
;
;PROCEDURE:
;	This is simply a call to GENPOS_GRID with hard-coded bounds
;	result=GENPOS_GRID(0.5,0.9,[nx,ny],xo=0.675,center=center)
;
;MODIFICATION HISTORY:
;	Written by:	ML Reinke - Summer 2006
;	1-20-07:	ML Reinke - added the zrange and rrange optional inputs
;-

FUNCTION grid_ves,nx=nx,ny=ny,center=center,div=div,zrange=zrange,rrange=rrange
	IF NOT keyword_set(nx) THEN nx=10
	IF NOT keyword_set(ny) THEN ny=10
	IF keyword_set(zrange) THEN BEGIN
		zo=mean(zrange)
		delz=zrange[1]-zrange[0]
	ENDIF ELSE BEGIN
		zo=0.0
		delz=0.9
	ENDELSE
	IF keyword_set(rrange) THEN BEGIN
		ro=mean(rrange)
		delr=rrange[1]-rrange[0]
	ENDIF ELSE BEGIN
		ro=0.675
		delr=0.5
	ENDELSE

	IF NOT keyword_set(div) THEN ves_grid=genpos_grid(delr,delz,[nx,ny],xo=ro,yo=zo,center=center) ELSE $
		ves_grid=genpos_grid(0.32,0.3,[nx,ny],xo=0.6,yo=-0.35,center=center)
	RETURN, ves_grid
END

;+
;NAME:
;	GENPOS_GPV2PSITH
;
;PURPOSE:
;	This procedure is used to visualize the voxel weighting of a multidetector GPV array
;	in (psi,theta) space.
;
;CALLING SEQUENCE:
;	GENPOS_GPV2PSITH,gpv,ves_cent
;
;INPUTS:
;	gpv:		FLTARR [n_det,n_pix] of voxel weightings generated using GENPOS_VOL_COEFS
;	ves_cent	STRUC of the voxel points created using GENPOS_GRID or GRID_VES with /center
;	
;OPTIONAL INPUTS:
;	det:		INT of the dector number to select gpv[det,*] DEFAULT: 10
;	shot:		LONG shot number for EFIT reconstructions DEFAULT: 1050426022
;	time:		FLT time point for EFIT reconstructions DEFFAULT: 1.0
;	nlevels:	INT number of levels for the contour plot DEFAULT: 20
;	cct:		INT color table for the contour plot DEFAULT: 39
;
;KEYWORD PARAMETERS:
;	ps:		/ps should be used if a postscript is to be generated to adjust plot formatting
;
;OUTPUTS:
;	All outputs are to the currently selected graphics device.
;
;MODIFICATION HISTORY:
;	Written by:	ML Reinke - Spring 2007
;
;-


PRO genpos_gpv2psith,gpv,ves_cent,pos=pos,det=det,shot=shot,time=time,nlevels=nlevels,cct=cct,ps=ps
	IF NOT keyword_set(shot) THEN shot=1050426022
	IF NOT keyword_set(time) THEN time=1.0
	IF NOT keyword_set(det) THEN det=10
	IF NOT keyword_set(nlevels) THEN nlevels=20
	IF NOT keyword_set(cct) THEN cct=39

	ves_psith=genpos_grid_rz2psith(ves_cent,shot,time)
	IF keyword_set(ps) THEN BEGIN
		title='!9(y,q)!3 of GPV Data'
		xtit='!9y!3 (norm)'
		ytit='!9q!3 (rad)'
	ENDIF ELSE BEGIN
		title='!7(w,h)!3 of GPV Data'
		xtit='!7w!3 (norm)'
		ytit='!7h!3 (rad)'
	ENDELSE
	plot,[0],[0],yr=[0,!pi],xr=[0,1],/xsty,/ysty,title=title,xtit=xtit,ytit=ytit,chars=1.5
	IF keyword_set(cct) THEN loadct,cct,/silent
	norm=max(gpv[det,*])
	nonzero=where(gpv[det,*] NE 0 AND ves_psith.psi GT 0)
	contour,reform(gpv[det,nonzero])/norm,ves_psith.psi[nonzero],ves_psith.th[nonzero],nlevels=nlevels,/irregular,/fill,/overplot,$
		levels=make(0.0,1.0,nlevels)
	loadct,12,/silent
	IF keyword_set(pos) THEN BEGIN
		psith=line_pos2psith(pos,shot,time)  
		oplot, psith.psi,psith.th,psym=8,color=200
	ENDIF

	oplot,[-1,2],[!pi/2.0,!pi/2.0],color=100,linestyle=2
	oplot,[-1,2],[3*!pi/2.0,3*!pi/2.0],color=100,linestyle=2
END

;+
;NAME:
;	GENPOS_GPV2CONTOUR
;
;PURPOSE:
;	This procedure is used to visualize a voxel weighting by making a contour plot in (R,Z) space.
;
;CALLING SEQUENCE:
;	GENPOS_GPV2CONTOUR,gpv,ves_cent
;
;INPUTS:
;;	gpv:		FLTARR [n_det,n_vox] of voxel weightings generated using GENPOS_VOL_COEFS
;	ves_cent	STRUC of the voxel points created using GENPOS_GRID or GRID_VES with /center
;	
;OPTIONAL INPUTS:
;	det:		INTARR of the dector number(s) to select gpv[det,*] DEFAULT: 0
;	shot:		LONG shot number for EFIT reconstructions DEFAULT: 1050426022
;	time:		FLT time point for EFIT reconstructions DEFFAULT: 1.0
;	nlevels:	INT number of levels for the contour plot DEFAULT: 20
;	cct:		INT color table for the contour plot DEFAULT: 39
;	pos:		FLTARR [Ro,Zo,Rt,Psi] of the detector(s) that is being visualized.  If included a line of sight 
;			will be plotted over the volume contour plot. For multiple detectors insert entire 2D POS array
;	win:		INT of the window number to send plot.  DEFAULT: see VESSEL_PLOT in MLR_FUNCTIONS
;	
;KEYWORD PARAMETERS:
;	invessel:	/invessel will prevent the program from plotting points that are outside the PFC boundary
;			by calling GENPOS_GRID_INVESSEL
;	ps:		/ps makes formatting changes and should be used when plotting to a postscript
;	edge:		/edge is sent to VESSEL_PLOT
;	div:		/div is sent to VESSEL_PLOT
;	fsplot:		/fsplot will plot the flux surfaces for the shot and time using FS_PLOT in MLR_FUNCTIONS
;	
;OUTPUTS:
;	All outputs are made to the currently selected graphics device.
;
;RESTRICTIONS:
;	GENPOS_GPV2CONTOUR was intially intended to check the contour plots of the output of GENPOS_PLANAR2GPV.  This means the
;	data input likes to be formatted in a [n_det,n_vox].  The output of GENPOS_VOL_COEFS will be a 1D array of length n_vox.
;	To use the contour plotting, make sure the optional input det=0 is used.
;	
;MODIFICATION HISTORY:
;	Written by: 	ML Reinke - Spring 2007
;       6-20-07:	ML Reinke - added the ability to plot multiple detectors  on a single cross-secion
;	7-14-07:	ML Reinke - added the ability to plot multiple POS lines when multiple detectors are plotted.
;
;-

PRO genpos_gpv2contour,gpv,ves_cent,det=det,shot=shot,time=time,nlevels=nlevels,cct=cct,pos=pos,invessel=invessel,win=win,ps=ps,$
		edge=edge,div=div,fsplot=fsplot,points=points,pntcol=pntcol,fscol=fscol
	
	IF NOT keyword_set(cct) THEN cct=39
        IF NOT keyword_set(pntcol) THEN pntcol=255
	tmp_gpv=gpv
	tmp_ves=ves_cent
	IF keyword_set(invessel) THEN BEGIN
		tmp_in=genpos_grid_invessel(ves_cent)
		in_ves=where(tmp_in EQ 1)
                ;stop
	ENDIF ELSE BEGIN
        	in_ves=lindgen(n(ves_cent.pnts[0,*])+1)
                tmp_in=fltarr(n(in_ves)+1)+1
        ENDELSE
	x=size(gpv)
	IF x[0] EQ 1 THEN gpv=transpose([[gpv],[gpv]])
	loadct,12,/silent
	IF NOT keyword_set(nlevels) THEN nlevels=50
	IF NOT keyword_set(shot) THEN shot=1050426022
	IF NOT keyword_set(time) THEN time=1.0
	IF NOT keyword_set(det) THEN det=0.0

	;if plots are set, get wall and limiter traces
	mdsopen, 'analysis',-1
	IF shot GT 1020000000 THEN  node_year='2002' else node_year='1994'
	r_lim=mdsvalue('.LIMITERS.WALL.WALL_'+node_year+':RLIM')
	z_lim=mdsvalue('.LIMITERS.WALL.WALL_'+node_year+':ZLIM')
	r_rf=mdsvalue('.LIMITERS.RF_LIMITER:R')
	z_rf=mdsvalue('.LIMITERS.RF_LIMITER:Z')
	mdsclose
	
	;plot the limiting structures
	;wreset,19
	;plot, r_lim,z_lim,title='SHOT: '+num2str(shot,1)+' TIME: '+num2str(time,dp=3),chars=1.3
	;oplot, r_rf,z_rf
	vessel_plot,title='SHOT: '+num2str(shot,1)+' TIME: '+num2str(time,dp=3),shot=shot,n=win,ps=ps,d_old=d_old,div=div,edge=edge
	
		
	;plot the contours
        FOR i=0,n(det) DO BEGIN
        	norm=max(gpv[det[i],in_ves])
                loadct,cct,/silent
                nonzero=where(gpv[det[i],*] NE 0)
                ;nonzero=where(gpv[det[i],*] NE 0 AND tmp_in EQ 1)
                contour,reform(gpv[det[i],nonzero])/norm,ves_cent.pnts[0,nonzero],ves_cent.pnts[1,nonzero],nlevels=nlevels,$
			/irregular,/fill,/overplot,levels=make(0.0,1.0,nlevels)
                loadct,12,/silent
                IF keyword_set(points) THEN oplot, ves_cent.pnts[0,nonzero],ves_cent.pnts[1,nonzero],color=pntcol,psym=3
        ENDFOR

	;plot the LCFS
;	efit_lcfs=line_getlcfs(shot)
;	efit_time=line_gettimes(shot)
;	efit_i=ipt(efit_time,time)
;	rbdry=efit_lcfs[0,where(efit_lcfs[0,*,efit_i] NE 0 AND efit_lcfs[1,*,efit_i] NE 0),efit_i]
;	zbdry=efit_lcfs[1,where(efit_lcfs[0,*,efit_i] NE 0 AND efit_lcfs[1,*,efit_i] NE 0),efit_i]
;	rbdry=rbdry[0:n(rbdry)-1]
;	zbdry=zbdry[0:n(zbdry)-1]
;	rb_plt=[rbdry,rbdry[0]]
;	zb_plt=[zbdry,zbdry[0]]
;	oplot,rb_plt,zb_plt,color=200,thick=2.5

	;plot the line path
	IF keyword_set(pos) THEN BEGIN
            	FOR i=0,n(det) DO BEGIN
			IF pos[2,det[i]] LT min(ves_cent.pnts[0,*]) THEN l=make(0,1,100) ELSE l=make(0,2,100)
			r=sqrt(pos[0,det[i]]^2+l*(l-2)*(pos[0,det[i]]^2-pos[2,det[i]]^2))
			z=pos[1,det[i]]-l*sqrt(pos[0,det[i]]^2-pos[2,det[i]]^2)*tan(pos[3,det[i]])
			oplot, r, z, color=30,thick=2.0,linestyle=2
                ENDFOR
	ENDIF
	vessel_plot,/oplot,shot=shot,ps=ps
	IF keyword_set(fsplot) THEN fs_plot,shot,time,color=fscol
	gpv=tmp_gpv
	ves_cent=tmp_ves
	IF keyword_set(ps) THEN device, xsize=float(d_old.x_size)/d_old.x_px_cm,ysize=float(d_old.y_size)/d_old.y_px_cm
END


FUNCTION grid_profile,gridpts,func,r,t,shot,cfunc=cfunc,sfunc=sfunc,tpts=tpts,verb=verb,sol=sol,rho=rho,debug=debug

	IF NOT keyword_set(tpts) THEN tpts=t
	rpts=gridpts.pnts[0,*]
	zpts=gridpts.pnts[1,*]
        
		
	rmid=efit_rz2rmid(rpts,zpts,tpts,shot=shot,rho=rho)
        efit_rmid=line_getrmid(shot)
        efit_time=line_gettimes(shot)
        efit_axis=line_getaxis(shot)
        
	;rmid[*,i] are the rmid values for all the (r,z) points for tpts[i]
	;rmid[i,*] are all rmid values for (r,z)_i for all time points

	num_r=n(rpts)+1
	num_t=n(tpts)+1

	IF n(t) GT 0 THEN BEGIN
		func_int=fltarr(num_r,num_t)
                cfunc_int=fltarr(num_r,num_t) 	;default = 0.0
                sfunc_int=fltarr(num_r,num_t) 	;default = 0.0
                costh=fltarr(num_r,num_t)
                sinth=fltarr(num_r,num_t)
                func_grid=fltarr(num_r,num_t)
		FOR i=0,num_t-1 DO BEGIN
                        index=ipt(efit_time,tpts[i])
                        ro=efit_axis[index,0]
                        zo=efit_axis[index,1]
                        rad=sqrt((rpts-ro)^2+(zpts-zo)^2)
                        costh[*,i]=(rpts-ro)/rad^2
                        sinth[*,i]=(zpts-zo)/rad^2
			r_interp=interp_vec_reform(r,rmid[*,i])
			t_interp=interp_vec_reform(t,fltarr(num_r)+tpts[i])
				IF keyword_set(verb) THEN BEGIN
					print, 'Interpolation for Temperature'
					print, 'R_INTERP'
					print, r_interp
					print, 'T_INTERP'
					print, t_interp
				ENDIF
			func_int[*,i]=interpolate(func,r_interp,t_interp, missing=1.0e-4)
                        IF keyword_set(cfunc) THEN cfunc_int[*,i]=interpolate(cfunc,r_interp,t_interp, missing=1.0e-4)
                        IF keyword_set(sfunc) THEN sfunc_int[*,i]=interpolate(sfunc,r_interp,t_interp, missing=1.0e-4)

                        IF NOT keyword_set(sol) THEN BEGIN
                        	efit_i=ipt(efit_time,tpts[i])
                        	tmp=where(rmid GT max(efit_rmid[efit_i,*]))
	                        IF tmp[0] NE -1 THEN func_interp[tmp,i]=0
                        ENDIF
                        func_grid[*,i]=func_int[*,i]+cfunc_int[*,i]*costh[*,i]+sfunc_int[*,i]*sinth[*,i]
		ENDFOR
	ENDIF ELSE BEGIN
            	func_interp=interpol(func,r, rmid)
                index=ipt(efit_time,t)
                ro=efit_axis[index,0]
                zo=efit_axis[index,1]
                rad=sqrt((rpts-ro)^2+(zpts-zo)^2)
                costh=(rpts-ro)/rad
                sinth=(zpts-zo)/rad
           
        	func_int=interpol(func,r, rmid)
                IF keyword_set(cfunc) THEN cfunc_int=interpol(cfunc,r,rmid) ELSE cfunc_int=fltarr(num_r)
                IF keyword_set(sfunc) THEN sfunc_int=interpol(sfunc,r,rmid) ELSE sfunc_int=fltarr(num_r)

               	IF NOT keyword_set(sol) THEN BEGIN
                	efit_i=ipt(efit_time,tpts[0])
                        tmp=where(rmid GT max(efit_rmid[efit_i,*]))
                	IF tmp[0] NE -1 THEN func_interp[tmp]=0
                ENDIF
                func_grid=func_int+cfunc_int*costh+sfunc_int*sinth
       ENDELSE

        IF keyword_set(debug) THEN stop
	RETURN, func_grid
END

;+
;NAME:
;	GENPOS_POS2RMIDTANG
;
;PURPOSE:
;	This function turns an array of pos vectors into an array of RMID values that correspond to the
;	flux surfaces of closest approach.
;
;CALLING SEQUENCE:
;	result=GENPOS_POS2RMIDTANG(pos,shot,time)
;
;INPUTS:
;	pos:	FLTARR 	[4,n_pos] of [Ro,Zo,Rt,PSI] x m position vectors
;	shot:	LONG	shot number for EFIT reference
;	time	FLTARR 	[n_time] of arbitrary length of the time points on which to perform operation
;	
;OPTIONAL INPUTS:
;	n_s:	INT	number of points to subdivide line before call to EFIT_RZ2RMID, DEFAULT: 40
;	rzbnd		FLTARR 	[r_min,r_max,|z|] of the rectangular bounding box for the lines DEFAULT: [0.44,1.0,0.45]
;	tree:		STRING	of the EFIT tree to use for calculating weighting matrix [DEFAULT: 'analysis']\
;
;KEYWORD PARAMETERS:
;	rho:	/rho is carred through to EFIT_RZ2RMID to output rho values instead of rmid
;	psin:	/psin runs EFIT_RZ2RHO instead of EFIT_RZ2RMID
;	efit:	/efit runs RZ2RMID on the EFIT time base instead of the time input.  The output data interpolates this result onto
;			the time base.  This is useful for diode data that has better time resolution than EFIT
;
;OUTPUTS:
;	result:	FLTARR 	[n_pos, n_time] of outboard major radii [m] of the flux surface of closest approach
;			for each POS vector (or r/a if /rho used or normalized poloidal flux if psinorm is used
;
;PROCEDURE:
;	The line of sight is subdivided into n_s parts and LINE_R and LINE_Z are used to find (R,Z) points
;	in the chamber.  One call, for all channels at all times, is sent to EFIT_RZ2RMID to save on computation.
;
;
;MODIFICATION HISTORY:
;	Written by:	ML Reinke - 7/22/07
;	1-17-08:	ML Reinke - added the keyword parameter EFIT which runs EFIT_RZ2RMID on the
;                                   EFIT time base and then interpolates results to the input time base
;	1-17-08:	ML Reinke - updated the s_min,s_max finder to include vertical bounderies and
;                                   added rzbnd optional input
;	1-30-09:	ML Reinke - fixed a bug by limiting the minimum RMID to those |z_pos| < rzbnd[2]
;	9-23-10:	ML Reinke - added the psinorm keyword 
;	1-05-11:	ML Reinke - if user doesn't confirm that time is inside valid EFITs, then RZ2XXX
;                                   will return NANs and crash.  Fixed it so rhotang = 0 will be returned
;                                   instead in these cases
;	6-13-11:	ML Reinke - added the tree optional input for use with EFIT_RZ2XXX codes.
;
;-
	
FUNCTION genpos_pos2rmidtang,pos,shot,time,rho=rho,n_s=n_s,efit=efit,rzbnd=rzbnd,debug=debug,psin=psin,tree=tree
	IF NOT keyword_set(n_s) THEN n_s=40
        IF NOT keyword_set(rzbnd) THEN rzbnd=[0.44,1.0,0.45] ;set path's [r_min,r_max,abs(z_max)] limit
        IF keyword_set(efit) THEN BEGIN
        	tau=time		 ;backup time 
                n_tau=n(time)+1
                time=line_gettimes(shot) ;run RZ2RMID at EFIT times
        ENDIF
	n_ch=n(pos[0,*])+1
	r_pos=fltarr(n_s*n_ch)
	z_pos=fltarr(n_s*n_ch)
	FOR i=0L,n_ch-1 DO BEGIN
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
		IF keyword_set(rz_ap) THEN s_ap=line_s(pos[*,i],z=rz_ap[1]) > 0 ELSE s_ap=0			
		s=make(s_ap,smin,n_s)	;n_s points going from the aperture to the boundary intersecion
		r_pos[i*n_s:(i+1)*n_s-1]=line_r(s,pos[*,i])
		z_pos[i*n_s:(i+1)*n_s-1]=line_z(s,pos[*,i])
        ENDFOR
        IF keyword_set(debug) THEN stop
	IF NOT keyword_set(psin) THEN rmid=efit_rz2rmid(r_pos,z_pos,time,shot=shot,rho=rho,tree=tree) ELSE rmid=efit_rz2rho(r_pos,z_pos,time,shot=shot,tree=tree)
        
	rmid_min=fltarr(n_ch,n(time)+1)
	FOR i=0L,n_ch-1 DO BEGIN
		FOR j=0,n(time) DO BEGIN
                    	check=where(finite(rmid[i*n_s:(i+1)*n_s-1,j]) EQ 1)
                        IF n(check)+1 EQ n_s THEN BEGIN
                        	rmid_ch=rmid[i*n_s:(i+1)*n_s-1,j]
                                tmp=where(rmid_ch GT 0 AND abs(z_pos[i*n_s:(i+1)*n_s-1]) LT rzbnd[2])
                                rmid_min[i,j]=min(rmid_ch[tmp])
                        ENDIF
		ENDFOR
	ENDFOR
	IF keyword_set(efit) THEN BEGIN
            	rmid_min_int=fltarr(n_ch,n_tau)
                FOR i=0,n_ch-1 DO rmid_min_int[i,*]=interpol(rmid_min[i,*],time,tau)
        	time=tau
                rmid_min=rmid_min_int
        ENDIF

	output=rmid_min
        IF keyword_set(debug) THEN stop
	RETURN,output
END
	
;+
;NAME:
;	GENPOS_RZ2PSITH
;
;PURPOSE:
;	This procedure uses EFIT data to covert a series of (r,z) points to (psi,th) points
;	where psi is normalized
;
;CALLING SEQUENCE:
;	GENPOS_RZ2PSITH,r,z,shot,time,psi,th
;
;INPUTS:
;	r:	FLTARR	of major radius points [meters] of length n_r
;	z:	FLTARR	of height points [meters] (must be same length as r)
;	shot:	LONG	shot number
;	time:	FLTARR	of time points at which to do computation [sec] of length n_t
;
;OUTPUTS
;	psi:	FLTARR [n_r,n_t] of normalized psi values
;	th:	FLTARR [n_r,n_z] of theta values
;
;PROCEDURE:
;	Normalized psi values are computed using EFIT_RZ2RHO with /psinorm optional input
;	Theta values are computed by taking a tangent of (z-zo)/(r-ro) which is not a strict
;	defination of theta, but it's what I'm working with for now.
;
;MODIFICATION HISTORY:
;	Written by:	ML Reinke 12/13/07
;
;-

PRO genpos_rz2psith,r,z,shot,time,psi,th,debug=debug
	psi=-1
        th=-1

	n_r=n(r)+1
        n_z=n(z)+1
        n_time=n(time)+1
        IF n_r NE n_z THEN RETURN
        
        ;initiate array
        th=fltarr(n_r,n_time)
        ;load EFIT axis data
        efit_t=line_gettimes(shot)
        axis=line_getaxis(shot)
        ;fill th values
        FOR i=0,n_time-1 DO BEGIN
        	FOR j=0,n_r-1 DO BEGIN
			pt=[r[j],z[j]]
			pt-=axis[ipt(efit_t,time[i]),*]
			ang=atan(pt[1]/pt[0])
			IF pt[1] GT 0 AND pt[0] LT 0 THEN ang=ang+!pi
			IF pt[1] LT 0 AND pt[0] LT 0 THEN ang=ang+!pi
			IF pt[1] LT 0 AND pt[0] GT 0 THEN ang=ang+2.0*!pi
			th[j,i]=ang
                ENDFOR
	ENDFOR

        ;fill psinorm data
	psi=efit_rz2rho(r,z,time,shot=shot,/psinorm)
	 
	IF keyword_set(debug) THEN stop
END

;+
;NAME:
;	GENPOS_POS_REFORM
;	
;PURPOSE:
;	This function reformulates POS vectors so that their (R,Z) starting points are on a perscribed boundary
;	
;CALLING SEQUENCE:
;	GENPOS_POS_REFORM,pos,rzbnd
;
;INPUTS:
;	pos		FLTARR	[4,n_pos] of the input POS vector
;	rzbnd:		FLTARR	[4] [R_in,R_out,Z_top,Z_bot] in [m] of the rectangular boundary
;
;OUTPUTS:
;	pos:		FLTARR	[4,n_pos] of the reformulated POS vector (copies over input)
;	
;OPTIONAL OUTPUTS:
;	old_pos:	FLTARR	[4,n_pos] a copy of the input POS vectors
;
;MODIFICATION HISTORY:
;	Written by:	ML Reinke - 1/30/08
;
;-

PRO genpos_pos_reform,pos,rzbnd,old_pos=old_pos
	
	n_pos=n(pos[0,*])+1
	IF n(rzbnd) NE 3 THEN RETURN
	new_pos=fltarr(4,n_pos)
	FOR i=0,n_pos-1 DO BEGIN
		;determine where the line of sight terminates in s due to radial boundries
		svec=line_s(pos[*,i],r=rzbnd[1])
		smin_r=svec[1]	 ;choose the smallest  where r=rmax

		;determine where the line of sight terminates in s due to z boundries
		IF pos[3,i] NE 0 THEN BEGIN
			s_top=line_s(pos[*,i],z=rzbnd[2])
			s_bot=line_s(pos[*,i],z=rzbnd[3])
		ENDIF ELSE smin_z=2.0*smin_r	;if no inclination, smin_z set > smin_r	
		s_bnd=[smin_r,s_top,s_bot]
		s=min(s_bnd[where(s_bnd GT 0)])
		r=line_r(s,pos[*,i])
		z=line_z(s,pos[*,i])
		new_pos[*,i]=[r,z,pos[2,i],pos[3,i]]
	ENDFOR
	old_pos=pos
	pos=new_pos
END

;+
;NAME:
;	GENPOS_POLOIDALPOS2RT
;
;PURPOSE:
;	This procedure caclulates tangency radii for poloial views, assuming that the flux surfaces
;	are circles with ROUT, the center of the LCFS as the center.
;
;CALLING SEQUENCE:
;	GENPOS_POLOIDALPOS2RT,shot,pos,rt,time
;
;INPUTS:
;	shot:		LONG 	shot number
;	pos:		FLTARR	[4,n_ch] of the view information
;	
;OUTPUTS:
;	rt:		FLTARR	[n_ch,n_time] of the tangency radii [m]
;	time:		FLTARR	[n_time] of the time points (EFIT TIME GRID) [s]	
;
;OPTIONAL OUTPUTS:
;	rout:		FLTARR	[n_time] of the center of the LCFS, assumed to be radial center of circular flux surfaces [m]
;	zmag:		FLTARR	[n_time] height of the magnetic axis, assumed to be vertical center of circular flux surfaces [m]
;
;MODIFICATION HISTORY:
;	Written by:	ML Reinke - 2/20/08
;
;-

PRO genpos_poloidalpos2rt,shot,pos,rt,time,debug=debug,rout=rout,zmag=zmag
	mdsopen,'analysis',shot
	rout=mdsvalue('\ANALYSIS::TOP:EFIT.RESULTS.A_EQDSK:ROUT')/100.0
	zmag=mdsvalue('\ANALYSIS::TOP:EFIT.RESULTS.A_EQDSK:ZMAGX')/100.0
	time=mdsvalue('dim_of(\ANALYSIS::TOP:EFIT.RESULTS.A_EQDSK:ZMAGX)')
	mdsclose,'analysis',shot
	n_ch=n(pos[0,*])+1
	n_time=n(time)+1

	rt=fltarr(n_ch,n_time)
	FOR i=0,n_ch-1 DO BEGIN
		delr=pos[0,i]-rout
		delz=pos[1,i]-zmag
		rt[i,*]=abs(sin(pos[3,i])*delr-delz*cos(pos[3,i]))
	ENDFOR
	IF keyword_set(debug) THEN stop
END
