; This contains various line integration functions for the genie code
; The main function in this code is LINE_BR

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;line integration tools 
;	see ML Reinke Personal
;	Logbook - 2 pg 10-12
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

;The functionality of LINE_MVEC, LINE_Z, LINE_R and LINE_S
;have been replaced with better functions that use the
;analytic expresstions derived in 'General Equations for 
;Pinhole Radiometry in Tokamaks' by Matt Reinke. (Note: title
;may have changed slightly).  See the file genie_line.bak
;for the older versions of these functions.


;compute height, z, at position,s, along line defined by POS
FUNCTION line_z,s,pos
	;assume pos=[r1,z1,rt,psi]
	output=pos[1]-s*sqrt(pos[0]^2-pos[2]^2)*tan(pos[3])
	RETURN, output
END

;compute major radius,r, at position,s, along line defined by POS
FUNCTION line_r,s,pos
	;assume pos=[r1,z1,rt,psi]
	output=sqrt(pos[0]^2+s*(s-2.0)*(pos[0]^2-pos[2]^2))
	RETURN, output
END

;comput angle, th, at position s, assuming th1=0 or th1 given
FUNCTION line_th,s,pos,th1=th1,cont=cont,neg=neg
	IF NOT keyword_set(th1) THEN th1=0.0
	IF keyword_set(neg) THEN ff=-1.0 ELSE ff=1.0
	th2=th1+ff*acos(pos[2]/pos[0])
	IF pos[2] EQ 0 THEN th2=th1
	a=(sin(th1)+s*(pos[2]/pos[0]*sin(th2)-sin(th1)))
	b=(cos(th1)+s*(pos[2]/pos[0]*cos(th2)-cos(th1)))
	output=atan(a/b)
	IF keyword_set(cont) THEN BEGIN
		tmp=where(output LT 0)
		IF tmp[0] NE -1 THEN output[tmp]+=!pi
	ENDIF
	RETURN,output
END

;computes the angle between toroidal angle
;and the line of sight projected into a constant z plane
;which is needed for rotation measurements.
FUNCTION line_torphi,s,pos
	phi=acos(pos[2]/line_r(s,pos))
	RETURN,phi
END

;scaling factor to turn parameter into physical line length
;the scaling of l is arbitrary with no physical signifigance except that
;l=0 is pt1 and l=1 is pt2 if the l values are multiplied by the distance
;between 1 and 2 then the qunatitative units are applied to the scaling.

FUNCTION line_scale,pos
	scale=fltarr(n(pos[0,*])+1)
	FOR i=0,n(scale) DO BEGIN
		ipos=reform(pos[*,i])
		scale[i]=sqrt((ipos[0]^2-ipos[2]^2)*(1.0+(tan(ipos[3]))^2))	
	ENDFOR
	RETURN,scale
END

;compute position,s, along line defined by POS relative to (ro, zo) in POS
;with s > 0 going inward and s < 0 going outward along line
FUNCTION line_s,pos,r=r,z=z, mvec=mvec
	
	IF NOT keyword_set(r) AND NOT keyword_set(z) THEN RETURN, 0
	IF keyword_set(r) AND keyword_set(z) THEN RETURN,0
	IF keyword_set(r) THEN IF r LT pos[2] THEN RETURN, -1 	;else will result in sqrt(<0)	

	IF keyword_set(r) THEN BEGIN
		s1=1+sqrt((r^2-pos[2]^2)/(pos[0]^2-pos[2]^2))
		s2=1-sqrt((r^2-pos[2]^2)/(pos[0]^2-pos[2]^2))
		s=[s1,s2]
	ENDIF
	
	IF keyword_set(z) THEN BEGIN
		s=(pos[1]-z)/(tan(pos[3])*sqrt(pos[0]^2-pos[2]^2))
	ENDIF

	output=s
	RETURN, output
END

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;last closed flux surface tools 
;	see ML Reinke Personal
;	Logbook - 2 pg 15
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

;+
;NAME:
;	LINE_VESSEL
;
;PURPOSE:
;	This function can be used to load (R,Z) points of the C-Mod wall and limiter structure from the
;	tree.  This function can also be used to create a more finely interpolated map of the outer wall
;	that (R,Z) grid points can be determined to be in/out of vessel by using LINE_INLCFS.
;
;CALLING SEQUENCE:
;	result=LINE_VESSEL()
;	
;OPTIONAL INPUTS:
;	shot:		LONINT shot number if prior to January 2002 vessel is desired
;	del:		FLT distance [meters] between grid points when /new is used DEFAULT = 0.01
;	path:		STR of path to file when /load is used DEFAULT: /home/username/idl/genie/data/line_vessel.dat
;
;KEYWORD PARAMETERS:
;	new:		/new generates a new set of (R,Z) points by interpolating the points from the tree down
;				to an interval set by del.  Small values of del are needed for LINE_INLCFS to work.
;	load:		/load saves on computation time by restoring the saveset pointed to by path.  This IDL saveset
;				must containt the variables r_wall, z_wall, r_lim and z_lim.  Running LINE_VESSEL with
;				a /new and /debug then saving these values is the best way to generate the saveset
;	debug:		/debug stops before the RETURN statement but after the output structure has been assembled.
;
;OUTPUTS:
;	result:		STRUC of containg (R,Z) points [meters] that define the outer limiter and wall of C-Mod
;				*.wall	- inner/outer wall points and divertor
;				*.*.r	- points in major radius [meters]
;				*.*.z	- points in z [meters]
;				*.lim	- outer (RF) limiter - not presently interpolated like the wall points
;				*.*.r	- points in major radius [meters]
;				*.*.z	- points in z [meters]
;
;PROCEDURE:
;	This function loads the data from the nodes in the ANALYSIS tree:
;	.LIMITERS.WALL.WALL_XXXX:RLIM		node year (XXXX) either 2002 or 1994
;	.LIMITERS.WALL.WALL_XXXX:ZLIM		with default being 2002
;	.LIMITERS.RF_LIMITER:R
;	.LIMITERS.RF_LIMITER:Z
;
;
;
;MODIFICATION HISTORY:
;	Written by:	ML Reinke: 8-04-06
;
;-


FUNCTION line_vessel,npts=npts,new=new,load=load,shot=shot,del=del,path=path,debug=debug
	IF NOT keyword_set(del) THEN del=0.01
	IF NOT keyword_set(shot) THEN shot=1020000001
	IF NOT keyword_set(path) THEN path='/home/'+logname()+'/idl/genie/data/line_vessel.dat'

	IF keyword_set(new) THEN BEGIN
			mdsopen, 'analysis',-1
			IF shot LT 1020000000 THEN node_year='1994'
			IF shot GT 1020000000 AND shot LT 1070300000 THEN node_year='2002'
			IF shot GT 1070300000 THEN node_year='2007'
			;IF shot GT 1090600000 THEN node_year='2009'
			r_wall=mdsvalue('.LIMITERS.WALL.WALL_'+node_year+':RLIM')
			z_wall=mdsvalue('.LIMITERS.WALL.WALL_'+node_year+':ZLIM')
			r_lim=mdsvalue('.LIMITERS.RF_LIMITER:R')
			z_lim=mdsvalue('.LIMITERS.RF_LIMITER:Z')
			mdsclose
			num_wall=n(r_wall)
			;fix some wall points (hardcoded to 2002 wall)
			IF node_year EQ '2002' THEN BEGIN
				z_wall[3]+=0.05
				r_wall[num_wall-5]=r_wall[num_wall-4]
				z_wall[num_wall-5]=z_wall[num_wall-4]-0.005
			ENDIF
			wall=fltarr(2,10000)
			cntr=0
			n_wall=n(r_wall)+1
			FOR i=0,n_wall-2 DO BEGIN
				length=dist([r_wall[i],z_wall[i]],[r_wall[i+1],z_wall[i+1]])
				IF length GT 2.0*del THEN BEGIN
					num=floor(length/del)
					l=make(0.0,1.0-del/length,num)
					rpts=r_wall[i]+l*(r_wall[i+1]-r_wall[i])
					zpts=z_wall[i]+l*(z_wall[i+1]-z_wall[i])
					wall[0,cntr:cntr+num-1]=rpts
					wall[1,cntr:cntr+num-1]=zpts
					cntr+=num					
				ENDIF ELSE BEGIN
					wall[*,cntr]=[r_wall[i],z_wall[i]]
					cntr+=1
				ENDELSE
			ENDFOR
			tmp=where(wall[0,*] NE 0)
			r_wall=reform(wall[0,tmp])
			z_wall=reform(wall[1,tmp])
		
	ENDIF ELSE BEGIN
		IF keyword_set(load) THEN BEGIN
			restore,path
		ENDIF ELSE BEGIN
			mdsopen, 'analysis',-1
			IF shot LT 1020000000 THEN node_year='1994'
			IF shot GT 1020000000 AND shot LT 1070300000 THEN node_year='2002'
			IF shot GT 1070300000 THEN node_year='2007'
			r_wall=mdsvalue('.LIMITERS.WALL.WALL_'+node_year+':RLIM')
			z_wall=mdsvalue('.LIMITERS.WALL.WALL_'+node_year+':ZLIM')
			r_lim=mdsvalue('.LIMITERS.RF_LIMITER:R')
			z_lim=mdsvalue('.LIMITERS.RF_LIMITER:Z')
			mdsclose
		ENDELSE
	ENDELSE
	
	wall={r:r_wall,z:z_wall}
	lim={r:r_lim, z:z_lim}
	output=create_struct('wall',wall,'lim',lim)
	IF keyword_set(debug) THEN stop
	RETURN,output
	
END


;find the times for which EFIT has been run
FUNCTION line_gettimes,shot
	mdsopen, 'mhd', (shot)
	t=mdsvalue('dim_of(\MHD::TOP.ANALYSIS.EFIT.RESULTS.G_EQDSK:RMAXIS)',$
		/quiet,status=status)
	IF status THEN output=t ELSE output=-1
	mdsclose, 'mhd', (shot)
	RETURN, t
END

FUNCTION line_getrmid,shot,psin=psin
	mdsopen, 'analysis', shot
	rmid=mdsvalue('\efit_rmid')
	psin=mdsvalue('dim_of(\efit_rmid,1)')
	mdsclose, 'analysis', shot
	output=rmid
	RETURN, output
END


;get EFIT output array that defines the LCFS
;[2,n_sp, n_time], so plot,[0,*,i],[1,*,i] will show you the
;LCFS it time point, i (reference t using line_gettimes)
FUNCTION line_getlcfs,shot
	mdsopen, 'mhd', (shot)
	bdry=mdsvalue('\MHD::TOP.ANALYSIS.EFIT.RESULTS.G_EQDSK:BDRY',$
		/quiet,status=status)
	IF NOT status THEN BEGIN
		rbbs=mdsvalue('\ANALYSIS::TOP:EFIT.RESULTS.G_EQDSK:RBBBS',/quiet, status=rstat)
		zbbs=mdsvalue('\ANALYSIS::TOP:EFIT.RESULTS.G_EQDSK:ZBBBS',/quiet,status=zstat)
		IF rstat AND zstat THEN BEGIN
			nbbs=n(rbbs[0,*])+1
			ntime=n(rbbs[*,0])+1
			bdry=fltarr(2,nbbs,ntime)
			FOR i=0,ntime-1 DO BEGIN
				bdry[0,*,i] = rbbs[i,*]
				bdry[1,*,i] = zbbs[i,*]
			ENDFOR
			output=bdry
		ENDIF ELSE output = -1			
	ENDIF ELSE output=bdry
	mdsclose,'mhd', (shot)
	RETURN, output
END

;get the [r,z] of the magnetic axis from EFIT, indices are time
FUNCTION line_getaxis,shot
	mdsopen, 'mhd', (shot)
	raxis=mdsvalue('\MHD::TOP.ANALYSIS.EFIT.RESULTS.G_EQDSK:RMAXIS',$
		/quiet,status=rstatus)
	zaxis=mdsvalue('\MHD::TOP.ANALYSIS.EFIT.RESULTS.G_EQDSK:ZMAXIS',$
		/quiet, status=zstatus)
	mdsclose, 'mhd', (shot)
	IF rstatus AND zstatus THEN output=[[raxis], [zaxis]] ELSE output = -1
	RETURN, output
END

;1/24/11 - MLR - added /psinorm capability 
FUNCTION line_getrmaj,shot,rho,time,psinorm=psinorm,debug=debug
	axis=line_getaxis(shot)
	rmid=line_getrmid(shot,psin=psin)
	efit_times=line_gettimes(shot)
	r=fltarr(n(rho)+1,n(time)+1)
	FOR i=0L,n(time) DO BEGIN
		efit_i=ipt(efit_times,time[i])
		IF NOT keyword_set(psinorm) THEN BEGIN
			Ro=axis[efit_i,0]
			a=rmid[efit_i,n(rmid[0,*])]-rmid[efit_i,0]
			r[*,i]=rho*a+Ro		
		ENDIF ELSE BEGIN
			r[*,i]=interpol(rmid[efit_i,*],psin,rho)
		ENDELSE
	ENDFOR
	output=r
	IF keyword_set(debug) THEN stop
	RETURN,output
END

FUNCTION line_getrho,shot,rmaj,time,debug=debug
	axis=line_getaxis(shot)
	rmid=line_getrmid(shot)
	efit_times=line_gettimes(shot)
	rho=fltarr(n(rmaj)+1,n(time)+1)
	FOR i=0,n(time) DO BEGIN
		efit_i=ipt(efit_times,time[i])
		Ro=axis[efit_i,0]
		a=rmid[efit_i,n(rmid[0,*])]-rmid[efit_i,0]
		rho[*,i]=(rmaj-Ro)/a		
	ENDFOR
	output=rho
	IF keyword_set(debug) THEN stop
	RETURN,output
END
	
;given a shot, time, and outboard midplane radial point return the r,z points for the
;flux surface that intersects the closest rmid at the point, rpt
FUNCTION line_getfs, shot, tpt, rpt,debug=debug
	efit_times=line_gettimes(shot)
	efit_rmid=line_getrmid(shot)
	time_i=ipt(efit_times, tpt)
	rad_i=ipt(efit_rmid[time_i,*], rpt)
	IF keyword_set(debug) THEn print, time_i, rad_i
	IF time_i EQ -1 OR rad_i EQ -1 THEN RETURN, -1
	efit_psicont,psi_r,psi_z,psi_n,psi_t,shot=shot,time=efit_times[time_i]
	out_r=reform(psi_r[*,rad_i])
	out_z=reform(psi_z[*,rad_i])
	tmp=where(out_r NE 0)
	output={r:out_r[tmp], z:out_z[tmp]}
	IF keyword_set(debug) THEN stop
	RETURN, output
END

;this function works in conjuction with line_inlcfs to generate
;the correct format from the EFIT last closed flux surface

FUNCTION line_genbdry,shot,time

	;load EFIT data
	efit_time=line_gettimes(shot)
	efit_lcfs=line_getlcfs(shot)
	efit_axis=line_getaxis(shot)

	efit_i=ipt(efit_time,time)
	efit_i=efit_i[0]

	;load relevant LCFS boundary and prepare for format needed for LINE_INCLFS
	rbdry=efit_lcfs[0,where(efit_lcfs[0,*,efit_i] NE 0 AND efit_lcfs[1,*,efit_i] NE 0),efit_i]
	zbdry=efit_lcfs[1,where(efit_lcfs[0,*,efit_i] NE 0 AND efit_lcfs[1,*,efit_i] NE 0),efit_i]
	rbdry=rbdry[0:n(rbdry)-1]
	zbdry=zbdry[0:n(zbdry)-1]
	axis=[efit_axis[efit_i,0],efit_axis[efit_i,1]]

	out={rbdry:rbdry, zbdry:zbdry, axis:axis}
	RETURN, out
END

;given the rbdry and zbdry points (from line_getlcfs) and the axis
;(from line_getaxis) determine if the pt=(r,z) is inside the LCFS
FUNCTION line_inlcfs, rbdry,zbdry,axis,point,debug=debug,shot=shot,time=time
	IF keyword_set(shot) and keyword_set(time) THEN BEGIN
		out=line_genbdry(shot,time)
		rbdry=out.rbdry
		zbdry=out.zbdry
		axis=out.axis
	ENDIF

	x=size(point)
	IF x[0] EQ 1 THEN output=intarr(1) ELSE output=intarr(x[2])

	FOR i=0,n(output) DO BEGIN
		pt=point[*,i]		
		d=(rbdry-pt[0])^2+(zbdry-pt[1])^2
		tmp=min(d)
		m=!C
		d[m]=max(d)
		tmp=min(d)
		n=!C
		pt_a=[rbdry[m], zbdry[m]]
		pt_b=[rbdry[n], zbdry[n]]
		pt_c=[0.5*(pt_a[0]+pt_b[0]),0.5*(pt_a[1]+pt_b[1])]
		IF keyword_set(debug) THEN BEGIN
			print, m,n
			print, 'Near Point #1', pt_a
			print, 'Near Point #2', pt_b
		ENDIF	

		vec_ab=pt_a-pt_b
		n1=[vec_ab[1],-vec_ab[0]]/sqrt(vec_ab[0]^2+vec_ab[1]^2)
		n2=[-vec_ab[1],vec_ab[0]]/sqrt(vec_ab[0]^2+vec_ab[1]^2)
	
		vec_c=pt_c-axis
		n_c=vec_c/sqrt(vec_c[0]^2+vec_c[1]^2)
	
		vec_pt=pt-pt_c
		n_p=vec_pt/sqrt(vec_pt[0]^2+vec_pt[1]^2)
	
		IF total(n_c*n1) GT 0 THEN nout=n1 ELSE nout=n2
		IF total(nout*n_p) GE 0 THEN output[i]=0 ELSE output[i]=1
		IF keyword_set(debug) THEN BEGIN
			print, 'Outward Normal', nout
			print, 'Vector to Point', n_p
		ENDIF	
	ENDFOR

	RETURN, output
END

FUNCTION line_invessel,r,z,del=del,plot=plot,debug=debug,shot=shot
	IF keyword_set(del) OR keyword_set(shot) THEN ves=line_vessel(/new, del=del,shot=shot) ELSE ves=line_vessel(/load)
	num_pts=n(r)+1
	good_pts=intarr(num_pts)


	;taken from GENPOS_GRID_INVESSEL on 8-1-07

	FOR i=0L,num_pts-1 DO BEGIN
		IF z[i] LT -0.35 AND r[i] LT 0.75 THEN BEGIN
                        CASE 1 OF
                            	r[i] GT 0.525 AND r[i] LT 0.6 AND z[i] GT -0.48 AND z[i] LE -0.35 : good_pts[i] = 1
                                r[i] GT 0.575 AND r[i] LT 0.615 AND z[i] GE -0.575 AND z[i] LE -0.48 : good_pts[i] = 1
                        	r[i] LE 0.525 AND z[i] GE -0.425 : BEGIN
                                	axis=[0.45,-0.4]
                                       IF line_inlcfs(ves.wall.r,ves.wall.z,axis,[r[i],z[i]]) EQ 0 $
                                        	THEN good_pts[i] = 1
                                END
                        	r[i] LE 0.525 AND z[i] LT -0.425 AND z[i] GT -0.48: BEGIN
                                	axis=[0.4,-0.5]
                                        IF line_inlcfs(ves.wall.r,ves.wall.z,axis,[r[i],z[i]]) EQ 0 $
                                        	THEN good_pts[i] = 1
                                END
                                r[i] GE 0.6 AND z[i] LE -0.35 AND z[i] GT -0.575 : BEGIN
                                	axis=[0.7,-0.45]
                                        IF line_inlcfs(ves.wall.r,ves.wall.z,axis,[r[i],z[i]]) EQ 0 $
                                        	THEN good_pts[i] = 1
                                END
                                r[i] LT 0.6 AND z[i] LE -0.48 AND z[i] GT -0.575 : BEGIN
                                	axis=[0.525,-0.525]
                                        IF line_inlcfs(ves.wall.r,ves.wall.z,axis,[r[i],z[i]]) EQ 0 $
                                       		THEN good_pts[i] = 1
                                END
                                ELSE:
                        ENDCASE
                ENDIF ELSE BEGIN 
                	axis=[0.68,0.2]
			IF line_inlcfs(ves.wall.r,ves.wall.z,axis,[r[i],z[i]]) EQ 1 THEN IF $
                          (z[i] LT -0.4 AND r[i] GT 0.65) EQ 0 THEN good_pts[i] = 1
                ENDELSE
                
        ENDFOR
        tmp=where(r GT 0.7 AND z LT -0.5)
        IF tmp[0] NE -1 THEN good_pts[tmp] = 0.0
        tmp=where(r LT 0.455 AND z LT -0.45)
        IF tmp[0] NE -1 THEN good_pts[tmp] = 0.0      
	
	IF keyword_set(plot) THEN BEGIN
            	vessel_plot,title='LINE_INVESSEL'
		makesym,9
		IF num_pts LE 3600 THEN oplot, r,z,color=200,psym=8
		makesym,10
		IF num_pts GT 3600 THEN use_sym=3 ELSE use_sym=8
		FOR i=0L,num_pts-1 DO IF good_pts[i] EQ 1 THEN oplot, [r[i]],[z[i]],color=200,psym=use_sym
		oplot, ves.wall.r,ves.wall.z
	ENDIF

	output=good_pts
	IF keyword_set(debug) THEN stop
	RETURN,output
END

		




;+
;NAME:
;	LINE_EMCONT
;
;PURPOSE:
;	The purpose of this function is to preprocess a midplane emissivity profile
;	so that it can be used to make an emissivity contour plot on a polodial
;	cross-section.
;
;CALLING_SEQUENCE:
;	output=LINE_EMCONT(emiss,radius,shot,tpt)
;
;INPUTS:
;	emiss:	FLTARR 	1D array of midplane emissivities 
;	radius:	FLTARR	1D array the radial scaling of emiss [m]
;	shot:	LONINT	shot number for EFIT to reference
;	tpt:	FLT	the time at which emiss is given computed [sec]
;
;OPTIONAL INPUTS:
;	ngrid:	ngrid=INT number of points in both R,Z to form emissivity [DEFAULT = 40]
;	rot:	STRUC containing
;		*.mach	1D array of ion Mach number at points of radius (M=v_tor/sqrt(2*Ti/mi))
;		*.z	INT value of the impurity ion
;		*.zave	1D array of average ionization of the impurity ion
;
;KEYWORD PARAMETERS:
;	debug:	/debug to stop code before the RETURN statement
;
;OUTPUTS:
;	output:	STRUC containing 3 arrays
;		*.em 1D array of emissivities in units of input emissivity
;		*.r  1D array of major radius points [m]
;		*.z  1D array of z-location points [m]
;
;RESTRICTIONS:
;	This code uses MLR_FUNCTIONS and EFIT_RZ2RMID as well as other functions in
;	GENIE_LINE.  It's best to use @genie_ini.bat to insure the latest codes are
;	compiled in the correct order.
;
;PROCEDURE:
;	This program uses the EFIT last closed flux surface (LCFS) at the time point
;	nearest to tpt, so to insure the most accurate information the input
;	emissivity should be computed at an EFIT time if possible.  The GENIE
;	function LINE_INLCFS is used to restrict the plotting grid to values
;	that are inside the LCFS.  This makes the arrays in the output structure have
;	an unpredictable length.  In order to make a contoure plot the following
;	should be used:
;
;	CONTOUR, emcon.em, emcon.r, emcon.z,/irregular [other contour keywords]
;
;MODIFICATION HISTORY:
;	Written by:	ML Reinke, 8-05-05
;
;-

FUNCTION line_emcont,emiss,r,shot,t_pt,ngrid=ngrid,debug=debug,sol=sol
	efit_time=line_gettimes(shot)
	efit_lcfs=line_getlcfs(shot)
	efit_axis=line_getaxis(shot)
	efit_i=ipt(efit_time,t_pt)
	IF efit_i[0] NE -1 THEN BEGIN
		rbdry=efit_lcfs[0,where(efit_lcfs[0,*,efit_i] NE 0 AND efit_lcfs[1,*,efit_i NE 0]),efit_i]
		zbdry=efit_lcfs[1,where(efit_lcfs[0,*,efit_i] NE 0 AND efit_lcfs[1,*,efit_i NE 0]),efit_i]
		rbdry=rbdry[0:n(rbdry)-1]
		zbdry=zbdry[0:n(zbdry)-1]
		axis=[efit_axis[efit_i,0],efit_axis[efit_i,1]]
	ENDIF ELSE RETURN, -1
	IF NOT keyword_set(ngrid) THEN ngrid=40
	IF NOT keyword_set(sol) THEN BEGIN
		rvec=make(min(rbdry),max(rbdry),ngrid)
		zvec=make(min(zbdry),max(zbdry),ngrid)
	ENDIF ELSE BEGIN
		rvec=make(0.43,0.93,ngrid)
		zvec=make(-0.5,0.5,ngrid)
	ENDELSE
	r2d=fltarr(ngrid*ngrid)
	z2d=fltarr(ngrid*ngrid)
	em2d=fltarr(ngrid*ngrid)
	cntr=0
	FOR i=0,ngrid-1 DO BEGIN
		FOR j=0,ngrid-1 DO BEGIN
			pt=[rvec[i],zvec[j]]
			IF NOT keyword_set(sol) THEN BEGIN
				IF line_inlcfs(rbdry,zbdry,axis,pt) THEN em2d[cntr]=1 ELSE em2d[cntr]=-1
			ENDIF			
			r2d[cntr]=pt[0]
			z2d[cntr]=pt[1]
			cntr+=1
		ENDFOR
	ENDFOR
	IF keyword_set(sol) THEN BEGIN
		invess=line_invessel(r2d,z2d,shot=shot)
		tmp=where(invess EQ 1)
	ENDIF ELSE tmp=where(em2d NE -1)
	rpts=r2d[tmp]
	zpts=z2d[tmp]
	;stop
	rmid=efit_rz2rmid(rpts,zpts,t_pt,shot=shot)
	empts=interpol(emiss,r,rmid)
	tmp=where(rmid LE max(r))
	IF keyword_set(debug) THEN stop
	emcon={em:empts[tmp], r:rpts[tmp],z:zpts[tmp]}
	RETURN, emcon
END

FUNCTION line_bfield,pos,s,shot,t_pts=t_pts
	
	n_s=n(s)+1
	rpts=line_r(s,pos)
	zpts=line_z(s,pos)
	
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
	psi=efit_rz2psi(rpts,zpts,t_pts,bz,br,shot=shot)
	
	bphi=fltarr(n_s,n_tpts)
	psi_norm=fltarr(n_s,n_tpts)
	FOR i=0, n(t_pts) DO BEGIN
		psi_norm[*,i]=(psi[*,i]-psi0[i])/(psibndry[i]-psi0[i])
		bphi[*,i]=Bt0[i]*R0[i]/rpts
	ENDFOR

	output={br:br, bphi:bphi, bz:bz, psi_norm:psi_norm, t_pts:t_pts, shot:shot}

	RETURN, output
END



FUNCTION line_pos2psith,pos,shot,time,debug=debug

	n=n(pos[0,*])
	tang=fltarr(2,n+1)-1.0
	psith=fltarr(2,n+1)-1.0
	pts=[pos[2,*],pos[1,*]-tan(pos[3,*])*sqrt(pos[0,*]^2-pos[2,*]^2)]
	inout=line_inlcfs(-1,-1,axis,pts,shot=shot,time=time)
	
	FOR i=0,n DO IF inout [i] THEN BEGIN
			pt=pts[*,i]
			pt-=axis
			th=atan(pt[1]/pt[0])
			IF pt[1] GT 0 AND pt[0] LT 0 THEN th=th+!pi
			IF pt[1] LT 0 AND pt[0] LT 0 THEN th=th+!pi
			IF pt[1] LT 0 AND pt[0] GT 0 THEN th=th+2.0*!pi
			psith[1,i]=th
	ENDIF

	tmp=where(inout NE 0)
	IF tmp[0] EQ -1 THEN RETURN,-1
	psith[0,tmp]=efit_rz2rho(reform(pts[0,tmp]),reform(pts[1,tmp]),time,shot=shot,/psinorm)
	 
	IF keyword_set(debug) THEN stop

	output={psi:reform(psith[0,*]),th:reform(psith[1,*])}
	RETURN,output
END				


;+
; NAME:
;	LINE_BR
;
; PURPOSE:
;	This function calculates the line integrated brightness from a given
;	emissivity profile and EFIT reconstruction for a given viewing geometry.
;	The point is to forward model theoretical emissivities (line or total) to
;	compare to the brightness of a diagnostic measurement.
;
; CALLING SEQUENCE:
;	BR = LINE_BR(pos, emiss, r, t, shot, t_pts)
;
; INPUTS:
;	POS: 	FLTARR 2D (m,4) defines the diagnostic view(s) (see below)
;	EMISS:	FLTARR 2D array in (nR,nt) of the emissivity (line or total) [(ph/s) or (W)/m^3]
;	R:	FLTARR 1D array of the radial scaling of emiss [meters]
;	T:	FLTARR 1D array of the time scaling of emiss [sec]
;	SHOT:	shot number for the EFIT reconstruction
;	T_PTS:	time point(s) to perform the line integration [sec]
;
; OPTIONAL INPUTS:
;	N_S:		INT  number of points along the line to calculate emiss(s) 
;				from emiss(r,z) the [DEFAULT = 75]
;	RZBND:		FLTARR boundary [r_in, r_out, abs(z)] for line of sight [DEFAULT = [0.44,1.0,0.5]]
;	PLOTLAB: 	STRUC *.line, *.wl are the line name (MO XXXI) and wavelength (115.99)
;				to be used then /plots is invoked. 
;
; KEYWORD PARAMETERS:
;	DEBUG:	/debug stops the code at various locations for debugging
;	VERB:	/verb prints terminal updates (very verbose terminal updates)
;	PLOTS:	/plots to shot a plot of C-Mod polodial cross section with views defined by pos
;	PS:	/ps (/plots must be set) to supress window calls so that a PS DEVICE can be used
;
; OUTPUTS:
;	The output is a 2D FLTARR (m,nT_PTS) of line integrated brightness, in the 
;	units of [(ph/s) or (W)/m^2], at the time points defined by T_PTS and the m-views defined
;	by the 2D POS array.  Values = 1.0e-4 are places where some part of the code failed.
;
; OPTIONAL OUTPUTS:
;	If a /plots is invoked, window 19 will be opened/cleared.  Each index of t_pts 
;	and each view will have a seperate plot and the user will need to enter a keystroke 
;	to cycle through them.  On each plot a cross-section of C-Mod with limiters will be 
;	displayed along with the LCFS and the line of sight defined by POS.  The function
;	LINE_EMCONT is used to make an (R,Z) contour plot of the emissivity assuming it is
;	constant on a flux surface.  A normalized emissivity as a function of pathlength is 
;	plotted and the peak value is displayed.  Color tables are forced to the values
;	I want them to be, cause I'm a color tyrant and you'll have to live with that. IF
;	a /ps is used in conjection with /plots user input is disabled and IDL.PS is filled
;	with all the plots.  The device is left open so that other programs can build the
;	postscript with multiple emissivity profiles. 
;
; PROCEDURE:
;	The line through the plasma is completely described by POS assuming toroidal
;	invariance in the plasma. Each row of POS = [ro, zo, rt, thz] where (ro, zo)
;	define a point on a line of sight outside the plasma (usually a pinhole location),
;	(rt) is the tangency radius of the view and (thz) is the angle the view makes w/r/t a 
;	z=constant surface (with inc < 0, dec > 0).
;
; RESTRICTIONS:
;	This function assumes that the plasma emissivity is toroidally symmetric and that
;	emissivity is constant on a flux surface (not generally true).  Use with caution
;	or some sort of verification methods.
;
;	This function calls EFIT_RZ2MID.PRO
;	This function needs MLR_FUNCTIONS
;	@GENIE_INI.BAT should be used to compile and load code in the correct order
;
; MODIFICATION HISTORY
;	Written by:	ML Reinke, June 2005
;	8-03-05:	ML Reinke, added some bug fixes if emiss was [n,1] in size (1 time point)
;	8-05-05:	ML Reinke, added a plot feature to show views on a C-Mod cross section
;	8-07-05:	ML Reinke, improved plot features and updated documentation
;      12-06-05:	ML Reinke, added a torodial plot if pos[2] (rtang) greater than zero
;-


FUNCTION line_br,pos,emiss,r,t,shot,t_pts, n_s=n_s, verb=verb,debug=debug,plots=plots,ps=ps,plotlab=plotlab,rmid=rmid,efit=efit

	IF NOT keyword_set(n_s) THEN n_s=75	;set number of segments to define line of sight
	IF NOT keyword_set(rzbnd) THEN rzbnd=[0.44,1.0,0.5]	;set path's [r_min,r_max,abs(z_max)] limit
	
	;determine the number of views
	xpos=size(pos)	
	IF xpos[0] EQ 1 THEN n_view=1 ELSE n_view=xpos[2]

	;determine the number of time points
	n_tpts=n(t_pts)+1
		
	;intialize output brightness array (number of views x number of time points)
	br=fltarr(n_view, n_tpts)

	;load EFIT data
	IF NOT keyword_set(efit) THEN BEGIN
		efit_time=line_gettimes(shot)
		efit_lcfs=line_getlcfs(shot)
		efit_axis=line_getaxis(shot)
		efit={time:efit_time,lcfs:efit_lcfs,axis:efit_axis}
	ENDIF ELSE BEGIN
		efit_time=efit.time
		efit_lcfs=efit.lcfs
		efit_axis=efit.axis
	ENDELSE
	

	;if plots are set, get wall and limiter traces
	IF keyword_set(plots) THEN BEGIN
		mdsopen, 'analysis',-1
		IF shot GT 1020000000 THEN  node_year='2002' else node_year='1994'
		r_lim=mdsvalue('.LIMITERS.WALL.WALL_'+node_year+':RLIM')
		z_lim=mdsvalue('.LIMITERS.WALL.WALL_'+node_year+':ZLIM')
		r_rf=mdsvalue('.LIMITERS.RF_LIMITER:R')
		z_rf=mdsvalue('.LIMITERS.RF_LIMITER:Z')
		mdsclose
	ENDIF

	;initialize the PS device if chosen
	IF keyword_set(plots) AND keyword_set(ps) THEN BEGIN
		set_plot,'ps'
		!p.thick=2.75
		!x.thick=2.75
		!y.thick=2.75
		device, /color
		device, /portrait
		device, encapsulated=0
		device, preview=0
		device, xsi=20.0, ysi=32.0
	ENDIF
	
	;debug information
	IF keyword_set(verb) THEN print, ' n_view = '+num2str(n_view)+'  n_tpts = '+num2str(n_tpts)
	IF keyword_set(debug) THEN STOP

	;start filling the brightness array cycling through each view
	FOR i=0,n_view-1 DO BEGIN
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
		;there is now array of points that go from the pinhole to the closest limiting
		;surface and is divided into n_s parts.  

		;cycle through all the time points
		FOR j=0,n_tpts-1 DO BEGIN
			IF keyword_set(verb) THEN print, ' j = '+num2str(j,1)
			emiss_s=fltarr(n_s)	;initialize emissivity as a function of path length array
			efit_i=ipt(efit_time,t_pts[j])
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

					;find their equivilent midplane radii
					IF NOT keyword_set(rmid) THEN r_on_mid=efit_rz2rmid(good_r, good_z,efit_time[efit_i], shot=shot) ELSE r_on_mid=rmid
					IF n_view EQ 1 AND n_tpts EQ 1 THEN rmid=r_on_mid
					IF keyword_set(verb) THEN BEGIN
						print, 'Good RMID'
						print, r_on_mid
					ENDIF

					;if any, set values of r_on_mid that are < min(r) to min(r)
					tmp=where(r_on_mid LT min(r))
					IF tmp[0] NE -1 THEN r_on_mid[tmp]=min(r)
		
					;use IDL function INTERPOLATE to finde emiss_s values
					;special case is isolated of emiss is a 1D array
					r_interp=interp_vec_reform(r, r_on_mid)
					t_interp=interp_vec_reform(t,fltarr(n(r_on_mid)+1)+t_pts[j])
					IF n(t) EQ 0 THEN t_interp=fltarr(n(r_on_mid)+1)	;if only 1 time slice, default to it
					IF keyword_set(verb) THEN BEGIN
						print, 'R_INTERP'
						print, r_interp
						print, 'T_INTERP'
						print, t_interp
					ENDIF
					IF n(t) NE 0 THEN emiss_good=interpolate(emiss,r_interp,t_interp, missing=1.0e-4) ELSE $
						emiss_good=interpolate(emiss,r_interp,missing=1.0e-4)
					emiss_s[goodpts]=emiss_good
					IF keyword_set(verb) THEN BEGIN
						print, 'Good EMISS_S'
						print, emiss_good
					ENDIF

					ipos=reform(pos[*,i])
					;the scaling of l is arbitrary with no physical signifigance except that
					;l=0 is pt1 and l=1 is pt2 if the l values are multiplied by the distance
					;between 1 and 2 then the qunatitative units are applied to the scaling.
					lscale=sqrt((ipos[0]^2-ipos[2]^2)*(1.0+(tan(ipos[3]))^2))
					br[i,j]=int_tabulated(s, emiss_s*lscale)

					;done with the math, but if PLOTS are on we have a ways to go
					IF keyword_set(plots) THEN BEGIN
						;if plotting to X-window setup a window w/ correct aspect ratio of
						;Alcator C-Mod (this is done by eye not by quantitative process)
						IF NOT keyword_set(ps) THEN BEGIN
							device, window_state=var
							IF var[19] EQ 0 THEN window,19,xsize=500,ysize=690,xpos=0,ypos=670,$
								title='vessel cx,19' ELSE wset,19
						ENDIF
			
						;plot the limiting structures
						plot, r_lim,z_lim,title='SHOT: '+num2str(shot,1)+' TIME: '+num2str(t_pts[j],dp=3),$
							chars=1.3
						oplot, r_rf,z_rf
						
						;find the 2D emissivity countours
						emcon=line_emcont(emiss[*,ipt(t,t_pts[j])],r,shot,efit_time[efit_i],ngrid=80)
						loadct,1,/silent
						contour,emcon.em/max(emcon.em),emcon.r,emcon.z,/irregular,nlevels=200,$
							/fill,/overplot,levels=make(0.0,1.0,200)
						loadct,12,/silent

						IF keyword_set(plotlab) THEN BEGIN
							xyouts,0.75,0.45,'Line: '+plotlab.line+' '+plotlab.wl,charsize=1.3
						ENDIF
	
						;plot the LCFS boundary
						loadct,12,/silent						
						rb_plt=[rbdry,rbdry[0]]
						zb_plt=[zbdry,zbdry[0]]
						oplot,rb_plt,zb_plt,color=100,thick=1.5

						;plot the line of sight on the polodial plane
						oplot,line_r(s, pos[*,i]),line_z(s,pos[*,i]),color=200,thick=2.0
						
						;label information about the view
						xyouts,0.65,0.4, 'POS = ['+num2str(ipos[0],dp=2)+','+num2str(ipos[1],dp=2)+','+$
							num2str(ipos[2],dp=2)+','+num2str(ipos[3],dp=2)+']',charsize=1.3
						
						;label and plot a normalized emissivity as along the path.
						xyouts,0.65,0.5, 'Max EMISS_S = '+num2str(max(emiss_s)),color=30,charsize=1.3
						oplot,reverse(s*lscale+0.44),(emiss_s/max(emiss_s)*0.25-0.6),$
							color=30,thick=2.0

						;plot a dotted line at the flux surface of the maximum of emiss
						fs=line_getfs(shot,efit_time[efit_i],r[maxloc(emiss)])
						fs_type=size(fs, /type)
						IF fs_type EQ 8 THEN oplot, fs.r,fs.z, linestyle=2, color=100, thick=3.0

						;if there is a finite tangency radius, plot a Z=const view of the line.	
						IF pos[2] NE 0 THEN BEGIN
							;if plotting to x-windows then setup top-down view
							IF NOT keyword_set(ps) THEN BEGIN
								device, window_state=var
								IF var[20] EQ 0 THEN window,20,xsize=800,ysize=400,xpos=500,ypos=670,$
									title='vessel cx,20' ELSE wset,20
							ENDIF
							alpha=make(0,2*!pi,100)
							;plot inner wall
							x=0.44*cos(alpha)
							y=0.44*sin(alpha)
							plot,x,y,xr=[-1.0,1.0],yr=[-1.0,0.0],title='SHOT: '+num2str(shot,1)+$
								' TIME: '+num2str(t_pts[j],dp=3),chars=1.3
							IF keyword_set(plotlab) THEN BEGIN
								xyouts,0.5,-0.25,'Line: '+plotlab.line+' '+plotlab.wl,charsize=1.3
							ENDIF
							xyouts,0.4,-0.2, 'Max EMISS_S = '+num2str(max(emiss_s)),color=30,charsize=1.3
							xyouts,0.4,-0.3, 'POS = ['+num2str(ipos[0],dp=2)+','+num2str(ipos[1],dp=2)+','+$
								num2str(ipos[2],dp=2)+','+num2str(ipos[3],dp=2)+']',charsize=1.3
							;plot outerlimiter
							rlim=max(r_rf)
							x=rlim*cos(alpha)
							y=rlim*sin(alpha)
							oplot,x,y
							;plot inner LCFS
							y=min(rb_plt)*sin(alpha)
							x=min(rb_plt)*cos(alpha)
							oplot,x,y,color=100,thick=1.5

							;plot outer LCFS
							y=max(rb_plt)*sin(alpha)
							x=max(rb_plt)*cos(alpha)
							oplot,x,y,color=100,thick=1.5

							;plot line of sight on z=const plane
							oplot,[-1.5,1.5],[-pos[2],-pos[2]],color=200,thick=2.0

							goodpts=[goodpts[0]-1,goodpts,goodpts[n(goodpts)]+1]
							good_s=s[goodpts]
							good_em=emiss_s[goodpts]
							delR=sqrt(ipos[0]^2-ipos[2]^2)	;same as lscale but for projection onto Z=const
							oplot,(good_s-good_s[0])*delR-sqrt(line_r(good_s[0],ipos)^2-ipos[2]^2),$
								good_em/max(good_em)*0.25*ipos[2]-ipos[2],color=30,thick=3.0
						ENDIF
						IF NOT keyword_set(ps) THEN BEGIN
							IF  i LT n_view-1 OR j LT n_tpts-1  THEN BEGIN
								print, 'Press ENTER key to continue'
								inputcommand=''
								read, inputcommand		
							ENDIF
						ENDIF
					ENDIF
					
				ENDIF ELSE br[i,j]=0.0
			ENDIF ELSE br[i,j]=-1	;value out of efit time range	
		ENDFOR
		
	ENDFOR	
	output=br
	IF keyword_set(debug) THEN stop
	RETURN, output
END

;+
;NAME:
;	LINE_PATH_PLOTS
;	
;PURPOSE:
;	This procedure plots the line of sight of an arbitrary number of POS vectors
;	in both poloidal and toroidal projections.
;
;CALLING SEQUENCE:
;	LINE_PATH_PLOTS,pos
;
;INPUTS:
;	pos: 	FLTARR 2D (m,4) defines the diagnostic view(s)
;
;OPTIONAL INPUTS:
;	shot:	LONINT of the shot number to display flux surface information DEFAULT=1050426022 (everyone's favorite)
;	tpt:	FLT of the time point to display flux surface information DEFAULT=1.0
;	thick:	FLT of the line thickness DEFAULT=2.0
;	fovpos:	FLTARR of pos vectors for the field of view error bars (see GENPOS)
;	rzbnd:	FLTARR boundary [r_in, r_out, abs(z)] for line of sight [DEFAULT = [0.44,1.0,0.6]]
;
;KEYWORD PARAMETERS:
;	div:	/div will make the poloidal view zoomed into the divertor region
;	tang:	/tang will display tangency points views that are constrained to a single Z-plane
;	verb:	/verb display line of sight info (same as LINE_BR)
;	debug:	/debug stops at various locations
;	ps:	/ps supresses WSET command.  Postscript device must be setup prior to call
;
;OUTPUTS:
;	Outputs to the X-windows plotting device will be made on windows 19 and 20.  If currently being
;	used, close so that LINE_PATH_PLOTS resizes them appropriately.  Note that switching between divertor
;	and full view requires the window to be closed for resizing.
;
;MODIFICATION HISTORY:
;	Written by:	M.L. Reinke - Spring 2006
;	8-24-06:	ML Reinke - added the rzbnd optional input 
;	7-13-07:	ML Reinke - made this program use the vessel cross-section from VESSEL_PLOT in MLR_FUNCTIONS
;-


PRO line_path_plots,pos,shot=shot,tpt=tpt,ps=ps,verb=verb,div=div,debug=debug,thick=thick,tang=tang,fovpos=fovpos,rzbnd=rzbnd,nsol=nsol,$
		sol=sol,tor=tor,edge=edge,bspline=bspline,n_s=n_s,inv=inv,col=col

	IF NOT keyword_set(tor) THEN notor=1 ELSE notor=0
	IF NOT keyword_set(col) THEN col=200
	IF NOT keyword_set(shot) THEN shot=1050426022	;shot for flux surfraces
	IF NOT keyword_set(tpt) THEN tpt =1.0		;time point for flux surfaces
	IF NOT keyword_set(thick) THEN thick=2.0
	IF NOT keyword_set(n_s) THEN n_s=100	;set number of segments to define line of sight
	IF NOT keyword_set(rzbnd) THEN rzbnd=[0.44,1.0,0.6]	;set path's [r_min,r_max,abs(z_max)] limit
	IF NOT keyword_set(nsol) THEN nsol=4
	
	;determine the number of views
	xpos=size(pos)	
	IF xpos[0] EQ 1 THEN n_view=1 ELSE n_view=xpos[2]

	;determine the number of time points
	n_tpts=1

	;load EFIT data
	efit_time=line_gettimes(shot)
	efit_lcfs=line_getlcfs(shot)
	efit_axis=line_getaxis(shot)
	efit_i=ipt(efit_time,tpt)
	efit_rmid=line_getrmid(shot)

	;load limiter data
	mdsopen, 'analysis',-1
	IF shot GT 1020000000 THEN  node_year='2002' else node_year='1994'
	r_lim=mdsvalue('.LIMITERS.WALL.WALL_'+node_year+':RLIM')
	z_lim=mdsvalue('.LIMITERS.WALL.WALL_'+node_year+':ZLIM')
	r_rf=mdsvalue('.LIMITERS.RF_LIMITER:R')
	z_rf=mdsvalue('.LIMITERS.RF_LIMITER:Z')
	mdsclose

	;plot the limiting structures
	IF keyword_set(div) THEN BEGIN
		yrange=[-0.6,-0.2]
		xrange=[0.4,0.8]
	ENDIF ELSE BEGIN
		yrange=[-0.6,0.6]
		xrange=[0.4,1.2]
	ENDELSE

	IF keyword_set(bspline) THEN BEGIN
		bspline_fs,shot,tpt,div=div,edge=edge
	ENDIF ELSE BEGIN
		vessel_plot,title='SHOT: '+num2str(shot,1)+' TIME: '+num2str(tpt,dp=3),ps=ps,shot=shot,d_old=d_old,div=div,edge=edge

		;plot the flux surfaces
		efit_psicont,psi_r,psi_z,psi_n,psi_t,shot=shot,time=efit_time[efit_i]
		FOR i=0,floor(n(efit_rmid[efit_i,*])/3.0)-1 DO BEGIN
			out_r=reform(psi_r[*,3*i+1])
			out_z=reform(psi_z[*,3*i+1])
			tmp=where(out_r NE 0)
			oplot, out_r[tmp], out_z[tmp], color=100,linestyle=2,thick=1.0
		ENDFOR
	ENDELSE

	IF keyword_set(debug) THEN stop

	;load relevant LCFS boundary and prepare for format needed for LINE_INCLFS
	rbdry=efit_lcfs[0,where(efit_lcfs[0,*,efit_i] NE 0 AND efit_lcfs[1,*,efit_i] NE 0),efit_i]
	zbdry=efit_lcfs[1,where(efit_lcfs[0,*,efit_i] NE 0 AND efit_lcfs[1,*,efit_i] NE 0),efit_i]
	rbdry=rbdry[0:n(rbdry)-1]
	zbdry=zbdry[0:n(zbdry)-1]
	rb_plt=[rbdry,rbdry[0]]
	zb_plt=[zbdry,zbdry[0]]
	IF NOT keyword_set(bspline) THEN oplot,rb_plt,zb_plt,color=100,thick=2.0

	IF keyword_set(sol) THEN BEGIN
		z_sol=zbdry[maxloc(rbdry)]
		sol_pts=make(max(rbdry),0.91,nsol+2)
		sol_pts=sol_pts[1:nsol]
		FOR i=0,nsol-1 DO BEGIN
			sol_tracet,shot,efit_time[efit_i],sol_pts[i],zbdry,xcontr,ycontr,ncontr,psivl,direction=1,/accurate,/nolim
			oplot, xcontr,ycontr,color=30
			sol_tracet,shot,efit_time[efit_i],sol_pts[i],zbdry,xcontr,ycontr,ncontr,psivl,direction=2, /accurate,/nolim
			oplot, xcontr,ycontr,color=30
		ENDFOR
	ENDIF
	

	;debug information
	IF keyword_set(debug) THEN STOP

	;start plotting each view
	smin=fltarr(n_view)
	FOR i=0L,n_view-1 DO BEGIN
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
		IF smin_r GT smin_z THEN smin[i] = smin_z ELSE smin[i]=smin_r		;set minimum s
		s=make(0,smin[i],n_s)
		IF keyword_set(verb) THEN print, s		
		;there is now array of points that go from the pinhole to the closest limiting
		;surface and is divided into n_s parts.  

		;plot the line of sight on the polodial plane
		
		IF total(abs(pos[3,*])) NE 0 THEN BEGIN
			rplot=line_r(s,pos[*,i])
			zplot=line_z(s,pos[*,i])
			IF keyword_set(inv) THEN BEGIN
				good=line_invessel(rplot,zplot,shot=shot)
				tmp=where(good EQ 0)
				oplot,rplot[0:tmp[0]],zplot[0:tmp[0]],color=col,thick=thick
			ENDIF ELSE oplot,rplot,zplot,color=col,thick=thick
		ENDIF
		IF keyword_set(tang) AND NOT keyword_set(fovpos) THEN BEGIN
			makesym,10
			oplot, [pos[2,i]],[pos[1,i]-tan(pos[3,i])*sqrt(pos[0,i]^2-pos[2,i]^2)],psym=8,color=col
		ENDIF
		IF keyword_set(fovpos) THEN BEGIN
			out=pos2width(pos,fovpos)
			oploterror,[pos[2,i]],[pos[1,i]-tan(pos[3,i])*sqrt(pos[0,i]^2-pos[2,i]^2)],out.r[i],$
				out.z[i],psym=8,color=col,errcolor=col
		ENDIF

	ENDFOR

	IF NOT keyword_set(notor) THEN BEGIN
		;if plotting to x-windows then setup top-down view
		IF NOT keyword_set(ps) THEN BEGIN
			device, window_state=var
			IF var[20] EQ 0 THEN window,20,xsize=585,ysize=785,xpos=500,ypos=670,title='vessel cx,20' ELSE wset,20
		ENDIF ELSE BEGIN
			device, xsize=float(d_old.x_size)/d_old.x_px_cm,ysize=float(d_old.y_size)/d_old.y_px_cm
		ENDELSE
			
		alpha=make(0,2*!pi,100)

		;plot inner wall
		x=0.44*cos(alpha)
		y=0.44*sin(alpha)	
		plot,x,y,xr=[-0.1,1.0],yr=[-1.0,0.6],title='SHOT: '+num2str(shot,1)+' TIME: '+num2str(tpt,dp=3),chars=1.3,/xsty,/ysty

		;plot outerlimiter
		rlim=max(r_rf)
		x=rlim*cos(alpha)
		y=rlim*sin(alpha)	
		oplot,x,y	

		;plot inner LCFS
		y=min(rb_plt)*sin(alpha)
		x=min(rb_plt)*cos(alpha)
		oplot,x,y,color=100,thick=2.0

		;plot outer LCFS
		y=max(rb_plt)*sin(alpha)
		x=max(rb_plt)*cos(alpha)
		oplot,x,y,color=100,thick=2.0

		FOR i=0,n_view-1 DO BEGIN
			ipos=reform(pos[*,i])
			s=make(0,smin[i],n_s)
			xpts=line_r(s,ipos)*cos(line_th(s,ipos,/cont)-!pi/2.0)
			ypts=line_r(s,ipos)*sin(line_th(s,ipos,/cont)-!pi/2.0)
					

			;plot line of sight on z=const plane
			oplot,xpts,ypts,color=200,thick=2.0
	
		ENDFOR
	ENDIF 
	IF keyword_set(debug) THEN stop

END







;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;scripts to show the above works
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;






;a script to show that the AXUV brightness profile can be reformed from the emissivity using LINE_BR
PRO test_linebr,shot,tpt,n_s=n_s,verb=verb,ape=ape

	load_emdiode_data,shot,emiss,r,t,/tree	
	load_axuvbr_data,shot,br,rtang,t_br,/tree

	IF NOT keyword_set(ape) THEN ape=[0.0,0.0]

	ro=fltarr(16)+0.935740+ape[0]
	zo=fltarr(16)+ape[1]
	psi=fltarr(16)
	pos=fltarr(4,16)
	pos[0,*]=ro
	pos[1,*]=zo
	pos[2,*]=rtang
	pos[3,*]=psi
	print, pos
	calc_br=line_br(pos,emiss,r,t,shot,tpt,n_s=n_s,verb=verb)

	raw_br=br[*,ipt(t_br,tpt)]

	plot,raw_br,psym=-6,xtitle='AXUV Channel', ytitle='Brightness [W/m^2-str]', $
		title='Raw_BR (BLACK) vs Calc_BR (RED)',chars=1.3
	oplot,calc_br,color=200
END

;just a little script to show line_inclfs works
PRO test_inlcfs,npts,seed
	zbdry=[0,1,1,1,0,-1,-1,-1]
	rbdry=[1.0,1,0,-1,-1,-1,0,1]
	axis=[0,0]
	pt_r=randomu(1.0*seed,npts,/normal)
	pt_z=randomu(100.0*seed,npts,/normal)
	plot, [rbdry,1.0],[zbdry,0],xrange=[-2.5,2.5],yrange=[-2.5,2.5] 
	FOR i=0,npts-1 DO BEGIN
		result=line_inlcfs(rbdry,zbdry,axis,[pt_r[i],pt_z[i]])
		IF result EQ 0 THEN oplot,[pt_r[i]],[pt_z[i]],psym=7
		IF result EQ 1 THEN oplot,[pt_r[i]],[pt_z[i]],psym=6
	ENDFOR
END

;just a little script to show line_inlcfs works with an actualy plasma boundary
PRO plasma_inlcfs,npts,seed
	;zbdry=[0,1,1,1,0,-1,-1,-1]
	;rbdry=[1,1,0,-1,-1,-1,0,1]
	;axis=[0,0]
	shot=1050426022
	efit_i=20
	lcfs=line_getlcfs(shot)
	axis=line_getaxis(shot)
	;remove any bad points which are set to (0,0)
	rbdry=lcfs[0,where(lcfs[1,*,efit_i] NE 0 AND lcfs[0,*,efit_i] NE 0),efit_i]
	zbdry=lcfs[1,where(lcfs[1,*,efit_i] NE 0 AND lcfs[0,*,efit_i] NE 0),efit_i]
	;remove last duplicated point
	rbdry=rbdry[0:n(rbdry)-1]
	zbdry=zbdry[0:n(zbdry)-1]
	pt_r=randomu(1.0*seed,npts,/normal)
	pt_z=randomu(100.0*seed,npts,/normal)
	plot, [rbdry],[zbdry],xrange=[0,1.0],yrange=[-1.0,1.0]
	print, axis[efit_i, *]
	FOR i=0,npts-1 DO BEGIN
		result=line_inlcfs(rbdry,zbdry,axis[efit_i,*],[pt_r[i],pt_z[i]])
		IF result EQ 0 THEN oplot,[pt_r[i]],[pt_z[i]],psym=7
		IF result EQ 1 THEN oplot,[pt_r[i]],[pt_z[i]],psym=6
	ENDFOR
END




