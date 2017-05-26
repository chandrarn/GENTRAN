;+
;NAME:
;	VESSEL_PLOT
;
;PURPOSE:
;	This procedures sends a properly scaled plot of a poloidal C-Mod vacuum vessel/tile cross-section 
;	to either the postscript or x-windows device.
;
;CALLING SEQUENCE:
;	VESSEL_PLOT
;
;OPTIONAL INPUTS:
;	n:		INT:	The window number to plot to
;	color:		INT:	Color of the vessel plot lines (not color of plot) DEFAULT: 0
;	title:		STR:	String of the plot title DEFAULT: "NSTX"
;	thick:		INT:	Thickness of the vessel plot lines DEFAULT: 0
;	shot:		LON:	Shot number to determine which vessel cross-section to use DEFAULT: 135440
;
;KEYWORD PARAMETERS:
;	force:		/force will force window number n to reopen and resize it accordingly
;	oplot:		/oplot will oplot the vessel structure on the currently selected plotting device
;	div:		/div will create a properly scaled zoomed plot of the lower divertor
;	edge:		/edge will create a properly scaled zoomed plot of the outer edge
;
;OUTPUTS:
;	All outputs are sent to the currently selected plotting device.  VESSEL_PLOTS will do nothing if a device other than
;	the X or PS is selected.
;
;PROCEDURE:
;	The vessel cross-section plots are taken from Brian Labombard's data files:
;			/home/labombard/minicad/vv_tiles_cryo_2007_s.vctr	2007-current
;			/home/labombard/minicad/vv_and_tiles_2002_s.vctr   	2002-2007
;			/home/labombard/minicad/vacuumvessel.vctr
;			/home/labombard/minicad/tiles_2002_s.vctr		1993-2002
;
;MODIFICATION HISTORY:
;	Written by:	ML Reinke, Fall 2006
;	5-10-07:	ML Reinke - added the new cyropump vessel cross-section and made it the default
;	7-02-15:	ML Reinke - modified for use with NSTX
;-

PRO vessel_plot,ps=ps,n=n, force=force,title=title,x_offset=x_offset,y_offset=y_offset,color=color,thick=thick,oplot=oplot,$
		shot=shot,d_old=d_old,div=div,edge=edge,pedestal=pedestal,wall=wall,isnstxu=isnstxu
	IF NOT keyword_set(shot) THEN shot=135440
	IF NOT keyword_set(title) THEN BEGIN
		title='NSTX'
		IF keyword_set(isnstxu) THEN title='NSTX-U'
	ENDIF
	IF NOT keyword_set(x_offset) THEN x_offset=0.0
	IF NOT keyword_set(y_offset) THEN y_offset=0.0
	
	CASE !d.name OF
		'PS': 	ps=1
		'X':	ps=0
		ELSE:	RETURN
	ENDCASE


	IF NOT keyword_set(div) AND NOT keyword_set(edge) AND NOT keyword_set(pedestal) AND NOT keyword_set(wall)  THEN BEGIN
		yrange=[-1.75,1.75]
		xrange=[0.0,1.75]
		IF keyword_set(ps) THEN BEGIN
			xsize=5.0
			ysize=8.15
		ENDIF ELSE BEGIN
			xsize=500.0
			ysize=880.0
		ENDELSE
	ENDIF ELSE BEGIN
		IF keyword_set(div) THEN BEGIN
			yrange=[-1.75,-0.7]
			xrange=[0.1,1.1]
			IF keyword_set(ps) THEN BEGIN
				xsize=7.0
				ysize=6.85
			ENDIF ELSE BEGIN
				xsize=700.0
				ysize=710.0
			ENDELSE
		ENDIF
; 		IF keyword_set(edge) THEN BEGIN
; 			yrange=[-0.35,0.35]
; 			xrange=[0.574,0.95]
; 			IF keyword_set(ps) THEN BEGIN
; 				xsize=6.0
; 				ysize=6.0*810.0/500.0
; 			ENDIF ELSE BEGIN
; 				xsize=500
; 				ysize=818
; 			ENDELSE	
; 		ENDIF
; 		IF keyword_set(wall) THEN BEGIN
; 			yrange=[-0.35,0.35]
; 			xrange=[0.574,0.95]-0.174
; 			IF keyword_set(ps) THEN BEGIN
; 				xsize=6.0
; 				ysize=6.0*810.0/500.0
; 			ENDIF ELSE BEGIN
; 				xsize=500
; 				ysize=818
; 			ENDELSE	
; 		ENDIF
; 		IF keyword_set(pedestal) THEN BEGIN
; 			yrange=[-0.1,0.1]
; 			xrange=[0.822,0.93]
; 			IF keyword_set(ps) THEN BEGIN
; 				xsize=6.0
; 				ysize=6.0*810.0/500.0
; 			ENDIF ELSE BEGIN
; 				xsize=500
; 				ysize=818
; 			ENDELSE	
; 		ENDIF
	ENDELSE

	IF NOT keyword_set(oplot) THEN BEGIN
		IF NOT keyword_set(ps) THEN BEGIN
			IF NOT keyword_set(n) THEN n=19
			device, window_state=var
			IF var[n] EQ 0 OR keyword_set(force) EQ 1 THEN window,n,xsize=xsize,ysize=ysize,xpos=0,ypos=670,title='vessel cx,'+num2str(n) $
				ELSE wset,n
		ENDIF ELSE BEGIN
			d_old=!d
			device, xsize=xsize, ysize=ysize, /inches
		ENDELSE
		plot, [0],[0],title=title,chars=1.3,xrange=xrange,yrange=yrange,xtit='R (m)',ytit='Z (m)',/xsty, /ysty
	ENDIF
	
	;rlimu=[581.02,414.94,414.94,341.74,341.74,341.74,414.94,414.94,581.02]/1.0e3	;from CAD drawing
	;zlimu=[1623.39,1623.39,1269.91,1050.06,0,-1050.06,-1269.91,-1623.39,-1623.39]/1.0e3
	IF keyword_set(isnstxu) THEN BEGIN
		;rlim=[reverse(rlimu[0:4]),rlim[6:30,0],reverse(rlimu[4:*])]
		;zlim=[reverse(zlimu[0:4]),zlim[6:30,0],reverse(zlimu[4:*])]
		mdsopen,'efit01',200184
		rlim=mdsvalue('\EFIT01::TOP.RESULTS.GEQDSK:RLIM')
		zlim=mdsvalue('\EFIT01::TOP.RESULTS.GEQDSK:ZLIM')
		mdsclose,'efit01',200184
  		rlim=rlim[*,0]
		zlim=zlim[*,0]
      	ENDIF ELSE BEGIN
		mdsopen,'efit01',shot
		rlim=mdsvalue('\EFIT01::TOP.RESULTS.GEQDSK:RLIM')
		zlim=mdsvalue('\EFIT01::TOP.RESULTS.GEQDSK:ZLIM')
		mdsclose,'efit01',shot	

		rlim=rlim[*,0]
		zlim=zlim[*,0]
        ENDELSE
	oplot,rlim,zlim


	;IF NOT keyword_set(oplot) THEN IF keyword_set(ps) THEN device, xsize=float(d_old.x_size)/d_old.x_px_cm,ysize=float(d_old.y_size)/d_old.y_px_cm
	
END
;+
;NAME:
;	FS_PLOT
;
;PURPOSE:
;	This procedure OPLOTS flux surfaces, presumabley onto a window where VESSEL_PLOT has been called.
;
;CALLING SEQUENCE
;	FS_PLOT,shot,time
;
;INPUTS:
;	shot:	LONG	shot number
;	time:	FLT	time point [sec]
;
;MODIFICATION HISTORY:
;	6/18/11		M.L. Reinke - modified to allow non-ANALYSIS EFIT's to be used
;	7/2/15		M.L. Reinke - modified to work with NSTX EFIT system
;-

PRO fs_plot,shot,time,num=num,color=color,style=style,thick=thick,sol=sol,nsol=nsol,rad_sol=rad_sol,invssel=invessel,tree=tree
	IF NOT keyword_set(color) THEN color=100
	IF NOT keyword_set(style) THEN style=2 ELSE IF style EQ -1 THEN style=0
	IF NOT keyword_set(thick) THEN thick=1.0
	IF NOT keyword_set(num) THEN num=10
	IF NOT keyword_set(nsol) THEN nsol=4
	IF NOT keyword_set(rad_sol) THEN rad_sol=0.905
		
	getefit,shot,time,a,g,runid=tree
	levels=make(0.0,1.0,num+1)
	levels=levels[0:num-1]		
	psin=(g.psirz-g.ssimag)/(g.ssibry-g.ssimag)
	contour,psin,g.r,g.z,levels=levels,color=color,/overplot, c_linestyle=style,thick=thick
	efit_lcfs=g.bdry
	oplot,[g.rmaxis],[g.zmaxis],psym=1,color=color,symsize=3.0

	rbdry=efit_lcfs[0,where(efit_lcfs[0,*] NE 0 AND efit_lcfs[1,*] NE 0)]
	zbdry=efit_lcfs[1,where(efit_lcfs[0,*] NE 0 AND efit_lcfs[1,*] NE 0)]
	rbdry=rbdry[0:n(rbdry)-1]
	zbdry=zbdry[0:n(zbdry)-1]
	rb_plt=[rbdry,rbdry[0]]
	zb_plt=[zbdry,zbdry[0]]
	oplot,rb_plt,zb_plt,color=color,thick=2.0


	IF keyword_set(sol) THEN BEGIN

	ENDIF

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
		sol=sol,tor=tor,edge=edge,bspline=bspline,n_s=n_s,inv=inv,col=col,tree=tree,win=win,oplot=oplot

	IF NOT keyword_set(tor) THEN notor=1 ELSE notor=0
	IF NOT keyword_set(col) THEN col=200
	IF NOT keyword_set(shot) THEN shot=135440	;shot for flux surfraces
	IF NOT keyword_set(tpt) THEN tpt =1.0		;time point for flux surfaces
	IF NOT keyword_set(thick) THEN thick=2.0
	IF NOT keyword_set(n_s) THEN n_s=100	;set number of segments to define line of sight
	IF NOT keyword_set(rzbnd) THEN rzbnd=[0.15,1.5,1.7]	;set path's [r_min,r_max,abs(z_max)] limit
	IF NOT keyword_set(nsol) THEN nsol=4
	IF NOT keyword_set(tree) THEN tree='efit01'
	
	;determine the number of views
	xpos=size(pos)	
	IF xpos[0] EQ 1 THEN n_view=1 ELSE n_view=xpos[2]

	IF NOT keyword_set(oplot) THEN BEGIN
		IF shot GT 0 THEN BEGIN
			vessel_plot,shot=shot,div=div,n=win
			fs_plot,shot,tpt,tree=tree
	        ENDIF ELSE BEGIN
			vessel_plot,div=div,/isnstxu,n=win
                    ENDELSE	
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

;find the times for which EFIT has been run
FUNCTION line_gettimes,shot,tree=tree
	IF NOT keyword_set(tree) THEN tree='efit01'
	mdsopen, tree, shot
	t=mdsvalue('\'+tree+'::TOP.RESULTS.AEQDSK:ATIME',/quiet,status=status)
	IF status THEN output=t ELSE output=-1
	mdsclose, tree, shot
	RETURN, t
END


;+
; USAGE: rmid=efit_rmid(shot,tree=tree,g=g)
;
; PURPOSE: 
;	Returns the vector of outboard midplane RMAJ (RMID) values corresponding to standard normalized psi
;	values. This routine is intended to be used with eqdsk outputs
;	or load from the tree, since NSTX doesn't have \EFIT_RMID :-(
;
; INPUTS:
;	shot	LONG	shot number
;
; OPTIONAL INPUTS:
;	tree	STRING	EFIT tree to use DEFAULT:'efit01'
;	g	STRUC	standard EFIT g-file (i.e. g=readg(shot,time,runid=tree)
;	
; OUTPUT:	FLTARR(nw,nt)	Values of RMID, the outer midplane radius for
;				which normalized psi = findgen(nw)/(nw-1).
;
; PROCEDURE: 
;	Adapted (really straight-up stolen) from EFIT_RMID from . 
;	/usr/local/cmod/codes/efit/idl/efit_rmid.pro on C-Mod
;
; Restrictions:
;	Only one time-slice is handled 
;
; MODIFICATION HISTORY:
;	M.L. Reinke	7/8/15 - stolen from the Wolfe 
;
;-

FUNCTION efit_rmid,shot,tree=tree,g=g,psin=psin,time=time
	IF NOT keyword_set(tree) THEN tree='efit01'
	IF NOT keyword_set(g) THEN BEGIN
		mdsopen,tree,shot
		psirz=mdsvalue('\'+tree+'::TOP.RESULTS.GEQDSK:PSIRZ')
		rmaxis=mdsvalue('\'+tree+'::TOP.RESULTS.GEQDSK:RMAXIS')
		zmaxis=mdsvalue('\'+tree+'::TOP.RESULTS.GEQDSK:ZMAXIS')
		ssimag=mdsvalue('\'+tree+'::TOP.RESULTS.GEQDSK:SSIMAG')
		ssibry=mdsvalue('\'+tree+'::TOP.RESULTS.GEQDSK:SSIBRY')
		xdim=mdsvalue('\'+tree+'::TOP.RESULTS.GEQDSK:XDIM')
		rgrid1=mdsvalue('\'+tree+'::TOP.RESULTS.GEQDSK:RGRID1')
		zdim=mdsvalue('\'+tree+'::TOP.RESULTS.GEQDSK:ZDIM')
		zmid=mdsvalue('\'+tree+'::TOP.RESULTS.GEQDSK:ZMID')
		time=mdsvalue('\'+tree+'::TOP.RESULTS.AEQDSK:ATIME')
		mdsclose,tree,shot
	ENDIF ELSE BEGIN
		psirz=g.psirz
		rmaxis=g.rmaxis
		zmaxis=g.zmaxis
		ssimag=g.ssimag
		ssibry=g.ssibry
		xdim=g.xdim
		rgrid1=g.rgrid1
		zdim=g.zdim
		zmid=g.zmid
		time=g.time
	ENDELSE
	ss = size(psirz)
	nw = ss(1) & nh = ss(2)
	IF ss[0] EQ 3 THEN nt=ss[3] ELSE nt=1
	rmid = fltarr(nw,nt)
	rgrid=rgrid1+xdim*findgen(nw)/(nw-1)
	zgrid=(zmid-zdim/2)+zdim*findgen(nh)/(nh-1)
	FOR it=0,nt-1 DO BEGIN
		IF keyword_set(norm) THEN psig=psirz[*,*,it] ELSE psig=(psirz[*,*,it]-ssimag[it])/(ssibry[it]-ssimag[it])
		; Set up the bicubic spline, using IMSLS interface
		wk = fltarr(2*(nw+1)*nh)
		ier=0l
		pds=0.
		c=fltarr(2,nw,2*nh)
		ibcccw,psig,rgrid,zgrid,c,bkx,bky,/idl
		jz = (where(zgrid GE zmaxis[it]))(0)
		xmax=rgrid((where(psig[*,jz] GE 1.0 AND rgrid GT rmaxis[it]))(0))
		rxx = rmaxis[it] + (xmax-rmaxis[it])*findgen(nw)/(nw-1)
		pxx = fltarr(nw)
		FOR i=1,nw-1 DO BEGIN
			ibcevw,rgrid,zgrid,c,bkx,bky,rxx[i],zmaxis[it],pds
			pxx(i)=pds
		ENDFOR
		ord = sort(pxx)
		rmid[*,it] = spline(pxx(ord),rxx(ord),findgen(nw)/(nw-1))
	ENDFOR
	psin=findgen(nw)/(nw-1)
	RETURN,rmid
END


FUNCTION line_getrmid,shot,psin=psin,tree=tree
	IF NOT keyword_set(tree) THEN tree='efit01'
	rmid=efit_rmid(shot,tree=tree,psin=psin)
	output=rotate(rmid,4)	;puts it in C-Mod format of [time,space]
	RETURN, output
END


;get EFIT output array that defines the LCFS
;[2,n_sp, n_time], so plot,[0,*,i],[1,*,i] will show you the
;LCFS it time point, i (reference t using line_gettimes)
FUNCTION line_getlcfs,shot,tree=tree
	IF NOT keyword_set(tree) THEN tree='efit01'
	mdsopen, tree, shot
	;bdry=mdsvalue('\MHD::TOP.ANALYSIS.EFIT.RESULTS.G_EQDSK:BDRY',/quiet,status=status); C-Mod node
	bdry=mdsvalue('\'+tree+'::TOP.RESULTS.GEQDSK:BDRY',/quiet,status=status)
        IF NOT status THEN BEGIN
              	;rbbs=mdsvalue('\ANALYSIS::TOP:EFIT.RESULTS.G_EQDSK:RBBBS',/quiet,status=rstat) ; C-Mod node
               	;zbbs=mdsvalue('\ANALYSIS::TOP:EFIT.RESULTS.G_EQDSK:ZBBBS',/quiet,status=zstat) ; C-Mod node
		rbbs=mdsvalue('\'+tree+'::TOP.RESULTS.GEQDSK:RBBBS',/quiet,status=rstat)
		zbbs=mdsvalue('\'+tree+'::TOP.RESULTS.GEQDSK:ZBBBS',/quiet,status=zstat)
		IF rstat AND zstat THEN BEGIN
			nbbs=n(rbbs[*,0])+1
			ntime=n(rbbs[0,*])+1
			bdry=fltarr(2,nbbs,ntime)
			FOR i=0,ntime-1 DO BEGIN
				bdry[0,*,i] = rbbs[*,i]
				bdry[1,*,i] = zbbs[*,i]
			ENDFOR
			output=bdry
		ENDIF ELSE output = -1		
	ENDIF ELSE output=bdry
	mdsclose,tree, shot
	RETURN, output
END

;get the [r,z] of the magnetic axis from EFIT, indices are time
FUNCTION line_getaxis,shot,tree=tree
	IF NOT keyword_set(tree) THEN tree='efit01'
	mdsopen, tree, shot
	;raxis=mdsvalue('\MHD::TOP.ANALYSIS.EFIT.RESULTS.G_EQDSK:RMAXIS',/quiet,status=rstatus) ; C-Mod node
	;zaxis=mdsvalue('\MHD::TOP.ANALYSIS.EFIT.RESULTS.G_EQDSK:ZMAXIS',/quiet,status=zstatus) ; C-Mod node
	raxis=mdsvalue('\'+tree+'::TOP.RESULTS.GEQDSK:RMAXIS',/quiet,status=rstatus) 
	zaxis=mdsvalue('\'+tree+'::TOP.RESULTS.GEQDSK:ZMAXIS',/quiet,status=zstatus)
	mdsclose, tree, (shot)
	IF rstatus AND zstatus THEN output=[[raxis], [zaxis]] ELSE output = -1
	RETURN, output
END

;holy total hack!  
FUNCTION efit_rz2rmid,r,z,t,shot=shot,rho=rho,tree=tree
	psin_rz=efit_rz2rho(r,z,t,shot=shot,tree=tree,/psin)
	rmid_rz=psin_rz*0.0
	rmid=efit_rmid(shot,psin=psin,time=time,tree=tree)
	FOR i=0,n(t) DO BEGIN
		dum=min(abs(time-t[i]),index)
		IF keyword_set(rho) THEN BEGIN
			irho=(rmid[*,index]-rmid[0,index])/(last(rmid[*,index])-rmid[0,index])
			rmid_rz[*,i]=interpol(irho,psin,psin_rz[*,i])
                ENDIF ELSE rmid_rz[*,i]=interpol(rmid[*,index],psin,psin_rz[*,i])
	ENDFOR
	RETURN,rmid_rz

END

PRO show_genpos_nstx,shot=shot,time=time,eps=eps,sigma=sigma,rtang=rtang,dpsi=dpsi,flat=flat,tree=tree
	IF NOT keyword_set(shot) THEN shot=135440
	IF NOT keyword_set(time) THEN time=1.0
	IF NOT keyword_set(tree) THEN tree='efit01'
	IF NOT keyword_set(eps) THEN eps=1.0
	IF NOT keyword_set(sigma) THEN sigma=0.0
	IF NOT keyword_set(rtang) THEN rtang=0.175
	IF NOT keyword_set(dpsi) THEN dpsi=1.0
	
	efit_time=line_gettimes(shot,tree=tree)
	efit_rmid=line_getrmid(shot,tree=tree)
	efit_axis=line_getaxis(shot,tree=tree)
	index=ipt(efit_time,time)
        irmid=reform(efit_rmid[index,*])
        rho=(irmid-irmid[0])/(last(irmid)-irmid[0])
	IF keyword_set(flat) THEN emiss=rho*0.0+1.0 ELSE emiss=(1-rho^2)+3*exp(-(rho-0.7)^2/(2*0.2^2))*(1-rho^1.5)
	nch=30
	pos=fltarr(4,nch)
	chnum=indgen(nch)+1
	pos[0,*]=1.7
	pos[1,*]=0.2
	pos[2,*]=rtang
	pos[3,*]=-1.0*make(-1.0*dpsi,1.0*dpsi,nch)
	;genpos_pos_reform,pos,[1.836,3.891,-2,2]
	rzbnd=[0.18,2.0,1.7]
	line_path_plots,pos,shot=shot,tpt=time,tree=tree,rzbnd=rzbnd
	br=genpos_line_br(pos,emiss,irmid,time,shot,time,rzbnd=rzbnd)
	err=randomn(seed,nch)*sigma*br
	br+=err
	vox=genpos_pos2voxel_matrix(pos,(fltarr(nch)+1.0)*(4.0*!pi),shot,tpts=time,rzbnd=rzbnd,n_rho=int(nch*1.5),rho_vec=rho_vec)
	eminv=genpos_profile_invert(br,vox,eps*max(vox),brchk=brchk,err=sigma*br)

	IF !d.name EQ 'PS' THEN BEGIN
		xsize=6.0
		ysize=5.0
		device, xsize=xsize, ysize=ysize, /inches
	ENDIF
	openwin,0
	ymax=max(brchk.mom) > max(br)
	plot,[0],[0],xr=[0,nch+1],yr=[0,ymax]*1.05,/xsty,/ysty,xtit='CH#',ytit='Brightness'
	oploterror,chnum,br,br*sigma,psym=8
	oplot,chnum,brchk.mom,color=30
	xyouts, 22	,1.5,'MEASURED'
	xyouts, 22,1.3,'BRCHK',color=30

	openwin,1
	ymax=max(eminv) > max(emiss)
	plot,[0],[0],xr=[0,1],yr=[0,ymax]*1.05,/xsty,/ysty,xtit='RHO',ytit='Emissivity'
	oploterror,rho_vec,eminv,brchk.inverr,color=30,errcolor=30
	oplot,rho,emiss
	xyouts, 0.05,2.0,'DEFINED'
	xyouts, 0.05,1.8,'RECONSTRUCTED',color=30

	stop
END


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
;	rrange:	FLTARR of the radial range DEFAULT [0.3, 1.6]
;	zrange: FLTARR of the z range DEFAULT [-1.65,1.65]
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
;       8-03-15:	M.L. Reinke - modified for defaults to be used for NSTX-U
;-

FUNCTION grid_ves,nx=nx,ny=ny,center=center,div=div,bot=bot,zrange=zrange,rrange=rrange
	IF NOT keyword_set(nx) THEN nx=10
	IF NOT keyword_set(ny) THEN ny=10
	IF keyword_set(zrange) THEN BEGIN
		zo=mean(zrange)
		delz=zrange[1]-zrange[0]
	ENDIF ELSE BEGIN
		zo=0.0
		delz=1.65*2.0
	ENDELSE
	IF keyword_set(rrange) THEN BEGIN
		ro=mean(rrange)
		delr=rrange[1]-rrange[0]
	ENDIF ELSE BEGIN
		ro=0.95
		delr=1.3
	ENDELSE

	ves_grid=genpos_grid(delr,delz,[nx,ny],xo=ro,yo=zo,center=center)
	IF keyword_set(div) THEN ves_grid=genpos_grid(0.75,0.8,[nx,ny],xo=0.47,yo=-1.5,center=center)
	IF keyword_set(bot) THEN ves_grid=genpos_grid(0.8,1.2,[nx,ny],xo=0.7,yo=-1.1,center=center)
	RETURN, ves_grid
END
