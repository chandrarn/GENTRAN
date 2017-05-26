;+
;NAME:
;	SPEC_BR
;
;PURPOSE:
;	This function is similar LINE_BR but instead of the line integrated brightness, the line-integrated
;	spectral brightness is calculated assuming the plasma emission is thermal and affected only by
;	Doppler broading and line shifts.
;
;CALLING SEQUENCE
;	result=SPEC_BR(pos,emiss,ti,vel,line,shot,t_pts)
;
;
;INPUTS:
;	pos:		FLTARR [4,m] where m is the number of views through the plasma (see PROCEDURE for 
;			references on generating a pos vector)
;	emiss:		STRUC containing the line emissivity and its spatial/temporal information
;			*.emiss	FLTARR [#r,#t] of the line emissivity (see MIST_ZZ_PROFILES)
;			*.r	FLTARR of the midplane major radii of emiss
;			*.t	FLTARR of the time points (can be a single point if *.emiss is 1D array)
;	ti:		STRUC containing the ion temperature and its spatial/temporal information
;			*.ti	FLTARR [#r,#t] of the ion temperature to be used in broadening calcs
;			*.r	FLTARR of the midplane major radii of ti
;			*.t	FLTARR of the time points (can be a single point if *.ti is 1D array)
;	vel:		STRUC containing the velocity and its spatial/temporal information
;			*.u	FLTARR [#r,#t] of the u function in incompressibility [m/s/T]
;			*.w	FLTARR [#r,#t] of the w function in incompressibility [1/s]
;			*.r	FLTARR of the midplane major radii of vel
;			*.t	FLTARR of the time points (can be a single point if *.vel is 1D array)	
;	line:		STRUC containing information about the line to be calculated
;			*.lam_o		FLT the rest wavelength of the line
;			*.del_lam	FLT used in the  wavelength interval [lam_o-del_lam, lam_o+del_lam]
;			*.mass		FLT the mass of ion in AMU (see READ_ATOMIC_MASS('XX')
;	shot:		LONINT shot number to reference flux surface information
;	t_pts:		FLTARR of length n of the desired time points where output is desired.  If input profiles
;			are at a single time point, make t_pts this time point as well.
;
;OPTIONAL INPUTS:
;	alpha:		FLT Zeeman splitting coefficient (see BLOM PPCF v 44 pg 1229)	
;	n_s:		FLT the number of points along the line of sight DEFAULT: 75
;	n_lam:		INT the number of spectral bins	DEFAULT: 100
;	rzbnd:		FLTARR boundary [r_in, r_out, abs(z)] for line of sight DEFAULT = [0.44,1.0,0.5]
;	ct_main:	INT the primary color table DEFAULT: 12
;	ct_pol:		INT the color table of the poloidal contour plots DEFAULT: 39
;	cct:		INT color table for the contour on GENPLT DEFAULT: 39
;	ngrid:		INT the number of points in (R,Z) for emissivity contours (see LINE_EMCONT) DEFAULT: 50
;	nlevels:	INT the number of levels for poloidal contour plots DEFAULT: 100
;	plotlab: 	STRUC of labels to be used for emissivity when /plots is invoked
;				*.line 	STR line label
;				*.wl	FLT wavelength 
;
;	***most GENPLT optional inputs are carried through.  See help, /rout, names='SPEC_BR' for list
;
;KEYWORD PARAMETERS:
;	lcfs:		/lcfs cuts all profiles off at the lcfs.  Better for core dominated profiles but invalid for
;			edge emission.
;	debug:		/debug stops at various points in SPEC_BR
;	gdebug:		/gdebug sends a debug command to GENPLT
;	plots:		/plots displays a large number of informative plots (see OPTIONAL OUTPUTS below)
;
;OUTPUTS:
;	result:		STRUC of the line integrated spectral brightness for each position and time point
;			*.br FLTARR [m,n_lam,n] for POS[*,i] at time point j the spectra is br[i,*,j]
;			*.lam FLTARR of length n_lam that is the wavelength scaling
;
;
;OPTIONAL OUTPUTS
;	emlam:		STRUC containing the spectral emissivity as a function of parameter for each
;			position and time point
;			*.emlam FLTARR [n_s, n_lam, m, n] of the spectral emissivity 
;			*.s FLTARR of length n_s that is the parameter scaling
;			*.lam FLTARR of length n_lam that is the wavelength scaling
;
;	If /plots is used a number of graphics will be output.
;		window      -           description
;	----------------------------------------------------------------------------------
;		 13		A contour plot of the spectral emissivity along the line of sight. (IF WIN=1 is
;				invoked for GENPLT then this will be window 16)
;		 14		Slices of the spectral emissivity vs parameter - l for various wavelengths.
;				(window,17 for WIN=1)
;		 15		Slices of the spectral emissivity vs. wavelength for various values of parameter
;				l. (window,18 for WIN=1)
;		 21		A contour plot of the line emissivity with the poloidal projection of the line
;				of sight.
;		 22		A contour plot of the ion temperature with the poloidal projection of the line
;				of sight.
;		 23		Line emissivity, line shift delta=(v dot b)/c and and line width
;				are plotted along the line of sight vs parameter l.
;		 24		The output of SPEC_BR, the line integrated spectral brightness is plotted
;				versus the wavelength scale
;
;	IF /plots and /ps are both used then windows 1,2,19,20 will be in /home/username/idl/plots/spec_br.ps
;	and the contour plots will be in /home/username/idl/plots/spec_br_*.ps where * means the typical
;	GENPLT naming (see HELPME -> DATAPLOTS)
;
;PROCEDURE
;	For information on specfying POS vectors, see the documentation for LINE_BR or see ML Reinke for
;	a copy of 'General Equations for Radioumetry in Tokamak Plasmas'.
;	
;	This function assumes a Gaussian line shape controlled by parameters that are flux functions.  In
;	general, MIST will be used to generate the line emissivity (see ZZ_UNCALIBRATED_MIST) and if
;	Te is assumed equal to Ti, the Te profile returned from MIST can be used.  The rotation
;	is specified seperately and will not be consistent with the impurity reconstructions.  This makes
;	this analysis invalid for arbitrarily high rotation speeds that will affect the impurity density
;
;	Currently only one line can be simulated at a time and the spectral brightness for close
;	lines can found, inefficiently, through multiple calls to SPEC_BR.	
;
;	Further upgrades to *_BR suite of codes will move away from assuming flux functions as well as
;	allow for non-thermal emission via a generalized spectral input instead of assuming Gaussians.
;
;MODIFICATION HISTORY:
;	Written by:	ML Reinke - June 15, 2006
;	6-16-06:	ML Reinke - tweaked the /ps output and redid the EMLAM output
;	8-01-06:	ML Reinke - updated to use new velocity structure which uses w(psi) and u(psi).
;				    removed normalized output plots and made a triple plot of them
;				    adjusted 
;-



FUNCTION spec_br,pos,emiss,ti,vel,line,shot,t_pts,n_lam=n_lam,n_s=n_s,ngrid=ngrid,nlevels=nlevels,verb=verb,$
		debug=debug,plots=plots,plotlab=plotlab,ct_main=ct_main,ct_pol=ct_pol,emlam=emlam,$
		io=io,jo=jo,ir=ir,jr=jr,maxpt=maxpt,minpt=minpt,ncntrs=ncntrs,path=path,suffix=suffix,prefix=prefix,$
		gdebug=gdebug,cct=cct,pct=pct,win=win,lcfs=lcfs,alpha=alpha,emcos=emcos,emsin=emsin

	IF NOT keyword_set(n_s) THEN n_s=75			;set number of segments to define line of sight
	IF NOT keyword_set(n_lam) THEN n_lam=100 		;set number of spectral bins
	IF NOT keyword_set(rzbnd) THEN rzbnd=[0.44,1.0,0.5]	;set path's [r_min,r_max,abs(z_max)] limit
	IF NOT keyword_set(ct_main) THEN ct_main=12		;set main color table
	IF NOT keyword_set(ct_pol) THEN ct_pol=39		;set contour color table for poloidal plot of EMISS
	IF NOT keyword_set(ngrid) THEN ngrid=80			;griding parameter for LINE_EMCONT
	IF NOT keyword_set(nlevels) THEN nlevels=200		;number of levels in poloidal contour plots
	IF NOT keyword_set(prefix) THEN prefix='spec_br_'	;set default prefix to GENPLT
	IF NOT keyword_set(cct) THEN cct=39
	IF NOT keyword_set(alpha) THEN alpha=0.0
	
	dev=!d.name
	IF !d.name EQ 'PS' THEN ps=1 ELSE ps=0

	;define physical constants
	c=3.0e8  		;speed of light
	angst=1.0e-10 		;conversion of lambda to meters
	e=1.60e-19		;conversion for eV -> J
	mconv=1.66e-27		;conversion for amu -> kg
		
	;determine the number of views
	xpos=size(pos)	
	IF xpos[0] EQ 1 THEN n_view=1 ELSE n_view=xpos[2]

	;set default asymmtery to be zero
	IF NOT keyword_set(emcos) THEN emcos=0.0*emiss.emiss
	IF NOT keyword_set(emsin) THEN emsin=0.0*emiss.emiss

	;determine the number of time points
	xtpts=size(t_pts)
	IF xtpts[0] EQ 1 OR xtpts[0] EQ 0 THEN n_tpts=1 ELSE n_tpts=xtpts[1]

	;determine the minimum radius point
	min_r=min(emiss.r) < min(ti.r) < min(vel.r)

	;determine the spectral range
	lam=make(line.lam_o-line.del_lam,line.lam_o+line.del_lam,n_lam)	

	;intialize output spectral brightness array (number of views * number of spectral bins *  number of time points)
	br=fltarr(n_view, n_lam, n_tpts)
	emlam=fltarr(n_s,n_lam,n_view,n_tpts)

	;load EFIT data
	efit_time=line_gettimes(shot)
	efit_lcfs=line_getlcfs(shot)
	efit_axis=line_getaxis(shot)
	
	;debug information
	IF keyword_set(verb) THEN print, 'n_s = '+num2str(n_s)+'  n_view = '+num2str(n_view)+'  n_lam = '$
			+num2str(n_lam)+'  n_tpts = '+num2str(n_tpts)
	IF keyword_set(debug) THEN stop

	;start filling the spectral brightness array cycling through each view
	FOR i=0,n_view-1 DO BEGIN
		IF keyword_set(verb) THEN print, ' pos - i = '+num2str(i,1)

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
		;there is now array of points (s) that go from point 1 [R1,Z1] to the closest limiting
		;surface and is divided into n_s parts.

		;form unit vectors
		lhat_th=atan((1.0-s)/pos[2,i]*sqrt(pos[0,i]^2-pos[2,i]^2))
		lhat_phi=cos(lhat_th)*cos(pos[3,i])
		lhat_r=-sin(lhat_th)*cos(pos[3,i])
		lhat_z=fltarr(n_s)-sin(pos[3,i])

		bfield=line_bfield(pos[*,i],s,shot,t_pts=t_pts)		

		;cycle through all the time points
		FOR j=0,n_tpts-1 DO BEGIN
			IF keyword_set(verb) THEN print, ' j = '+num2str(j,1)
			;intialize spectral emissivity array (number of points along view * number of spectral bins)
			emlam_s=fltarr(n_s,n_lam)
			efit_i=ipt(efit_time,t_pts[j])
			IF efit_i NE -1 THEN BEGIN
				axis=[efit_axis[efit_i,0],efit_axis[efit_i,1]]
				IF keyword_set(lcfs) THEN BEGIN
					;load relevant LCFS boundary and prepare for format needed for LINE_INCLFS
					rbdry=efit_lcfs[0,where(efit_lcfs[0,*,efit_i] NE 0 AND efit_lcfs[1,*,efit_i] NE 0),efit_i]
					zbdry=efit_lcfs[1,where(efit_lcfs[0,*,efit_i] NE 0 AND efit_lcfs[1,*,efit_i] NE 0),efit_i]
					rbdry=rbdry[0:n(rbdry)-1]
					zbdry=zbdry[0:n(zbdry)-1]
				
					;define an array of 1's (inside) and 0's (outside) where s is in/out LCFS
					s_inlcfs=intarr(n_s)
					FOR k=0,n_s-1 DO s_inlcfs[k]=line_inlcfs(rbdry,zbdry,axis,$
						[line_r(s[k],pos[*,i]), line_z(s[k],pos[*,i])])
				ENDIF ELSE s_inlcfs=intarr(n_s)+1
				goodpts=where(s_inlcfs NE 0)
				IF keyword_set(verb) THEN print, goodpts					
				
				;if any points are in the LCFS fill emlam_s
				IF goodpts[0] NE -1 THEN BEGIN
					;print verb info
					IF keyword_set(verb) THEN BEGIN
						print, 'Good R'
						print, line_r(s[goodpts], pos[*,i])
						print, 'Good Z'
						print, line_z(s[goodpts], pos[*,i])
						print, 'Time Point'
						print, efit_time[efit_i]
					ENDIF
					;determine the set of (R,Z) points for the view
					good_r=line_r(s[goodpts], pos[*,i])
					good_z=line_z(s[goodpts], pos[*,i])
					good_s=s[goodpts]

					;find their equivilent midplane radii
					r_on_mid=efit_rz2rmid(good_r, good_z,efit_time[efit_i], shot=shot)
					cth=(good_r-axis[0])/sqrt((good_r-axis[0])^2+(good_z-axis[1])^2)
					sth=(good_z-axis[1])/sqrt((good_r-axis[0])^2+(good_z-axis[1])^2)		
					IF keyword_set(verb) THEN BEGIN
						print, 'Good RMID'
						print, r_on_mid
					ENDIF

					;if any, set values of r_on_mid that are < min_r  to min_r
					;this is generally necessary in the core when min_r is like 0.69-0.70
					tmp=where(r_on_mid LT min_r)
					IF tmp[0] NE -1 THEN r_on_mid[tmp]=min_r
					
					;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
					;GENERATE SPECTRAL EMISSIVITY;
					;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

					;find total line emissivity as a function of parameter
					;-------------------------------------------------------------------------------------------
			
					;use IDL function INTERPOLATE
					IF n(emiss.t) EQ 0 THEN BEGIN
						emiss_rad=interpol(emiss.emiss,emiss.r,r_on_mid)
						emiss_cos=interpol(emcos,emiss.r,r_on_mid)
						emiss_sin=interpol(emsin,emiss.r,r_on_mid)
						emiss_good=emiss_rad+emiss_cos*cth+emiss_sin*sth
					ENDIF ELSE BEGIN
						r_interp=interp_vec_reform(emiss.r, r_on_mid)
						t_interp=interp_vec_reform(emiss.t,fltarr(n(r_on_mid)+1)+efit_time[efit_i])
						emiss_good=interpolate(emiss.emiss,r_interp,t_interp, missing=1.0e-4)
					ENDELSE						
					line_emiss=fltarr(n_s)
					line_emiss[goodpts]=emiss_good
					IF keyword_set(verb) THEN BEGIN
						print, 'Good LINE_EMISS'
						print, line_emiss
					ENDIF
					;--------------------------------------------------------------------------------------------

					;find line width as a function of parameter
					;--------------------------------------------------------------------------------------------

					;use IDL function INTERPOLATE
					IF n(ti.t) EQ 0 THEN BEGIN
						ti_good=interpol(ti.ti,ti.r,r_on_mid)
					ENDIF ELSE BEGIN
						r_interp=interp_vec_reform(ti.r, r_on_mid)
						t_interp=interp_vec_reform(ti.t,fltarr(n(r_on_mid)+1)+efit_time[efit_i])
						ti_good=interpolate(ti.ti,r_interp,t_interp, missing=1.0e-4) 
					ENDELSE
					sigma=fltarr(n_s)+1.0e-10 ;to prevent divide by zero error in gaussian
					sigma[goodpts]=line.lam_o*sqrt(ti_good*1.0e3*e/(line.mass*mconv))/c   ;in units of lam_o
					IF keyword_set(verb) THEN BEGIN
						print, 'Good SIGMA'
						print, sigma
					ENDIF
					;--------------------------------------------------------------------------------------------


					;find line shift as a function of parameter
					;--------------------------------------------------------------------------------------------
					;use IDL function INTERPOLATE
					IF n(vel.t) EQ 0 THEN BEGIN
						w_good=interpol(vel.w,vel.r,r_on_mid)
						u_good=interpol(vel.u,vel.r,r_on_mid)
					ENDIF ELSE BEGIN
						r_interp=interp_vec_reform(vel.r, r_on_mid)
						t_interp=interp_vec_reform(vel.t,fltarr(n(r_on_mid)+1)+efit_time[efit_i])
						w_good=interpolate(vel.w,r_interp,t_interp, missing=1.0e-4)
						u_good=interpolate(vel.u,r_interp,t_interp, missing=1.0e-4)
					ENDELSE
					delta=fltarr(n_s)
				        delta[goodpts]=1.0/c*(u_good*(lhat_r[goodpts]*bfield.br[goodpts,j]+lhat_z[goodpts]*bfield.bz[goodpts,j]$
						+lhat_phi[goodpts]*bfield.bphi[goodpts,j])+w_good*lhat_phi[goodpts]*good_r)	
					IF keyword_set(verb) THEN BEGIN
						print, 'Good DELTA'
						print, delta
					ENDIF
					;--------------------------------------------------------------------------------------------
				
					
					;the scaling of s is arbitrary with no physical signifigance except that
					;s=0 is pt1 and s=1 is pt2.  If the s values are multiplied by the distance
					;between 1 and 2 then the qunatitative units are applied to the scaling.
					ipos=reform(pos[*,i])
					scale=line_scale(ipos)

					mag_b=sqrt(bfield.br[*,j]^2+bfield.bz[*,j]^2+bfield.bphi[*,j]^2)
					lhat_dot_bhat=(lhat_r*bfield.br[*,j]+lhat_z*bfield.bz[*,j]+lhat_phi*bfield.bphi[*,j])/mag_b
					;form the spectral emissivity as a function of parameter
					FOR k=0,n_lam-1 DO BEGIN
						emlam_sigma_plus=(1.0+lhat_dot_bhat^2)*line_emiss/(4.0*sqrt(2.0*!pi)*sigma)*$
							exp(-(lam[k]-line.lam_o*(1.0+delta)+alpha*mag_b/2.0)^2/(2.0*sigma^2))
						emlam_sigma_minus=(1.0+lhat_dot_bhat^2)*line_emiss/(4.0*sqrt(2.0*!pi)*sigma)*$
							exp(-(lam[k]-line.lam_o*(1.0+delta)-alpha*mag_b/2.0)^2/(2.0*sigma^2))
						emlam_pi=(1.0-lhat_dot_bhat^2)*line_emiss/(2.0*sqrt(2.0*!pi)*sigma)*$
							exp(-(lam[k]-line.lam_o*(1.0+delta))^2/(2.0*sigma^2))
						emlam_s[*,k]=emlam_sigma_plus+emlam_sigma_minus+emlam_pi
					ENDFOR
					emlam[*,*,i,j]=emlam_s

					;integrate over parameter to get the line-integrated spectral brightness (include scaling)
					FOR k=0,n_lam-1 DO br[i,k,j]=scale*int_tabulated(s, reform(emlam_s[*,k]))
		
					IF keyword_set(debug) THEN stop

					;done with the math, but if PLOTS are on we have a ways to go
					IF keyword_set(plots) THEN BEGIN
						IF keyword_set(lcfs) THEN sol=0 ELSE sol=1
						;plot the poloidal projection of a line of sight over emissivity
						;-------------------------------------------------------------------------------------
						vessel_plot,edge=edge,n=21,d_old=d_old,ps=ps
					

						;find and plot normalized 2D line emissivity countours
						loadct,ct_pol,/silent
						emcon=line_emcont(emiss.emiss,emiss.r,shot,efit_time[efit_i],ngrid=50,sol=sol)
						contour,emcon.em/max(emcon.em),emcon.r,emcon.z,/irregular,nlevels=nlevels,$
							/fill,/overplot,levels=make(0.0,1.0,nlevels)
						loadct,ct_main,/silent
						fs_plot,shot,t_pts[j]
						IF keyword_set(plotlab) THEN BEGIN
							xyouts,0.8,0.40,'Line: '+plotlab.line+' '+plotlab.wl,charsize=1.3
						ENDIF
	
						;plot the line of sight on the polodial plane
						oplot,line_r(s, pos[*,i]),line_z(s,pos[*,i]),color=200,thick=2.0
						
						;label information about the view
						xyouts,0.8,0.5, 'POS = ['+num2str(ipos[0],dp=2)+','+num2str(ipos[1],dp=2)+','+$
							num2str(ipos[2],dp=2)+','+num2str(ipos[3],dp=2)+']',charsize=1.3
						xyouts, 0.75,0.45, 'LINE EMISSIVITY',chars=2.0
						;plot a dotted line at the flux surface of the maximum of emiss
						;fs=line_getfs(shot,efit_time[efit_i],emiss.r[maxloc(emiss.emiss)])
						;fs_type=size(fs, /type)
						;IF fs_type EQ 8 THEN oplot, fs.r,fs.z, linestyle=2, color=100, thick=3.0
						;-------------------------------------------------------------------------------------


						;plot the poloidal projection of a line of sight over ion temperature
						;-------------------------------------------------------------------------------------
						vessel_plot,edge=edge,n=22,ps=ps

						;find and plot normalized 2D Ti countours
						loadct,ct_pol,/silent
						ticon=line_emcont(ti.ti,ti.r,shot,efit_time[efit_i],ngrid=50,sol=sol)
						contour,ticon.em/max(ticon.em),ticon.r,ticon.z,/irregular,nlevels=nlevels,/fill,/overplot,$
							levels=make(0.0,1.0,nlevels)
						loadct,ct_main,/silent
						fs_plot,shot,t_pts[j]

						;plot the line of sight on the polodial plane
						oplot,line_r(s, pos[*,i]),line_z(s,pos[*,i]),color=200,thick=2.0
						
						;label information about the view
						xyouts,0.8,0.5, 'POS = ['+num2str(ipos[0],dp=2)+','+num2str(ipos[1],dp=2)+','+$
							num2str(ipos[2],dp=2)+','+num2str(ipos[3],dp=2)+']',charsize=1.3
						xyouts, 0.75,0.45, 'ION TEMPERATURE',chars=2.0
						;-------------------------------------------------------------------------------------

						;Plot normalized line emissivity,ion temp,velocity*dotprod along the path.
						;-------------------------------------------------------------------------------------
						IF keyword_set(ps) THEN BEGIN
							xsize=6.0
							ysize=6.0*1000/850.0
							ls=0.5
						ENDIF ELSE BEGIN
							xsize=700.0
							ysize=950.0
							ls=1.0
						ENDELSE
						!p.multi=[0,0,3]
						openwin,23,xsize=xsize,ysize=ysize
						IF keyword_set(ps) THEN device, xsize=xsize, ysize=ysize, /inches
						xtit='Parameter - l'
						tit='Profiles Along Line of Sight'
						xr=[min(good_s),max(good_s)]
						plot,good_s,line_emiss[goodpts],xr=xr,/xsty,tit=tit,chars=3.0*ls,ytit='Line Emiss',$
							thick=2.0*ls
						plot,good_s,delta[goodpts],xr=xr,/xsty,chars=3.0*ls,ytit='v/c dot lhat ('+n2g('delta')+')',$
							thick=2.0*ls
						plot,good_s,sigma[goodpts],xr=xr,/xsty,chars=3.0*ls,ytit='Line Width ('+n2g('sigma')+')',$
							thick=2.0*ls,xtit=xtit
						!p.multi=0
						IF keyword_set(ps) THEN device, xsize=float(d_old.x_size)/d_old.x_px_cm,ysize=$
							float(d_old.y_size)/d_old.y_px_cm
						;-------------------------------------------------------------------------------------

						;Plot line integrated spectral brightness.
						;-------------------------------------------------------------------------------------
						openwin,24
						plot,lam,br[i,*,j],thick=2.0,chars=1.3,xtit='Wavelength',ytit='Spec. Brightness',$
							tit='Line-Integrated Spectral Brightness'
						;-------------------------------------------------------------------------------------

						IF keyword_set(ps) THEN BEGIN
							device, /close
							spawn, 'cp idl.ps /home/'+logname()+'/idl/plots/spec_br.ps'
						ENDIF
							
						;plot contour plot of the spectral emissivity along the line of sight
						;-------------------------------------------------------------------------------------
						labels={ilab:'Parameter - l', jlab:'Wavelength', klab:'Line Emiss',$
							ctit:'Spectral Emissivity Contours',jtit:'',itit:'Spectra Along Line of Sight'}
						genplt,emlam_s,s,lam,io=io,jo=jo,ir=ir,jr=jr,maxpt=maxpt,minpt=minpt,ncntrs=ncntrs,$
							labels=labels,path=path,suffix=suffix,prefix=prefix,dp=dp,debug=gdebug,cct=cct,$
							pct=pct,ps=ps,win=win
						;-------------------------------------------------------------------------------------


						IF NOT keyword_set(ps) AND i LT n_view-1 THEN BEGIN
							print, 'Press ENTER key to continue'
							inputcommand=''
							read, inputcommand		
						ENDIF
					ENDIF
					
				ENDIF 
			ENDIF ELSE br[i,*,j]=-1	;value out of efit time range	
		ENDFOR
		
	ENDFOR	
	output={br:br, lam:lam}
	emlam={emlam:emlam,s:s,lam:lam}

	IF keyword_set(debug) THEN stop
	RETURN, output
END


PRO test_spec_br,pos=pos,n_s=n_s,out=out,debug=debug
	IF NOT keyword_set(pos) THEN pos=[1.0,0.0,0.68,0.0]
	;mist_zz_profiles,18,1051202015,1.0,/h_mode,pro_out=out,/local,/fit
	zz_uncalibrated_mist,18,1051202015,1.0,0.1,lam,q_lam,radpts,em,csden,/h_mode,/fit,te_profile=te,/local,psf_te=1.25
	t=[1.0]
	ti={ti:te/1.25, r:radpts, t:t}
	emiss={emiss:reform(em[0,*]), r:radpts, t:t}
	vel={u:fltarr(51)+0.0e5,w:fltarr(51)+0.0e5, r:radpts, t:t}
	line={lam_o:lam[0],del_lam:0.02,mass:read_atomic_mass('Ar')}
	
	out=spec_br(pos,emiss,ti,vel,line,1051202015,t,/plots,n_s=n_s)
	IF keyword_set(debug) THEN stop
END



FUNCTION lineint_spectrum,pos,emiss,r,shot,time,n_s=n_s, verb=verb,debug=debug,plots=plots,ps=ps,plotlab=plotlab

	;emiss=[nlam,nr]

	IF NOT keyword_set(n_s) THEN n_s=75			;set number of segments to define line of sight
	IF NOT keyword_set(rzbnd) THEN rzbnd=[0.44,1.0,0.5]	;set path's [r_min,r_max,abs(z_max)] limit
	
	;determine the number of views
	xpos=size(pos)	
	IF xpos[0] EQ 1 THEN n_view=1 ELSE n_view=xpos[2]
	
	;determin the number of wavelengths
	xem=size(emiss)
	n_lam=xem[1]
	
	;intialize output brightness array (number of views x number of time points)
	br=fltarr(n_lam, n_view)
	br_tot=fltarr(n_view)

	;load EFIT data
	efit_time=line_gettimes(shot)
	efit_lcfs=line_getlcfs(shot)
	efit_axis=line_getaxis(shot)

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

		IF keyword_set(verb) THEN print, ' j = '+num2str(j,1)
		efit_i=ipt(efit_time,time)
		IF efit_i[0] EQ -1 THEN RETURN,-1

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
		goodpts=where(s_inlcfs NE 0)
		IF keyword_set(verb) THEN print, goodpts
		emiss_s=fltarr(n_s)	;initialize emissivity as a function of path length array
		
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
			r_on_mid=efit_rz2rmid(good_r,good_z,efit_time[efit_i],shot=shot)
			IF keyword_set(verb) THEN BEGIN
				print, 'Good RMID'
				print, r_on_mid
			ENDIF

			ipos=reform(pos[*,i])	;the scaling of l is arbitrary with no physical signifigance except that
						;l=0 is pt1 and l=1 is pt2 if the l values are multiplied by the distance
						;between 1 and 2 then the qunatitative units are applied to the scaling.

			lscale=sqrt((ipos[0]^2-ipos[2]^2)*(1.0+(tan(ipos[3]))^2))
			FOR j=0,n_lam-1 DO BEGIN
				emiss_s[goodpts]=interpol(reform(emiss[j,*]),r,r_on_mid)
				br[j,i]=int_tabulated(s, emiss_s)*lscale
			ENDFOR
		ENDIF ELSE br[*,i]=0.0		
	ENDFOR	
	output=br
	IF keyword_set(debug) THEN stop
	RETURN, output
END