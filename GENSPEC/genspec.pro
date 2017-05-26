;+
;NAME:
;	_GENSPEC
;
;PURPOSE:
;	The GENSPEC_* suite of procedures and functions are the
;	application of formulas and methods in 'General Equations for
;	Radiometry in Tokamak Plasmas' for spectral emissivity.  There
;	is a substantial draw on the lower level upon GENPOS codes.
;
;	Version 0.5 - testing of inversions
;	Version 0.9 - updated to use incompressibility in velocity (06/07)
;	Version 1.0 - added the ability to use a matrix inversion
;                     technique (07/07)	
;	Version 1.1 - upgraded GENSPEC_GPV2SPECTRA to use new input
;                     formats and multiple lines.  Added the ability
;                     to determine the u-profile and use edge zeros in
;                     GENSPEC_MATRIX_INVERT (08/07)
;	Version 1.2 - various tweaks and the addition of
;                     GENSPEC_SHELLINVERT
;.
;MODIFICATION HISTORY:
;	Version 1.2
;	ML Reinke 
;	February 22nd, 2008
;-


;+
;NAME:
;	GENSPEC_GPV2SPECTRA
;
;PURPOSE:
;	This function takes a view described by a GPV and LHAT from GENPOS_VOL_COEFS and forms the line integrated spectral power.  
;	The Doppler broadened/shifted spectral emissivity at anypoint in the plasma is formed by assuming input Ti, vel 
;	and total emissivity profiles are flux functions.
;
;CALLING SEQUENCE:
;	result=GENSPEC_GPV2SPECTRA(gridpts,gpv,lhat,emiss,ti,vel,line,shot)
;
;INPUTS:
;	gridpts:	STRUC of points describing the grid.  This is a GENPOS_GRID output using /center
;	gpv: 		STRUC containing gpv for each upos for each detector only for non-zero elements
;			*.d0	gpv and spatial information for detector 0
;				*.tmp	LONARR 	[#non_zero] of array elements (gridpts.pnts[*,gpv.tmp])
;				*.gpv	FLTARR 	[#non_zero,l] of gpv values [1/m^3]
;			*.d1	gpv and spatial information for detector 1
;			  .			 .
;			  .			
;			*.dn	gpv and spatial information for detector n
;	lhat: 		STRUC containing lhat structures for each UPOS for each detector only for non-zero elements
;			*.d0	lhat and spatial information for detector 0	
;				*.tmp	LONARR 	[#non_zero] of array elements (gridpts.pnts[*,gpv.tmp])
;				*.lhat	STRUC	of lhat data at each point in *.tmp
;					*.lr	STRUC of unit vectors in the radial (major radius) direction (in and out)
;						*.in	FLTARR [#non_zero] of the average radial unit vector going in
;						*.dvin	FLTARR [#non_zero] of the volume weightings going in
;						*.out	FLTARR [#non_zero] of the average radial unit vector coming out
;						*.dvout	FLTARR [#non_zero] of the volume weightings coming out
;					*.lphi 	FLTARR [#non_zero] average unit vector in the polar angle (toroidal) direction
;					*.lz	FLTARR [#non_zero] average unit vector in the vertical (vertical) direction
;	emiss:		STRUC containing the line emissivity and its spatial/temporal information
;			*.emiss	FLTARR [#r,#lines,#t] of the line emissivity [X/vol] (see MIST_ZZ_PROFILES)
;			*.r	FLTARR of the midplane major radii of emiss
;			*.t	FLTARR of the time points (can be a single point if *.emiss is 1D array)
;	ti:		STRUC containing the ion temperature and its spatial/temporal information
;			*.ti	FLTARR [#r,#t] of the ion temperature [keV] to be used in broadening calcs
;			*.r	FLTARR of the midplane major radii of ti
;			*.t	FLTARR of the time points (can be a single point if *.ti is 1D array)
;	vel:		STRUC containing the velocity and its spatial/temporal information
;			*.u	FLTARR [#r,#t] of the u function in incompressibility [m/s/T]
;			*.w	FLTARR [#r,#t] of the w function in incompressibility [1/s]
;			*.r	FLTARR of the midplane major radii of vel
;			*.t	FLTARR of the time points (can be a single point if *.vel is 1D array)	
;	line:		STRUC containing information about the line to be calculated
;			*.lam_o		FLTARR [#lines] the rest wavelength of the line [units set wavelength scale]
;			*.del_lam	FLT used in the  wavelength interval [min(lam_o)-del_lam, max(lam_o)+del_lam]
;			*.mass		FLT the mass of ion in AMU (see READ_ATOMIC_MASS('XX')
;	shot:		LONINT shot number
;
;OPTIONAL INPUTS:
;	gpv:		FLATRR of the gpv value at each (R,Z) from					***/NORAW FORMAT***
;	lhat:		STRUC of the average unit vectors in (r,phi,z)					***/NORAW FORMAT***
;			*.lr	STRUC sub structure for the radial unit vector data		
;				*.in	FLTARR [n, m] of the average radial unit vector going in
;				*.dvin	FLTARR [n, m] of the volume weightings going in
;				*.out	FLTARR [n, m] of the average radial unit vector coming out
;				*.dvout	FLTARR [n, m] of the volume weightings coming out
;			*.lphi	FLTARR [n, m] of the average toroidal unit vector
;			*.lz	FLTARR [n, m] of the average vertical unit vector

;	t_pts:		FLTARR of length t time points for which output data is calculated DEFAULT: emiss.t
;	n_lam:		INT the number of spectral bins DEFAULT: 100
;	grid2rmid:	FLTARR of the output of GENPOS_GRID2RMID(gridpts, shot, tpts=t_pts).  This can be input to save
;				time if running GENPOS_GPV2SPECTRA for multiple detectors
;	missing:	FLT value emiss,ti and vel for points outside the LCFS DEFAULT: 1.0e-4 (non-zero because of division error)
;	interp:		STRUC of ti, emiss and vel interpolated points.  Again, this data would be the same for multiple
;				detectors.  See OPTIONAL OUTPUTS.
;
;KEYWORD PARAMETERS:
;	debug:		/debug stops the code before the RETURN statement
;	noraw:		/noraw will process the gpv and lhat inputs assuming not "raw" format to increase speed.
;
;OUTPUTS:
;	result:		STRUC of line integrated spectral brightness
;			*.plam 		FLTARR of size n_lam x t that is the total line integrated spectral power summed over all upos
;						at the time points given by t_pts.  Units are of [emiss]*[raw_gpv]/[lam_o] such as W/nm 
;			*.lam 		FLTARR of length n_lam that is the wavelength scaling 
;
;			Error infromation is encoded in the output if GENPOS_GPV2SPECTRA fails.
;			result 	= -1	error in Ti interpolation (usually demand a time point out of bounds of the inputs)
;				= -2	error in emissivity interpolation 
;				= -3	error in velocity interpolations 	
;
;OPTIONAL OUTPUTS
;	interp:		STRUC of the ti, emiss and vel data at each of the (R,Z) locations on the plasma grid.  If this struc is not
;				input then these will be calculated from the ti, emiss, vel inputs.  
;				*.ti 	FLTARR of size n x t of Ti on each pixel point in ves_grid for each time point in t_pts
;				*.em 	FLTARR of size n x t of emissivity
;				*.u	FLTARR of size n x t of the u function in the velocity equation
;				*.w	FLTARR of size n x t of the u function in the velocity equation				
;	gpvspec:	FLTARR of size n x n_lam x t of the spectral emissivity for each plasma pixel as seen by the detector
;				described by upos.
;
;PROCEDURE:
;	The depencedy of the spectral emissivity on line of sight due to Doppler shift, delta=dot(v,l_hat), prevents the
;	spectral emissivity to be expanded onto a grid for general usage like GRID_PROFILE does for the total emissivity.  Instead, the
;	spectral power (W/nm or photons/s/nm for example) deposited on the detector must be calculated by using information for each
;	sub-pixel in upos before being summed.  This is why the raw output of GENPOS_VOL_COEFS must be input.  The equations
;	are listed in 'General Equations for Radiometery in Tokamak Plasmas' and more discussion is presented.
;
;	The emiss,ti,vel and line input STRUC can be easily generated using MIST_GENPOS_INPUT_DATA.  If this function is called
;	multiple times, it should be run first to fill the optional outputs of grid2rmid and interp.  These can then be put
;	back into subsequent calls for a large speed increase.
;
;	After 6-05-07, the way GENPOS_GPV2SPECTRA calculates the doppler shift has been changed.  Now the incompressible
;	velocity relation v = u*B + w*R is used.  See section XX in 'General Equations for Radiometery in Tokamak Plasmas' for more discussion
;
;MODIFICATION HISTORY:
;	Written by:	ML Reinke - 8-12-06
;	8-18-06:	ML Reinke - overhauled the loop to fill the data to increase the speed (about x10 increase).  I also noticed
;				that the profile interpolations will not work correctly if more than one time point is called.  Since
;				this program is linked with MIST its not a concern right now.
;	8-20-06:	ML Reinke - fixed the time dependency of the profile interpolations and added the missing optional input
;	8-21-06:	ML Reinke - renamed to GENSPEC_GPV2SPECTRA
;	6-05-07:	ML Reinke - reformatted the velocity input structure and the spectra calculation to use u and w from
;                                   incompressibility equation (see procedure)
;	8-14-07:	ML Reinke - overhauled the code to use the new format of raw GPV and LHAT from GENPOS_VOL_COEFS.  Changed input upos to lhat.
;	8-15-07:	ML Reinke - allowed code to accept emiss and line structures that include inputs for multiple lines.  
;	8-15-07:	ML Reinke - allowed code to accept GPV and LHAT data for a single view that is not raw (keyword noraw)
;-

FUNCTION genspec_gpv2spectra,gridpts,gpv,lhat,emiss,ti,vel,line,shot,t_pts=t_pts,n_lam=n_lam,grid2rmid=grid2rmid,interp=interp,$
		quiet=quiet,debug=debug,gpvspec=gpvspec,missing=missing,noraw=noraw

	IF NOT keyword_set(n_lam) THEN n_lam=100 							;set number of spectral bins
	IF NOT keyword_set(t_pts) THEN t_pts=emiss.t							;set default time points
	IF NOT keyword_set(missing) THEN missing=1.0e-4
	IF NOT keyword_set(grid2rmid) THEN grid2rmid=genpos_grid2rmid(gridpts,shot,tpts=t_pts)		;get grid2rmid if not supplied
        
	;determine the number of time points
	xtpts=size(t_pts)
	IF xtpts[0] EQ 1 OR xtpts[0] EQ 0 THEN n_tpts=1 ELSE n_tpts=xtpts[1]

        ;determine number of lines
        n_lines=n(line.lam_o)+1

	;determine number of grid points
	n_grid=gridpts.n[0]*gridpts.n[1]
	       
        ;determine the Bphi, Br, Bz at the grid points
        bfields=genpos_grid_bfield(gridpts,shot,t_pts=t_pts)
        br=bfields.br
        bphi=bfields.bphi
        bz=bfields.bz

	;determine number of upos
	IF NOT keyword_set(noraw) THEN n_upos=n(gpv.gpv[0,*])+1 ELSE upos=1

	;define physical constants
	c=3.0e8  		;speed of light
	angst=1.0e-10 		;conversion of lambda to meters
	e=1.60e-19		;conversion for eV -> J
	mconv=1.66e-27		;conversion for amu -> kg

	;determine the spectral range
	lam=make(min(line.lam_o)-line.del_lam,max(line.lam_o)+line.del_lam,n_lam)	

	;initialize output arrays
	gpvspec=fltarr(n_grid,n_lam,n_tpts)
	plam=fltarr(n_lam,n_tpts)

	IF NOT keyword_set(interp) THEN BEGIN
		;initialize interpolation variables
		ti_interp=fltarr(n_grid,n_tpts)+missing
		em_interp=fltarr(n_grid,n_lines,n_tpts)+missing
		u_interp=fltarr(n_grid,n_tpts)+missing
                w_interp=fltarr(n_grid,n_tpts)+missing

		;interpolate ion temperature
		IF NOT keyword_set(quiet) THEN print, 'Interpolating Ti'
		FOR i=0,n_tpts-1 DO BEGIN
			tmp=where(ti.t EQ t_pts[i])
			IF tmp[0] EQ -1 THEN BEGIN
				bnd=ibound(ti.t,t_pts[i])
				IF bnd[0] NE -1 THEN BEGIN
					a=bnd[0]
					b=bnd[1]
					ti_tmp=ti.ti[*,a]+(ti.ti[*,b]-ti.ti[*,a])/(ti.t[b]-ti.t[a])*(t_pts[i]-ti.t[a])
				ENDIF ELSE RETURN, -1
			ENDIF ELSE BEGIN
				ti_tmp=ti.ti[*,tmp[0]]
			ENDELSE
			tmp=where(grid2rmid[*,i] GE min(ti.r) AND grid2rmid[*,i] LE max(ti.r)) 
			ti_interp[tmp,i]=interpol(ti_tmp,ti.r,grid2rmid[tmp,i])
			
		ENDFOR

		;interpolate line emissivity
		IF NOT keyword_set(quiet) THEN print, 'Interpolating Emissivity'
		FOR i=0,n_tpts-1 DO BEGIN
                    	FOR j=0,n_lines-1 DO BEGIN
                            	tmp=where(emiss.t EQ t_pts[i])
				IF tmp[0] EQ -1 THEN BEGIN
					bnd=ibound(emiss.t,t_pts[i])
					IF bnd[0] NE -1 THEN BEGIN
						a=bnd[0]
						b=bnd[1]
						em_tmp=emiss.emiss[*,j,a]+(emiss.emiss[*,j,b]-emiss.emiss[*,j,a])/(emiss.t[b]-emiss.t[a])*(t_pts[i]-emiss.t[a])
					ENDIF ELSE RETURN, -2
				ENDIF ELSE BEGIN
					em_tmp=emiss.emiss[*,j,tmp[0]]
				ENDELSE
				tmp=where(grid2rmid[*,i] GE min(emiss.r) AND grid2rmid[*,i] LE max(emiss.r))
				em_interp[tmp,j,i]=interpol(em_tmp,emiss.r,grid2rmid[tmp,i])
                        ENDFOR
		ENDFOR

		;interpolate velocity
		IF NOT keyword_set(quiet) THEN print, 'Interpolating Velocity'
		FOR i=0,n_tpts-1 DO BEGIN
			tmp=where(vel.t EQ t_pts[i])
			IF tmp[0] EQ -1 THEN BEGIN
				bnd=ibound(vel.t,t_pts[i])
				IF bnd[0] NE -1 THEN BEGIN
					a=bnd[0]
					b=bnd[1]
					u_tmp=vel.u[*,a]+(vel.u[*,b]-vel.u[*,a])/(vel.t[b]-vel.t[a])*(t_pts[i]-vel.t[a])
                                        w_tmp=vel.w[*,a]+(vel.w[*,b]-vel.w[*,a])/(vel.t[b]-vel.t[a])*(t_pts[i]-vel.t[a])
				ENDIF ELSE RETURN, -3
			ENDIF ELSE BEGIN
				u_tmp=vel.u[*,tmp[0]]
                                w_tmp=vel.w[*,tmp[0]]
			ENDELSE
			tmp=where(grid2rmid[*,i] GE min(vel.r) AND grid2rmid[*,i] LE max(vel.r))
			u_interp[tmp,i]=interpol(u_tmp,vel.r,grid2rmid[tmp,i])
                        w_interp[tmp,i]=interpol(w_tmp,vel.r,grid2rmid[tmp,i])
		ENDFOR
		
		;setup optional output if function will be run for multiple detectors
		interp={ti:ti_interp,em:em_interp,u:u_interp,w:w_interp}
	ENDIF ELSE BEGIN
		ti_interp=interp.ti
		em_interp=interp.em
                u_interp=interp.u
		w_interp=interp.w
	ENDELSE

	IF keyword_set(debug) THEN stop
                                                    
	;define constants that need not be calculated during loops
	sq2pi=sqrt(2.0*!pi)					;
	r_major=gridpts.pnts[0,*]	;rename for clarity
        z_major=gridpts.pnts[1,*]	
        IF NOT keyword_set(noraw) THEN BEGIN
        	lz=lhat.lhat.lz
	        lphi=lhat.lhat.lphi
        	lr_in=lhat.lhat.lr.in		;decompile the lhat structure for easier coding
	        lr_out=lhat.lhat.lr.out
        	vin=lhat.lhat.lr.vin
	        vout=lhat.lhat.lr.vout
        ENDIF ELSE BEGIN
  		lz=lhat.lz
                lphi=lhat.lphi
                lr_in=lhat.lr.in
                lr_out=lhat.lr.out
                vin=lhat.lr.vin
                vout=lhat.lr.vout
        ENDELSE

        FOR h=0,n_lines-1 DO BEGIN
            	const0=line.lam_o[h]*sqrt(1.0e3*e/(line.mass*mconv))/c	;multiplier for Ti
                IF NOT keyword_set(quiet) AND NOT keyword_set(noraw) THEN print, 'Adding Line '+num2str(h,1)+' of '+num2str(n_lines-1,1)
            	FOR i=0,n_tpts-1 DO BEGIN                    	
			IF NOT keyword_set(noraw) THEN BEGIN
                        	FOR j=0L,n_upos-1 DO BEGIN
					tmp=gpv.tmp						;isolate pixels that matter to save time (it saves a lot!)	
                                        sigma=const0*sqrt(ti_interp[tmp,i]) ;gaussian width from Ti
                                        vec0=em_interp[tmp,h,i]/(sq2pi*sigma) ;peak spectral emissivity			
                                        delta_in=1.0/c*(u_interp[tmp,i]*(lr_in[*,j]*br[tmp]+lz[*,j]*bz[tmp]+lphi[*,j]*bphi[tmp])$
	       	                 		+w_interp[tmp,i]*lphi[*,j]*r_major[tmp]) ;1/c*dot(v,l_hat)
                                        vec1_in=line.lam_o[h]*(1.0+delta_in) ;line center adjusted for velocity
                                        delta_out=1.0/c*(u_interp[tmp,i]*(lr_out[*,j]*br[tmp]+lz[*,j]*bz[tmp]+lphi[*,j]*bphi[tmp])$
                        			+w_interp[tmp,i]*lphi[*,j]*r_major[tmp]) ;1/c*dot(v,l_hat)
                                        vec1_out=line.lam_o[h]*(1.0+delta_out) ;line center adjusted for velocity
					FOR k=0,n_lam-1 DO BEGIN
        	        	       		vec3=lam[k]-vec1_in
						spec=vec0*exp(-(vec3)^2/(2.0*sigma^2))
						gpvspec[tmp,k,i]+=vin[*,j]*spec		;add at each lambda add the volume weighting going in
						vec3=lam[k]-vec1_out
						spec=vec0*exp(-(vec3)^2/(2.0*sigma^2))
						gpvspec[tmp,k,i]+=vout[*,j]*spec	;add at each lambda add the volume weighting going out              
					ENDFOR							
        	        	ENDFOR
				FOR k=0,n_lam-1 DO plam[k,i]=total(gpvspec[*,k,i])	
                        ENDIF ELSE BEGIN
                        	tmp=where(gpv NE 0)
                                sigma=const0*sqrt(ti_interp[tmp,i])
                                vec0=em_interp[tmp,h,i]/(sq2pi*sigma)
                                delta_in=1.0/c*(u_interp[tmp,i]*(lr_in[tmp]*br[tmp]+lz[tmp]*bz[tmp]+lphi[tmp]*bphi[tmp])$
	       	                	+w_interp[tmp,i]*lphi[tmp]*r_major[tmp]) ;1/c*dot(v,l_hat)
                                vec1_in=line.lam_o[h]*(1.0+delta_in) ;line center adjusted for velocity
                                delta_out=1.0/c*(u_interp[tmp,i]*(lr_out[tmp]*br[tmp]+lz[tmp]*bz[tmp]+lphi[tmp]*bphi[tmp])$
                        		+w_interp[tmp,i]*lphi[tmp]*r_major[tmp]) ;1/c*dot(v,l_hat)
                                vec1_out=line.lam_o[h]*(1.0+delta_out) ;line center adjusted for velocity
				FOR k=0,n_lam-1 DO BEGIN
        	        	   	vec3=lam[k]-vec1_in
					spec=vec0*exp(-(vec3)^2/(2.0*sigma^2))
					plam[k,i]+=total(vin[tmp]*spec)		;add at each lambda add the volume weighting going in
					vec3=lam[k]-vec1_out
					spec=vec0*exp(-(vec3)^2/(2.0*sigma^2))
					plam[k,i]+=total(vout[tmp]*spec)	;add at each lambda add the volume weighting going out    
				ENDFOR
			ENDELSE
            	ENDFOR
	ENDFOR
	
	;sigma=sqrt(ti_interp[tmp,i])*line.lam_o*sqrt(1.0e3*e/(line.mass*mconv))/c
	;delta=vel_interp[tmp,i]/c*upos[2,j]/r_major[tmp]
	;spec = em_interp[tmp,i]/(sqrt(2.0*!pi)*sigma)*exp(-(line.lam_o*(1.0+delta)-lam)^2/(2.0*sigma^2))


	output={plam:plam,lam:lam}
	IF keyword_set(debug) THEN stop
	RETURN,output
 END


;+
;NAME:
;	GENSPEC_MOMENTS
;
;PURPOSE:
;	This function forms the zeroth, first and second moment vectors from an array of spectra.  This data is
;	needed as an input to GENPOS_INVERT_MOMENTS.
;
;CALLING SEQUENCE:
;	result=GENSPEC_MOMENTS(spec,lam,lam_o)
;
;INPUTS:
;	spec:	FLTARR of size n_lam x n_ch x n_time of the spectra units [X/[lam_o units] - like ph/s/ang or mW/ang]
;	lam:	FLTARR of length n_lam that is the wavelength scale for spec.  If the wavelength scale is different for
;			different channels then lam can be input as a n_lam x n_ch array.
;	lam_o:	FLT of the rest wavelength of the line.
;	
;KEYWORD PARMAETERS:
;	debug:		/debug stops the codes before the RETURN statement
;
;OUTPUTS:
;	result:	FLTARR of size 3 x n_ch,n_time where [i,*,k] is the ith-moment for all the channels at the kth time point
;		The units are those of X, X*lam_o and X*lam_o^2 for each moment.
;
;		An output of result = -1 indicate that the lam input is improperly formatted compared to spec
;PROCEDURE:
;	The moments are formulated using INT_TABULATED in the following manner (j-det k-time):
;		mom_i=int_tabulated(lam,spec[*,j,k]*(lam-lam_o)^i)
;
;MODIFICATION HISTORY:
;	Written by:	ML Reinke: 8-21-06
;
;-


FUNCTION genspec_moments,spec,lam,lam_o,debug=debug
	
	num_ch=n(spec[0,*,0])+1
	num_time=n(spec[0,0,*])+1
	num_lam=n(lam[*,0])+1
	
	;check sizing of lam and spec
	IF n(spec[*,0,0])+1 NE num_lam THEN RETURN,-1	
	x=size(lam)
	IF x[0] EQ 2 AND x[2] NE num_ch THEN RETURN,-1

	mom=fltarr(3,num_ch,num_time) 		;initialize ouptut

	;choose between lam being the same for all detectors and lam being a function of each detector
	IF x[0] EQ 1 THEN BEGIN
 		del=lam-(fltarr(num_lam)+lam_o)	
		FOR i=0,num_ch-1 DO BEGIN
			FOR j=0,num_time-1 DO BEGIN
				mom[0,i,j]=int_tabulated(lam,spec[*,i,j])
				mom[1,i,j]=int_tabulated(lam, spec[*,i,j]*del)
				mom[2,i,j]=int_tabulated(lam, spec[*,i,j]*del^2)
			ENDFOR
		ENDFOR
	ENDIF ELSE BEGIN
		FOR i=0,num_ch-1 DO BEGIN
			del=lam[*,i]-(fltarr(num_lam)+lam_o)	
			FOR j=0,num_time-1 DO BEGIN
				mom[0,i,j]=int_tabulated(lam[*,i],spec[*,i,j])
				mom[1,i,j]=int_tabulated(lam[*,i], spec[*,i,j]*del)
				mom[2,i,j]=int_tabulated(lam[*,i], spec[*,i,j]*del^2)
			ENDFOR
		ENDFOR	
        ENDELSE

	output=mom
	IF keyword_set(debug) THEN stop
	RETURN,mom
	
END

;+
;NAME:
;	GENSPEC_LI_MOMENTS
;
;PURPOSE:
;	This function generates the spectral moments used in GENSPEC inversion from sets of (lambda, intensity)
;	pairs instead of ordered data like GENSPEC_MOMENTS needs.
;
;CALLING SEQUENCE
;	result=GENSPEC_LI_MOMENTS(lam,int,lam_o)	
;
;INPUTS:
;	lam:	FLTARR [n_pts,n_ch,n_time] of the wavelength values [in lam_o units] for each channel at each time
;			OPTIONALLY you can make lam a FLTARR [n_pts]
;	int:	FLTARR [n_pts,n_ch,n_time] of the intensity values [photons/lam] at those wavelength values
;	lam_o	FLT	of the unshifted wavelength of the line
;
;OPTIONAL INPUTS:
;	lr:	FLTARR [2] (lam_low, lam_high) of the inclusive wavelength region to truncate dataset.  See
;		GENSPEC_MOMENT_RANGE for generating values of this number.
;	ikplot:	FLTARR [2] of the (i,k) (ch, time) point that you only want to plot.  Still use /plot to invoke.
;
;KEYWORD PARAMETERS:
;	backsub:	/backsub will use the first two and last two lambda points to form an average for background subtraction
;	debug:		/debug stops the code before the return statement
;	plot:		/plot will plot int vs lam for each time point for each channel
;	sort:		/sort invokes /sort in INT_TABULATED
;	double:		/double invokes /double in INT_TABULATED
;	
;OUTPUTS:
;	result:	FLTARR 	[3, n_ch, n_time] of the spectral moments
;			[0, *, *] is the zeroth moment [photons]
;			[1, *, *] is the first moment  "wavelength*counts" [angstroms*photons]  consistent units with lam_o
;			[2, *, *] is the second moment  "wavelength^2*counts" [angstroms^2*photons]  consistent units with lam_o
;
;OPTIONAL OUTPUTS:
;	err:	FLTARR  [3, n_ch, n_time] of statisical uncertainty (sigma)
;			[0, *, *] is the zeroth moment error
;			[1, *, *] is the first moment error
;			[2, *, *] is the second moment error
;
;PROCEDURE
;	This function uses INT_TABULATED to form the numerical integral.  Thus values of lam can be irregularly gridded and even in a
;	random order (use /sort then) but they must be unqiue.
;	
;	Use GENSPEC_MOMENT_RANGE to inform the decesion of an lr interval.  This will help to reduce photon and background subtraction noise
;	from affecting the first and second moment values.
;
;MODIFICATION HISTORY;
;	Written by: 	ML Reinke - 7/19/07
;       8-22-07:	ML Reinke - modified so that lam could be a 1D as well as 3D
;       10-30-08:	ML Reinke - added the optional output err, the statistical errors of the moment profiles
;	12-11-08	ML Reinke - made the moments error absolute value since when doing background subtraction negative brightness
;                                   can result.  Need to modify to include bsub error
;
;-

FUNCTION genspec_li_moments,lam,int,lam_o,debug=debug,plot=plot,ikplot=ikplot,lr=lr,backsub=backsub,sort=sort,double=double,err=err

	size=size(int)
        IF size[0] EQ 3 THEN BEGIN
        	num_lam=size[1]
                num_ch=size[2]
                num_time=size[3]
        ENDIF ELSE BEGIN
        	num_lam=size[1]
                num_ch=size[2]
                num_time=1
        ENDELSE
        
        ;allow for 1D and 3D formatting of lam input
        size=size(lam)
        IF size[0] EQ 3 THEN lam1d=0 ELSE lam1d=1
        a=0
        b=0

	del=lam-(fltarr(num_lam)+lam_o)

	mom=fltarr(3,num_ch, num_time)
	err=fltarr(3,num_ch, num_time)
	FOR i=0,num_ch-1 DO BEGIN
            	FOR j=0,num_time-1 DO BEGIN
                       ;del=reform(lam[*,i,j])-(fltarr(num_lam)+lam_o)
                    	IF NOT lam1d THEN BEGIN
                        	a=i
                                b=j
                        ENDIF
			IF total(int[*,i,j]) NE 0 THEN BEGIN
				IF keyword_set(backsub) THEN BEGIN
                                    	baseline=mean([int[0:1,i,j],int[num_lam-2:num_lam-1,i,j]])
                                	int[*,i,j]-=baseline
                                        ;print, baseline
                                ENDIF
				IF keyword_set(lr) THEN BEGIN
					tmp=where(lam[*,a,b] GE lr[0] AND lam[*,a,b] LE lr[1])
					x=reform(lam[tmp,a,b])
					y=reform(int[tmp,i,j])
				ENDIF ELSE BEGIN
					x=reform(lam[*,a,b])
					y=reform(int[*,i,j])
				ENDELSE
				del=x-(fltarr(n(x)+1)+lam_o)
				mom[0,i,j]=int_tabulated(x,y,sort=sort,double=double)
				mom[1,i,j]=int_tabulated(x,y*del,sort=sort,double=double)
				mom[2,i,j]=int_tabulated(x,y*del^2,sort=sort,double=double)
                                err[0,i,j]=sqrt(abs(mom[0,i,j]))	;use absolute error since when doing background can get < 0
                                err[1,i,j]=sqrt(abs(mom[2,i,j]))	;use absolute error since when doing background can get < 0
                                err[2,i,j]=sqrt(abs(3.0*mom[0,i,j]))*(mom[2,i,j]/mom[0,i,j]-(mom[1,i,j]/mom[0,i,j])^2) ;use absolute error since when doing background can get < 0
				IF keyword_set(plot) THEN BEGIN
					IF keyword_set(ikplot) THEN BEGIN
						IF i EQ ikplot[0] AND j EQ ikplot[1] THEN BEGIN
                                                	plot, x,y,xtit='Wavelength',ytit='Intensity'
							print, 'CH :'+num2str(i,1)+' TIME: '+num2str(j,1)
							print, 'MOM0 '+num2str(mom[0,i,j])
							print, 'MOM1 '+num2str(mom[1,i,j])
							print, 'MOM2 '+num2str(mom[2,i,j])
							stop		
                                                ENDIF
					ENDIF ELSE BEGIN
						plot, x,y,xtit='Wavelength',ytit='Intensity'
						print, 'CH :'+num2str(i,1)+' TIME: '+num2str(j,1)
						print, 'MOM0 '+num2str(mom[0,i,j])
						print, 'MOM1 '+num2str(mom[1,i,j])
						print, 'MOM2 '+num2str(mom[2,i,j])
						stop
					ENDELSE
				ENDIF
			ENDIF
		ENDFOR
	ENDFOR	
	IF keyword_set(debug) THEN stop
	RETURN, mom
END

;+
;NAME:
;	GENSPEC_VEL_COEF_MATRIX
;
;PURPOSE:
;	This function generates the coefficent matrix used to invert the v*em profile.  The velocity is assumed to
;	be purely toroidal and a flux function.  The radial profile is expanded using polynomials or orthogonal bessel functions.
;
;CALLING SEQUENCE:
;	result=GENSPEC_VEL_COEF_MATRIX(gpv,pos,line,rhopts,ves_grid,order)
;
;INPUTS:
;	gpv:		FLTARR	n x m where n is the number of channels and m is the number of pixels.  These are the volume coefficents
;				generated using GENPOS_VOL_COEFS [meters^3]
;	pos:		FLTARR 	4 x n of the central (or weighted) POS vector for each channels.
;	bfield:		STRUC   of data describing the magnetic field on the grid (See GENPOS_GRID_BFIELD for method)
;			*.br		FLTARR [n_grid, n_time] radial component of the magnetic field [T]
;			*.bphi 		FLTARR [n_grid, n_time] toroidal component of the magnetic field [T]
;			*.bz		FLTARR [n_grid, n_time] veritical component of the magnetic field [T]
;			*.psi_norm	FLTARR [n_grid, n_time] normalized psi values ***not used***
;			*.t_pts		FLTARR [n_time] time points [s] ***not used***
;			*.shot		LONG shot number ***not used***
;	line:		STRUC 	containing information about the line to be calculated
;			*.lam_o		FLT the rest wavelength of the line [units set wavelength scale]        
;			*.del_lam	FLT used in the  wavelength interval [lam_o-del_lam, lam_o+del_lam]    ***not used***
;			*.mass		FLT the mass of ion in AMU (see READ_ATOMIC_MASS('XX')
;	rhopts:		FLTARR	m x t where t is the number of time points.  These are normalized radius (rho) values that correspond 
;				to the m pixels and should be generated using GENPOS_GRID2RMID.
;	ves_grid:	STRUC	the vessel grid generated from GENPOS_GRID using /center
;	order:		INT	order of the expansion (see PROCEDURE for details)
;
;KEYWORD PARAMETERS:
;	bessel:		/bessel will use an orthogonal bessel function expansion instead of polynomials
;	solidbody:	/solidbody will generate a coef matrix using only w(rho) which is pure toroidal rotation
;
;OUTPUTS:
;	result:		FLTARR n x order x t where [*,*,i] is the matrix, A, to be used in the least squares inversion of (A*x=b) where
;			b is the vector of length n and x is the vector of length order and are the coeffcients of the
;			See http://en.wikipedia.org/wiki/Linear_least_squares for a brief overview and Numerical Recipes for
;			something more indepth.
;
;PROCEDUE:
;	There will be a note in 'General Equations for Radiometry in Tokamak Plasmas' regarding how the expansion works but
;	for now just talk to ML Reinke directly if you've got questions on the method.
;
;MODIFICATION HISTORY:
;	Written by:	ML Reinke: 8-21-06
;	6-7-08:	 	ML Reinke - modified to use u(rho) w(rho) and added solidbody keyword to just use w(rho)
;				    Also changed the WHERE that total is taken over to be both the non-zero values of GPV and rhopts
;-


FUNCTION genspec_vel_coef_matrix,gpv,pos,bfield,line,rhopts,ves_grid,order,debug=debug,bessel=bessel,solidbody=solidbody
	
	c=3.0e8  		;speed of light

	num_det=n(gpv[*,0])+1
	num_time=n(rhopts[0,*])+1
       	IF keyword_set(bessel) THEN  bessel_zeros=bessel_zeros()

	IF keyword_set(solidbody) THEN BEGIN
		coef_matrix=dblarr(num_det,order,num_time)
		FOR i=0,num_det-1 DO BEGIN
			FOR j=0,order-1 DO BEGIN
				FOR k=0,num_time-1 DO BEGIN
					tmp=where(rhopts[*,k] GE 0 AND gpv[i,*] GT 0)
					lhat=genpos_lhat(pos[*,i],ves_grid,good=tmp)
                                        Bi=lhat.lphi*ves_grid.pnts[0,tmp]
					IF NOT keyword_set(bessel) THEN BEGIN 
                                    		IF j GT 0 THEN jplus=1 ELSE jplus=0
	                                	coef_matrix[i,j,k]=total(gpv[i,tmp]*Bi*rhopts[tmp,k]^(j+jplus))
        	                        ENDIF ELSE BEGIN                                      
						coef_matrix[i,j,k]=total(gpv[i,tmp]*Bi*beselj(rhopts[tmp,k]*bessel_zeros[2*j],2*j))
                                	ENDELSE
				ENDFOR
			ENDFOR
		ENDFOR
	ENDIF ELSE BEGIN
		coef_matrix=dblarr(num_det,order*2,num_time)
		FOR i=0,num_det-1 DO BEGIN
			FOR j=0,order-1 DO BEGIN
				FOR k=0,num_time-1 DO BEGIN
					tmp=where(rhopts[*,k] GE 0 AND gpv[i,*] GT 0)
					lhat=genpos_lhat(pos[*,i],ves_grid,good=tmp)
					Ai=lhat.lr*bfield.br[tmp,k]+lhat.lz*bfield.bz[tmp,k]+lhat.lphi*bfield.bphi[tmp,k]
					Bi=lhat.lphi*ves_grid.pnts[0,tmp]
					IF NOT keyword_set(bessel) THEN BEGIN 
                                    		IF j GT 0 THEN jplus=1 ELSE jplus=0
	                                	coef_matrix[i,j,k]=total(gpv[i,tmp]*Bi*rhopts[tmp,k]^(j+jplus))
						coef_matrix[i,j+order,k]=total(gpv[i,tmp]*Ai*rhopts[tmp,k]^(j+jplus))
        	                        ENDIF ELSE BEGIN                                      
						coef_matrix[i,j,k]=total(gpv[i,tmp]*Bi*beselj(rhopts[tmp,k]*bessel_zeros[2*j],2*j))
						coef_matrix[i,j+order,k]=total(gpv[i,tmp]*Ai*beselj(rhopts[tmp,k]*bessel_zeros[2*(j+1)],2*(j+1)))
                                	ENDELSE
				ENDFOR
			ENDFOR

		ENDFOR
	ENDELSE
	

	output=line.lam_o/c*coef_matrix
	IF keyword_set(debug) THEN stop
	RETURN,output
END

;+
;NAME:
;	GENSPEC_GPV2VOXEL_VEL_MATRIX
;
;PURPOSE
;	This function will take a set of GPVs, the rho values of their locations,the pos vector, ves_cent and bfield  and
;	turn it into a weighting matrix for velocity weighted emissivity profile inversion.
;
;CALLING SEQUENCE
;	result=GENSPEC_GPV2VOXEL_VEL_MATRIX(gpv, rhopts,pos,ves_cent,bfield)
;
;INPUTS:
;	gpv		FLTARR [n_ch, n_pts] of the volume weightings generated using GENPOS_VOL_COEFS
;	rhophts		FLTARR [n_pts,n_time] of the rho locations of the (R,Z) points that locate each GPV at each time point
;	pos:		FLTARR 	4 x n of the central (or weighted) POS vector for each channel.
;	ves_cent:	STRUC	the vessel grid generated from GENPOS_GRID or GRID_VES using /center
;	bfield:		STRUC   of data describing the magnetic field on the grid (See GENPOS_GRID_BFIELD for method)
;			*.br		FLTARR [n_pts] radial component of the magnetic field [T]
;			*.bphi 		FLTARR [n_pts] toroidal component of the magnetic field [T]
;			*.bz		FLTARR [n_pts] veritical component of the magnetic field [T]
;			
;			****BFIELD NOT BEING USED YET****
;
;OPTIONAL INPUTS:
;	n_rho:		INT number of points from 0.0 -> 1.0 (inclusive) for the 
;	rho_vec:	FLTARR [n_rho] of the rho points that each column of the result
;	posrev:		INTARR [n_ch] of 0's and 1's to indicate for each pos to reverse or not in GENPOS_LHAT
;
;KEYWORD PARAMETERS:
;	solidbody:	/solidbody forces a matrix to be generated that will be used to find w profile (just toroidal rotation DEFAULT)
;	diff:		/diff is passed to GENPOS_LHAT
;
;OUTPUTS
;	result:		FLTARR [n_ch, n_rho, n_time] of the voxel weightings for each channel at each rho to be used in profile inversion
;
;OPTIONAL OUTPUTS:
;	rho_vec:	FLTARR if set to an unnamed variable it will be filled with the rho.
;
;MODIFIATION HISTORY:
;	Written by:	ML Reinke 6/11/07 - only allowing toroidal rotation
;	7-12-07:	ML Reinke - fixed a WHERE bug that would crash if no weighting in the plasma was foun
;       7-17-07:	ML Reinke - allowed 2D rhopts so that multiple time slices can be calculated
;	7-19-07:	ML Reinke - moved into GENSPEC.PRO changed name from GENPOS_*** to GENSPEC_***
;	2-18-08:	ML Reinke - updated the u and w voxel matrices to interpolate between rho points
;                                   instead of locking to nearest one.
;	2-18-08:	ML Reinke - added posrev optional input and diff keyword parameter to be used with GENPOS_LHAT
;-

			
FUNCTION genspec_gpv2voxel_vel_matrix,gpv,rhopts,pos,ves_cent,bfield,rho_vec=rho_vec,n_rho=n_rho,solidbody=solidbody,posrev=posrev,diff=diff

	IF NOT keyword_set(n_rho) THEN n_rho=20
	IF NOT keyword_set(rho_vec) THEN rho_vec=make(0.0,1.0,n_rho)

	n_ch=n(gpv[*,0])+1
	n_grid=n(rhopts[*,0])+1
        n_time=n(rhopts[0,*])+1

	w_voxel=fltarr(n_ch,n_rho,n_time)
        u_voxel=fltarr(n_ch,n_rho,n_time)
       	FOR i=0,n_ch-1 DO BEGIN
		good_gpv=where(gpv[i,*] GT 0)
	        IF good_gpv[0] NE -1 THEN BEGIN
                    	IF keyword_set(posrev) THEN rev=posrev[i] ELSE rev=0
        		lhat=genpos_lhat(pos[*,i],ves_cent,good=good_gpv,rev=rev,diff=diff)	;find lhat at points where gpv > 0
	                Bi=lhat.lphi*ves_cent.pnts[0,good_gpv]			;calculate lhat*R
        	        FOR k=0,n_time-1 DO BEGIN
                		tmp=where(rhopts[good_gpv,k] GT 0)		;find good rho values (rho < 1 for SOL if chosen)
                                IF NOT keyword_set(solidbody) THEN Ai=lhat.lr[tmp]*bfield.br[good_gpv[tmp],k]+$			;calc parallel weighting
                                      	lhat.lz[tmp]*bfield.bz[good_gpv[tmp],k]+lhat.lphi[tmp]*bfield.bphi[good_gpv[tmp],k]
                        	IF tmp[0] NE -1 THEN FOR j=0,n(tmp) DO BEGIN	;fill for each rho point
                                    	rho_ch=rhopts[good_gpv[tmp[j]],k]
                                        ibnd=ibound(rho_vec,rho_ch)		;find bounding rhos 
	                                IF ibnd[0] NE -1 THEN BEGIN
                                            	del_rho=rho_vec[ibnd[1]]-rho_vec[ibnd[0]]
                                                del_high=rho_vec[ibnd[1]]-rho_ch
                                                del_low=rho_ch-rho_vec[ibnd[0]]
                                                w_vox=gpv[i,good_gpv[tmp[j]]]*Bi[tmp[j]]
                                                IF NOT keyword_set(solidbody) THEN u_vox=gpv[i,good_gpv[tmp[j]]]*Ai[j] ELSE u_vox=0.0
                                                IF ibnd[0] NE ibnd[1] THEN BEGIN	;fill bounding points, splitting the weighting
                                                	w_voxel[i,ibnd[0],k]+=w_vox*(1.0-del_low/del_rho)
                                                        w_voxel[i,ibnd[1],k]+=w_vox*(1.0-del_high/del_rho)
                                                     	u_voxel[i,ibnd[0],k]+=u_vox*(1.0-del_low/del_rho)   
                                                        u_voxel[i,ibnd[1],k]+=u_vox*(1.0-del_high/del_rho)
                                                ENDIF ELSE BEGIN
                                                    	w_voxel[i,ibnd[0],k]+=w_vox
                                                        u_voxel[i,ibnd[0],k]+=u_vox
                                                ENDELSE
                                        ENDIF
                                ENDFOR
                        ENDFOR
                ENDIF
        ENDFOR

	output={w:w_voxel, u:u_voxel}
	RETURN,output
END

;+
;NAME:
;	GENSPEC_POS2VOXEL_VEL_MATRIX
;
;PURPOSE:
;	This function creates a velocity voxel matrix using a pos vector and an etendue vector using an EFIT reconstruction
;
;CALLING SEQUENCE:
;	result=genspec_pos2voxel_matrix(pos,u,shot)
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
;	m:		INT	of the poloidal m number of the sine/cosine weighting matrix (if selected) DEFAULT: 1
;	tree:		STRING	of the EFIT tree to use for calculating weighting matrix [DEFAULT: 'analysis']
;
;KEYWORD PARAMETERS:
;	debug:		/debug stops the function before the RETURN statement
;	sine:		/sine multiples each voxel weighting by sine (using (R,Z) pt and magnetic axis at time point)
;	cosine:		/cosine multiples each voxel weighting by cosine (using (R,Z) pt and magnetic axis at time point)
;
;OUTPUTS
;	result:		FLTARR [n_ch, n_rho, n_time,2] of the voxel weightings for each channel at each rho to be used in profile inversion
;				at each time point.  [*,*,*,0] are for the w inversion while [*,*,*,1] are for the u profile.
;
;OPTIONAL OUTPUTS:
;	rho_vec:	FLTARR if set to an unnamed variable it will be filled with the rho vector used if it is not input through rho_vec.
;	rhopts:		FLTARR [n_ch*n_s,n_time] of the rho values for all the channels along the lines of sight for each time
;				slice.  This output is tied to a specific pos,tpts and n_sinput and is useful for iterating on rho_vec.
;	br:		FLTARR [n_ch*n_s,n_time] of the radial B-field along the line of sight for each time slice
;	bz:		FLTARR [n_ch*n_s,n_time] of the vertical B-field along the line of sight for each time slice
;
;
;MODFICATION HISTORY:
;	Written by:	ML Reinke 8/13/10 (Adapted from GENPOS_POS2VOXEL_MATRIX and GENSPEC_GPV2VOXEL_VEL_MATRIX)
;	4/4/11		ML Reinke - added ability to get arb. m cosine/sine velvoxel matrices
;	4/5/11		ML Reinke - noticed an error in the l_hat_r term, also realized that this lhat
;                                   parameterization will NOT work for RT > Rmin where l_hat_r changes sign
;	6-13-11:	ML Reinke - added the tree optional input for use with EFIT_RZ2XXX codes.	
;
;-


FUNCTION genspec_pos2voxel_vel_matrix,pos,u,shot,tpts=tpts,n_s=n_s,rzbnd=rzbnd,rho_vec=rho_vec,n_rho=n_rho,rho_min=rho_min,rz_ap=rz_ap,r_ap=r_ap,exp_rho=exp_rho,$
		kdebug=kdebug,debug=debug,verb=verb,rhopts=rhopts,br=br,bz=bz,psinorm=psinorm,sine=sine,cosine=cosine,m=m,tree=tree

	IF NOT keyword_set(n_rho) THEN IF keyword_set(rho_vec) THEN n_rho=n(rho_vec)+1 ELSE n_rho=20
        IF NOT keyword_set(rho_min) THEN rho_min=0.0
	IF NOT keyword_set(rho_vec) THEN rho_vec=make(rho_min,1.0,n_rho)
	IF NOT keyword_set(rzbnd) THEN rzbnd=[0.44,1.0,0.45]	;set path's [r_min,r_max,abs(z_max)] limit
	IF NOT keyword_set(n_s) THEN n_s=300
	IF NOT keyword_set(tpts) THEN tpts=line_gettimes(shot)
	IF NOT keyword_set(exp_rho) THEN exp_rho=42.5*n_rho/23.0
	n_ch=n(pos[0,*])+1
	n_time=n(tpts)+1
	n_s=long(n_s)

	;load the toroidal field
	mdsopen, 'magnetics', shot
        bt=abs(mdsvalue("\magnetics::BTOR",/quiet,status=status1))
  	bt_time=mdsvalue("dim_of(\magnetics::BTOR)",/quiet,status=status2)
	bt_int=interpol(bt,bt_time,tpts)	;interpolate the toroidal field onto the tpts grid
	mdsclose,'magnetics',shot
	
	;load the magnetic axis
	mdsopen,'analysis',shot
	raxis=mdsvalue('\ANALYSIS::TOP:EFIT.RESULTS.G_EQDSK:RMAXIS')
	zaxis=mdsvalue('\ANALYSIS::TOP:EFIT.RESULTS.G_EQDSK:ZMAXIS')
	taxis=mdsvalue('dim_of(\ANALYSIS::TOP:EFIT.RESULTS.G_EQDSK:RMAXIS)')
	raxis_int=interpol(raxis,taxis,tpts)	;interpolate the magnetic axis onto the tpts grid
	zaxis_int=interpol(zaxis,taxis,tpts)	;interpolate the magnetic axis onto the tpts grid
	mdsclose,'analysis',shot
        
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
        IF keyword_set(verb) THEN print, '(R,Z) points found, calling RZ2RMID and RZ2PSI'
	IF NOT keyword_set(rhopts) THEN BEGIN
        	IF keyword_set(psinorm) THEN rhopts=efit_rz2rho(pos_r,pos_z,tpts,shot=shot,/psinorm,tree=tree) ELSE rhopts=efit_rz2rmid(pos_r,pos_z,tpts,shot=shot,/rho,tree=tree)
        ENDIF
	IF NOT keyword_set(br) AND NOT keyword_set(bz) THEN psi=efit_rz2psi(pos_r,pos_z,tpts,bz,br,shot=shot)
        IF keyword_set(verb) THEN print, 'RZ2RMID done, filling voxel matrix'
	voxel=fltarr(n_ch,n_rho,n_time,2)	;[*,*,*,0] is the w matrix and [*,*,*,1] is the u matrix
        FOR k=0,n_time-1 DO BEGIN
		FOR i=0,n_ch-1 DO BEGIN
			rho_ch=rhopts[i*n_s:(i+1)*n_s-1,k]
			s_ch=pos_s[i*n_s:(i+1)*n_s-1]
			r_ch=pos_r[i*n_s:(i+1)*n_s-1]
			z_ch=pos_z[i*n_s:(i+1)*n_s-1]
			br_ch=br[i*n_s:(i+1)*n_s-1]				;br and bz at least (R,Z) point from efit_rz2psi
			bz_ch=bz[i*n_s:(i+1)*n_s-1]	
			bt_ch=bt_int[k]*raxis_int[k]/r_ch			;bt at each (R,Z) point from 1/R (ignores diamagnetic effect)
			w_vox=pos[2,i]*cos(pos[3,i])				;the R*l_phi=Rt*cos(psi) from the POS vector and doesn't depend on (R,Z)
			u_vox_r=-cos(pos[3,i])*sqrt(1.0-pos[2,i]^2/r_ch^2)*br_ch
			u_vox_phi=cos(pos[3,i])*pos[2,i]/r_ch*bt_ch
			u_vox_z=-sin(pos[3,i])*bz_ch		
			u_vox=u_vox_r+u_vox_phi+u_vox_z				;add all three terms to get the l*B weighting
			IF keyword_set(kdebug) THEN if i EQ kdebug THEN stop
			FOR j=1,n_s-2 DO BEGIN			;using +/- points to determine delta_s 
                                ibnd=ibound(rho_vec,rho_ch[j])
                                del_s=0.5*(s_ch[j+1]-s_ch[j-1])	
				vol=u[i]/(4.0*!pi)*sqrt(((pos[0,i])^2-(pos[2,i])^2)*(1.0+(tan(pos[3,i]))^2))*del_s
                                IF keyword_set(cosine) OR keyword_set(sine) THEN BEGIN
                                        deltaR=r_ch[j]-raxis_int[k]
                                        deltaZ=z_ch[j]-zaxis_int[k]
                                        IF deltaR GE 0 THEN BEGIN
                                            	IF deltaZ GE 0 THEN sign=1.0 ELSE sign=-1.0
                                                IF deltaZ GE 0 THEN dth=0 ELSE dth=2.0*!pi
                                        ENDIF
                                        IF deltaR LT 0 THEN BEGIN
                                            	IF deltaZ GE 0 THEN sign=-1.0 ELSE sign=1.0
                                                IF deltaR LT 0 THEN dth=!pi
                                        ENDIF
                                        th=sign*atan(abs(deltaZ/deltaR))+dth
                                        IF keyword_set(sine) THEN asym_term=sin(m*th)
                                        IF keyword_set(cosine) THEN asym_term=cos(m*th)
                                ENDIF ELSE asym_term=1.0
                        	IF ibnd[0] NE -1 THEN BEGIN
					del_rho=rho_vec[ibnd[1]]-rho_vec[ibnd[0]]
					del_high=rho_vec[ibnd[1]]-rho_ch[j]
					del_low=rho_ch[j]-rho_vec[ibnd[0]]

                                        IF ibnd[0] NE ibnd[1] THEN BEGIN	;rho[j] between rho vec values
						voxel[i,ibnd[0],k,0]+=vol*(1.0-del_low/del_rho)*w_vox*asym_term
						voxel[i,ibnd[1],k,0]+=vol*(1.0-del_high/del_rho)*w_vox*asym_term
						voxel[i,ibnd[0],k,1]+=vol*(1.0-del_low/del_rho)*u_vox[j]*asym_term
						voxel[i,ibnd[1],k,1]+=vol*(1.0-del_high/del_rho)*u_vox[j]*asym_term
                                        ENDIF ELSE BEGIN			;rho exactly on a rho_vec value
                                        	voxel[i,ibnd[0],k,0]+=vol*w_vox*asym_term
                                                voxel[i,ibnd[0],k,1]+=vol*u_vox[j]*asym_term
                                        ENDELSE
                                ENDIF ELSE IF rho_ch[j] GE rho_vec[n_rho-1] THEN BEGIN	;rho outside the rho_vec adds decaying amount to last channel (I think this is dumb)
					voxel[i,n_rho-1,k,0]+=vol*exp(-(rho_ch[j]-rho_vec[n_rho-1])*exp_rho)*w_vox*asym_term
					voxel[i,n_rho-1,k,1]+=vol*exp(-(rho_ch[j]-rho_vec[n_rho-1])*exp_rho)*u_vox[j]*asym_term
				ENDIF
                        ENDFOR
                        IF keyword_set(debug) THEN stop
                        IF i EQ 0 AND k EQ 0 AND keyword_set(verb) THEN print, 'vel_voxel matrix filled for ch '+num2str(i,1)+' of '+num2str(n_ch-1,1)
                        IF i EQ 0 AND k EQ 1 AND keyword_set(verb) THEN print, 'vel_voxel matrix filled for time '+num2str(k,1)+' of '+num2str(n_time-1,1)
		ENDFOR
        ENDFOR
	output=voxel
	IF keyword_set(debug) THEN stop
	RETURN,output
END




;+
;NAME:
;	GENSPEC_SUB_VECTOR
;
;PURPOSE:
;	This function generates the vector to be subtracted from the 2nd moment vector in the inversion of
;	a Doppler spectra.  Data from the emissivity and velocity inversion are necessary to generate the output.
;
;CALLING SEQUENCE:
;	result=GENSPEC_SUB_VECTOR(gpv,pos,line,rhopts,ves_grid,coefs_0,order_0,coefs_1,order_1)
;
;INPUTS:
;	gpv:		FLTARR	n x m where n is the number of channels and m is the number of pixels.  These are the volume coefficents
;				generated using GENPOS_VOL_COEFS [meters^3]
;	pos:		FLTARR 	4 x n of the central (or weighted) POS vector for each channel.
;	bfield:		STRUC   of data describing the magnetic field on the grid (See GENPOS_GRID_BFIELD for method)
;			*.br		FLTARR [n_grid, n_time] radial component of the magnetic field [T]
;			*.bphi 		FLTARR [n_grid, n_time] toroidal component of the magnetic field [T]
;			*.bz		FLTARR [n_grid, n_time] veritical component of the magnetic field [T]
;			*.psi_norm	FLTARR [n_grid, n_time] normalized psi values ***not used***
;			*.t_pts		FLTARR [n_time] time points [s] ***not used***
;			*.shot		LONG shot number ***not used***
;	line:		STRUC 	containing information about the line to be calculated
;			*.lam_o		FLT the rest wavelength of the line [units set wavelength scale]        
;			*.del_lam	FLT used in the  wavelength interval [lam_o-del_lam, lam_o+del_lam]    ***not used***
;			*.mass		FLT the mass of ion in AMU (see READ_ATOMIC_MASS('XX')
;	rhopts:		FLTARR	m x t where t is the number of time points.  These are normalized radius (rho) values that correspond 
;				to the m pixels and should be generated using GENPOS_GRID2RMID.
;	ves_grid:	STRUC	the vessel grid generated from GENPOS_GRID using /center
;	coefs_0:	FLTARR	of size order_0 x t of the emissivity expansion coefficents
;	order_0:	INT	order of the emissivity expansion
;	coefs_1:	FLTARR	of size order_1 x t of the velocity expansion coefficents 
;	order_1:	INT	order of the velocity expansion
;
;KEYWORD PARAMETERS:
;	debug:		/debug stops the code before the RETURN statement
;	bessel:		/bessel uses an orthogonal bessel expansion (see GENPOS_COEFS2PROFILE for specifics)
;	solidbody:	/solidbody will generate a sub vector using only w(rho) which is pure toroidal rotation
;
;OUTPUTS:
;	result:		FLTARR of size n x t.  This data must be subtracted from the 2nd moment vector to adjust for
;			velocity effects on the 2nd momemnt prior to solving the matrix equation for Ti.  Specifics are 
;			(or will be soon) in 'General Equations for Radiometry in Tokamak Plasmas'
;
;MODIFICATION HISTORY:
;	Written by:	ML Reinke: 8-21-06
;	6-7-07:		ML Reinke - adjusted program to use incompressible velocity togenerate the subtraction vector.
;                                   also made all totals sum over points where gpv > 0 and rhopts > 0
;
;-

FUNCTION genspec_sub_vector,gpv,pos,bfield,line,rhopts,ves_grid,coefs_0,order_0,coefs_1,order_1,debug=debug,bessel=bessel,solidbody=solidbody

	c=3.0e8  		;speed of light

	num_det=n(gpv[*,0])+1
	num_time=n(rhopts[0,*])+1
	num_grid=ves_grid.n[0]*ves_grid.n[1]

	sub_vector=dblarr(num_det,num_time)
        bessel_zeros=bessel_zeros()

	;form emissivity on the pixels
	em=fltarr(num_grid,num_time)
	FOR k=0,num_time-1 DO BEGIN
		FOR j=0,order_0-1 DO BEGIN
			tmp=where(rhopts[*,k] GE 0) 
			IF NOT keyword_set(bessel) THEN BEGIN 
                        	IF j GT 0 THEN jplus=1 ELSE jplus=0
                               	em[tmp,k]+=coefs_0[j]*rhopts[tmp,k]^(j+jplus)
                        ENDIF ELSE BEGIN                                      
				em[tmp,k]+=coefs_0[j]*beselj(rhopts[tmp,k]*bessel_zeros[2*j],2*j)
                        ENDELSE
		ENDFOR
	ENDFOR
	
        FOR i=0,num_det-1 DO BEGIN
 		;form em*vel on the pixels
		em_vel=fltarr(num_grid,num_time)
       		IF keyword_set(solidbody) THEN BEGIN
       		 	FOR k=0,num_time-1 DO BEGIN
               		 	FOR j=0,order_1-1 DO BEGIN
					tmp=where(rhopts[*,k] GE 0 AND gpv[i,*] GT 0)
                        	        lhat=genpos_lhat(pos[*,i],ves_grid,good=tmp)
                                	Bi=lhat.lphi*ves_grid.pnts[0,tmp]
					IF NOT keyword_set(bessel) THEN BEGIN 
       	 		                	IF j GT 0 THEN jplus=1 ELSE jplus=0
        	        	               	em_vel[tmp,k]+=coefs_1[j]*Bi*rhopts[tmp,k]^(j+jplus)
                        		ENDIF ELSE BEGIN                                      
						em_vel[tmp,k]+=coefs_1[j]*Bi*beselj(rhopts[tmp,k]*bessel_zeros[2*j],2*j)
		                        ENDELSE
				ENDFOR
                	ENDFOR
	        ENDIF ELSE BEGIN
        		FOR k=0,num_time-1 DO BEGIN
                		FOR j=0,order_1-1 DO BEGIN
					tmp=where(rhopts[*,k] GE 0 AND gpv[i,*] GT 0)
                                	lhat=genpos_lhat(pos[*,i],ves_grid,good=tmp)
					Ai=lhat.lr*bfield.br[tmp,k]+lhat.lz*bfield.bz[tmp,k]+lhat.lphi*bfield.bphi[tmp,k]
        	                        Bi=lhat.lphi*ves_grid.pnts[0,tmp]
					IF NOT keyword_set(bessel) THEN BEGIN 
        	        	        	IF j GT 0 THEN jplus=1 ELSE jplus=0
                	        	       	em_vel[tmp,k]+=coefs_1[j]*Bi*rhopts[tmp,k]^(j+jplus)
	                                        em_vel[tmp,k]+=coefs_1[j+order_1]*Ai*rhopts[tmp,k]^(j+jplus)
        	                	ENDIF ELSE BEGIN                                      
						em_vel[tmp,k]+=coefs_1[j]*Bi*beselj(rhopts[tmp,k]*bessel_zeros[2*j],2*j)
                        	                em_vel[tmp,k]+=coefs_1[j+order_1]*Ai*beselj(rhopts[tmp,k]*bessel_zeros[2*(j+1)],2*(j+1))
	                        	ENDELSE
				ENDFOR
        	        ENDFOR
                ENDELSE
                ;form the full subtraction vector using em*vel^2=(em*vel)^2/em
    		FOR k=0,num_time-1 DO BEGIN
			tmp=where(rhopts[*,k] GE 0 AND em[*,k] NE 0 AND gpv[i,*] GT 0)
                        sub_vector[i,k]=total(gpv[i,tmp]*(line.lam_o/c)^2*double(em_vel[tmp,k])^2/(em[tmp,k]))
                ENDFOR
        ENDFOR

	

	output=sub_vector
	IF keyword_set(debug) THEN stop
	RETURN,output
END

FUNCTION genspec_matrix_sub_vector,voxel,vel_voxel,vel_inv,emiss
	
	x=size(voxel)     
	sub_matrix=dblarr(x[1],x[2])
        tmp=where(voxel NE 0)
        sub_matrix[tmp]=vel_voxel[tmp]^2/voxel[tmp]
        sub_vector=sub_matrix#(vel_inv^2/emiss)

        output=sub_vector
        RETURN,output
END
;+
;NAME:
;	GENSPEC_INVERT_MOMENTS
;
;PURPOSE:
;	This is a higer level function that calls other GENPOS and GENSPEC functions to create matrices and then 
;	does a least squares inversion on each moment vector to find the coefficents of expansion for the line emissivity,
;	velocity and ion temperature profiles.
;
;CALLING SEQUENCE:
;	result=GENSPEC_INVERT_MOMENTS(m0,m1,m2,gpv,pos,bfield,line,ves_grid,rhopts,order)
;
;INPUTS:
;	m0:		FLTARR 	of length n, where n is the number of channels, of the zeroth moment of the spectra     ]
;	m1:		FLTARR 	of length n of the first moment								](see GENSPEC_MOMENTS)
;	m2:		FLTARR 	of lenght n of the second moment							]
;	gpv:		FLTARR	n x m  where m is the number of pixels.  These are the volume coefficents
;				generated using GENPOS_VOL_COEFS [meters^3]
;	pos:		FLTARR 	4 x n of the central (or weighted) POS vector for each channels.
;	bfield:		STRUC   of data describing the magnetic field on the grid (See GENPOS_GRID_BFIELD for method)
;			*.br		FLTARR [n_grid, n_time] radial component of the magnetic field [T]
;			*.bphi 		FLTARR [n_grid, n_time] toroidal component of the magnetic field [T]
;			*.bz		FLTARR [n_grid, n_time] veritical component of the magnetic field [T]
;			*.psi_norm	FLTARR [n_grid, n_time] normalized psi values ***not used***
;			*.t_pts		FLTARR [n_time] time points [s] ***not used***
;			*.shot		LONG shot number ***not used***
;	line:		STRUC 	containing information about the line to be calculated
;			*.lam_o		FLT the rest wavelength of the line [units set wavelength scale]        
;			*.del_lam	FLT used in the  wavelength interval [lam_o-del_lam, lam_o+del_lam]    ***not used***
;			*.mass		FLT the mass of ion in AMU (see READ_ATOMIC_MASS('XX')
;	rhopts:		FLTARR	m x t where t is the number of time points.  These are normalized radius (rho) values that correspond 
;				to the m pixels and should be generated using GENPOS_GRID2RMID.
;	ves_grid:	STRUC	the vessel grid generated from GENPOS_GRID using /center
;
;KEYWORD PARAMETERS:
;	solidbody:	/solidbody will generate a sub vector and velocity coef matrix using only w(rho) which is pure toroidal rotation
;
;OUTPUTS:
;	result:		STRUC of the least squares fit of the profile coefficients to the moments;
;				*.c0	FLTARR of length order[0] for the total emissivity profile coeffcients
;				*.c1	FLTARR of length 2*order[1] for the emissivity weighted "velocity" profile coeffcients
;					c1[0:order[1]-1] are the em(rho)*w(rho) coefs
;					c1[order[1]:2*order[1]-1] are the em(rho)*u(rho) coefs (if /solidbody not set)
;				*.c2	FLTARR of length order[2] for the emissivity weighted temperature profile coefficients
;
;PROCEDUE:
;	The matrix equation Ax=b is solved where x are the expansion coefficents, b are the moments (see GENSPEC_MOMENTS) 
;	and A is necessary geometry information.  This is solved using x=LA_INVERT(transpose(A)#A,/double)#transpos(A)#b 
;	which is the least squares inversion.
;	
;	More detailed equations can be found (in the near future) in 'General Equations for Radiometry in Tokamak Plasmas'.
;
;RESTRICTIONS
;	This function calls GENPOS_GPV2COEF_MATRIX, GENSPEC_VEL_COEF_MATRIX and GENSPEC_SUB_VECTOR to calculate the
;	necessary matrices.  All inversions are done using LA_INVERT.
;
;MODIFICATION HISTORY:
;	Written by:	ML Reinke: 8-21-06
;	6-7-07:		ML Reinke: added solidbody optional input to carry thru to GENPSPEC_VEL_COEF_MATRIX and GENSPEC_SUB_VECTOR
;				   added the bfield input for use in these function and updated help file
;	6-8-07:		ML Reinke: fixed a bug that prevented the 2nd moment from having a different order than the 0th moment
;-


FUNCTION genspec_invert_moments,m0,m1,m2,gpv,pos,bfield,line,ves_grid,rhopts,order,bessel=bessel,solidbody=solidbody,debug=debug

	c=3.0e8  		;speed of light
	e=1.60e-19		;conversion for eV -> J
	mconv=1.66e-27		;conversion for amu -> kg

	;generate coeffcient matrices
	coef_0_matrix=genpos_gpv2coef_matrix(gpv,rhopts,order[0],bessel=bessel)
	coef_1_matrix=genspec_vel_coef_matrix(gpv,pos,bfield,line,rhopts,ves_grid,order[1],bessel=bessel,solidbody=solidbody)
	IF order[0] EQ order [2] THEN coef_2_matrix=coef_0_matrix*(line.lam_o/c)^2*(e*1.0e3/(line.mass*mconv))	ELSE $
        	coef_2_matrix=genpos_gpv2coef_matrix(gpv,rhopts,order[2],bessel=bessel)*(line.lam_o/c)^2*(e*1.0e3/(line.mass*mconv))

	coefs_0=la_invert(transpose(coef_0_matrix)#coef_0_matrix,/double)#transpose(coef_0_matrix)#m0
	IF total(pos[2,*]) NE 0 THEN coefs_1=la_invert(transpose(coef_1_matrix)#coef_1_matrix,/double)#transpose(coef_1_matrix)#m1 $
		ELSE coefs_1=fltarr(order[1])
	sub_vector=genspec_sub_vector(gpv,pos,bfield,line,rhopts,ves_grid,coefs_0,order[0],coefs_1,order[1],bessel=bessel,solidbody=solidbody)
	coefs_2=la_invert(transpose(coef_2_matrix)#coef_2_matrix,/double)#transpose(coef_2_matrix)#(m2-sub_vector)

	output={c0:coefs_0, c1:coefs_1, c2:coefs_2}
	IF keyword_set(debug) THEN stop
	RETURN, output
END

;+
;NAME:
;	GENSPEC_MATRIX_INVERT
;
;PURPOSE;
;	This is an upper level function that takes spectral moments and geometry weighting information and calculates
;	emissivity, velocity and ion temperatures for a given impurity.  Some might call this Unleashing The Fury (TM) 
;	but I feel that's being a bit melodramatic.
;
;CALLING SEQUENCE:
;	result=GENSPEC_MATRIX_INVERT(moments,gpv,pos,shot,time,ves_cent,lam_o,z)
;
;INPUTS:
;	moments:	FLTARR 	[3, n_ch, n_time] of the spectral moments
;			       	[0, *, *] is the zeroth moment [photons/s] or other "power" unit
;			       	[1, *, *] is the first moment  [angstroms*photons/s]  consistent units with lam_o
;			       	[2, *, *] is the first moment  [angstroms^2*photons/s]  consistent units with lam_o
;	gpv:		FLTARR 	[n_ch, n_grid] of the GENPOS volume coefficients for each channel.  See GENPOS_VOL_COEFS
;	pos:		FLTARR 	[4,n_ch] of the average POS vectors for each channel.  Used to calculate projection of Doppler shift
;	shot:		LONG 	shot number
;	time:		FLTARR 	[n_time] of the time points.  Note that all geometry info will be calculated at EFIT time points
;			       	so either run EFIT at your time points or calculate the moments at EFIT time points for most accurate
;			       	results
;	ves_cent:	STRUCT 	the vessel grid structure generated from GENPOS_GRID or GRID_VES using /center, consistent with gpv
;	lam_o:		FLT	the unshifted wavelength of the line.  Units consistent with how the moments are calculated.
;	z:		INT	the z of the element that's being observed
;
;OPTIONAL INPUTS:
;	good:		INTARR	[n_ch] with 1's (use) or 0's (do not use) indicating which channels to use during the inversion
;	rhopts:		FLTARR	[n_grid, n_time] of the rho values corresponding to the (R,Z) grid points.  This will be calculated
;				by GENSPEC_MATRIX_INVERT but is an optional output that can be reused on different spectral lines to
;				save computation time.
;	bfield:		STRUC	see GENPOS_GRID_BFIELD for structure but like rhopts this needs only be calculated once per shot
;				and can be taken as an optional output and used in later inversions for different lines.
;	n_rho:		INT	number of points from minimum < rho < 1.0 (inclusive) to be used in inversion DEFAULT: 25
;	rho_vec:	FLTARR 	the actually rho points to be used DEFAULT: not used
;	u_rho_lim:	FLTARR  of the rho value below which the u-profile is forced to zero (helps convergance)  DEFAULT: 0.35
;	eta:		FLT	edge zero weighting factor DEFAULT: 0.0
;	eps_em:		FLT	emissivity profile smoothing factor. DEFAULT: 1.0
;	eps_w:		FLT	toroidal rotation profile smoothing factor. DEFAULT: 1.0
;	eps_u:		FLT	field aligned rotation profile smoothing factor. DEFAULT: 1.0
;	eps_ti:		FLT	ion temperature profile smoothing factor. DEFAULT: 1.0
;	n_iter:		INT 	number of iterations (if finding u-profile) to perform to find u from first moment residual DEFAULT: 50
;	posrev:		INTARR	[n_ch] of 0's and 1's for reversing the pos in GENPOS_LHAT through GENSPEC_GPV2VOXEL_VEL_MATRIX
;	err:		FLTARR	[3,n_ch,n_time] of the error in the moment profiles.  Use output from GENSPEC_LI_MOMENTS
;	n_err_inv:	INT	number of inversion to complete to calculate error.  DEFAULT set in GENPOS_PROFILE_INVERT
;
;KEYWORD PARAMETERS:
;	solidbody:	/solidbody does the first moment inversion forcing the u-profile = 0
;	parallel:	/parallel does the first moment inversion forcing the w-profile = 0 [NOT TESTED YET]
;	nofirst:	/nofirst prevents the matrix inversion from weighting the profile so that it has zero derivative at rho=0
;	debug:		/debug stops the code in various places and just before the RETURN statement
;	quiet:		/quiet suppresses terminal messages displaying computation times of various procedures
;	diff:		/diff is passed onto GENPOS_LHAT through GENSPEC_GPV2VOXEL_VEL_MATRIX
;
;OUTPUTS:
;	result:		STRUC	containing the kinetic profiles and their moment profile checks
;			*.rho		FLTARR 	[n_rho] of the rho values at which each profile has been calculated
;			*.time		FLTARR 	[n_time] of the time points (copy of time input for convience)
;			*.emiss 	FLTARR 	[n_rho, n_time] of the line emissivities [photons/s/m^3] or "power"/m^3
;			*.w		FLTARR	[n_rho, n_time] of the "w/2pi"-profile or toroidal rotation frequency [Hz]
;			*.u		FLTARR  [n_rho, n_time] of the u-profile [m/s/T]
;			*.ti		FLTARR  [n_rho, n_time] of the impurity ion temperature [keV]
;			*.ch		INTARR 	[n_ch] of only the GOOD channels used in the inversion
;			*.brchk		FLTARR 	[3, n_ch, n_time] of the moments calculated from the calculated profiles.  This can be
;						be used to check the quality of the output versus the input data.
;			*.inverr	FLTARR	[3,n_rho,n_time] of the absolute inversion error in the kinetic profiles.  Will 
;						be zeros if optional input err is not specified.  Units for each are in the same as above.
;
;OPTIONAL OUTPUTS:
;	rhopts:		FLTARR [n_grid, n_time] the output of GENPOS_GRID2RMID(ves_cent,shot,tpts=time,/rho) for future use
;	bfield:		STRUC the output of GENPOS_GRID_BFIELD(ves_cent,shot,t_pts=time) for future use
;	moment1_conv	FLTARR [n_iter, n_time] showing the convergence of the STDEV of the first moment residual over the
;				mean of the first moment * 100
;
;PROCEDURE:
;	This program calculates the spatial weighting matrices using GENPOS_GPV2VOXEL_MATRIX and GENPOS_GPV2VOXEL_VEL_MATRIX and then
;	inverts them using GENPOS_PROFILE_INVERT.  The eps values for each inversion, at each time point are weighted to the max value of
;	the voxel matrix.  For example, the GENPOS_PROFILE_INVERT eps input for the emissivity profile, at each time point i is: 
;		
;		double(eps_em*(max(voxel[*,*,i]))^2)
;
;	To determine the u-profile and w-profile from the same data set, the assumption is made that it is small compared to the w-profile.
;	The inversion is done assuming all of moment[1,*,i] is due to the w-profile.  Then the residual is found from brchk and the u-profile is
;	determined from this residual.  This brchk is subtracted from  moment[1,*,i] and new w-profile is calculated.  This
;	process is iterated n_iter number of times.  This type of analysis assumes that your GPV/POS are setup such that you can distinguish
;	the solidbody and parallel parts (ie a full ~poloidal cross section like HIREX-SR).
;
;RESTRICTIONS:
;	Use GENIE_INI.bat to compile the GENPOS and GENSPEC functions in the correct order.  You should know what you're doing before
;	using this function.  It's the third rail of spectroscopy.
;
;MODFICATION HISTORY:
;	Written by:	ML Reinke - 7/18/07
;	7/23/07:	ML Reinke - changed default n_rho to 25 from n_ch
;	8/02/07:	ML Reinke - added the nofirst keyword and  made the function calculate the
;				    minimum rho and use as optional input in GENPOS_GPV2VOXEL_MATRIX
;	8/03/07:	ML Reinke - fixed a bug that was performing the sub_vector operation incorrectly.  Now 
;				    GENSPEC_MATRIX_SUB_VECTOR is called at each time interval.
;	8/03/07:	ML Reinke - updated to be able to calculate a u-profile using an iterative method on the first moment residual
;				    added keywords: solidbody and parallel, optional inputs: eps_u, u_rho_lim and n_iter,
;				    outputs: *.u, optional output: moment1_conv
;	8/15/07:	ML Reinke - added the eta optional input for edge zero
;	2/01/08:	ML Reinke - fixed voxel indexing bugs in the Ti profile that caused time range dependant profiles effects.
;	2/12/08:	ML Reinke - fixed indexing bug for the u profile brightness check.
;	2/18/08:	ML Reinke - added diff and posrev for GENPOS_LHAT called in GENSPEC_GPV2VOXEL_VEL_MATRIX
;	10/29/08:	ML Reinke - added the iterate keyword to switch between methods for velocity inverstion
;	10/30/08:	ML Reinke - added the err,n_err_inv optional inputs for calculating the inversion error profiles
;-	


FUNCTION genspec_matrix_invert,moments,gpv,pos,shot,time,ves_cent,lam_o,z,good=good,rhopts=rhopts,bfield=bfield,n_rho=n_rho,rho_vec=rho_vec,eta=eta,$
		eps_em=eps_em,eps_w=eps_w,eps_u=eps_u,u_rho_lim=u_rho_lim,eps_ti=eps_ti,n_iter=n_iter,nofirst=nofirst,solidbody=solidbody,$
                parallel=parallel,moment1_conv=moment1_conv,debug=debug,quiet=quiet,rhosol=rhosol,diff=diff,posrev=posrev,iterate=iterate,err=err,n_err_inv=n_err_inv

	c=3.0e8  			;speed of light
	e=1.60e-19			;conversion for eV -> J
	mconv=1.66e-27			;conversion for amu -> kg
	mass=read_atomic_mass(z)

	n_ch=n(pos[0,*])+1
	n_time=n(time)+1
	ch=indgen(n_ch)
        ;debug=1
	IF NOT keyword_set(eps_em) THEN eps_em=1.0
	IF NOT keyword_set(eps_w) THEN eps_w=1.0
        IF NOT keyword_set(eps_u) THEN eps_u=1.0
	IF NOT keyword_set(eps_ti) THEN eps_ti=1.0
        IF NOT keyword_set(eta) THEN eta=0.0
        IF NOT keyword_set(u_rho_lim) THEN u_rho_lim=0.1
        IF NOT keyword_set(n_iter) THEN n_iter=50

	IF NOT keyword_set(n_rho) THEN IF NOT keyword_set(rho_vec) THEN n_rho=22 ELSE n_rho=n(rho_vec)+1
	IF NOT keyword_set(good) THEN good=intarr(n_ch)+1
	start_time=systime(/seconds)
	IF NOT keyword_set(rhopts) THEN rhopts=genpos_grid2rmid(ves_cent,shot,tpts=time,/rho,sol=rhosol)
	ctime=systime(/seconds)-start_time
	IF NOT keyword_set(quiet) THEN print, 'VES_CENT converted to RHOPTS: '+num2str(ctime,dp=2)
	start_time=systime(/seconds)
	IF NOT keyword_set(bfield) THEN bfield=genpos_grid_bfield(ves_cent,shot,t_pts=time)
	ctime=systime(/seconds)-start_time
	IF NOT keyword_set(quiet) THEN print, 'BFIELD calculated at VES_CENT points: '+num2str(ctime,dp=2)
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
	start_time=systime(/seconds)
	voxel=genpos_gpv2voxel_matrix(gpv[tmp,*],rhopts,rho_vec=rho_vec,n_rho=n_rho,rho_min=rho_min)
        IF NOT keyword_set(quiet) THEN print, 'RHO values from '+num2str(min(rho_vec),dp=3)+' < rho < '+num2str(max(rho_vec),dp=3)+' selected'
	ctime=systime(/seconds)-start_time
	IF NOT keyword_set(quiet) THEN print, 'VOXEL arrays calculated: '+num2str(ctime,dp=2)
	start_time=systime(/seconds)
	vel_voxel=genspec_gpv2voxel_vel_matrix(gpv[tmp,*],rhopts,pos,ves_cent,bfield,rho_vec=rho_vec,n_rho=n_rho,solidbody=solidbody,posrev=posrev,diff=diff)
	ctime=systime(/seconds)-start_time
	IF NOT keyword_set(quiet) THEN print, 'VEL_VOXEL arrays calculated: '+num2str(ctime,dp=2)
	IF keyword_set(debug) THEN stop

	emiss=fltarr(n_rho,n_time)
	w=fltarr(n_rho,n_time)
        u=fltarr(n_rho,n_time)
        moment1_conv=fltarr(n_iter,n_time)
	ti=fltarr(n_rho,n_time)
	brchk=fltarr(4,n_good,n_time)
        sub_vector=fltarr(n_good,n_time)
        inverr=fltarr(4,n_rho,n_time)
        bp=fltarr(n_rho,n_time)

	start_time=systime(/seconds)
	FOR i=0,n_time-1 DO BEGIN
            	;emissivity inversions
		a=max(voxel[*,*,i])
                IF keyword_set(err) THEN moment_err = err[0,tmp,i] ELSE moment_err=0
		emiss[*,i]=genpos_profile_invert(reform(moments[0,tmp,i]),voxel[*,*,i],double(eps_em*a^2),brchk=br,nofirst=nofirst,eta=double(eta*a^2),$
			err=moment_err,n_err_inv=n_err_inv)
		brchk[0,*,i]=br.mom
                inverr[0,*,i]=br.inverr

                ;velocity inversions
		IF NOT keyword_set(parallel) THEN BEGIN
                    	a=max(vel_voxel.w[*,*,i])
                        IF keyword_set(err) THEN moment_err = err[1,tmp,i] ELSE moment_err=0
                	w_vel_inv=genpos_profile_invert(reform(moments[1,tmp,i]),vel_voxel.w[*,*,i],double(eps_w*a^2),brchk=br_w,nofirst=nofirst,$
                        	eta=double(eta*a^2),err=moment_err,n_err_inv=n_err_inv)
                        brchk[1,*,i]=br_w.mom
                        sub_vector[*,i]=genspec_matrix_sub_vector(voxel[*,*,i],vel_voxel.w[*,*,i],w_vel_inv,emiss[*,i])
                        u_vel_inv=fltarr(n_rho)
                ENDIF ELSE BEGIN
                        a=max(abs(vel_voxel.u[*,*,i]))
                        u_vel_inv=genpos_profile_invert(reform(moments[1,tmp,i]),vel_voxel.u[*,*,i],double(eps_u*a^2),brchk=br_u,nofirst=nofirst,$
                        	eta=double(eta*a^2))
                        brchk[1,*,i]=br_u.mom
                        sub_vector[*,i]=genspec_matrix_sub_vector(voxel[*,*,i],vel_voxel.u[*,*,i],u_vel_inv,emiss[*,i])
                        w_vel_inv=fltarr(n_rho)
                ENDELSE
                IF NOT keyword_set(solidbody) AND NOT keyword_set(parallel) THEN BEGIN
                	IF keyword_set(iterate) THEN BEGIN
	                     	i_rho=ipt(rho_vec, u_rho_lim)
                                moment1_u=reform(moments[1,tmp,i])-br_w.mom
                                a=max(vel_voxel.u[*,*,i])
                                u_vel_inv=genpos_profile_invert(moment1_u,vel_voxel.u[*,*,i],double(eps_u*a^2),brchk=br_u,nofirst=nofirst,err=moment_err,n_err_inv=n_err_inv)
                                IF i_rho NE 0 AND i_rho NE -1 THEN u_vel_inv[0:i_rho]=0.0
                                FOR j=0,n_iter-1 DO BEGIN
                         		moment1_w=reform(moments[1,tmp,i])-vel_voxel.u[*,*,i]#u_vel_inv
                                        a=max(vel_voxel.w[*,*,i])
                                        w_vel_inv=genpos_profile_invert(moment1_w,vel_voxel.w[*,*,i],double(eps_w*a^2),brchk=br_w,nofirst=nofirst,$
                                                     eta=double(eta*a^2),err=moment_err,n_err_inv=n_err_inv)
                                        moment1_u=reform(moments[1,tmp,i])-br_w.mom
                                        a=max(vel_voxel.u[*,*,i])
                                        u_vel_inv=genpos_profile_invert(moment1_u,vel_voxel.u[*,*,i],double(eps_u*a^2),brchk=br_u,nofirst=nofirst,$
                                                     eta=double(eta*a^2),err=moment_err,n_err_inv=n_err_inv)
                                        moment1_conv[j,i]=stdev(moments[1,*,i]-br_w.mom-br_u.mom)/mean(moments[1,*,i])*100.0
                                        IF i_rho NE 0 AND i_rho NE -1 THEN u_vel_inv[0:i_rho]=0.0
                                ENDFOR
                                brchk[3,*,i]=br_w.mom+br_u.mom
                         ENDIF ELSE BEGIN
  	                 	vel_vox_tot=[[vel_voxel.w[*,*,i]],[vel_voxel.u[*,*,i]]]
                                a=max(abs(vel_vox_tot))
                                vel_inv=genpos_profile_invert(reform(moments[1,tmp,i]),vel_vox_tot,double(eps_w*a^2),brchk=br,nofirst=nofirst,eta=double(eta*a^2),$
                                	err=moment_err,n_err_inv=n_err_inv)
                                w_vel_inv=vel_inv[0:n_rho-1]
                                u_vel_inv=vel_inv[n_rho:*]
                                brchk[3,*,i]=br.mom                                
                                sub_vector[*,i]=genspec_matrix_sub_vector(voxel[*,*,i],vel_voxel.w[*,*,i],w_vel_inv,emiss[*,i])
                         ENDELSE 
                ENDIF
                u[*,i]=u_vel_inv/(lam_o/c)/emiss[*,i]
                w[*,i]=w_vel_inv/(lam_o/c)/emiss[*,i]/(2.0*!pi)
                IF keyword_set(solidbody) THEN BEGIN
                	inverr[1,*,i]=w[*,i]*sqrt(br_w.inverr^2/w_vel_inv^2+(inverr[0,*,i])^2/(emiss[*,i])^2)
                ENDIF ELSE BEGIN
                    	IF keyword_set(iterate) THEN BEGIN
                        	inverr[1,*,i]=w[*,i]*sqrt(br_w.inverr^2/w_vel_inv^2+(inverr[0,*,i])^2/(emiss[*,i])^2)
                                inverr[3,*,i]=u[*,i]*sqrt(br_u.inverr^2/u_vel_inv^2+(inverr[0,*,i])^2/(emiss[*,i])^2)
                        ENDIF ELSE BEGIN
                           	inverr[1,*,i]=w[*,i]*sqrt(br.inverr[0:n_rho-1]^2/w_vel_inv^2+(inverr[0,*,i])^2/(emiss[*,i])^2)
                                inverr[3,*,i]=u[*,i]*sqrt(br.inverr[n_rho:*]^2/u_vel_inv^2+(inverr[0,*,i])^2/(emiss[*,i])^2)                     	
                        ENDELSE
                ENDELSE 

                ;ion temperature inversions
		a=max(voxel[*,*,i])  
                IF keyword_set(err) THEN moment_err = err[2,tmp,i] ELSE moment_err=0 
		ti_inv=genpos_profile_invert(reform(moments[2,tmp,i])-sub_vector[*,i],voxel[*,*,i],double(eps_ti*a^2),brchk=br,nofirst=nofirst,$
                	eta=double(eta*a^2),err=moment_err,n_err_inv=n_err_inv)
		brchk[2,*,i]=br.mom+sub_vector[*,i]
		ti[*,i]=ti_inv/((lam_o/c)^2*(e*1.0e3/(mass*mconv)))/emiss[*,i]
                inverr[2,*,i]=ti[*,i]*sqrt(br.inverr^2/ti_inv^2+(inverr[0,*,i])^2/(emiss[*,i])^2)
	ENDFOR
	ctime=systime(/seconds)-start_time
	IF NOT keyword_set(quiet) THEN print, 'inversions completed for '+num2str(n_time)+' time points: '+num2str(ctime,dp=2)

	output={rho:rho_vec,time:time,emiss:emiss,w:w,u:u,ti:ti,ch:good_ch,brchk:brchk,inverr:inverr}
	IF keyword_set(debug) THEN stop
	RETURN,output
	
END

;+
;NAME:
;	GENSPEC_SHELLINVERT
;
;PURPOSE;
;	This is an upper level function that takes spectral moments and geometry weighting information and calculates
;	emissivity, velocity and ion temperatures for a given impurity using a shell inversion technique.
;
;CALLING SEQUENCE:
;	result=GENSPEC_SHELLINVERT(moments,pos,shot,time,lam_o,z)
;
;INPUTS:
;	moments:	FLTARR 	[3, n_ch, n_time] of the spectral moments
;			       	[0, *, *] is the zeroth moment [photons/s/m^2/str] or other "brightness" unit
;			       	[1, *, *] is the first moment  [angstroms*"bright"] consistent units with lam_o  ***will get altered***
;			       	[2, *, *] is the first moment  [angstroms^2*"bright"] consistent units with lam_o
;	pos:		FLTARR 	[4,n_ch] of the POS vectors for each channel.
;	shot:		LONG 	shot number
;	time:		FLTARR 	[n_time] of the time points.  Note that all geometry info will be calculated at EFIT time points
;				and linearly interpolated between them
;	lam_o:		FLT	the unshifted wavelength of the line.  Units consistent with how the moments are calculated.
;	z:		INT	the z of the element that's being observed
;
;OPTIONAL INPUTS:
;	good:		INTARR	[n_ch] with 1's (use) or 0's (do not use) indicating which channels to use during the inversion
;	npts:		INT	number of points to be used in SHELLINVERT DEFAULT: 25
;	redge:		FLT	edge of the emissivity in major radius (will be adjusted for shell inversion) DEFAULT: calculated for outermost choord
;	eps_em:		FLT	emissivity profile smoothing factor. DEFAULT: 0.1
;	eps_w:		FLT	toroidal rotation profile smoothing factor. DEFAULT: 0.1
;	eps_u:		FLT	field aligned rotation profile smoothing factor. DEFAULT: 0.1
;	eps_ti:		FLT	ion temperature profile smoothing factor. DEFAULT: 0.1
;	posrev:		INTARR	[n_ch] of 0's and 1's indicating whether or not the POS has been reveresed (for velocity information)
;
;KEYWORD PARAMETERS:
;	solidbody:	/solidbody does the first moment inversion forcing the u-profile = 0
;	parallel:	/parallel does the first moment inversion forcing the w-profile = 0
;	debug:		/debug stops the code in various places and just before the RETURN statement
;	quiet:		/quiet suppresses terminal messages displaying computation times of various procedures
;
;OUTPUTS:
;	result:		STRUC	containing the kinetic profiles and their moment profile checks
;			*.r		FLTARR 	[npts,n_time] of the midplane major radius values [m]
;			*.time		FLTARR 	[n_time] of the time points (copy of time input for convience)
;			*.emiss 	FLTARR 	[npts, n_time] of the line emissivities [photons/s/m^3] or "power"/m^3
;			*.w		FLTARR	[npts, n_time] of the "w/2pi"-profile or toroidal rotation frequency [Hz]
;			*.u		FLTARR  [npts, n_time] of the u*gradPsi-profile m^2/s (divide by R to get poloidal velocity)
;			*.ti		FLTARR  [npts, n_time] of the impurity ion temperature [keV]
;			*.ch		INTARR 	[n_good+1] of only the GOOD channels used in the inversion plus the forced zero at redge
;			*.brchk		FLTARR 	[3, n_good+1, n_time] of the moments calculated from the calculated profiles.  This can be
;						be used to check the quality of the output versus the input data.
;			*.rt		FLTARR	[n_ch,n_time] of the tangency radii [m] (time evolving for /parallel only)
;
;PROCEDURE:
;	This program calculates the spatial weighting matrices using CALC_GEOM in SHELLINVERT.PRO and then
;	inverts them using SHELLINVERT.  For a poloidal view, the "tangency radii" are calculated from the POS vector using
;	GENPOS_POLOIDALPOS2RT while for toroidal views, the tangency radii are in the POS vector.
;
;	The first moment profile will be altered so that it can be compared with the BRCHK.  See Section 7.1 in General Radiometry
;
;RESTRICTIONS:
;	Use GENIE_INI.bat to compile the GENPOS and GENSPEC functions in the correct order.  This also requires SHELLINVERT.
;
;MODFICATION HISTORY:
;	Written by:	ML Reinke - 2/22/08 - form adapted from GENSPEC_MATRIX_INVERT
;
;-

FUNCTION genspec_shellinvert,moments,pos,shot,time,lam_o,z,good=good,eta=eta,eps_em=eps_em,eps_w=eps_w,eps_u=eps_u,eps_ti=eps_ti,$
		npts=npts,redge=redge,solidbody=solidbody,parallel=parallel,debug=debug,quiet=quiet,posrev=posrev

	c=3.0e8  			;speed of light
	e=1.60e-19			;conversion for eV -> J
	mconv=1.66e-27			;conversion for amu -> kg
	mass=read_atomic_mass(z)

	n_ch=n(pos[0,*])+1
	IF NOT keyword_set(eps_em) THEN eps_em=0.1
	IF NOT keyword_set(eps_w) THEN eps_w=0.1
        IF NOT keyword_set(eps_u) THEN eps_u=0.1
	IF NOT keyword_set(eps_ti) THEN eps_ti=0.1
	IF NOT keyword_set(good) THEN good=intarr(n_ch)+1
	IF NOT keyword_set(npts) THEN npts=25
	IF NOT keyword_set(posrev) THEN posrev=fltarr(n_ch)+1.0

	n_ch=n(pos[0,*])+1
	n_time=n(time)+1
	ch=indgen(n_ch)
	n_good=total(good)
	good_ch=ch[where(good EQ 1)]
	IF keyword_set(parallel) THEN BEGIN
		genpos_poloidalpos2rt,shot,pos,rt,rt_time,debug=debug,rout=rout,zmag=zmag 
		rout_interp=interpol(rout,rt_time,time)
		zmag_interp=interpol(zmag,rt_time,time)
	ENDIF ELSE BEGIN
		rt=reform(pos[2,*])
		rt_time=time
		rout_interp=fltarr(n(time)+1)
	ENDELSE

	IF keyword_set(debug) THEN stop

	emiss=fltarr(npts,n_time)
	w=fltarr(npts,n_time)
        u=fltarr(npts,n_time)
	ti=fltarr(npts,n_time)
	rmaj=fltarr(npts,n_time)
	brchk=fltarr(3,n_good+1,n_time)
	brchk_rad=fltarr(n_good+1,n_time)
        sub_vector=fltarr(n_good,n_time)

	start_time=systime(/seconds)
	FOR i=0,n_time-1 DO BEGIN
		ipt=ipt(rt_time,time[i])
		bright_in=reform(moments[0,*,i])
		IF keyword_set(parallel) THEN BEGIN
			rtang_in=rt[*,ipt]+(rt[*,ipt+1]-rt[*,ipt])/(rt_time[ipt+1]-rt_time[ipt])*(time[i]-rt_time[ipt])
			IF keyword_set(redge) THEN redge_i=redge-rout_interp[i] ELSE redge_i=rtang_in[n_ch-1]+2.0*(rtang_in[n_ch-1]-rtang_in[n_ch-2])
		ENDIF ELSE BEGIN
			rtang_in=rt
			IF NOT keyword_set(redge) THEN redge_i=rtang_in[n_ch-1]+2.0*(rtang_in[n_ch-1]-rtang_in[n_ch-2])
		ENDELSE
		shellinvert,bright_in,rtang_in,emiss_i,radii,good=good,npts=npts,eps=eps_em,redge=redge_i,brchk=brchk_i
		rmaj[*,i]=radii+rout_interp[i]
		emiss[*,i]=emiss_i
		brchk[0,*,i]=brchk_i.br
		brchk_rad[*,i]=brchk_i.rad

		IF keyword_set(parallel) THEN BEGIN
			bright_in=reform(moments[1,*,i])*c/(lam_o*rtang_in)*(rout_interp[i]+rtang_in*sin(posrev*pos[3,*]))
			shellinvert,bright_in,rtang_in,emiss_i,radii,good=good,npts=npts,eps=eps_u,redge=redge_i,brchk=brchk_i
			u[*,i]=emiss_i/emiss[*,i]*radii
			brchk[1,*,i]=brchk_i.br
			moments[1,*,i]*=c/(lam_o*rtang_in)*(rout_interp[i]+rtang_in*sin(posrev*pos[3,*]))
		ENDIF
		IF keyword_set(solidbody) THEN BEGIN
			bright_in=reform(moments[1,*,i])*c/(lam_o*rtang_in)
			shellinvert,bright_in,rtang_in,emiss_i,radii,good=good,npts=npts,eps=eps_ti,redge=redge_i,brchk=brchk_i
			w[*,i]=emiss_i/emiss[*,i]/(2.0*!pi)
			brchk[1,*,i]=brchk_i.br
			moments[1,*,i]*=c/(lam_o*rtang_in)
                ENDIF
		;currently no subtraction vector calculated
       		bright_in=reform(moments[2,*,i])
		shellinvert,bright_in,rtang_in,emiss_i,radii,good=good,npts=npts,eps=eps_em,redge=redge_i,brchk=brchk_i
		ti[*,i]=emiss_i/((lam_o/c)^2*(e*1.0e3/(mass*mconv)))/emiss[*,i]
		brchk[2,*,i]=brchk_i.br

	ENDFOR
	ctime=systime(/seconds)-start_time
	IF NOT keyword_set(quiet) THEN print, 'inversions completed for '+num2str(n_time)+' time points: '+num2str(ctime,dp=2)

	output={r:rmaj,time:time,emiss:emiss,w:w,u:u,ti:ti,ch:good_ch,brchk:brchk,rt:rt}
	IF keyword_set(debug) THEN stop
	RETURN,output
	
END

