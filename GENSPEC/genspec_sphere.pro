;+
;NAME:
;	XCIS_LINESHAPE
;
;-

PRO xcis_lineshape,x,dphi,th_o=th_o,nphi=nphi,delta=delta,ti=ti,mean=mean,width=width,pix=pix,pow=p,noplot=noplot
	xin=x	
	R=1.442
	z=18

	lam_o=3.9494					;w-line in Angstroms
	IF NOT keyword_set(ti) THEN ti=1.0e3		;ion temperature in eV
	c=2.998e8 					;speed of light
	e=1.602e-19					;conversion for eV -> J
	mconv=1.661e-27					;conversion for amu -> kg
	m=read_atomic_mass(z)*mconv

	w=sqrt(lam_o^2*e*ti/(m*c^2)) 			;linewidth in Angstroms
	twod=4.56216					;2d spacing in angstroms
	omega=w/twod  					;normalized linewidth

	th_o=59.961/180*!pi
	IF NOT keyword_set(delta) THEN delta=0.0		;shift of line center in mAngstroms
	th_o_new=asin((lam_o+delta*1.0e-3)/twod)
	delta=delta*1.0e-3/twod					;normalized line shift	

	x+=R*(cos(th_o_new)-cos(th_o))*0.7

	IF NOT keyword_set(nphi) THEN nphi=100
	phi=make(-dphi,dphi,nphi)
	n_x=n(x)+1
	tanbragg=fltarr(n_x,nphi)
	sinbragg=fltarr(n_x,nphi)
	cosbragg=fltarr(n_x,nphi)
	sinbeta=fltarr(n_x,nphi)
	lsq=fltarr(n_x,nphi)
	FOR i=0,n_x-1 DO BEGIN
		FOR j=0,nphi-1 DO BEGIN
			tanbragg[i,j]=(1.0/(cos(th_o)+x[i])-cos(th_o+phi[j]))/sin(th_o+phi[j])
			sinbragg[i,j]=tanbragg[i,j]/sqrt(1.0+tanbragg[i,j]^2)
			cosbragg[i,j]=1.0/sqrt(1.0+tanbragg[i,j]^2)
			sinbeta[i,j]=cosbragg[i,j]*cos(th_o+phi[j])+sinbragg[i,j]*sin(th_o+phi[j])
			lsq[i,j]=1+(cos(th_o)+x[i])^2-2*(cos(th_o)+x[i])*cos(th_o+phi[j])

		ENDFOR
	ENDFOR
	dp=sinbragg*sinbeta/lsq*exp(-(sinbragg-sin(th_o)-delta)^2/(2*omega^2))
	IF NOT keyword_set(noplot) THEN genplt,dp,x,phi

	n_x=n(x)+1
	p=fltarr(n_x)
	FOR i=0,n_x-1 DO p[i]=int_tabulated(phi,dp[i,*])

	pix=x*R/172.0e-6
	mean=int_tabulated(pix,pix*p)/int_tabulated(pix,p)
	width=int_tabulated(pix,(pix-mean)^2*p)/int_tabulated(pix,p)
	
	IF NOT keyword_set(noplot) THEN BEGIN
		openwin,0
		p/=max(p)
		plot,pix,p,xtit='Pilatus Pixel # (172 ['+n2g('mu')+'m])',ytit='Normalized Power Deposition Profile',chars=1.2,/xsty,$
			tit='d'+n2g('phi')+'='+num2str(dphi*180/!pi,dp=2)+' [deg]   T!li!n='+num2str(ti*1.0e-3,dp=1)+' [keV]'
		xyouts,pix[maxloc(p)]-15,0.8,' MEAN: '+num2str(mean,dp=4),chars=1.4
		xyouts,pix[maxloc(p)]-15,0.65,'WIDTH: '+num2str(width,dp=1),chars=1.4
	ENDIF
	x=xin
END

PRO plot_xcis_lineshape,load=load,noplot=noplot
	IF !d.name EQ 'PS' THEN ps=1 ELSE ps=0
	R=1.442
	ti=0.05e3			;ion temperature [eV]
	cry_width=[0.1,1,5,10,15,20,25,30,35,40,45,50,55,60]		;1/2 width of the crystal [mm]
	dlam=[0,5,10,15,20,25,30,35,40,45,50,55,60,65,70]	;shift [mAng] from focused line
	save_path='/home/mlreinke/hirex_sr/genspec_sphere_lineshape.dat'

	n_w=n(cry_width)+1
	n_l=n(dlam)+1
	
	mean=fltarr(n_w,n_l)
	sigma=fltarr(n_w,n_l)
	n_x=600
	x=make(-2.0e-3,2.0e-3,n_x)
	pixel=fltarr(n_w,n_l,n_x)
	line=fltarr(n_w,n_l,n_x)

	IF NOT keyword_set(load) THEN BEGIN
		FOR i=0,n_w-1 DO BEGIN
			print, i
			FOR j=0,n_l-1 DO BEGIN
				xcis_lineshape,x,(cry_width[i]*1.0e-3)/R,ti=ti,delta=dlam[j],mean=m,width=s,pix=pix,pow=pow,noplot=noplot
				mean[i,j]=m
				sigma[i,j]=s
				pixel[i,j,*]=pix
				line[i,j,*]=pow
			ENDFOR
		ENDFOR
		save,mean,sigma,pixel,line,filename=save_path
	ENDIF ELSE restore,save_path
	
	labels={ilab:'Crystal 1/2 Width [mm]',jlab:'Line Center Shift from Focus [mAng]',klab:'Approx. T!lINST!n [eV]',ctit:'',itit:'',jtit:''}
	genplt,((sigma/sigma[0,0])^2-1.0)*ti,cry_width,dlam,ncntrs=25,cct=39,lab=labels,ps=ps,ir=[0,40],jr=[0,40]
	
	openwin,0
	!p.multi=[0,0,2]
	plot,dlam,mean[0,*],/xsty,/ysty,chars=1.2,ytit='Line Center [pix]',yr=[min(mean),0]
	oplot,dlam,mean[n_w-1,*],linestyle=2.0

	m=(mean[0,1]-mean[0,0])/(dlam[1]-dlam[0])
	meanfit=m*(dlam-dlam[0])+mean[0,0]
	mu_small=mean[0,*]-meanfit
	mu_large=mean[n_w-1,*]-meanfit
	minplt=min(mu_small) < min(mu_large)
	maxplt=max(mu_small) > max(mu_large)

	plot,dlam,mu_small,yr=[minplt,maxplt]*1.03,/ysty,/xsty,xtit='Wavelength from Focus [mAng]',ytit=n2g('Delta')+' From Linear [pix]',chars=1.2
	oplot,dlam,mu_large,linestyle=2.0
	xyouts,10,-3,'Crystal 1/2 Width [mm]'
	xyouts,10,-4,'DASHED: '+num2str(cry_width[n_w-1],dp=1)
	xyouts,10,-5,' SOLID: '+num2str(cry_width[0],dp=1)
	!p.multi=0
	

	openwin,1
	pts=[0,3,5,7,9,11,13]
	colors=[0,30,70,100,130,150,200]
	plot,pixel[0,0,*]*0.172,line[0,0,*],xr=[-4,11]*0.172,/xsty,yr=[0,1.03],/ysty,xtit='Detector Position [mm]',ytit='Normalized Power on Detector',chars=1.2,$
		tit='Line Shape Varition @ Focus'
	xyouts,3.75*0.172,0.9,'Crystal 1/2 Width [mm]',chars=1.4
	xyouts,10*0.172,0.1,'T!li!n = '+num2str(ti*1.0e-3,dp=2)+ ' [kev]',orient=90

	FOR i=0,n(pts) DO BEGIN
		oplot,pixel[pts[i],0,*]*0.172,line[pts[i],0,*],color=colors[i]
		xyouts,5*0.172,0.85-i*0.075,num2str(cry_width[pts[i]],dp=1), color=colors[i]
	ENDFOR
	openwin,2
	pts=[0,3,5,7,9,11,13]
	colors=[0,30,70,100,130,150,200]
	lpt=10
	plot,pixel[0,lpt,*]*0.172,line[0,lpt,*],/xsty,xr=[-167,-152]*0.172,yr=[0,1.03],/ysty,xtit='Detector Position [mm]',ytit='Normalized Power on Detector',chars=1.2,$
		tit='Line Shape Varition @ '+n2g('delta')+n2g('lambda')+'='+num2str(dlam[lpt],1)+' [mAng]'
	xyouts,-159*0.172,0.9,'Crystal 1/2 Width [mm]',chars=1.4
	xyouts,-153*0.172,0.1,'T!li!n = '+num2str(ti*1.0e-3,dp=2)+ ' [kev]',orient=90
	FOR i=0,n(pts) DO BEGIN
		oplot,pixel[pts[i],lpt,*]*0.172,line[pts[i],lpt,*],color=colors[i]
		xyouts,-157*0.172,0.85-i*0.075,num2str(cry_width[pts[i]],dp=1), color=colors[i]
	ENDFOR
	
	stop
END


PRO genspec_sphere_power,info,ti=ti,lam_o=lam_o,nphi=nphi,npsi=npsi,path=path,pow=pow,dp=dp
	IF keyword_set(path) THEN info=genpos_spherical_info(path)
	IF NOT keyword_set(ti) THEN ti=[1.0e3]		;[keV]
	IF NOT keyword_set(lam_o) THEN lam_o=[3.946]	;Ar+16 w-line in Angstroms
	IF NOT keyword_set(nphi) THEN nphi=30		
	IF NOT keyword_set(npsi) THEN npsi=10
	n_lines=n(lam_o)+1

	;detector info
	xi=info.det.xi
	zeta=info.det.zeta
	n_pix=n(xi)+1
	det_pos=fltarr(2,n_pix)
	det_pos[0,*]=info.det.xi
	det_pos[1,*]=info.det.zeta
	x0=info.det.x0
	x1=info.det.x1
	x2=info.det.x2
        norm=crossp((x2-x0),(x1-x0))
        na_hat=norm/sqrt(total(norm*norm))
	xyz_pix=genpos_det2xyz(x0,x1,x2,det_pos)	;calculate xyz of all detector positions

	;mirror info
	mdy=info.m.size[0]		;dy of mirror [m]
	mdz=info.m.size[1]		;dz of mirror [m]
	R=info.m.rad
	twod=info.m.bragg.twod		;twod in Ang

	dphi=atan(mdy/2.0/R)
	phi=make(-dphi,dphi,nphi)
	dpsi=atan(mdz/2.0/R)
	psi=make(!pi/2.0-dpsi,!pi/2.0+dpsi,npsi)

	e=1.61e-19					;conversion for Joules/eV
	m=40.0*1.66e-27					;argon mass in kg
	c=2.998e8					;speed on light m/s

	w=sqrt(lam_o^2*e*ti/(m*c^2)) 			;linewidth in Angstroms
	omega=w/twod  					;normalized linewidth
	delta=lam_o/twod				;normalized line shift	

	;initialize arrays
	dp=fltarr(nphi,npsi,n_pix)
	pow=fltarr(n_pix)

	FOR i=0L,n_pix-1 DO BEGIN
		xa=xyz_pix[*,i]
		IF i MOD 1000 EQ 0 THEN print, 'Pixel #: '+num2str(i,1)+' of '+num2str(n_pix-1,1)
		FOR j=0,npsi-1 DO BEGIN
			rproj=R*sin(psi[j])
			s_z=R*cos(psi[j])
			diff_psi=psi[1]-psi[0]	
			FOR k=0,nphi-1 DO BEGIN
				xs=[R-rproj*cos(phi[k]),rproj*sin(phi[k]),s_z]
				norm=[R,0.0,0.0]-xs
				ns_hat=norm/sqrt(total(norm*norm))
				l=xa-xs
				l_hat=l/sqrt(total(l*l))
				dU=R^2*sin(psi[j])*total(ns_hat*l_hat)*(-1.0)*total(na_hat*l_hat)/(4.0*!pi*total(l*l))
				FOR m=0,n_lines-1 DO dp[k,j,i]+=dU*exp(-(total(ns_hat*l_hat)-delta[m])^2/(2.0*omega[m]^2))/(omega[m]*sqrt(2.0*!pi))
			ENDFOR
			pow[i]+=int_tabulated(phi,dp[*,j,i])*diff_psi
		ENDFOR
	ENDFOR
END


;+
;NAME:
;	GENSPEC_SPHERE_OBSEMISS
;
;PURPOSE:
;	This function calculates the "observed" emissivity profile for a Bragg reflector along a line of sight through the plasma
;
;CALLING SEQUENCE:
;	result=GENSPEC_SPHERE_OBSEMISS(lam_b,lam_o,em,vel,wid,lhat)
;
;INPUTS:
;	lam_b	FLT	the wavelength that satisfies the Bragg equation [Ang]
;	lam_o	FLTARR	[nlines] of the unshifted line centers for the spectrum of interest
;	em	FLTARR	[ns,nlines] of the emissivity for each line along the line of sight [AU]
;	vel	FLTARR	[ns,3] of the velocity [v_R, v_phi,v_Z] along the line of sight[m/s]
;	ti	FLTARR	[ns] of the line width along the line of sight [keV]
;	lhat	FLTARR	[ns,3] of the unit vector [R,phi,Z] along the line of sight
;	
;OUTPUTS:
;	result	FLTARR	[ns] of the observed emissivity along the line of sight [AU]
;
;PROCEDURE:
;	This function calculates the Doppler line shape at each spatial point in the plasma along the line of sight and 
;
;
;MODIFICATION HISTORY:
;	Written by:	M.L. Reinke - 8/20/10
;	M.L. Reinke	8/22/10 - modified to properly include the Rocking curve effect and use temperature rather than width
;	M.L. Reinke	8/23/10 - modified to asssume the ti/vel are the same for all lines (for now)
;
;-	

FUNCTION genspec_sphere_obsemiss,lam_b,lam_o,z,em,vel,ti,lhat,flat=flat
	nlines=n(lam_o)+1
	x=size(em)
	ns=x[1]
	c=2.998e8 			;speed of light
	e=1.602e-19			;conversion for eV -> J
	mconv=1.661e-27			;conversion for amu -> kg		

	IF keyword_set(flat) THEN BEGIN
		i=maxloc(em[*,0])
		mass=read_atomic_mass(z[0])
		wid=sqrt((lam_o[0]/c)^2*abs(ti[i])*e*1.0e3/(mass*mconv))
		obsem=fltarr(ns)+em[i,0]/(wid*sqrt(2.0*!pi))*exp(-1.0*(lam_b-lam_o[0])^2/(2.0*wid^2))
	ENDIF ELSE BEGIN
		obsem=fltarr(ns)
		FOR k=0,nlines-1 DO BEGIN
			mass=read_atomic_mass(z[k])
			;wid=sqrt((lam_o[k]/c)^2*abs(ti[*,k])*e*1.0e3/(mass*mconv))	;force Ti to be positive definite so that interpolation errors done lead to NaNs
			;lam_d=lam_o[k]+lam_o[k]/c*(lhat[*,0]*vel[*,0,k]+lhat[*,1]*vel[*,1,k]+lhat[*,2]*vel[*,2,k])
			wid=sqrt((lam_o[k]/c)^2*abs(ti)*e*1.0e3/(mass*mconv))
			lam_d=lam_o[k]+lam_o[k]/c*(lhat[*,0]*vel[*,0]+lhat[*,1]*vel[*,1]+lhat[*,2]*vel[*,2])
			obsem+=em[*,k]/(wid*sqrt(2.0*!pi))*exp(-1.0*(lam_b-lam_d)^2/(2.0*wid^2))
		ENDFOR
	ENDELSE

	RETURN,obsem
END

;+
;NAME:
;	GENSPEC_ROCKING_MAP
;
;PURPOSE:
;	This function calculates the deposition on the detector plane for a finite rocking curve width.  This
;	can then be multiplied by the line-integrated brightness
;
;CALLING SEQUENCE
;	result=GENSPEC_ROCKING_MAP(xyz_m,xyz_b,xyz_det,sig_r,int_r)
;
;INPUTS:
;	xyz_m	FLTARR	[3] of the <x,y,z> vector in the mirror coordinates of the submirror element
;	xyz_b	FLTARR	[3] of the <x,y,z> vector of the central ray that defines the bragg angle
;	xyz_det	FLTARR	[3,ndet] of the <x,y,z> positions of the pixels
;	sig_r	FLT	of the rocking curve width [arcsec]
;	
;OUTPUTS:
;	result	FLTARR	[ndet] of the power deposition fraction over the variety of pixels
;
;OPTIONAL OUTPUTS:
;	dth	FLTARR	[ndet] of the angular deviation of the ray from the primary ray
;
;PROCEDURE:
;	From the three points given (m,b,det) the angular deviation of the ray is found using law of cosines.
;	The falloff is assumed to be gaussian w/r/t the angular deviation and controlled by the rocking
;	curve width, sig_r.  This is normalized so that the total power to the detector is constant as sig_r is varied.
;
;	DO NOT USE IF SIG_R=0.  Check beforehand and just define your own map.
;
;MODIFICATION HISTORY:
;	Written by:	M.L. Reinke 8/20/10
;
;-

FUNCTION genspec_rocking_map,xyz_m,xyz_b,xyz_det,sig_r,dth=dth
	x_bragg=sqrt(total((xyz_b-xyz_m)*(xyz_b-xyz_m),/double))
	x=size(xyz_det)
	dth=dblarr(x[2])
	map=fltarr(x[2])
	FOR i=0,x[2]-1 DO BEGIN
		x_ch=sqrt(total((xyz_det[*,i]-xyz_m)*(xyz_det[*,i]-xyz_m),/double))
		dx=sqrt(total((xyz_det[*,i]-xyz_b)*(xyz_det[*,i]-xyz_b),/double))
		dth[i]=acos((x_bragg^2+x_ch^2-dx^2)/(2.0*x_bragg*x_ch)) 	;	cos(dth) is the dot product between the "bragg vector" and another point on the detector
		map[i]=exp(-dth[i]^2/(2.0*(sig_r*!pi/6.48e5)^2))
	ENDFOR
	map/=total(map)		;this way as sig_r increases the total power to the detector will remain constant
	RETURN,map
END
	
;+
;NAME:
;	GENSPEC_SPHERE_IMAGE
;
;PURPOSE:
;	This function does the raytracing for a spherical reflector, w/ finite size and rocking curve,
;	through the plasma to create synthetic spectral images for testing of analysis and instrumental effects
;
;CALLING SEQUENCE:
;	result=GENSPEC_SPHERE_IMAGE(info,shot,time,lam_o,z,emiss,w,u,ti,rho)
;
;INPUTS:
;	info	STRUC	stucture of an spherical info file (see GENPOS_SPHERICAL_INFO)
;	shot	LONG	shot number to use for EFIT map
;	time	FLT	time point to use in shot for EFIT map
;	lam_o	FLTARR	[nlines] of line centers to include in the spectra
;	z	INTARR	[nlines] of the atomic number of the transitions
;	emiss	FLTARR	[nrho,nlines] of the emissivity profile for each line
;	w	FLTARR	[nrho] of the toroidal rotation frequency [kHz]
;	u	FLTARR	[nrho] of the poloidal rotation flux function (vp=u*Bp) [km/s/T]
;	ti	FLTARR	[nrho] of the ion temperature [keV]
;	rho	FLTARR	[nrho] of the rho values for the above profiles (actually PSINORM)
;
;OPTIONAL INPUTS:
;	n_s	INT	number of points in the plasma to use to calculate the observed emissivity DEFAULT: 100
;	n_y	INT	number of mirror bins in the y-direction DEFAULT: 5
;	n_z	INT	number of mirror bins in the z-direction DEFAULT: 5
;	rzbnd	FLTARR	[3] of the [rmin,rmax,abs(z)] in [m] of the rectangular bin to terminate ray tracing DEFAULT: [0.44,1.0,0.45]
;	r_ap	FLT	of the major radius [m] to translate the POS vector to DEFAULT: not used, ray starts at (R,Z)=(pos[0],pos[1])
;	z_ap	FLT	of the height [m] to translate the POS vector to DEFAULT: not used, ray starts at (R,Z)=(pos[0],pos[1])
;
;KEYWORD PARAMETERS:
;	/verb turns on messages to the terminal
;	/uposload will load the UPOS data from a save file to increase speed for multiple calls
;	
;PROCEDURE:
;	This function operates on similar principles as LINE_BR, defining the emissivity along a line of sight and then integrating along the
;	path to get the brightness.  The observed emissivity profiles is calculated for each pixel/mirror combination using the Bragg angle
;	to select the wavelength.  This can be deposited over multiple pixels using the map from GENSPEC_ROCKING_MAP.
;	
;	A loop is done over the detector points defined by the INFO file and the impact on the image is saved after each pixel is traced over
;	all mirror elements so that the image can be viewed prior to completion due to the long run times for large numbers of
;	pixels.
;
;RESTRICTIONS:
;	This program will save data in your home directory at /idl/hirexsr/genspec/ so that file must be created otherwise things will crash.
;	The username is autodetected using LOGNAME in MLR_FUNCTIONS.
;
;MODIFICATION HISTORY:
;	Written by:	M.L. Reinke - 8/22/10
;
;-

FUNCTION genspec_sphere_image,info,shot,time,lam_o,z,emiss,w,u,ti,rho,n_s=n_s,ny=ny,nz=nz,uposload=uposload,rzbnd=rzbnd,verb=verb,r_ap=r_ap,z_ap=z_ap,flat=flat
	IF NOT keyword_set(rzbnd) THEN rzbnd=[0.44,1.0,0.45]	;set path's [r_min,r_max,abs(z_max)] limit
	IF NOT keyword_set(n_s) THEN n_s=100
	IF NOT keyword_set(ny) THEN ny=5
	IF NOT keyword_set(nz) THEN nz=5
	logname=logname()
	upos_save='/home/'+logname+'/idl/hirexsr/genspec/upos.dat'
	IF NOT keyword_set(uposload) THEN BEGIN
		upos=genpos_spherical2upos(info,ny=ny,nz=nz,du=du,thb=thb,xyz_m=xyz_m,xyz_det=xyz_det)
		save,upos,du,thb,xyz_m,xyz_det,ny,nz,filename=upos_save
	ENDIF ELSE restore, upos_save
	IF keyword_set(verb) THEN print,'UPOS and DU calculated'
	IF total(u) EQ 0 THEN u=-1.0
	int_r=info.m.bragg.iref
	sig_r=info.m.bragg.rwid
	lam_b=info.m.bragg.twod*sin(thb)	;bragg wavelength for pixel/mirror combos

	n_ch=n(upos[0,0,*])+1
	n_m=ny*nz
	n_lines=n(lam_o)+1
	n_s=long(n_s)
	image=fltarr(n_ch)			;create the base image and save to disk after every pixel
	image_save='/home/'+logname+'/idl/hirexsr/genspec/image.dat'

	;load the toroidal field
	mdsopen, 'magnetics', shot
        bt=abs(mdsvalue("\magnetics::BTOR",/quiet,status=status1))
  	bt_time=mdsvalue("dim_of(\magnetics::BTOR)",/quiet,status=status2)
	bt_int=interpol(bt,bt_time,time)	;interpolate the toroidal field onto the tpts grid
	mdsclose,'magnetics',shot
	
	;load the magnetic axis
	mdsopen,'analysis',shot
	raxis=mdsvalue('\ANALYSIS::TOP:EFIT.RESULTS.G_EQDSK:RMAXIS')
	taxis=mdsvalue('dim_of(\ANALYSIS::TOP:EFIT.RESULTS.G_EQDSK:RMAXIS)')
	raxis_int=interpol(raxis,taxis,time)	;interpolate the magnetic axis onto the tpts grid
	mdsclose,'analysis',shot
        
	;intiate (r,z) and s points as vectors for fast EFIT_RZ2RMID and determine r_z points
	ipos_r=fltarr(n_m*n_s)
	ipos_z=fltarr(n_m*n_s)
	ipos_s=fltarr(n_m*n_s)
	FOR i=0,n_ch-1 DO BEGIN
		ipos=upos[*,*,i]
		IF ipos[0] GT 1.0 THEN genpos_pos_reform,ipos,[0.44,1.0,-0.6,0.6]
		iu=du[*,i]
		ilam_b=sin(thb[*,i])*info.m.bragg.twod	;bragg wavelength
		IF keyword_set(verb) THEN IF i MOD 195 EQ 0 THEN print,'calculating power for pixel '+num2str(i,1)+' of '+num2str(n_ch-1,1)
		FOR j=0,n_m-1 DO BEGIN
			;determine where the line of sight terminates in s due to radial boundries
			IF ipos[2,j] LT rzbnd[0] THEN BEGIN	;if r_tang < r_min
				svec=line_s(ipos[*,j],r=rzbnd[0])
				smin_r=svec[1]	 ;choose the smallest s where r=rmin
			ENDIF ELSE BEGIN			;if r_tang > r_min
				svec=line_s(ipos[*,j],r=rzbnd[1])
				smin_r=svec[0]	 ;choose the largest (and non-negative) s where r=rmax
 			ENDELSE

			;determine where the line of sight terminates in s due to z boundries
			IF ipos[3,j] NE 0 THEN BEGIN
				IF line_s(ipos[*,j],z=rzbnd[2]) LT 0 THEN $
					smin_z=line_s(ipos[*,j],z=-1.0*rzbnd[2]) ELSE $
					smin_z=line_s(ipos[*,j],z=rzbnd[2])
			ENDIF ELSE smin_z=2.0*smin_r	;if no inclination, smin_z set > smin_r	
			smin=smin_r < smin_z 	;set minimum s
                	s_ap=0
			IF keyword_set(z_ap) THEN s_ap=line_s(ipos[*,j],z=z_ap) > 0
                	IF keyword_set(r_ap) THEN s_ap=min(line_s(ipos[*,j],r=r_ap)) > 0
			s=make(s_ap,smin,n_s)	;n_s points going from the aperture to the boundary intersecion
			ipos_s[j*n_s:(j+1)*n_s-1]=s
			ipos_r[j*n_s:(j+1)*n_s-1]=line_r(s,ipos[*,j])
			ipos_z[j*n_s:(j+1)*n_s-1]=line_z(s,ipos[*,j])
		ENDFOR
		rhopts=efit_rz2rho(ipos_r,ipos_z,time,shot=shot)
		;rhopts=efit_rz2rmid(ipos_r,ipos_z,time,shot=shot,/rho)		;slower
		
		;interpolate the profiles
		emint=fltarr(n_m*n_s,n_lines)						
		FOR k=0,n_lines-1 DO emint[*,k]=interpol(emiss[*,k],rho,rhopts)	;allow different emissivity for each line
		IF u[0] NE -1 THEN BEGIN
			psi=efit_rz2psi(ipos_r,ipos_z,time,bz,br,shot=shot)
			uint=interpol(u,rho,rhopts)*1.0e3		;convert from km/s/T to m/s/T
		ENDIF ELSE BEGIN
			uint=fltarr(n_m*n_s)				;set u and b-field data to zero if not needed
			bz=fltarr(n_m*n_s)
			br=fltarr(n_m*n_s)
		ENDELSE
		wint=interpol(w,rho,rhopts)*2.0*!pi*1.0e3		;assume rotation profiles (convert from kHz to radians/s) and ti are the same
		tiint=interpol(ti,rho,rhopts)
		FOR j=0,n_m-1 DO BEGIN
			ch=indgen(n_s)+j*n_s
			lhat=[[-cos(ipos[3,j])*(1.0-ipos[2,j]/ipos_r[ch])],[ipos[2,j]*cos(ipos[3,j])/ipos_r[ch]],[fltarr(n_s)-sin(ipos[3,j])]]
			tmp=where(ipos_s[ch] GT 1.0)	;this indicates the line is coming back OUT of the torus
			IF tmp[0] NE -1 THEN lhat[tmp,0]=-1.0*lhat[tmp,0]
			velint=[[uint[ch]*br[ch]],[uint[ch]*bt_int*raxis_int/ipos_r[ch]+wint[ch]*ipos_r[ch]],[uint[ch]*bz[ch]]]	;[R,phi,Z] velocity at each (R,Z) along the line of sight
			obsemiss=genspec_sphere_obsemiss(lam_b[j,i],lam_o,z,emint[ch,*],velint,tiint[ch],lhat,flat=flat)
			lscale=sqrt((ipos[0,j]^2-ipos[2,j]^2)*(1.0+(tan(ipos[3,j]))^2))
			IF sig_r EQ 0 THEN BEGIN
				map=fltarr(n_ch)
				map[i]=1.0
			ENDIF ELSE map=genspec_rocking_map(xyz_m[*,j],xyz_det[*,i],xyz_det,sig_r)
			image+=int_tabulated(ipos_s[ch],obsemiss*lscale)*map*du[j,*]/(4.0*!pi)			;pow=u/4pi*br
		ENDFOR
		save,image,filename=image_save
        ENDFOR
	RETURN,image
END

PRO run_gsi,image,uposload=uposload,verb=verb,ny=ny,nz=nz,info=info
	nrho=100
	rho=make(0,1.1,nrho)
	emiss=(1.05-rho^4)+10.0*(rho)^2*exp(-4.0*rho^2)-0.1
	u=fltarr(nrho)
	w=20.0e3/0.68*(1.05-rho^3)
	ti=2.0*(1.05-rho^2)
	tmp=where(rho GE 1.0)
	w[tmp]=w[tmp[0]]*exp(-(rho[tmp]-rho[tmp[0]])/0.01)
	emiss[tmp]=emiss[tmp[0]]*exp(-(rho[tmp]-rho[tmp[0]])/0.01)
	ti[tmp]=ti[tmp[0]]*exp(-(rho[tmp]-rho[tmp[0]])/0.01)

	;set flat emissivity and ti profiles
	ti=fltarr(nrho)+0.1
	ti[nrho-3:nrho-1]=0.001
	emiss=fltarr(nrho)+1.0
	emiss[nrho-3:nrho-1]=0.0
	w=fltarr(nrho)			;set velocity to zero for now

	infopath='/home/mlreinke/idl/genie/data/info/hirexsr/hirexsr_genspec_test.info'
	IF NOT keyword_set(info) THEN info=genpos_spherical_info(infopath)
	shot=1070830020
	time=1.25
	lam_o=[3.994]
	z=[18]
	emiss_o=emiss
	FOR i=0,n(lam_o)-1 DO emiss=[[emiss],[emiss_o]]
	image=genspec_sphere_image(info,shot,time,lam_o,z,emiss,w,u,ti,rho,uposload=uposload,r_ap=1.0,verb=verb,ny=ny,nz=nz)
END

PRO rr_gsi
	num=10
	image=fltarr(100,num)
	FOR i=1,num DO BEGIN
		print,i
		ny=i
		nz=i
		start_time=systime(/seconds)
		run_gsi,im,ny=nz,nz=nz
		image[*,i-1]=im
		ctime=systime(/seconds)-start_time
		print,ctime
	ENDFOR
	stop
END

PRO gsi_zfocus_convergance,load=load,z0=z0
	IF NOT keyword_set(z0) THEN z0=0.0
	infosave='/home/mlreinke/idl/hirexsr/genspec/gsi_info_f3994.dat'
	imagesave='/home/mlreinke/idl/hirexsr/genspec/gsi_zfocus_conv_z'+num2str(int(z0*100.0),1)+'cm.dat'
	restore,infosave

	;double the dxi resolution around center to handle smaller Ti
	nxi=200
	dxi=max(info.det.xi)-min(info.det.xi)
	xi0=mean(info.det.xi)-(0.07*(z0*100)^2)*1.0e-3
	xi=make(xi0-0.15*dxi,xi0+0.25*dxi,nxi)
	zeta=fltarr(nxi)
	info.m.rot[1]=0.001

	det={x0:info.det.x0,x1:info.det.x1,x2:info.det.x2,xi:xi,zeta:zeta+z0,size:info.det.size,n_xi:nxi,n_zeta:1}
	info={name:info.name,m:info.m,det:det,type:info.type,author:info.author}

	nrho=100
	rho=make(0,1.1,nrho)

	;set flat emissivity and ti profiles
	ti=fltarr(nrho)+0.002
	ti[nrho-3:nrho-1]=0.00001
	emiss=fltarr(nrho)+1.0
	emiss[nrho-3:nrho-1]=0.0
	w=fltarr(nrho)			;set velocity to zero for now
	u=w
	shot=1070830020
	time=1.25
	lam_o=[3.994]
	z=[18]

	grids=[1,2,3,4,6,8,10,12,15,20,25]
	;grids=[1,4]
	num=n(grids)+1
	numpix=info.det.n_xi
	image2d=fltarr(numpix,num)
	image1d=fltarr(numpix,num)
	IF NOT keyword_set(load) THEN BEGIN
		FOR i=0,num-1 DO BEGIN
			print,i
			ny=grids[i]
			nz=grids[i]
			start_time=systime(/seconds)
			image2d[*,i]=genspec_sphere_image(info,shot,time,lam_o,z,emiss,w,u,ti,rho,r_ap=1.0,verb=verb,ny=ny,nz=nz)
			image1d[*,i]=genspec_sphere_image(info,shot,time,lam_o,z,emiss,w,u,ti,rho,r_ap=1.0,verb=verb,ny=ny,nz=1)
			ctime=systime(/seconds)-start_time
			print,ctime
		ENDFOR
		save,image1d,image2d,info,ti,grids,filename=imagesave
	ENDIF ELSE restore, imagesave
	xi=info.det.xi*1.0e3
	image1d*=1.0e10
	image2d*=1.0e10
	ymax=max(image1d[*,0])*1.05
	xr=xi[maxloc(image1d[*,0])]+[-0.2,0.3]

	openwin,0
	plot,xi,image2d[*,0],xtit=n2g('xi')+' [mm]',ytit='Spectral Brightness [AU]',chars=1.2,xr=xr,/xsty,yr=[0,ymax],/ysty
	color=[0,30,70,100,120,150,200]
	pts=[0,2,4,6,8,9,10]
	;pts=[0,1]
	FOR i=0,n(pts) DO BEGIN
		oplot,xi,image2d[*,pts[i]],color=color[i]
		xyouts,5.4,(0.9-0.07*i)*ymax,'n!ly!n=n!lz!n='+num2str(grids[pts[i]],1),color=color[i]
	ENDFOR
	

	openwin,1
	plot,xi,image2d[*,0],xtit=n2g('xi')+' [mm]',ytit='Spectral Brightness [AU]',chars=1.2,xr=xr,/xsty,yr=[0,ymax],/ysty,$
		tit=n2g('lambda')+' = '+num2str(lam_o,dp=3)+'[Ang]   T = '+num2str(ti[0]*1.0e3,dp=1)+' [eV]'
	oplot,xi,image2d[*,num-1]
	oplot,xi,image1d[*,0],linestyle=2.0,color=200
	oplot,xi,image1d[*,num-1],linestyle=2.0,color=200
	out=gaussfit(xi,image2d[*,0],a0,nterms=3)
	oplot,a0[1]*[1,1],[0,ymax],color=30,linestyle=3.0
	db=deriv(xi,image2d[*,num-1])
	chk=where(xi GT a0[1] AND db LT 0)
	xfit=xi[indgen(4)-2+chk[0]]
	yfit=db[indgen(4)-2+chk[0]]
	coefs=poly_fit(xfit,yfit,1)
	peak=-1.0*coefs[0]/coefs[1]
	oplot,peak*[1,1],[0,ymax],color=30,linestyle=3.0
	delta=peak-a0[1]
	jshift=fltarr(num)
	jwidth=fltarr(num)
	FOR i=0,num-1 DO BEGIN
		jshift[i]=int_tabulated(xi,xi*image2d[*,i])/int_tabulated(xi,image2d[*,i])
		jwidth[i]=sqrt(int_tabulated(xi,(xi-jshift[i])^2*image2d[*,i])/int_tabulated(xi,image2d[*,i]))
	ENDFOR
	xyouts,a0[1]+0.1,0.9*ymax,n2g('Delta')+' PEAK '+num2str(delta*1.0e3,dp=1)+' ['+n2g('mu')+'m]',chars=1.5,color=30
	xyouts,a0[1]-0.17,0.9*ymax,'2-D Crystal'
	xyouts,a0[1]-0.17,0.8*ymax,'1-D Crystal',color=200
	johann=cos(asin(3.994/info.m.bragg.twod))*(info.m.size[0])^2/(8.0*info.m.rad)*1.0e3	;johann error in mm
	oplot,a0[1]+[0.1,0.1+johann],0.8*ymax*[1,1],thick=8
	oplot,a0[1]+[0.1,0.1],ymax*[0.76,0.84],thick=8
	oplot,a0[1]+0.1+[johann,johann],ymax*[0.76,0.84],thick=8
	xyouts,a0[1]+0.13,0.745*ymax,'Johann Error',chars=1.3
	openwin,2
	!p.multi=[0,2]
	plot,grids,(jshift-jshift[0])*1.0e3,psym=-7,yr=[-0.0,70],/ysty,xr=[1,26],/xsty,xtit='# Crystal Grids',ytit='Intsr. Shift: l-l11 [microns]'
	plot,grids,jwidth*1.0e3,psym=-2,yr=[0,max(jwidth)*1.05e3],/ysty,xr=[1,26],/xsty,xtit='# Crystal Grids',ytit='Instr. Width [microns]'
	!p.multi=0

	;stop 
END

PRO gsi_wfocus_convergance,load=load
	infosave='/home/mlreinke/idl/hirexsr/genspec/gsi_info_f3944.dat'
	imagesave='/home/mlreinke/idl/hirexsr/genspec/gsi_wfocus_conv.dat'
	restore,infosave

	nrho=100
	rho=make(0,1.1,nrho)

	;set flat emissivity and ti profiles
	ti=fltarr(nrho)+0.1
	ti[nrho-3:nrho-1]=0.001
	emiss=fltarr(nrho)+1.0
	emiss[nrho-3:nrho-1]=0.0
	w=fltarr(nrho)			;set velocity to zero for now
	u=w
	shot=1070830020
	time=1.25
	lam_o=[3.944]
	z=[18]

	num=10
	numpix=info.det.n_xi
	image2d=fltarr(numpix,num)
	image1d=fltarr(numpix,num)
	IF NOT keyword_set(load) THEN BEGIN
		FOR i=1,num DO BEGIN
			print,i
			ny=i
			nz=i
			start_time=systime(/seconds)
			image2d[*,i-1]=genspec_sphere_image(info,shot,time,lam_o,z,emiss,w,u,ti,rho,r_ap=1.0,verb=verb,ny=ny,nz=nz)
			image1d[*,i-1]=genspec_sphere_image(info,shot,time,lam_o,z,emiss,w,u,ti,rho,r_ap=1.0,verb=verb,ny=ny,nz=1)
			ctime=systime(/seconds)-start_time
			print,ctime
		ENDFOR
		save,image1d,image2d,filename=imagesave
	ENDIF ELSE restore, imagesave
	xi=info.det.xi*1.0e3
	image1d*=1.0e10
	image2d*=1.0e10

	openwin,0
	plot,xi,image2d[*,0],xtit=n2g('xi')+' [mm]',ytit='Spectral Brightness [AU]',chars=1.2,xr=[27.0,28.4],/xsty,yr=[0,2.25],/ysty
	color=[30,70,100,120,150,200]
	pts=[1,2,3,5,7,9]
	FOR i=0,n(pts) DO BEGIN
		oplot,xi,image2d[*,pts[i]],color=color[i]
		xyouts,5.6,2.0-0.2*i,num2str(pts[i],1),color=color[i]
	ENDFOR

	openwin,1
	plot,xi,image2d[*,0],xtit=n2g('xi')+' [mm]',ytit='Spectral Brightness [AU]',chars=1.2,/xsty,xr=[27.0,28.4],yr=[0,2.25],/ysty,$
		tit=n2g('lambda')+' = '+num2str(lam_o,dp=3)+'[Ang]   T = '+num2str(ti[0],dp=1)+' [keV]'
	oplot,xi,image2d[*,num-1]
	oplot,xi,image1d[*,0],linestyle=2.0,color=200
	oplot,xi,image1d[*,num-1],linestyle=2.0,color=200
	oplot,xi[maxloc(image2d[*,0])]*[1,1],[0,3],color=30,linestyle=3.0
	oplot,xi[maxloc(image2d[*,num-1])]*[1,1],[0,3],color=30,linestyle=3.0
	delta=xi[maxloc(image2d[*,num-1])]-xi[maxloc(image2d[*,0])]
	xyouts,27.9,2.0,n2g('Delta')+' PEAK '+num2str(delta*1.0e3,dp=1)+' ['+n2g('mu')+'m]',chars=1.5,color=30


	stop
END

PRO gsi_oozfocus_convergance,load=load
	infosave='/home/mlreinke/idl/hirexsr/genspec/gsi_info_f3994.dat'
	imagesave='/home/mlreinke/idl/hirexsr/genspec/gsi_oozfocus_conv.dat'
	restore,infosave

	nrho=100
	rho=make(0,1.1,nrho)

	;set flat emissivity and ti profiles
	ti=fltarr(nrho)+0.1
	ti[nrho-3:nrho-1]=0.001
	emiss=fltarr(nrho)+1.0
	emiss[nrho-3:nrho-1]=0.0
	w=fltarr(nrho)			;set velocity to zero for now
	u=w
	shot=1070830020
	time=1.25
	lam_o=[3.944]
	z=[18]
	info.det.xi+=28.3e-3

	num=10
	numpix=info.det.n_xi
	image2d=fltarr(numpix,num)
	image1d=fltarr(numpix,num)
	IF NOT keyword_set(load) THEN BEGIN
		FOR i=1,num DO BEGIN
			print,i
			ny=i
			nz=i
			start_time=systime(/seconds)
			image2d[*,i-1]=genspec_sphere_image(info,shot,time,lam_o,z,emiss,w,u,ti,rho,r_ap=1.0,verb=verb,ny=ny,nz=nz)
			image1d[*,i-1]=genspec_sphere_image(info,shot,time,lam_o,z,emiss,w,u,ti,rho,r_ap=1.0,verb=verb,ny=ny,nz=1)
			ctime=systime(/seconds)-start_time
			print,ctime
		ENDFOR
		save,image1d,image2d,filename=imagesave
	ENDIF ELSE restore, imagesave
	xi=info.det.xi*1.0e3
	image1d*=1.0e10
	image2d*=1.0e10		

	openwin,0
	plot,xi,image2d[*,0],xtit=n2g('xi')+' [mm]',ytit='Spectral Brightness [AU]',chars=1.2,/xsty,yr=[0,2.25],/ysty
	color=[30,70,100,120,150,200]
	pts=[1,2,3,5,7,9]
	FOR i=0,n(pts) DO BEGIN
		oplot,xi,image2d[*,pts[i]],color=color[i]
		xyouts,33.8,2.0-0.2*i,'n!ly!n=n!lz!n='+num2str(pts[i]+1,1),color=color[i]
	ENDFOR
	

	openwin,1
	plot,xi,image2d[*,0],xtit=n2g('xi')+' [mm]',ytit='Spectral Brightness [AU]',chars=1.2,/xsty,yr=[0,2.25],/ysty,$
		tit=n2g('lambda')+' = '+num2str(lam_o,dp=3)+'[Ang]   T = '+num2str(ti[0],dp=1)+' [keV]'
	oplot,xi,image2d[*,num-1]
	oplot,xi,image1d[*,0],linestyle=2.0,color=200
	oplot,xi,image1d[*,num-1],linestyle=2.0,color=200
	out=gaussfit(xi,image2d[*,0],a0,nterms=3)
	oplot,a0[1]*[1,1],[0,3],color=30,linestyle=3.0
	out=gaussfit(xi,image2d[*,num-1],a2,nterms=3)
	oplot,a2[1]*[1,1],[0,3],color=30,linestyle=3.0
	delta=a2[1]-a0[1]
	xyouts,5.5,2.0,n2g('Delta')+' PEAK '+num2str(delta*1.0e3,dp=1)+' ['+n2g('mu')+'m]',chars=1.5,color=30
	xyouts,33.8,1.5,'2-D Crystal'
	xyouts,33.8,1.25,'1-D Crystal',color=200
	stop
END


PRO gsi_oowfocus_convergance,load=load
	infosave='/home/mlreinke/idl/hirexsr/genspec/gsi_info_f3944.dat'
	imagesave='/home/mlreinke/idl/hirexsr/genspec/gsi_oowfocus_conv.dat'
	restore,infosave

	nrho=100
	rho=make(0,1.1,nrho)

	;set flat emissivity and ti profiles
	ti=fltarr(nrho)+0.1
	ti[nrho-3:nrho-1]=0.001
	emiss=fltarr(nrho)+1.0
	emiss[nrho-3:nrho-1]=0.0
	w=fltarr(nrho)			;set velocity to zero for now
	u=w
	shot=1070830020
	time=1.25
	lam_o=[3.994]
	z=[18]
	info.det.xi-=27.5e-3

	num=10
	numpix=info.det.n_xi
	image2d=fltarr(numpix,num)
	image1d=fltarr(numpix,num)
	IF NOT keyword_set(load) THEN BEGIN
		FOR i=1,num DO BEGIN
			print,i
			ny=i
			nz=i
			start_time=systime(/seconds)

			image2d[*,i-1]=genspec_sphere_image(info,shot,time,lam_o,z,emiss,w,u,ti,rho,r_ap=1.0,verb=verb,ny=ny,nz=nz)
			image1d[*,i-1]=genspec_sphere_image(info,shot,time,lam_o,z,emiss,w,u,ti,rho,r_ap=1.0,verb=verb,ny=ny,nz=1)
			ctime=systime(/seconds)-start_time
			print,ctime
		ENDFOR
		save,image1d,image2d,filename=imagesave
	ENDIF ELSE restore, imagesave
	xi=info.det.xi*1.0e3
	image1d*=1.0e10
	image2d*=1.0e10		

	openwin,0
	plot,xi,image2d[*,0],xtit=n2g('xi')+' [mm]',ytit='Spectral Brightness [AU]',chars=1.2,/xsty,yr=[0,2.25],/ysty
	color=[30,70,100,120,150,200]
	pts=[1,2,3,5,7,9]
	FOR i=0,n(pts) DO BEGIN
		oplot,xi,image2d[*,pts[i]],color=color[i]
		xyouts,5.6,2.0-0.2*i,num2str(pts[i],1),color=color[i]
	ENDFOR
	

	openwin,1
	plot,xi,image2d[*,0],xtit=n2g('xi')+' [mm]',ytit='Spectral Brightness [AU]',chars=1.2,/xsty,yr=[0,2.25],/ysty,$
		tit=n2g('lambda')+' = '+num2str(lam_o,dp=3)+'[Ang]   T = '+num2str(ti[0],dp=1)+' [keV]'
	oplot,xi,image2d[*,num-1]
	oplot,xi,image1d[*,0],linestyle=2.0,color=200
	oplot,xi,image1d[*,num-1],linestyle=2.0,color=200
	oplot,xi[maxloc(image2d[*,0])]*[1,1],[0,3],color=30,linestyle=3.0
	oplot,xi[maxloc(image2d[*,num-1])]*[1,1],[0,3],color=30,linestyle=3.0
	delta=xi[maxloc(image2d[*,num-1])]-xi[maxloc(image2d[*,0])]
	xyouts,5.5,2.0,n2g('Delta')+' PEAK '+num2str(delta*1.0e3,dp=1)+' ['+n2g('mu')+'m]',chars=1.5,color=30
	stop
END



PRO gsi_zfocus_gamma,load=load
	infosave='/home/mlreinke/idl/hirexsr/genspec/gsi_info_f3994.dat'
	imagesave='/home/mlreinke/idl/hirexsr/genspec/gsi_zfocus_gamma.dat'
	restore,infosave

	nrho=100
	rho=make(0,1.1,nrho)

	;set flat emissivity and ti profiles
	ti=fltarr(nrho)+0.1
	ti[nrho-3:nrho-1]=0.001
	emiss=fltarr(nrho)+1.0
	emiss[nrho-3:nrho-1]=0.0
	w=fltarr(nrho)			;set velocity to zero for now
	u=w
	shot=1070830020
	time=1.25
	lam_o=[3.994]
	z=[18]

	gamma=[2.03,2.05,2.07,2.09]
	num=n(gamma)+1
	numpix=info.det.n_xi
	image0d=fltarr(numpix,num)
	image1d=fltarr(numpix,num)
	IF NOT keyword_set(load) THEN BEGIN
		FOR i=0,num-1 DO BEGIN
			print,i
			info.m.rot[2]=gamma[i]
			start_time=systime(/seconds)
			image0d[*,i]=genspec_sphere_image(info,shot,time,lam_o,z,emiss,w,u,ti,rho,r_ap=1.0,verb=verb,ny=1,nz=1,/flat)
			image1d[*,i]=genspec_sphere_image(info,shot,time,lam_o,z,emiss,w,u,ti,rho,r_ap=1.0,verb=verb,ny=50,nz=1,/flat)
			stop
			ctime=systime(/seconds)-start_time
			print,ctime
		ENDFOR
		save,image0d,image1d,filename=imagesave
	ENDIF ELSE restore, imagesave
	xi=info.det.xi*1.0e3
	image0d*=1.0e10
	image1d*=1.0e10

	openwin,0
	


	stop
END

PRO gsi_zfocus_temp,load=load,radial=radial,z0=z0
	IF NOT keyword_set(z0) THEN z0=0.0
	infosave='/home/mlreinke/idl/hirexsr/genspec/gsi_info_f3994.dat'
	imagesave='/home/mlreinke/idl/hirexsr/genspec/gsi_zfocus_temp_z'+num2str(int(z0*100.0),1)+'cm.dat'
	restore,infosave
	info.m.rot[1]=0.001

	;temp=[0.1,0.2,0.5,0.7,1.0,2.0,3.0]
	temp=[0.1,0.8,1.8,3.5]
	;adjust to nominally radial
	;IF keyword_set(radial) THEN BEGIN
	;	info.m.rot[2]=2.075
	;	imagesave='/home/mlreinke/idl/hirexsr/genspec/gsi_zfocus_temp_radial.dat'
	;ENDIF ELSE imagesave='/home/mlreinke/idl/hirexsr/genspec/gsi_zfocus_temp.dat'

	;double the dxi range around center to handle larger Ti
	nxi=250
	dxi=max(info.det.xi)-min(info.det.xi)
	xi0=mean(info.det.xi)-(0.07*(z0*100)^2)*1.0e-3
	xi=make(xi0-1.2*dxi,xi0+1.2*dxi,nxi)
	zeta=fltarr(nxi)
	info.m.rot[1]=0.001

	det={x0:info.det.x0,x1:info.det.x1,x2:info.det.x2,xi:xi,zeta:zeta+z0,size:info.det.size,n_xi:nxi,n_zeta:1}
	info={name:info.name,m:info.m,det:det,type:info.type,author:info.author}

	nrho=100
	rho=make(0,1.1,nrho)

	;set flat emissivity and ti profiles
	emiss=fltarr(nrho)+1.0
	emiss[nrho-3:nrho-1]=0.0
	w=fltarr(nrho)			;set velocity to zero for now
	u=w
	shot=1070830020
	time=1.25
	lam_o=[3.994]
	z=[18]

	num=n(temp)+1
	numpix=info.det.n_xi
	image2d=fltarr(numpix,num)
	image0d=fltarr(numpix,num)
	IF NOT keyword_set(load) THEN BEGIN
		FOR i=0,num-1 DO BEGIN
			print,i
			start_time=systime(/seconds)
			ti=fltarr(nrho)+temp[i]
			ti[nrho-3:nrho-1]=0.001
			image2d[*,i]=genspec_sphere_image(info,shot,time,lam_o,z,emiss,w,u,ti,rho,r_ap=1.0,verb=verb,ny=3,nz=3)
			image0d[*,i]=genspec_sphere_image(info,shot,time,lam_o,z,emiss,w,u,ti,rho,r_ap=1.0,verb=verb,ny=1,nz=1)
			ctime=systime(/seconds)-start_time
			print,ctime
		ENDFOR
		save,image0d,image2d,info,ti,filename=imagesave
	ENDIF ELSE restore, imagesave
	xi=info.det.xi*1.0e3
	image0d*=1.0e10
	image2d*=1.0e10
	ymax=max(image0d[*,0])*1.05
	xr=xi[maxloc(image0d[*,0])]+[-2.0,2.0]

	openwin,0
	plot,xi,image0d[*,0],xtit=n2g('xi')+' [mm]',ytit='Spectral Brightness [AU]',chars=1.2,xr=xr,/xsty,yr=yr,/ysty
	color=[0,30,100,200]
	jshift=fltarr(num)
	jwidth=fltarr(num)
	;pts=[0,2,4,6]
	pts=[0,1,2,3]
	FOR i=0,n(temp) DO BEGIN
		out=gaussfit(xi,image0d[*,i],a0,nterms=3)
		out=gaussfit(xi,image2d[*,i],a2,nterms=3)
		jshift[i]=a2[1]-a0[1]
		jwidth[i]=a2[2]/a0[2]
	ENDFOR

	FOR k=0,n(pts) DO BEGIN
		i=pts[k]
		oplot,xi,image0d[*,i],color=color[k]
		oplot,xi,image2d[*,i],color=color[k],linestyle=3.0
		IF k EQ 0 THEN oplot,a0[1]*[1,1],[0,3],linestyle=1,color=color[k]
		oplot,jshift[i]*[1,1],[0,3],linestyle=1,color=color[k]
		;xyouts,5.6,2.0-0.2*i,'n!ly!n=n!lz!n='+num2str(pts[i]+1,1),color=color[k]
	ENDFOR
	
	stop

END