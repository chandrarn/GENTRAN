;+
;NAME:
;	PHA_POS
;
;PURPOSE:
;	Loads line of sight information for the PILATUS-based PHA installed at H-Port
;	(will be stored in tree eventually)
;
;CALLING SEQUENCE:
;	output=PHA_POS()
;
;OPTIONAL INPUTS:
;	beta	FLT	of the declination of the system DEFAULT (in .info file)
;
;KEYWORD PARAMETERS:
;	load	/load will load the data from an IDL saveset at /usr/local/cmod/codes/spectroscopy/pha/
;	save	/save will save the data to that saveset
;
;OUTPUTS:
;	output	FLTARR [4,487*195] of the POS vectors for each pixel on the PILATUS
;	
;OPTIONAL OUTPUTS:
;	u	FLTARR [487*195] of the etendue values [m^2-str] of each pixel
;	info	STRUC of the info file that describes the setup
;
;RESTRICTIONS:
;	This function loads the info data from:
;	/home/mlreinke/GENIE/usr/mlreinke/pilatus_pha.info
;
;MODIFICATION HISTORY:
;	Written by:	M.L. Reinke - 10/16/15
;
;-

FUNCTION pha_pos,u=u,info=info,beta=beta,load=load,save=save
	savepath='/usr/local/cmod/codes/spectroscopy/pha/pha_pos.sav'
	IF NOT keyword_set(load) THEN BEGIN
		path='/home/mlreinke/GENIE/usr/mlreinke/pilatus_pha.info'
		info=genpos_planar_info(path)
		IF keyword_set(beta) THEN info.ap.rot[1]=beta*!pi/180.0		;beta in degrees
		pos=genpos_planar2pos(info,etendue=u)
		IF keyword_set(save) THEN save,pos,u,info,filename=savepath
	ENDIF ELSE restore, savepath
	RETURN, pos
END

;+
;NAME:
;	INVERT_PHA
; 
;PURPOSE:
;	This procedure uses GENPOS to create spatial weighting matrices from given POS vector arrays and then invert given brightness profiles 
;
;CALLING SEQUENCE
;	INVERT_PHA,shot,time,br,pos,emiss,rho,brchk
;
;INPUTS:
;	shot	LONG 	shot number
;	time	FLTARR	[ntime] of time time points to do inversion
;	br	FLTARR	[nch,ntime] of the 'brightness' signal to be inverted
;	pos	FLTARR	[4,nch] of the line of sight parameterization
;
;OPTIONAL INPUTS:
;	brerr	FLTARR	[nch,ntime] of the uncertainties in the br array values	DEFAULT: 0.0
;	eps	FLT	of the regularization weighting DEFAULT: 1.0 (norm. to longest length element)
;	nrho	INT	[nrho] of the number of channels in the inversion DEFAULT: 1.25*nch
;
;KEYWORD PARAMETERS:
;	psinorm	/psinorm will make the radial grid be normalized psi instead of LFS midplane normalized minor radius.	
;
;OUTPUTS:
;	emiss	FLTARR	[nrho,ntime] of the 'emissivities'
;	rho	FLTARR	[nrho] of the spatial locations in r/a (unless /psinorm used)
;	brchk	FLTARR	[nch,ntime] of the 'brightness' values consistent with emissivity for checking against input
;	
;OPTIONAL OUTPUTS:
;	emerr	FLTARR	[nrho,ntime] of the uncertainty in the emissivity
;
;RESTRICTIONS:
;	This procedure uses GENPOS_POS2VOXEL_MATRIX for spatial weighting and GENPOS_PROFILE_INVERT for the inversion
;
;MODIFICATION HISTORY:
;	Written by: 	M.L. Reinek - 10/19/15
;
;-

PRO invert_pha,shot,time,br,pos,emiss,rho,brchk,eps=eps,nrho=nrho,brerr=brerr,emerr=emerr,psinorm=psinorm
	IF NOT keyword_set(eps) THEN eps=1.0
	nch=n(pos[0,*]+1)
	ntime=n(time)+1
	IF NOT keyword_set(nrho) THEN nrho=nch*1.25
	IF NOT keyword_set(brerr) THEN brerr=br*0.0		;set zero error if none set
	emiss=fltarr(nrho,ntime)
	emerr=emerr(nrho,ntime)
	brchk=fltarr(nch,ntimee)

	vox=genpos_pos2voxel_matrix(pos,(fltarr(nch)+1.0)*(4.0*!pi),shot,tpts=time,n_rho=nrho,rho_vec=rho,psinorm=psinorm)
	FOR i=0,n(time) DO BEGIN
		emiss[*,i]=genpos_profile_invert(br[*,i],vox[*,*,i],eps*max(vox[*,*,i]),brchk=ibrchk,err=brerr[*,i])
		brchk[*,i]=ibrchk.mom
		emerr[*,i]=ibrchk.inverr
	ENDFOR
END


FUNCTION example_pha_pos,etendue=etendue,xch=xch
	ipos=pha_pos(u=iu,/load)
	IF NOT keyword_set(xch) THEN xch=13
	pixels=findgen(487/xch)*195*xch+195/2	;every x pixels along the center
	pos=ipos[*,pixels]
	genpos_pos_reform,pos,[0.4,1.0,-0.8,0.8]
	tmp=where(pos[1,*] LT 0.3)	;rough filter to remove upper channels
	etendue=iu[pixels]
	pos=pos[*,tmp]
	etendue=etendue[tmp]
	RETURN,pos
END

PRO show_genpos_pha,shot=shot,time=time,eps=eps,sigma=sigma,flat=flat,xch=xch
	IF NOT keyword_set(shot) THEN shot=1140815012
	IF NOT keyword_set(time) THEN time=1.0
	IF NOT keyword_set(eps) THEN eps=1.0
	IF NOT keyword_set(sigma) THEN sigma=0.0
	
	efit_time=line_gettimes(shot)
	efit_rmid=line_getrmid(shot)
	efit_axis=line_getaxis(shot)
	index=ipt(efit_time,time)
        irmid=reform(efit_rmid[index,*])
        rho=(irmid-irmid[0])/(last(irmid)-irmid[0])
	IF keyword_set(flat) THEN emiss=rho*0.0+1.0 ELSE emiss=(1-rho^2)+3*exp(-(rho-0.7)^2/(2*0.2^2))*(1-rho^1.5)
	pos=example_pha_pos(etendue=u,xch=xch)
	nch=n(pos[0,*])+1
	chnum=indgen(nch)+1
	line_path_plots,pos,shot=shot,tpt=time,rzbnd=rzbnd
	br=genpos_line_br(pos,emiss,irmid,time,shot,time)
	err=randomn(seed,nch)*sigma*br
	pos[3,*]+=0.08
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
