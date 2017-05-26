; following GENTRAN/gentran_z/gentran_moly.pro

FUNCTION pos_loweus,shot
; Returns the pos vector (see Reinke Thesis appendix C) for the McPherson spectrometer
; as deduced from CAD images in ~mlreinke/xeus/xeus_loweus_POS_*.tif
	RETURN, [2.561, -0.215, 0.196,-0.1136] ; [R, Z, major radius of closest approach, azimuth]
END

FUNCTION pos_xeus,shot
; Returns the pos vector (see Reinke Thesis appendix C) for the McPherson spectrometer
; as deduced from CAD images in ~mlreinke/xeus/xeus_loweus_POS_*.tif
	RETURN, [2.561, 0.215, 0.196, 0.1136] ; [R, Z, major radius of closest approach, azimuth]
END

PRO load_loweus_line,shot,br,verb=verb
	mdsopen,'spectroscopy',shot
	  br0={d:mdsvalue('\SPECTROSCOPY::TOP.IMPSPEC.Mo.LINE0:BR'), $
	       t:mdsvalue('dim_of(\SPECTROSCOPY::TOP.IMPSPEC.Mo.LINE0:BR,0)'), $
	       sig:mdsvalue('dim_of(\SPECTROSCOPY::TOP.IMPSPEC.Mo.LINE0:BR,1)'), $
	       lab:(mdsvalue('\SPECTROSCOPY::TOP.IMPSPEC.Mo.LINE0:LABEL'))[0]+' '+num2str((mdsvalue('\SPECTROSCOPY::TOP.IMPSPEC.Mo.LINE0:LAM'))[0]), $
	       lam:(mdsvalue('\SPECTROSCOPY::TOP.IMPSPEC.Mo.LINE0:LAM'))[0] $
	      }
	  br1={d:mdsvalue('\SPECTROSCOPY::TOP.IMPSPEC.Mo.LINE1:BR'), $
	       t:mdsvalue('dim_of(\SPECTROSCOPY::TOP.IMPSPEC.Mo.LINE1:BR)'), $
	       lab:(mdsvalue('\SPECTROSCOPY::TOP.IMPSPEC.Mo.LINE1:LABEL'))[0]+' '+num2str((mdsvalue('\SPECTROSCOPY::TOP.IMPSPEC.Mo.LINE1:LAM'))[0]), $
	       lam:(mdsvalue('\SPECTROSCOPY::TOP.IMPSPEC.Mo.LINE1:LAM'))[0] $
	      }
	mdsclose,'spectroscopy',shot

	br=ptrarr(2,/allocate_heap)
	*br[0]=br0
	*br[1]=br1
END

PRO moly_density_loweus,shot,dff=dff,plot=plot,xtomo=xtomo
	; load time histories of lines from loweus
	load_loweus_line,shot,br
	 br_MoXXXI_loweus=*br[0]
	 br_MoXXXII_loweus=*br[1]
	 t=br_MoXXXI_loweus.t
	 dt=0.05 ; 50ms lowest timescale on which steady state transport makes sense

	; load profiles and midplane radius
	qfit_ng,shot,dens,temp,rmaj,time,data=profdata
	mdsopen,'analysis',shot
	 rmid=mdsvalue('\efit_rmid')
	 efit_t=mdsvalue('dim_of(\efit_rmid)')
	mdsclose,'analysis',shot
	profdata=create_struct('dens',dens,'temp',temp,'rmaj',rmaj,'time',time,'rmid',rmid,'efit_t',efit_t, profdata)

	; reject times where efits or thompson profiles don't exist
	tmp=where(t ge efit_t[0] and t le max(efit_t)  and  t ge time[0] and t le max(time))
	if tmp[0] ne -1 then t=t[tmp]

	t=make(min(t),max(t),2*ceil((max(t)-min(t))/dt))
	tmp=where(t ge t[0]+dt)
	if tmp[0] ne -1 then begin
		t1=t[tmp-tmp[0]]
		t2=t[tmp]
	endif else stop, 'not enough time slices for steady state time...'

	IF n(t1) NE n(t2) THEN STOP, 't1 and t2 must be the same size!'
	ntime=n(t1)+1

	IF NOT keyword_set(dff) THEN dff=0.5

	charge=findgen(42)
	nlam=3000
	lam0=100 ;[A]
	lam1=290 ;[A]
	dlam=0.27 ; line width
	lam=make(lam0,lam1,nlam)

	;setup view
	pos_loweus=pos_loweus(shot)
	genpos_pos_reform,pos_loweus,[0.44,1.0,0.6,-0.6]

	; empirical calibration factors, see Reinke APS 2012 poster	
	eta_MoXXXII=1 ;0.65e19
	eta_MoXXXI=1 ;1.05e19

	print,'calculating densities for shot '+num2str(shot,1)
	FOR j=0,ntime-1 DO BEGIN
		print, ' j='+num2str(j,1)+' '+num2str(t1[j],dp=3)+' < t < '+num2str(t2[j],dp=3)
		efit=0
		data=0

		; calculates the abundance fraction profile for various charge states of molybdenum
		;gentran_profiles,shot,42,t1[j],t2[j],/pt,dff=dff,data=data,profdata=profdata,plot=plot,/qfit,/noly
		gentran_profiles,shot,42,t1[j],t2[j],/pt,dff=dff,data=data,profdata=profdata,plot=plot,/qfit,/noly,/fq ; for profiles determined only by coronal equilibrium
		;gentran_profiles,1120927010,42,0.6,0.7,/pt,data=data,/plot,/qfit,/noly
		;gentran_profiles,1120613007,42,0.65,0.69,/pt,data=data,/plot,/qfit,/noly

		; loads MIST/HULLAC data for excitation rates and calculates emissivity
		sxr_vuv_molybdenum,data.csden,data.temp,data.dens,wave,chg,iemiss,molyatomic_dat=molyatomic_dat
		;iemiss=fltarr(size(iemiss,/dim))+1.0 ; to calculate relative brightness due only to transport changes
		;print, 'emissivities calculated'

		nlines=n(wave)+1
		nrad=n(data.rho)+1
		IF j EQ 0 THEN BEGIN
			emiss=fltarr(nlines,nrad,ntime)
			bright=fltarr(nlines,ntime)
			spec=fltarr(nlam,ntime)
			zeff=fltarr(nrad,ntime)
		ENDIF
		emiss[*,*,j]=iemiss

		; calculate line-integrated brightness
		FOR i=0,nlines-1 DO BEGIN
			bright[i,j]=genpos_line_br(pos_loweus,reform(emiss[i,*,j]),data.rmaj,[data.time],data.shot,data.time,efit=efit,verb=verb)
			;bright[i,j]=line_br(pos_mcp,reform(emiss[i,*,j]),data.rmaj,[data.time],data.shot,data.time,plots=plot,efit=efit,rmid=rmid,verb=verb)
		ENDFOR

		; calculate line-integrated spectra
		FOR i=0,nlines-1 DO spec[*,j] += bright[i,j]/(dlam*sqrt(2.0*!pi))*exp(-(lam-wave[i])^2/(2*dlam^2))
		FOR i=0,nrad-1 DO zeff[i,j] = total(data.csden[*,i]*charge^2/data.dens[i])

		IF keyword_set(xtomo) THEN xtomo_genrad_profiles,shot,t1[j],t2[j],$
			data.csden,data.cserr,data.temp,data.terr,data.dens,data.derr,data.rho,$
			plotwin=plotwin,t=t,nz=nz_MoXXXI[0],zeff=zeff_bck,out=outx
	ENDFOR

	l=br_MoXXXI_loweus.lam[0]
	tmp=min(abs(l-wave),i) & i_MoXXXI=i
	 br_MoXXXI_gentran=reform(bright[i,*])  ;115.999
	 print,'Loweus Be-like line at '+num2str(l)+'A closest to HULLAC line at '+num2str(wave[i])+'A'

	l=br_MoXXXII_loweus.lam[0]
	tmp=min(abs(l-wave),i) & i_MoXXXII=i
	 br_MoXXXII_gentran=reform(bright[i,*]) ;127.873
	 print,'Loweus Li-like line at '+num2str(l)+'A closest to HULLAC line at '+num2str(wave[i])+'A'

	nz_MoXXXII=fltarr(ntime)
	nz_MoXXXI=fltarr(ntime)
	FOR j=0,ntime-1 DO BEGIN
		tmp=where(br_MoXXXII_loweus.t ge t1[j] $
		      and br_MoXXXII_loweus.t le t2[j] ) 
		IF tmp[0] NE -1 $
		 THEN ave_br=mean(br_MoXXXII_loweus.d[tmp]) $
		 ELSE ave_br=br_MoXXXII_loweus.d[ipt((t1[j]+t2[j])/2,br_MoXXXII_loweus.t)]
		nz_MoXXXII[j]=eta_MoXXXII*ave_br/br_MoXXXII_gentran[j]
		IF tmp[0] NE -1 $
		 THEN ave_br=mean(br_MoXXXI_loweus.d[tmp]) $
		 ELSE ave_br=br_MoXXXI_loweus.d[ipt((t1[j]+t2[j])/2,br_MoXXXII_loweus.t)]
		nz_MoXXXI[j]=eta_MoXXXI*ave_br/br_MoXXXI_gentran[j]
	ENDFOR

	d={t1:t1, t2:t2, profdata:profdata, br_MoXXXI_loweus:br_MoXXXI_loweus, br_MoXXXII_loweus:br_MoXXXII_loweus, $
	   bright:bright, spec:spec, zeff:zeff, br_MoXXXI_gentran:br_MoXXXI_gentran, br_MoXXXII_gentran:br_MoXXXII_gentran, $
	   nz_MoXXXII:nz_MoXXXII, nz_MoXXXI:nz_MoXXXI, lam:lam, rho:data.rho, wave:wave, i_MoXXXI:i_MoXXXI, i_MoXXXII:i_MoXXXII}

	fnam=num2str(shot,1)+'-moly_density_loweus.dat'
	print,'saving results to '+fnam
	save,d,filename=fnam
	stop
END
