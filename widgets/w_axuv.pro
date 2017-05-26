;load diode emissivity
PRO load_emiss_data,shot,array,array_status,em,r,t,brchk,chkrad,rmid,z,em_err,rm_err,err_stat=err_stat,tree=tree
	mdsopen,'spectroscopy',shot
	em=[0]
	r=[0]
	t=[0]
	brchk=[0]
	chkrad=[0]
	rmid=[0]
	z=[0]
	tmp=mdsvalue('\SPECTROSCOPY::TOP.BOLOMETER.RESULTS.DIODE.'+array+':EPS',/quiet,status=status)
	IF status AND array_status THEN BEGIN
		path='\SPECTROSCOPY::TOP.BOLOMETER.RESULTS.DIODE.'+array+':EMISS'
		em=mdsvalue('_sig='+path)*1.0e-6
		r=mdsvalue('dim_of(_sig,0)')
		t=mdsvalue('dim_of(_sig,1)')
		rmid=mdsvalue('dim_of(_sig,2)',/quiet,status=rm_stat)
		em_err=mdsvalue('dim_of(_sig,3)',/quiet,status=err_stat)
		rmidp=mdsvalue('dim_of(_sig,4)',/quiet,status=err_stat)
		rmidm=mdsvalue('dim_of(_sig,5)',/quiet,status=err_stat)
		tree=mdsvalue('dim_of(_sig,6)',/quiet,status=tree_stat)
		IF n(rmidp) EQ 0 THEN err_stat=0
		IF err_stat THEN BEGIN 
			em_err*=1.0e-6
			rm_err=0.5*(rmidp-rmidm)
		ENDIF
		IF NOT tree_stat THEN tree='ANALYSIS'
		path='\SPECTROSCOPY::TOP.BOLOMETER.RESULTS.DIODE.'+array+':BRCHK'
		brchk=mdsvalue('_sig='+path,/quiet,status=br_stat)
		chkrad=mdsvalue('dim_of(_sig,0)',/quiet)
		z=mdsvalue('\SPECTROSCOPY::TOP.BOLOMETER.DIODE_CALIB.'+array+':Z_O',/quiet,status=zstat)
		IF NOT zstat THEN z=mdsvalue('\SPECTROSCOPY::TOP.BOLOMETER.DIODE_CALIB.'+array+':Z_0',/quiet)
		IF NOT br_stat THEN BEGIN
			brchk=[0]
			chkrad=[0]			
		ENDIF ELSE brchk*=1.0e-6
		IF NOT rm_stat THEN rmid=[0]
	ENDIF
	mdsclose,'spectroscopy',shot
END

;load diode brightness
PRO load_bright_data,shot,array,array_status,br,rt,t,good,br_err,rt_err,err_stat=err_stat
	mdsopen,'spectroscopy',shot
	br=[0]
	rt=[0]
	t=[0]
	good=[-1]
	tmp=mdsvalue('\SPECTROSCOPY::TOP.BOLOMETER.RESULTS.DIODE.'+array+':DEL_I',/quiet,status=status)
	IF status AND array_status THEN BEGIN
		path='\SPECTROSCOPY::TOP.BOLOMETER.RESULTS.DIODE.'+array+':BRIGHT'
		br=mdsvalue('_sig='+path)*1.0e-6
		rt=mdsvalue('dim_of(_sig,0)')
		t=mdsvalue('dim_of(_sig,1)')
		br_err=mdsvalue('dim_of(_sig,2)',/quiet,status=err_stat)
		rt_err=mdsvalue('\SPECTROSCOPY::TOP.BOLOMETER.DIODE_CALIB.'+array+':DRT',/quiet,status=err_stat)
		IF err_stat THEN BEGIN
			br_err*=1.0e-6
		ENDIF
		good=mdsvalue('\SPECTROSCOPY::TOP.BOLOMETER.DIODE_CALIB.'+array+':GOOD')
	ENDIF
	mdsclose,'spectroscopy',shot
END

PRO calc_cosine_profile,em,rmid,r,z,raxis,zaxis,emr,emc,rpts,err=err,coserr=coserr,thlfs=thlfs,thhfs=thhfs
	cent=minloc(rmid)
	rpts=rmid[0:cent-1]
	emhfs=em[0:cent-1]
	drhfs=(r[0:cent-1]-raxis)
	dzhfs=drhfs*0.0+(z-zaxis)	
	thhfs=acos(drhfs/sqrt((drhfs^2+dzhfs^2)))
	emlfs=interpol(em[cent:*],rmid[cent:*],rpts)
	drlfs=(r[cent:*]-raxis)
	dzlfs=drlfs*0.0+(z-zaxis)
	thlfs=acos(drlfs/sqrt((drlfs^2+dzlfs^2)))
	thlfs=interpol(thlfs,rmid[cent:*],rpts)	
	emc=(emlfs-emhfs)/(cos(thlfs)-cos(thhfs))
	emr=emlfs-emc*cos(thlfs)
	IF keyword_set(err) THEN BEGIN
		errhfs=err[0:cent-1]
		errlfs=interpol(err[cent:*],rmid[cent:*],rpts)
		coserr=0.5*sqrt(errlfs^2+errhfs^2)
	ENDIF ELSE coserr=-1
END

PRO plot_bright_time,u
	tau=u.plot.tau[where(u.plot.tau NE -1)]
	intime=u.stat.time
	FOR i=0,n(tau) DO BEGIN
		u.stat.time=tau[i]
		plot_bright,u
	ENDFOR
	u.stat.time=intime	
END

PRO plot_bright,u
	IF u.stat.ps THEN BEGIN
		xsize=7.0
		ysize=5.0
		ls=1.2
		col=u.stat.pscol 
		tit=num2str(u.shot,1)+' t='+num2str(u.stat.time,dp=3)
		d_old=!d
		device,xsize=xsize,ysize=ysize,/in
	ENDIF ELSE BEGIN
		ls=1.0
		col=u.stat.col
		widget_control,u.id.draw1,get_value=draw_win
		window,0,xsize=u.plot.bsize[0],ysize=u.plot.bsize[1],/pixmap
		tit=''
	ENDELSE
	xr=[u.plot.low[0],u.plot.up[0]]
	yr=[u.plot.low[1],u.plot.up[1]]

	plot,[0],[0],xr=xr,yr=yr,xtit='R!lTANG!n [m]',ytit='Brightness [MW/m!u2!n]',/xsty,/ysty,tit=tit,chars=ls
	narr=n(u.stat.arr)+1
	FOR i=0,narr-1 DO BEGIN
		IF u.stat.larr[i] AND u.plot.parr[i] THEN BEGIN
			index=ipt(u.dat.br.(i).t,u.stat.time)
			xplot=u.dat.br.(i).rt
			yplot=*u.dat.br.(i).br[index]
			xerr=u.dat.br.(i).rterr
			yerr=*u.dat.br.(i).err[index]
			IF u.plot.psym[i] EQ 8 THEN makesym,u.plot.msym[i]
			tmp=where(u.dat.br.(i).good EQ 1)
			IF u.stat.err[i] THEN oploterror,xplot[tmp],yplot[tmp],xerr[tmp],yerr[tmp],psym=u.plot.psym[i],color=col[i],symsize=0.5,errcolor=col[i] ELSE $
				oplot,xplot[tmp],yplot[tmp],psym=u.plot.psym[i],color=col[i],symsize=0.75
			xplot=u.dat.br.(i).chkrad
			yplot=*u.dat.br.(i).brchk[index]
			oplot,xplot,yplot,col=col[i],linestyle=u.plot.sty[i]	
		ENDIF
	ENDFOR

	IF NOT u.stat.ps THEN BEGIN
		wset,draw_win
		device,copy=[0,0,u.plot.bsize[0],u.plot.bsize[1],0,0,0]
	ENDIF ELSE device, xsize=float(d_old.x_size)/d_old.x_px_cm,ysize=float(d_old.y_size)/d_old.y_px_cm
END

PRO plot_paper_emiss,u,nt=nt,nm=nm,thick=thick
	IF u.stat.ps THEN BEGIN		
		ls=0.85
		IF keyword_set(thick) THEN BEGIN
			!p.thick=6
			!p.charthick=2.0
			!x.thick=5.0
			!y.thick=5.0	
			ls=0.85
		ENDIF
		xsize=7.5
		ysize=4.0
		col=u.stat.pscol 
		tit=num2str(u.shot,1)+'  t='+num2str(u.stat.time,dp=3)
		d_old=!d
		device,xsize=xsize,ysize=ysize,/in
	ENDIF ELSE BEGIN
		ls=1.0
		col=u.stat.col
		window,0,xsize=u.plot.csize[0],ysize=u.plot.csize[1]
		tit=' '
	ENDELSE
	IF NOT keyword_set(nt) THEN nt=3
	IF NOT keyword_set(nm) THEN nm=3
	lfscol=35
	hfscol=220

	xr=[u.plot.low[2],u.plot.up[2]]
	yr=[u.plot.low[4],u.plot.up[4]]
	plot,[0],[0],xr=xr,yr=yr,xtit='R [m]',ytit='Emissivity [MW/m!u3!n]',/xsty,/ysty,chars=ls,pos=[0.085,0.085,0.55,0.98],yticks=nt,yminor=nm
	narr=n(u.stat.arr)+1
	FOR i=0,narr-1 DO BEGIN
		IF u.stat.larr[i] AND u.plot.parr[i] THEN BEGIN
			index=ipt(u.dat.em.(i).t,u.stat.time)
			irmid=reform(u.dat.rmid[ipt(u.dat.t,u.stat.time),*])
			xplot=u.dat.em.(i).r
			yplot=*u.dat.em.(i).em[index]
			yerr=*u.dat.em.(i).err[index]
			IF u.plot.psym[i] EQ 8 THEN makesym,u.plot.msym[i]
			IF u.stat.arr[i] EQ 'axj' THEN BEGIN
				loadct,39,/silent
				tmp=where(xplot LT irmid[0])
				oploterror,xplot[tmp],yplot[tmp],yerr[tmp],color=hfscol,linestyle=u.plot.sty[i],errcolor=hfscol
				tmp=where(xplot GE irmid[0])
				oploterror,xplot[tmp],yplot[tmp],yerr[tmp],color=lfscol,linestyle=u.plot.sty[i],errcolor=lfscol
				loadct,12,/silent
			ENDIF ELSE oploterror,xplot,yplot,yerr,color=col[i],linestyle=u.plot.sty[i],errcolor=col[i] 

			IF NOT u.plot.rho THEN BEGIN
				oplot,irmid[0]*[1,1],yr,linestyle=last(u.plot.sty),color=last(col)
				oplot,last(irmid)*[1,1],yr,linestyle=last(u.plot.sty),color=last(col)	
			ENDIF
			oplot,[0.6,0.65],yr[0]+(i+1)*0.1*(yr[1]-yr[0])*[1,1],linestyle=u.plot.sty[i]
			xyouts,0.67,yr[0]+(i+1)*0.1*(yr[1]-yr[0])-0.01*(yr[1]-yr[0]),strupcase(u.stat.arr[i])
		ENDIF
	ENDFOR

	xr=[u.plot.low[3],u.plot.up[3]]
	yr=[u.plot.low[4],u.plot.up[4]]


	plot,[0],[0],xr=xr,yr=yr,xtit='r/a',/xsty,/ysty,chars=ls,/noerase,pos=[0.55,0.085,0.97,0.98],$
		yticks=nt,yminor=nm,ytickname=replicate(' ',nt+2)
	narr=n(u.stat.arr)+1
	u.plot.parr[0]=0
	FOR i=0,narr-1 DO BEGIN
		IF u.stat.larr[i] AND u.plot.parr[i] THEN BEGIN
			irmid=reform(u.dat.rmid[ipt(u.dat.t,u.stat.time),*])
			index=ipt(u.dat.em.(i).t,u.stat.time)
			xplot=*u.dat.em.(i).rmid[index]
			yplot=*u.dat.em.(i).em[index]
			yerr=*u.dat.em.(i).err[index]
			xerr=*u.dat.em.(i).rerr[index]
			IF u.plot.psym[i] EQ 8 THEN makesym,u.plot.msym[i]
			IF u.plot.rho THEN BEGIN
				xplot=(xplot-irmid[0])/(last(irmid)-irmid[0])
				xerr=xerr/(last(irmid)-irmid[0])
			ENDIF ELSE BEGIN
				oplot,irmid[0]*[1,1],yr,linestyle=last(u.plot.sty),color=last(col)
				oplot,last(irmid)*[1,1],yr,linestyle=last(u.plot.sty),color=last(col)
			ENDELSE
			k=minloc(xplot)
			loadct,39,/silent
			oploterror,xplot[0:k],yplot[0:k],xerr[0:k],yerr[0:k],color=hfscol,linestyle=u.plot.sty[i],errcolor=hfscol
			oploterror,xplot[k+1:*],yplot[k+1:*],xerr[k+1:*],yerr[k+1:*],color=lfscol,linestyle=u.plot.sty[i],errcolor=lfscol

			j=i
			oplot,[0.15,0.25],yr[0]+(j+1)*0.1*(yr[1]-yr[0])*[1,1],linestyle=u.plot.sty[i],color=lfscol
			xyouts,0.27,yr[0]+(j+1)*0.1*(yr[1]-yr[0])-0.01*(yr[1]-yr[0]),strupcase(u.stat.arr[i])+' - LFS',color=lfscol
			j=i-1
			oplot,[0.15,0.25],yr[0]+(j+1)*0.1*(yr[1]-yr[0])*[1,1],linestyle=u.plot.sty[i],color=hfscol
			xyouts,0.27,yr[0]+(j+1)*0.1*(yr[1]-yr[0])-0.01*(yr[1]-yr[0]),strupcase(u.stat.arr[i])+' - HFS',color=hfscol
			loadct,12,/silent
		ENDIF
	ENDFOR
	u.plot.parr[0]=1
	IF yr[0] LT 0 THEN oplot,xr,[0,0],linestyle=2,color=last(col)
	xyouts,xr[1]+0.04*(xr[1]-xr[0]),yr[0]+0.01*(yr[1]-yr[0]),tit,orient=90,chars=ls*0.8
	IF NOT u.stat.ps THEN BEGIN

	ENDIF ELSE device, xsize=float(d_old.x_size)/d_old.x_px_cm,ysize=float(d_old.y_size)/d_old.y_px_cm

END

PRO plot_rmaj_emiss,u
	IF u.stat.ps THEN BEGIN
		xsize=7.0
		ysize=5.0
		ls=1.2
		col=u.stat.pscol 
		tit=num2str(u.shot,1)+' t='+num2str(u.stat.time,dp=3)
		d_old=!d
		device,xsize=xsize,ysize=ysize,/in
	ENDIF ELSE BEGIN
		ls=1.0
		col=u.stat.col
		widget_control,u.id.draw3,get_value=draw_win
		window,0,xsize=u.plot.csize[0],ysize=u.plot.csize[1],/pixmap
	ENDELSE
	xr=[u.plot.low[2],u.plot.up[2]]
	yr=[u.plot.low[4],u.plot.up[4]]

	plot,[0],[0],xr=xr,yr=yr,xtit='R [m]',ytit='Emissivity [MW/m!u3!n]',/xsty,/ysty,tit=tit,chars=ls
	narr=n(u.stat.arr)+1
	FOR i=0,narr-1 DO BEGIN
		IF u.stat.larr[i] AND u.plot.parr[i] THEN BEGIN
			index=ipt(u.dat.em.(i).t,u.stat.time)
			irmid=reform(u.dat.rmid[ipt(u.dat.t,u.stat.time),*])
			xplot=u.dat.em.(i).r
			yplot=*u.dat.em.(i).em[index]
			yerr=*u.dat.em.(i).err[index]
			IF u.plot.psym[i] EQ 8 THEN makesym,u.plot.msym[i]
			IF u.stat.err[i] THEN oploterror,xplot,yplot,yerr,color=col[i],linestyle=u.plot.sty[i],errcolor=col[i] ELSE $
				oplot,xplot,yplot,color=col[i],linestyle=u.plot.sty[i]
			IF NOT u.plot.rho THEN BEGIN
				oplot,irmid[0]*[1,1],yr,linestyle=last(u.plot.sty),color=last(col)
				oplot,last(irmid)*[1,1],yr,linestyle=last(u.plot.sty),color=last(col)	
			ENDIF
		ENDIF
	ENDFOR

	IF NOT u.stat.ps THEN BEGIN
		wset,draw_win
		device,copy=[0,0,u.plot.csize[0],u.plot.csize[1],0,0,0]
	ENDIF ELSE device, xsize=float(d_old.x_size)/d_old.x_px_cm,ysize=float(d_old.y_size)/d_old.y_px_cm
END

PRO plot_rmaj_time,u
	tau=u.plot.tau[where(u.plot.tau NE -1)]
	intime=u.stat.time
	FOR i=0,n(tau) DO BEGIN
		u.stat.time=tau[i]
		plot_rmaj_emiss,u
	ENDFOR
	u.stat.time=intime	
END

PRO plot_rmid_emiss,u
	IF u.stat.ps THEN BEGIN
		xsize=6.0
		ysize=6.0
		ls=1.2
		col=u.stat.pscol
		tit=num2str(u.shot,1)+' t='+num2str(u.stat.time,dp=3)
		d_old=!d
		device,xsize=xsize,ysize=ysize,/in
	ENDIF ELSE BEGIN
		ls=1.0
		col=u.stat.col
		widget_control,u.id.draw2,get_value=draw_win
		window,0,xsize=u.plot.esize[0],ysize=u.plot.esize[1],/pixmap
	ENDELSE
	xr=[u.plot.low[3],u.plot.up[3]]
	yr=[u.plot.low[4],u.plot.up[4]]
	IF u.plot.rho THEN xtit='r/a' ELSE xtit='RMID [m]'

	plot,[0],[0],xr=xr,yr=yr,xtit=xtit,ytit='Emissivity [MW/m!u3!n]',/xsty,/ysty,tit=tit,chars=ls
	narr=n(u.stat.arr)+1
	FOR i=0,narr-1 DO BEGIN
		IF u.stat.larr[i] AND u.plot.parr[i] THEN BEGIN
			irmid=reform(u.dat.rmid[ipt(u.dat.t,u.stat.time),*])
			index=ipt(u.dat.em.(i).t,u.stat.time)
			xplot=*u.dat.em.(i).rmid[index]
			yplot=*u.dat.em.(i).em[index]
			yerr=*u.dat.em.(i).err[index]
			xerr=*u.dat.em.(i).rerr[index]
			IF u.plot.psym[i] EQ 8 THEN makesym,u.plot.msym[i]
			IF u.plot.rho THEN BEGIN
				xplot=(xplot-irmid[0])/(last(irmid)-irmid[0])
				xerr=xerr/(last(irmid)-irmid[0])
			ENDIF ELSE BEGIN
				oplot,irmid[0]*[1,1],yr,linestyle=last(u.plot.sty),color=last(col)
				oplot,last(irmid)*[1,1],yr,linestyle=last(u.plot.sty),color=last(col)
			ENDELSE
			IF u.stat.err[i] THEN oploterror,xplot,yplot,xerr,yerr,color=col[i],linestyle=u.plot.sty[i],errcolor=col[i] ELSE $
				oplot,xplot,yplot,color=col[i],linestyle=u.plot.sty[i]
		ENDIF
	ENDFOR
	IF yr[0] LT 0 THEN oplot,xr,[0,0],linestyle=2,color=last(col)

	IF NOT u.stat.ps THEN BEGIN
		wset,draw_win
		device,copy=[0,0,u.plot.esize[0],u.plot.esize[1],0,0,0]
	ENDIF ELSE device, xsize=float(d_old.x_size)/d_old.x_px_cm,ysize=float(d_old.y_size)/d_old.y_px_cm
END

PRO plot_paper_cos,u,z=z,ffmin=ffmin,pcol=pcol,nofric=nofric,notheory=notheory,tht=tht,zeff=zeff,xplotr=xplotr,pos=pos,noerase=noerase

	posd=[0.1,0.1,0.95,0.98]
	IF NOT keyword_set(xplotr) THEN xplotr=[0.1,0.87]
	IF NOT keyword_set(pos) THEN pos=posd
	IF u.stat.ps THEN BEGIN
		xsize=4.5
		ysize=4.0
		ls=0.95
		col=u.stat.pscol 
		tit=num2str(u.shot,1)+' t='+num2str(u.stat.time,dp=3)
		d_old=!d
		IF total(pos-posd) EQ 0 THEN  device,xsize=xsize,ysize=ysize,/in
	ENDIF ELSE BEGIN
		ls=1.0
		col=u.stat.col
		window,0,xsize=u.plot.esize[0],ysize=u.plot.esize[1]
	ENDELSE
	xr=[u.plot.low[3],u.plot.up[3]]
	IF u.plot.ratio THEN BEGIN
		yr=[u.plot.low[6],u.plot.up[6]]
		ytit='!S!A~!R!nn!iz!n/<n>'
	ENDIF ELSE BEGIN
		yr=[u.plot.low[4],u.plot.up[4]]
		ytit='Emissivity [MW/m!u3!n]
	ENDELSE	
	IF u.plot.rho THEN xtit='r/a' ELSE xtit='RMID [m]'	
	IF keyword_set(noerase) THEN ytit=''

	plot,[0],[0],xr=xr,yr=yr,xtit=xtit,ytit=ytit,/xsty,/ysty,chars=ls,pos=pos,noerase=noerase
	narr=n(u.stat.arr)+1
	u.plot.parr[0] = 0
	FOR i=0,narr-1 DO BEGIN
		IF u.stat.larr[i] AND u.plot.parr[i] THEN BEGIN
			irmid=reform(u.dat.rmid[ipt(u.dat.t,u.stat.time),*])
			index=ipt(u.dat.em.(i).t,u.stat.time)			
			xplot=*u.dat.emc.(i).r[index]
			emr=*u.dat.emc.(i).emr[index]
			emc=*u.dat.emc.(i).emc[index]
			err=*u.dat.emc.(i).err[index]

			IF u.plot.psym[i] EQ 8 THEN makesym,u.plot.msym[i]
			IF u.plot.rho THEN BEGIN
				xplot=(xplot-irmid[0])/(last(irmid)-irmid[0])
			ENDIF ELSE BEGIN
				oplot,irmid[0]*[1,1],yr,linestyle=last(u.plot.sty),color=last(col)
				oplot,last(irmid)*[1,1],yr,linestyle=last(u.plot.sty),color=last(col)
			ENDELSE
			IF u.plot.ratio THEN BEGIN
				IF NOT keyword_set(pcol) THEN pcol=col[i]
				yplot=emc/emr
				yerr=sqrt(yplot^2*err^2*(1.0/emr^2+1.0/emc^2))
				tmp=where(xplot GE xplotr[0] AND xplot LE xplotr[1])
				oploterror,xplot[tmp],yplot[tmp],yerr[tmp],color=pcol,errcolor=pcol
			ENDIF ELSE BEGIN
				IF u.stat.err[i] THEN BEGIN
					oploterror,xplot,emr,err,color=col[i],linestyle=0,errcolor=col[i]
					oploterror,xplot,emc,err,color=col[i],linestyle=3,errcolor=col[i]
				ENDIF ELSE BEGIN
					oplot,xplot,emr,color=col[i],linestyle=0
					oplot,xplot,emc,color=col[i],linestyle=3
				ENDELSE
			ENDELSE
			
		ENDIF
	ENDFOR
	xyouts,1.045,yr[0]+0.05*(yr[1]-yr[0]),num2str(u.shot,1)+'  t='+num2str(u.stat.time,dp=3),orient=90,chars=0.7*ls,color=pcol
	oplot,[0.03,0.09],(yr[0]+0.85*(yr[1]-yr[0]))*[1,1],color=pcol
	xyouts,0.12,yr[0]+0.84*(yr[1]-yr[0]),'EXPERIMENT',color=pcol,chars=0.9*ls

	IF NOT keyword_set(notheory) THEN BEGIN
		u.plot.parr[0] = 1
		IF yr[0] LT 0 THEN oplot,xr,[0,0],linestyle=2,color=last(col)
		IF NOT keyword_set(z) THEN z=42
		IF NOT keyword_set(ffmin) THEN ffmin=[1.0]
		efit=load_efit_struc(u.shot,time=u.stat.time)
		prof=load_profile_struc(efit,z,/inst,tht=tht)
		IF keyword_set(zeff) THEN prof.zeff=zeff
		linestyle=[2,3,1]
		FOR i=0,n(ffmin) DO BEGIN
			fh_numerical_asym,efit,prof,asym,/ban,ffmin=ffmin[i],nofric=nofric
			efitrho=(efit.rmid-efit.rmid[0])/(last(efit.rmid)-efit.rmid[0])
			tmp=where(efitrho GE xplotr[0] AND efitrho LE xplotr[1])
			oploterror,efitrho[tmp],asym.ncos[tmp],asym.cerr[tmp],linestyle=linestyle[i]
			oplot,[0.03,0.09],(yr[0]+(0.94-0.08*i+0.01)*(yr[1]-yr[0]))*[1,1],linestyle=linestyle[i]
			IF ffmin[i] LT 0.1 THEN thstr='THEORY (INERTIA)' ELSE thstr='THEORY (INERTIA+ICRH) '
			IF ffmin[i] GT 0.1 AND ffmin[i] LT 1.0 THEN thstr=thstr+'@'+num2str(ffmin[i]*100,dp=0)+'%'
			xyouts,0.12,yr[0]+(0.93-0.08*i)*(yr[1]-yr[0]),thstr,chars=0.9*ls
		ENDFOR
	ENDIF
	
	IF NOT u.stat.ps THEN BEGIN

	ENDIF ELSE device, xsize=float(d_old.x_size)/d_old.x_px_cm,ysize=float(d_old.y_size)/d_old.y_px_cm
END	

PRO plot_cos_emiss,u
	IF u.stat.ps THEN BEGIN
		xsize=6.0
		ysize=6.0
		ls=1.2
		col=u.stat.pscol 
		tit=num2str(u.shot,1)+' t='+num2str(u.stat.time,dp=3)
		d_old=!d
		device,xsize=xsize,ysize=ysize,/in
	ENDIF ELSE BEGIN
		ls=1.0
		col=u.stat.col
		widget_control,u.id.draw2,get_value=draw_win
		window,0,xsize=u.plot.esize[0],ysize=u.plot.esize[1],/pixmap
	ENDELSE
	xr=[u.plot.low[3],u.plot.up[3]]
	IF u.plot.ratio THEN BEGIN
		yr=[u.plot.low[6],u.plot.up[6]]
		ytit='EMC/EMR'
	ENDIF ELSE BEGIN
		yr=[u.plot.low[4],u.plot.up[4]]
		ytit='Emissivity [MW/m!u3!n]'
	ENDELSE	
	IF u.plot.rho THEN xtit='r/a' ELSE xtit='RMID [m]'	

	plot,[0],[0],xr=xr,yr=yr,xtit=xtit,ytit=ytit,/xsty,/ysty,tit=tit,chars=ls
	narr=n(u.stat.arr)+1
	FOR i=0,narr-1 DO BEGIN
		IF u.stat.larr[i] AND u.plot.parr[i] THEN BEGIN
			irmid=reform(u.dat.rmid[ipt(u.dat.t,u.stat.time),*])
			index=ipt(u.dat.em.(i).t,u.stat.time)			
			xplot=*u.dat.emc.(i).r[index]
			emr=*u.dat.emc.(i).emr[index]
			emc=*u.dat.emc.(i).emc[index]
			err=*u.dat.emc.(i).err[index]

			IF u.plot.psym[i] EQ 8 THEN makesym,u.plot.msym[i]
			IF u.plot.rho THEN BEGIN
				xplot=(xplot-irmid[0])/(last(irmid)-irmid[0])
			ENDIF ELSE BEGIN
				oplot,irmid[0]*[1,1],yr,linestyle=last(u.plot.sty),color=last(col)
				oplot,last(irmid)*[1,1],yr,linestyle=last(u.plot.sty),color=last(col)
			ENDELSE
			IF u.plot.ratio THEN BEGIN
				yplot=emc/emr
				yerr=sqrt(yplot^2*err^2*(1.0/emr^2+1.0/emc^2))
				IF u.stat.err[i] THEN oploterror,xplot,yplot,yerr,color=col[i],errcolor=col[i] ELSE oplot,xplot,yplot,color=col[i]
			ENDIF ELSE BEGIN
				IF u.stat.err[i] THEN BEGIN
					oploterror,xplot,emr,err,color=col[i],linestyle=0,errcolor=col[i]
					oploterror,xplot,emc,err,color=col[i],linestyle=3,errcolor=col[i]
				ENDIF ELSE BEGIN
					oplot,xplot,emr,color=col[i],linestyle=0
					oplot,xplot,emc,color=col[i],linestyle=3
				ENDELSE
			ENDELSE
			
		ENDIF
	ENDFOR
	IF yr[0] LT 0 THEN oplot,xr,[0,0],linestyle=2,color=last(col)

	IF NOT u.stat.ps THEN BEGIN
		wset,draw_win
		device,copy=[0,0,u.plot.esize[0],u.plot.esize[1],0,0,0]
	ENDIF ELSE device, xsize=float(d_old.x_size)/d_old.x_px_cm,ysize=float(d_old.y_size)/d_old.y_px_cm
END

PRO plot_emiss,u
	IF u.plot.cos THEN plot_cos_emiss,u ELSE plot_rmid_emiss,u
END

PRO plot_cos_time,u
	IF u.stat.ps THEN BEGIN
		xsize=6.0
		ysize=6.0
		ls=1.2
		col=u.stat.pscol 
		tit=num2str(u.shot,1)+'
		d_old=!d
		device,xsize=xsize,ysize=ysize,/in
	ENDIF ELSE BEGIN
		ls=1.0
		col=u.stat.col
		widget_control,u.id.draw2,get_value=draw_win
		window,0,xsize=u.plot.esize[0],ysize=u.plot.esize[1],/pixmap
	ENDELSE
	
	xr=[u.plot.low[3],u.plot.up[3]]
	dx=xr[1]-xr[0]
	yr=[u.plot.low[6],u.plot.up[6]]
	dy=yr[1]-yr[0]
	ytit='EMC/EMR'
	IF u.plot.rho THEN xtit='r/a' ELSE xtit='RMID [m]'	
	plot,[0],[0],xr=xr,yr=yr,xtit=xtit,ytit=ytit,/xsty,/ysty,tit=tit,chars=ls

	j=1			;hardcode AXJ
	tau=u.plot.tau[where(u.plot.tau NE -1)]
	col=colormap(tau)
	
	FOR i=0,n(tau) DO BEGIN
		index=ipt(u.dat.em.(j).t, tau[i])
		irmid=reform(u.dat.rmid[ipt(u.dat.t,tau[i]),*])
		xplot=*u.dat.emc.(j).r[index]
		emc=*u.dat.emc.(j).emc[index]
		emr=*u.dat.emc.(j).emr[index]
		err=*u.dat.emc.(j).err[index]
		IF u.plot.rho THEN xplot=(xplot-irmid[0])/(last(irmid)-irmid[0])	
		yplot=emc/emr
		yerr=sqrt(yplot^2*err^2*(1.0/emr^2+1.0/emc^2))
		oploterror,xplot,yplot,yerr,color=col[i],errcolor=col[i]
		xyouts,xr[0]+0.1*dx,yr[1]-(i+1)*0.05*dy,num2str(tau[i],dp=4),col=col[i]
	ENDFOR
	IF NOT u.stat.ps THEN BEGIN
		wset,draw_win
		device,copy=[0,0,u.plot.esize[0],u.plot.esize[1],0,0,0]
	ENDIF ELSE device, xsize=float(d_old.x_size)/d_old.x_px_cm,ysize=float(d_old.y_size)/d_old.y_px_cm		
END

PRO plot_legend,win
	;plot legend
	dy=yr[1]-yr[0]
	yo=yr[0]
	dx=xr[1]-xr[0]
	xo=xr[0]-pdx*dx
	oplot,xo+dx*[0.8,0.85],yo+dy*0.9*[1.0,1.0],color=col[0],linestyle=linestyle
	xyouts,xo+dx*0.9,yo+dy*0.9,'RADIAL ONLY'
	IF u.stat.m1stat THEN BEGIN
		oplot,xo+dx*[0.8,0.85],yo+dy*0.85*[1,1],color=col[1]
		xyouts,xo+dx*0.9,yo+dy*0.85,'m=0 RAD',color=col[1]
		oplot,xo+dx*[0.8,0.85],yo+dy*0.8*[1,1],color=col[2]
		xyouts,xo+dx*0.9,yo+dy*0.8,'m=1 COS',color=col[2]
		oplot,xo+dx*[0.8,0.85],yo+dy*0.75*[1,1],color=col[4]
		xyouts,xo+dx*0.9,yo+dy*0.75,'m=1 SIN',color=col[4]
	ENDIF
END

PRO plot_cont,u
	index=u.index
	IF u.plot.rho THEN BEGIN
		r=u.dat.rho 
		ylab='r/a'
	ENDIF ELSE BEGIN
		r=*u.dat.r[index]
		ylab='RMID [m]'
	ENDELSE
	CASE u.plot.cont OF
		0 : BEGIN
			labels={ilab:n2g('theta')+' ['+n2g('pi')+']',jlab:ylab,klab:'Emissivity [kW/m!u3!n]',ctit:'2D Emissivity Profiles',$
				itit:'Radial Profiles at Fixed '+n2g('theta'),jtit:'Polar Profiles at Fixed r/a'}
			emr=*u.dat.emr[index]-u.dat.emrdel
			emc1=*u.dat.emc1[index]-u.dat.emc1del
			emc2=*u.dat.emc2[index]-u.dat.emc2del
			ems1=*u.dat.ems1[index]-u.dat.ems1del
			th=make(0,2.0*!pi,u.stat.nth)
			array=fltarr(u.stat.nth,u.stat.nr)
			FOR i=0,u.stat.nth-1 DO BEGIN
				FOR j=0,u.stat.nr-1 DO BEGIN
					IF u.stat.m2stat THEN array[i,j]=emr[j]+emc1[j]*cos(th[i])+emc2[j]*cos(2.0*th[i])+ems1[j]*sin(th[i]) ELSE $
						array[i,j]=emr[j]+emc1[j]*cos(th[i])+ems1[j]*sin(th[i])
				ENDFOR
			ENDFOR
			ivec=th/!pi
			jvec=r
			IF NOT u.plot.auto[3] THEN jr=[u.plot.low[3],u.plot.up[3]]
		END
		1 : BEGIN
			labels={jlab:'t [sec]',ilab:ylab,klab:'Emissivity [kW/m!u3!n]',ctit:'Radial Emissivity Profiles',$
				jtit:'Radial Profiles at Fixed Time',itit:'Temporal Profiles at Fixed Radius'}
			array=*u.dat.xemr
			;FOR i=0,u.stat.ntime-1 DO array[*,i]-=u.dat.emrdel
			jvec=u.dat.t
			ivec=r
			jr=[u.plot.low[5],u.plot.up[5]]
			ir=[u.plot.low[3],u.plot.up[3]]
			
		END
		2 : BEGIN
			labels={jlab:'t [sec]',ilab:ylab,klab:'Emissivity [kW/m!u3!n]',ctit:'m=1Cosine Emissivity Profiles',$
				jtit:'Radial Profiles at Fixed Time',itit:'Temporal Profiles at Fixed Radius'}
			array=*u.dat.xemc1
			IF u.plot.ratio THEN BEGIN
				emr=*u.dat.xemr
				array=array/emr
			ENDIF
			;FOR i=0,u.stat.ntime-1 DO array-=u.dat.emc1del
			jvec=u.dat.t
			ivec=r
			jr=[u.plot.low[5],u.plot.up[5]]
			ir=[u.plot.low[3],u.plot.up[3]]
		END
		3 : BEGIN
			labels={jlab:'t [sec]',ilab:ylab,klab:'Emissivity [kW/m!u3!n]',ctit:'m=2 Cosine Emissivity Profiles',$
				jtit:'Radial Profiles at Fixed Time',itit:'Temporal Profiles at Fixed Radius'}
			array=*u.dat.xemc2
			;FOR i=0,u.stat.ntime-1 DO array[*,i]-=u.dat.emc2del
			jvec=u.dat.t
			ivec=r
			jr=[u.plot.low[5],u.plot.up[5]]
			ir=[u.plot.low[3],u.plot.up[3]]
		END
		4 : BEGIN
			labels={jlab:'t [sec]',ilab:ylab,klab:'Emissivity [kW/m!u3!n]',ctit:'m=1Sine Emissivity Profiles',$
				jtit:'Radial Profiles at Fixed Time',itit:'Temporal Profiles at Fixed Radius'}
			array=*u.dat.xems1
			;FOR i=0,u.stat.ntime-1 DO array[*,i]-=u.dat.ems1del
			jvec=u.dat.t
			ivec=r
			jr=[u.plot.low[5],u.plot.up[5]]
			ir=[u.plot.low[3],u.plot.up[3]]
		END
	ENDCASE

	widget_control,u.id.draw3,get_value=dw3
	widget_control,u.id.draw4,get_value=dw4
	widget_control,u.id.draw5,get_value=dw5
	io=*u.plot.io
	jo=*u.plot.jo
	IF NOT u.plot.auto[4] THEN BEGIN
		minpt=u.plot.low[4]
		maxpt=u.plot.up[4]
		IF u.plot.ratio THEN BEGIN
			minpt=-1.0
			maxpt=1.0
		ENDIF
	ENDIF ELSE BEGIN
		minpt=0.0
		maxpt=0.0
	ENDELSE
	genplt,array,ivec,jvec,cct=39,ncntrs=60,labels=labels,pixmap=u.plot.csize,set=set,jr=jr,ir=ir,io=io,jo=jo,minpt=minpt,maxpt=maxpt
	IF n(*u.plot.io) EQ 0 THEN BEGIN
		*u.plot.io=[io[0],io]
		iostr=num2str(io[0],dp=2)
		FOR i=1,n(io) DO iostr=iostr+', '+num2str(io[i],dp=2)
		u.plot.iostr=iostr
		widget_control,u.id.io_val,set_value=u.plot.iostr
	ENDIF
	IF n(*u.plot.jo) EQ 0 THEN BEGIN
		*u.plot.jo=[jo[0],jo]
		jostr=num2str(jo[0],dp=2)
		FOR i=1,n(jo) DO jostr=jostr+', '+num2str(jo[i],dp=2)
		u.plot.jostr=jostr
		widget_control,u.id.jo_val,set_value=u.plot.jostr
	ENDIF
	IF NOT u.stat.ps THEN BEGIN
		wset,dw3
		device,copy=[0,0,u.plot.csize[0],u.plot.csize[1],0,0,set[0]]
		wset,dw4
		device,copy=[0,0,u.plot.csize[2],u.plot.csize[3],0,0,set[1]]
		wset,dw5
		device,copy=[0,0,u.plot.csize[4],u.plot.csize[5],0,0,set[2]]
	ENDIF
END

PRO calc_delta,u
	;calculate the time-avereged profile to subtract off (good for showing injection evolution)
	tmp=where(u.dat.t GE u.plot.t1[0] AND u.dat.t LE u.plot.t1[1])
	IF tmp[0] NE -1 THEN BEGIN
		npts=n(tmp)+1.0
		brdel=u.dat.brdel*0.0
		emdel=u.dat.emdel*0.0
		emrdel=u.dat.emrdel*0.0
		emc1del=u.dat.emc1del*0.0
		emc2del=u.dat.emc2del*0.0
		ems1del=u.dat.ems1del*0.0
		FOR i=0,npts-1 DO BEGIN
			ibr=*u.dat.br[tmp[i]]
			brdel+=ibr
			iem=*u.dat.em[tmp[i]]
			emdel+=iem
			iemr=*u.dat.emr[tmp[i]]
			emrdel+=iemr
			iemc1=*u.dat.emc1[tmp[i]]
			emc1del+=iemc1
			iemc2=*u.dat.emc2[tmp[i]]
			emc2del+=iemc2
			iems1=*u.dat.ems1[tmp[i]]
			ems1del+=iems1
		ENDFOR
		u.dat.brdel=brdel/npts
		u.dat.emdel=emdel/npts
		u.dat.emrdel=emrdel/npts
		u.dat.emc1del=emc1del/npts
		u.dat.emc2del=emc2del/npts
		u.dat.ems1del=ems1del/npts
	ENDIF
END

PRO load_waxuv_data,u

	;load data from the tree
	shot=u.shot
	narr=n(u.stat.arr)+1

	;load EFIT rmid
	IF u.tree EQ 'ANALYSIS' THEN BEGIN
		mdsopen,'analysis',shot
		efit_rmid=mdsvalue('\efit_rmid')
		efit_time=mdsvalue('dim_of(\efit_rmid,0)')
		raxis=mdsvalue('\ANALYSIS::TOP:EFIT.RESULTS.A_EQDSK:RMAGX')/100.0
		zaxis=mdsvalue('\ANALYSIS::TOP:EFIT.RESULTS.A_EQDSK:ZMAGX')/100.0
		mdsclose,'analysis',shot
	ENDIF ELSE BEGIN
		mdsopen,u.tree,shot
		efit_rmid=mdsvalue('\efit_rmid')
		efit_time=mdsvalue('dim_of(\efit_rmid,0)')	
		raxis=mdsvalue('\'+u.tree+'::TOP.RESULTS.A_EQDSK:RMAGX')/100.0
		zaxis=mdsvalue('\'+u.tree+'::TOP.RESULTS.A_EQDSK:ZMAGX')/100.0
		mdsclose,u.tree,shot
	ENDELSE

	em={shot:shot,arr:u.stat.arr,z:fltarr(narr)}
	emc={shot:shot,arr:u.stat.arr}
	br={shot:shot,arr:u.stat.arr}
	FOR k=0,narr-1 DO BEGIN
		i=narr-1-k
		IF u.stat.larr[i] THEN BEGIN
			load_emiss_data,shot,u.stat.arr[i],1,iem,r,t,brchk,chkrad,rmid,z,em_err,rm_err,err_stat=err_stat,tree=loadtree
			load_bright_data,shot,u.stat.arr[i],1,ibr,rt,t,good,br_err,rt_err
			em.z[i]=z
			IF err_stat THEN BEGIN
				u.stat.err[i]=1
			ENDIF ELSE BEGIN
				rt_err=-1
				u.stat.err[i]=0
			ENDELSE
			widget_control,u.id.message,set_value='SHOT: '+num2str(shot,1)+' '+u.stat.arr[i]+'  LOADED - '+loadtree,/app
			ntime=n(t)+1
			xem=ptrarr(ntime,/allocate)
			xemerr=ptrarr(ntime,/allocate)			
			xrmid=ptrarr(ntime,/allocate)
			xrmerr=ptrarr(ntime,/allocate)
			xbr=ptrarr(ntime,/allocate)
			xbrerr=ptrarr(ntime,/allocate)
			xbrchk=ptrarr(ntime,/allocate)
			xemr=ptrarr(ntime,/allocate)
			xemc=ptrarr(ntime,/allocate)
			xerr2=ptrarr(ntime,/allocate)
			xrcos=ptrarr(ntime,/allocate)
			FOR j=0L,ntime-1 DO BEGIN
				;load profile data
				*xem[j]=iem[*,j]
				*xrmid[j]=rmid[*,j]
				*xbr[j]=ibr[*,j]
				*xbrchk[j]=brchk[*,j]
				IF u.stat.err[i] THEN BEGIN
					*xemerr[j]=em_err[*,j]
					*xrmerr[j]=rm_err[*,j]
					*xbrerr[j]=br_err[*,j]
				ENDIF ELSE BEGIN
					*xemerr[j]=-1
					*xrmerr[j]=-1
					*xbrerr[j]=-1
				ENDELSE
		
				;calculate cosine profiles
				cent=minloc(rmid[*,j])
				efit_index=ipt(efit_time,t[j])	
				IF cent NE 0 AND efit_index[0] NE -1 THEN BEGIN
					IF u.stat.err[i] THEN err=em_err[*,j] ELSE err=0
					calc_cosine_profile,iem[*,j],rmid[*,j],r,z,raxis[efit_index],zaxis[efit_index],yemr,yemc,yrpts,err=err,coserr=coserr
					*xemc[j]=yemc
					*xemr[j]=yemr
					*xrcos[j]=yrpts
					*xerr2[j]=coserr	
				ENDIF ELSE BEGIN
					*xemr[j]=-1
					*xemc[j]=-1
					*xerr2[j]=-1		
					*xrcos[j]=-1
				ENDELSE

			ENDFOR
			emstr={em:xem,r:r,t:t,rmid:xrmid,err:xemerr,rerr:xrmerr}
			emcstr={emc:xemc,emr:xemr,r:xrcos,err:xerr2}
			brstr={br:xbr,rt:rt,t:t,brchk:xbrchk,chkrad:chkrad,err:xbrerr,rterr:rt_err,good:good}
			em=create_struct(u.stat.arr[i],emstr,em)
			emc=create_struct(u.stat.arr[i],emcstr,emc)
			br=create_struct(u.stat.arr[i],brstr,br)		
		ENDIF ELSE BEGIN
			em=create_struct(u.stat.arr[i],-1,em)
			emc=create_struct(u.stat.arr[i],-1,emc)
			br=create_struct(u.stat.arr[i],-1,br)
		ENDELSE

	ENDFOR

	u.plot.low[5]=min(efit_time)
	widget_control,u.id.tmin,set_value=num2str(u.plot.low[5],dp=4)
	u.plot.up[5]=max(efit_time)
	widget_control,u.id.tmax,set_value=num2str(u.plot.up[5],dp=4)
	dat={em:em,emc:emc,br:br,rmid:efit_rmid,raxis:raxis,zaxis:zaxis,t:efit_time}

	u.stat.dat=1
	u={id:u.id,shot:u.shot,stat:u.stat,plot:u.plot,dat:dat,tree:u.tree}
	widget_control,u.id.base, set_uvalue=u	
END


PRO reset_cont_button,u
	id=[u.id.emr,u.id.emc]
	widget_control,id[u.plot.cont],set_button=0
END

PRO w_axuv_event,event
	widget_control,event.top,get_uvalue=u
	id = u.id
	tag = tag_names(event,/st)
	button=' '
	idtags=tag_names(id)
	FOR i=0,n(idtags) DO IF id.(i) EQ event.id THEN ename=idtags[i]
	CASE tag OF
		"WIDGET_BASE" : BEGIN

		END
		"WIDGET_BUTTON": BEGIN
			widget_control,event.id,get_value=button,get_uvalue=uvalue
			IF NOT keyword_set(uvalue) THEN uvalue='none'
			CASE ename OF
				"LOAD" : BEGIN
					WIDGET_CONTROL,/hourglass
					IF u.stat.dat THEN BEGIN
						heap_free,u.dat.em
						heap_free,u.dat.emc
						heap_free,u.dat.br
						heap_gc
					ENDIF
					load_waxuv_data,u
					plot_bright,u
					plot_rmaj_emiss,u
					plot_emiss,u
				END
				"QUIT" : BEGIN
					widget_control,event.top,/destroy
					heap_free,u
					heap_gc
					!except=1
				END
				"PRINT" : BEGIN
					u.stat.ps=1
					psplot
					plot_bright,u
					plot_rmaj_emiss,u
					plot_emiss,u
					psc
					xwplot
					set_plot,'x'
					u.stat.ps=0
				END
				"STOP" : BEGIN
					stop
				END
				"EMR" : BEGIN
					WIDGET_CONTROL,/hourglass
					reset_cont_button,u
					widget_control,u.id.emr,set_button=1
					u.plot.cont=1
					IF u.stat.m1stat THEN plot_cont,u
				END
				"EMC" : BEGIN
					WIDGET_CONTROL,/hourglass
					reset_cont_button,u
					widget_control,u.id.emc1,set_button=1
					u.plot.cont=2
					IF u.stat.m1stat THEN plot_cont,u
				END
				"RTAUTO" : IF event.select EQ 1 THEN BEGIN
						u.plot.auto[0]=1 
					   ENDIF ELSE u.plot.auto[0]=0
				"BRAUTO" : IF event.select EQ 1 THEN BEGIN
						u.plot.auto[1]=1 
					   ENDIF ELSE u.plot.auto[1]=0
				"RAUTO" : IF event.select EQ 1 THEN BEGIN
						u.plot.auto[2]=1 
					    ENDIF ELSE u.plot.auto[2]=0
				"RDAUTO" : IF event.select EQ 1 THEN BEGIN
				     		u.plot.auto[3]=1 
					  ENDIF ELSE u.plot.auto[3]=0
				"EMAUTO" : IF event.select EQ 1 THEN BEGIN
						u.plot.auto[4]=1 
					   ENDIF ELSE u.plot.auto[4]=0
				"TAUTO" : IF event.select EQ 1 THEN BEGIN
						u.plot.auto[5]=1 
					  ENDIF ELSE u.plot.auto[5]=0

				"XRMAJ" : BEGIN
					IF event.select EQ 1 THEN BEGIN
						widget_control,u.id.xrho,set_button=0
						u.plot.rho=0	
						u.plot.low[3]=0.67
						u.plot.up[3]=0.90
						widget_control,u.id.rdmax,set_value=num2str(u.plot.up[3],dp=2)
						widget_control,u.id.rdmin,set_value=num2str(u.plot.low[3],dp=2)
						plot_emiss,u
						*u.plot.io=0
						*u.plot.jo=0
					ENDIF ELSE widget_control,u.id.xrmaj,set_button=1
				END
				"XRHO" : BEGIN
					IF event.select EQ 1 THEN BEGIN
						widget_control,u.id.xrmaj,set_button=0
						u.plot.rho=1
						u.plot.low[3]=0
						u.plot.up[3]=1.0
						widget_control,u.id.rdmax,set_value=num2str(u.plot.up[3],dp=2)
						widget_control,u.id.rdmin,set_value=num2str(u.plot.low[3],dp=2)	
						plot_emiss,u
						*u.plot.io=0
						*u.plot.jo=0
					ENDIF ELSE widget_control,u.id.xrho,set_button=1
				END
				"PCOS" : BEGIN
					IF event.select EQ 1 THEN u.plot.cos=1 ELSE u.plot.cos=0
					plot_emiss,u
				END
				"PRAT" : BEGIN
					IF event.select EQ 1 THEN u.plot.ratio=1 ELSE u.plot.ratio=0
					plot_emiss,u
				END
				"LAXA" : IF event.select EQ 1 THEN u.stat.larr[0]=1 ELSE u.stat.larr[0]=0
				"LAXJ" : IF event.select EQ 1 THEN u.stat.larr[1]=1 ELSE u.stat.larr[1]=0
				"LWB1" : IF event.select EQ 1 THEN u.stat.larr[2]=1 ELSE u.stat.larr[2]=0
				"LWB2" : IF event.select EQ 1 THEN u.stat.larr[3]=1 ELSE u.stat.larr[3]=0
				"LWB3" : IF event.select EQ 1 THEN u.stat.larr[4]=1 ELSE u.stat.larr[4]=0
				"LWB4" : IF event.select EQ 1 THEN u.stat.larr[5]=1 ELSE u.stat.larr[5]=0
				"PAXA" : BEGIN
					IF event.select EQ 1 THEN u.plot.parr[0]=1 ELSE u.plot.parr[0]=0
					plot_bright,u
					plot_rmaj_emiss,u
					plot_emiss,u
				END
				"PAXJ" : BEGIN
					IF event.select EQ 1 THEN u.plot.parr[1]=1 ELSE u.plot.parr[1]=0
					plot_bright,u
					plot_rmaj_emiss,u
					plot_emiss,u
				END
				"PWB1" : BEGIN
					IF event.select EQ 1 THEN u.plot.parr[2]=1 ELSE u.plot.parr[2]=0
					plot_bright,u
					plot_rmaj_emiss,u
					plot_emiss,u
				END
				"PWB2" : BEGIN
					IF event.select EQ 1 THEN u.plot.parr[3]=1 ELSE u.plot.parr[3]=0
					plot_bright,u
					plot_rmaj_emiss,u
					plot_emiss,u
				END
				"PWB3" : BEGIN
					IF event.select EQ 1 THEN u.plot.parr[4]=1 ELSE u.plot.parr[4]=0
					plot_bright,u
					plot_rmaj_emiss,u
					plot_emiss,u
				END
				"PWB4" : BEGIN
					IF event.select EQ 1 THEN u.plot.parr[5]=1 ELSE u.plot.parr[5]=0
					plot_bright,u
					plot_rmaj_emiss,u
					plot_emiss,u
				END
				"PTIME" : plot_cos_time,u
				ELSE:
			ENDCASE
		END
		"WIDGET_SLIDER": BEGIN
			widget_control,event.id,get_value=slider
			IF NOT keyword_set(uvalue) THEN uvalue='none'
			CASE ename OF 
				'T_SLIDER' : BEGIN
					IF u.stat.dat THEN BEGIN
						newtime=slider/5.0e3
						IF newtime GE max(u.dat.t) THEN newtime=max(u.dat.t)
						IF newtime LE min(u.dat.t) THEN newtime=min(u.dat.t)
						IF u.stat.time NE newtime THEN BEGIN
							u.stat.time=newtime
							plot_bright,u
							plot_rmaj_emiss,u
							plot_emiss,u
							widget_control,u.id.t_text,set_value=num2str(u.stat.time,dp=4)
						ENDIF
					ENDIF
				END
				ELSE:
			ENDCASE
		END
   		"WIDGET_TEXT_CH": BEGIN
			widget_control,event.id,get_value=data
			CASE event.id OF
			u.id.shotid : u.shot=long(data)
			u.id.rmin : BEGIN
				u.plot.low[3]=float(data)
				plot_emiss,u
				plot_cont,u
			END		
			u.id.rmax : BEGIN
				u.plot.up[3]=float(data)
				plot_emiss,u
				plot_cont,u
			END
			u.id.tmin : BEGIN
				u.plot.low[5]=float(data)
				IF u.plot.cont NE 0 THEN plot_cont,u
			END		
			u.id.tmax : BEGIN
				u.plot.up[5]=float(data)
				IF u.plot.cont NE 0 THEN plot_cont,u
			END
			u.id.ratmin : BEGIN
				u.plot.low[6]=float(data)
				IF u.plot.ratio THEN plot_emiss,u
			END		
			u.id.ratmax : BEGIN
				u.plot.up[6]=float(data)
				IF u.plot.ratio THEN plot_emiss,u
			END
			u.id.brmin : BEGIN
				u.plot.low[1]=float(data)
				plot_bright,u
			END		
			u.id.brmax : BEGIN
				u.plot.up[1]=float(data)
				plot_bright,u
			END
			u.id.rdmin : BEGIN
				u.plot.low[3]=float(data)
				plot_emiss,u
			END		
			u.id.rdmax : BEGIN
				u.plot.up[3]=float(data)
				plot_emiss,u
			END
			u.id.emmin : BEGIN
				u.plot.low[4]=float(data)
				plot_emiss,u
				IF u.plot.cont NE 0 THEN plot_cont,u
			END		
			u.id.emmax : BEGIN
				u.plot.up[4]=float(data)
				plot_emiss,u
				IF u.plot.cont NE 0 THEN plot_cont,u
			END

			u.id.t_text : BEGIN
				u.stat.time=float(data)
				u.index=ipt(u.dat.t,u.stat.time)
				widget_control,u.id.t_slider,set_value=u.stat.time*5.0e4
				plot_bright,u
				plot_emiss,u
				IF u.stat.m1stat AND u.plot.cont EQ 0 THEN plot_cont,u
			END
			u.id.io_val : BEGIN
				widget_control,u.id.io_val,get_value=iostr
				io=float(strsplit(iostr,',',/extract))
				IF u.stat.dat THEN BEGIN	
					*u.plot.io=io
					plot_cont,u
				ENDIF

			END
			u.id.jo_val : BEGIN
				widget_control,u.id.jo_val,get_value=jostr
				jo=float(strsplit(jostr,',',/extract))
				IF u.stat.dat THEN BEGIN	
					*u.plot.jo=jo
					plot_cont,u
				ENDIF

			END
			u.id.t1_val : BEGIN
				widget_control,u.id.t1_val,get_value=t1str
				t1=float(strsplit(t1str,',',/extract))
				IF u.stat.dat THEN BEGIN	
					u.plot.t1=t1
				ENDIF
			END
			u.id.t2_val : BEGIN
				widget_control,u.id.t2_val,get_value=t2str
				t2=float(strsplit(t2str,',',/extract))
				IF u.stat.dat THEN BEGIN	
					u.plot.t2=t2
				ENDIF
			END
			ELSE: 
			ENDCASE
		END
		ELSE:
	ENDCASE
	IF button NE 'QUIT' THEN widget_control,event.top,set_uvalue=u		
END

;+
;NAME:
;	W_AXUV
;
;PURPOSE:
;	This procedure launches the W_AXUV widget which is used to
;	visualize the results of the SHELLINVERT inversions of the multiple
;	horizontally viewing AXUV diode arrays. 
;
;CALLING SEQUENCE:
;	@run_waxuv.pro 		will compile the proper files, set the color
;				tables and launch the widget
;PROCEUDRE:
;	See the C-Mod Wiki page on the W_AXUV widget for instructions
;	(listed under Shot Analysis Tools)
;
;RESTRICTIONS:
;	Due to the limited namespace, many of the GENIE widgets have
;	overlapping function/procedure names so these widgets may not run
;	well in contiguous IDL sessions.
;
;MODIFICATION HISTORY:
;	Written by:	M.L. Reinke - 2010
;
;-
PRO w_axuv,shot=shot,time=time,ptau=ptau,tree=tree
	IF NOT keyword_set(shot) THEN shot=1110105031
	IF NOT keyword_set(time) THEN time=1.0
	user=logname()
	loadct,12,/silent
	mdim=get_screen_size()
	IF mdim[0] LT 1600 OR  mdim[1] LT 1100 THEN base=widget_base(title='AXUV ASYMMETRY ANALYSIS',/row,tlb_size_events=1,/scroll,$
		x_scroll_size=mdim[0]*0.85,y_scroll_size=mdim[1]*0.85) ELSE base=widget_base(title='AXUV ASYMMETRY ANALYSIS',/row,tlb_size_events=1)
	A=widget_base(base,/column)
	B=widget_base(base,/column)

	xsA=600
	dum=widget_label(A,value='BRIGHTNESS PROFILES')
	A1=widget_base(A,frame=2)
	draw1=widget_draw(A1,xsize=xsA,ysize=570)
	


	AX=widget_base(A,/row)
	A2=widget_base(AX,/column,xsize=xsA*0.55,ysize=460,/frame)
	dum=widget_label(A2,value='SETUP')
	A2p1=widget_base(A2,/row)
	dum = widget_label(A2p1,value='SHOT: ')
	shotid = widget_text(A2p1,xsize=10,ysize=1,/edit)
	dum = widget_label(A2p1,value='')
	load= widget_button(A2p1,value='LOAD')
	dum = widget_label(A2p1,value='')
	quit= widget_button(A2p1,value='QUIT')
	dum = widget_label(A2p1,value='')
	print= widget_button(A2p1,value='PRINT')
	dum = widget_label(A2p1,value='')
	stop= widget_button(A2p1,value='STOP')
	A2p1x=widget_base(A2,/row)
	dum = widget_label(A2p1x,value='LOAD: ')	
	A2l=widget_base(A2p1x,/row, /nonexcl)	
	laxa=widget_button(A2l,value='AXA')
	laxj=widget_button(A2l,value='AXJ')
	lwb1=widget_button(A2l,value='WB1')
	lwb2=widget_button(A2l,value='WB2')
	lwb3=widget_button(A2l,value='WB3')
	lwb4=widget_button(A2l,value='WB4')
	A2p2=widget_base(A2,/row)
	message = widget_text(A2p2,xsize=40,ysize=5,/scroll)

	dum=widget_label(A2,value='CONTOUR OPTIONS')
	A2p3=widget_base(A2,/row)
	dum=widget_label(A2p3,value='CONTOUR: ')
	A2p3x=widget_base(A2p3,/row, /nonexcl)
	emr=widget_button(A2p3x,value='EMR')
	emc=widget_button(A2p3x,value='EMC')
	A2p4=widget_base(A2,/row)
	dum=widget_label(A2p4,value='io: ')
	io_val=widget_text(A2p4,xsize=35,ysize=1,/edit)
	A2p5=widget_base(A2,/row)
	dum=widget_label(A2p5,value='jo: ')
	jo_val=widget_text(A2p5,xsize=35,ysize=1,/edit)
	A2p6=widget_base(A2,/row)
	dum=widget_label(A2p6,value='t1: ')
	t1_val=widget_text(A2p6,xsize=10,ysize=1,/edit)
	dum=widget_label(A2p6,value=' t2: ')
	t2_val=widget_text(A2p6,xsize=10,ysize=1,/edit)


	
	;plotting options
	Px=widget_base(AX,/column,xsize=xsA*0.4,ysize=320,/frame)
	dum=widget_label(Px, value='PLOTTING OPTIONS')
	P1=widget_base(Px, /row)
	P1a=widget_base(P1,/row,/nonexcl)
	rtauto=widget_button(P1a,value='AUTO ')
	rtmin=widget_text(P1,xsize=5,ysize=1,/edit)
	dum=widget_label(P1,value=' < RT < ')
	rtmax=widget_text(P1,xsize=5,ysize=1,/edit)
	P2=widget_base(Px, /row)
	P2a=widget_base(P2,/row,/nonexcl)
	brauto=widget_button(P2a,value='AUTO ')
	brmin=widget_text(P2,xsize=5,ysize=1,/edit)
	dum=widget_label(P2,value=' < BR < ')
	brmax=widget_text(P2,xsize=5,ysize=1,/edit)
	P3=widget_base(Px, /row)
	P3a=widget_base(P3,/row,/nonexcl)
	rauto=widget_button(P3a,value='AUTO ')
	rmin=widget_text(P3,xsize=5,ysize=1,/edit)
	dum=widget_label(P3,value=' < R < ')
	rmax=widget_text(P3,xsize=5,ysize=1,/edit)
	P4=widget_base(Px, /row)
	P4a=widget_base(P4,/row,/nonexcl)
	rdauto=widget_button(P4a,value='AUTO ')
	rdmin=widget_text(P4,xsize=5,ysize=1,/edit)
	dum=widget_label(P4,value=' < RMID < ')
	rdmax=widget_text(P4,xsize=5,ysize=1,/edit)
	P5=widget_base(Px, /row)
	P5a=widget_base(P5,/row,/nonexcl)
	emauto=widget_button(P5a,value='AUTO ')
	emmin=widget_text(P5,xsize=5,ysize=1,/edit)
	dum=widget_label(P5,value=' < EM < ')
	emmax=widget_text(P5,xsize=5,ysize=1,/edit)
	P6=widget_base(Px, /row)
	P6a=widget_base(P6,/row,/nonexcl)
	tauto=widget_button(P6a,value='AUTO ')
	tmin=widget_text(P6,xsize=5,ysize=1,/edit)
	dum=widget_label(P6,value=' < t < ')
	tmax=widget_text(P6,xsize=5,ysize=1,/edit)
	P11=widget_base(Px, /row)
	P11a=widget_base(P11,/row,/nonexcl)
	ratauto=widget_button(P11a,value='AUTO ')
	ratmin=widget_text(P11,xsize=5,ysize=1,/edit)
	dum=widget_label(P11,value=' < RAT < ')
	ratmax=widget_text(P11,xsize=5,ysize=1,/edit)
	P9=widget_base(Px, /row)
	dum=widget_label(P9,value='ARRAY:  ')
	P9a=widget_base(P9,/row,/nonexcl)
	paxa=widget_button(P9a,value='AXA ')
	paxj=widget_button(P9a,value='AXJ ')
	pwb1=widget_button(P9a,value='WB1 ')
	P10=widget_base(Px, /row)
	dum=widget_label(P10,value='ARRAY:  ')
	P10a=widget_base(P10,/row,/nonexcl)
	pwb2=widget_button(P10a,value='WB2 ')
	pwb3=widget_button(P10a,value='WB3 ')
	pwb4=widget_button(P10a,value='WB4 ')
	P7=widget_base(Px, /row)
	dum=widget_label(P7,value='RADIAL:  ')
	P7a=widget_base(P7,/row,/nonexcl)
	xrmaj=widget_button(P7a,value='RMAJ ')
	xrho=widget_button(P7a,value='r/a ')
	P8=widget_base(Px, /row)
	dum=widget_label(P8,value='COSINE:  ')
	P8a=widget_base(P8,/row,/nonexcl)
	pcos=widget_button(P8a,value=' ')
	dum=widget_label(P8,value='RATIO:  ')
	P8b=widget_base(P8,/row,/nonexcl)
	prat=widget_button(P8b,value=' ')
	ptime=widget_button(P8,value='TIME')

	xsB=900
	dum=widget_label(B,value='RADIAL PROFILES')
	B1=widget_base(B,frame=2)
	draw2=widget_draw(B1,xsize=xsB+50,ysize=550)
	Bt=widget_base(B,/row)
	dum=widget_label(Bt, value='TIME: ')
	t_slider=widget_slider(Bt,xsize=xsB-75,min=0,max=10000,value=5000,/drag,/suppress)
	t_text=widget_text(Bt,xsize=8,ysize=1,/edit)
	B2=widget_base(B,/row)
	draw3=widget_draw(B2,xsize=550,ysize=350)
	B2a=widget_base(B2,/column)
	draw4=widget_draw(B2a,xsize=xsb-510,ysize=175)
	draw5=widget_draw(B2a,xsize=xsb-510,ysize=175)
	
	id={base:base,draw1:draw1,draw2:draw2,draw3:draw3,draw4:draw4,draw5:draw5,t_slider:t_slider,t_text:t_text,$
		shotid:shotid,load:load,quit:quit,print:print,stop:stop,message:message,$
		laxa:laxa,laxj:laxj,lwb1:lwb1,lwb2:lwb2,lwb3:lwb3,lwb4:lwb4,$
		paxa:paxa,paxj:paxj,pwb1:pwb1,pwb2:pwb2,pwb3:pwb3,pwb4:pwb4,$
		emr:emr,emc:emc,io_val:io_val,jo_val:jo_val,t1_val:t1_val,t2_val:t2_val,$
		rtauto:rtauto,brauto:brauto,rauto:rauto,rdauto:rdauto,emauto:emauto,tauto:tauto,ratauto:ratauto,$
		rtmax:rtmax,brmax:brmax,rmax:rmax,rdmax:rdmax,emmax:emmax,tmax:tmax,ratmax:ratmax,$
		rtmin:rtmin,brmin:brmin,rmin:rmin,rdmin:rdmin,emmin:emmin,tmin:tmin,ratmin:ratmin,$
		xrmaj:xrmaj,xrho:xrho,pcos:pcos,prat:prat,ptime:ptime}


	;plots ch,br,r,rmid,em,t,rat
	auto=[0,0,0,0,0,1,0]
	up=float([0.92,1.5,0.92,0.90,1.5,0,0.3])
	low=float([0.44,0,0.44,0.67,0,0,-0.1])
	csize=[550,350,390,175,390,175]
	bsize=[600,570]
	esize=[950,550]
	arr=['axa','axj','wb1ax','wb2ax','wb3ax','wb4ax']
	larr=[1,1,0,0,0,0]
	parr=[1,1,0,0,0,0]
	err=[0,0,0,0,0,0]	;keep track of if error is to be plotted
	IF NOT keyword_set(ptau) THEN ptau=fltarr(20)-1.0
	IF NOT keyword_set(tree) THEN tree='ANALYSIS'
	plot={auto:auto,up:up,low:low,cos:0,rho:0,cont:0,csize:csize,esize:esize,bsize:bsize,ratio:0,io:ptr_new([0],/allocate_heap),jo:ptr_new([0],/allocate_heap),t1:[0.0,0.0],t2:[0.0,0.0],$
		iostr:'',jostr:'',t1str:'0.0,0.0',t2str:'0.0,0.0',parr:parr,msym:[10,9,12,13,0,0],psym:[8,8,8,8,6,5],sty:[0,3,0,0,3,3,1],tau:ptau}
	stat={ps:0,dat:0,time:time,col:[255,255,50,100,200,120,255],pscol:[0,0,30,100,200,120,0],nr:0,larr:larr,arr:arr,err:err}
	u={id:id,shot:shot,stat:stat,plot:plot,tree:tree}
	widget_control,base,set_uvalue=u
	widget_control,u.id.shotid,set_value=num2str(u.shot,1)
	
	autoid=[u.id.rtauto,u.id.brauto,u.id.rauto,u.id.rdauto,u.id.emauto,u.id.tauto,u.id.ratauto]
	upid=[u.id.rtmax,u.id.brmax,u.id.rmax,u.id.rdmax,u.id.emmax,u.id.tmax,u.id.ratmax]
	lowid=[u.id.rtmin,u.id.brmin,u.id.rmin,u.id.rdmin,u.id.emmin,u.id.tmin,u.id.ratmin]
	FOR i=0,n(autoid) DO BEGIN
		IF u.plot.auto[i] THEN widget_control,autoid[i],set_button=1 ELSE BEGIN
			widget_control,upid[i],set_value=num2str(u.plot.up[i],dp=2)
			widget_control,lowid[i],set_value=num2str(u.plot.low[i],dp=2)
		ENDELSE
	ENDFOR
	
	IF u.plot.rho THEN widget_control,u.id.xrho,set_button=1 ELSE widget_control,u.id.xrmaj,set_button=1
	CASE u.plot.cont OF
		1 : widget_control,u.id.emr,set_button=1
		2 : widget_control,u.id.emc,set_button=1
		ELSE : 
	ENDCASE
	ids=[u.id.laxa,u.id.laxj,u.id.lwb1,u.id.lwb2,u.id.lwb3,u.id.lwb4]
	FOR i=0,n(ids) DO IF u.stat.larr[i] THEN widget_control,ids[i],set_button=1
	ids=[u.id.paxa,u.id.paxj,u.id.pwb1,u.id.pwb2,u.id.pwb3,u.id.pwb4]
	FOR i=0,n(ids) DO IF u.stat.larr[i] THEN widget_control,ids[i],set_button=1

	widget_control,u.id.t1_val,set_value=u.plot.t1str		
	widget_control,u.id.t2_val,set_value=u.plot.t2str		


	!except=0
	widget_control,base,/realize
	xmanager,'w_axuv',base
END
